#!/usr/bin/env python3
"""Search reference proteins directly against the assembly (a tblastn equivalent)
and report where each gene sits, so a gene the annotation missed, a frameshifted
pseudogene, and a true absence can be told apart.

How:
  1. translate every contig in all six frames and split at stop codons
  2. DIAMOND blastp: reference proteins -> translated segments
  3. group the hits of one query that lie close together on the same contig and
     strand into a locus (so a frameshifted gene, which hits in two frames, is one
     locus with n_hsps > 1 and frames listed)
  4. report the best loci per query with identity, query coverage and the Bakta
     feature each locus overlaps (needs --gff3)
  5. optionally list every annotated feature around the loci of --window-genes,
     flagging mobile-element words (transposase, IS, insertion, integrase, traT ...)

Standard bacterial code; stdlib + the diamond binary only.

Outputs:
  --out-loci          one row per (query, locus), top --top-n per query:
                      query, rank, contig, strand, nt_start, nt_end, n_hsps, frames,
                      pident, qcov, bitscore, evalue, contig_len, nt_to_contig_end,
                      overlapping_feature (nt_to_contig_end = distance from the locus to the
                      nearer contig end: a gene "missing" next to a contig end may just be
                      on the next contig)
  --out-neighbourhood one row per annotated feature in the window (optional)
"""
import argparse
import bisect
import collections
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AAS[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
MOBILE = re.compile(r"(?i:transposase|insertion|integrase|recombinase|phage|mobile)|\bIS\d+|\bIS[A-Z]|\btra[A-Z]\b|\bTn\d+")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def translate(dna):
    dna = dna.upper()
    return "".join(CODON.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def segments(contig, dna, min_aa):
    """yield (id, protein, meta) for every stop-free stretch >= min_aa in six frames."""
    length = len(dna)
    rev = dna.translate(COMP)[::-1]
    for strand, seq in (("+", dna), ("-", rev)):
        for frame in range(3):
            prot = translate(seq[frame:])
            pos = 0
            for piece in prot.split("*"):
                if len(piece) >= min_aa:
                    # aa i of the piece is codon (pos + i) of this frame
                    yield (f"{contig}|{strand}{frame}|{pos}", piece,
                           (contig, strand, frame, pos, length))
                pos += len(piece) + 1


def aa_to_nt(meta, aa_start, aa_end):
    """1-based aa coordinates inside a segment -> 1-based nt coordinates (forward strand)."""
    contig, strand, frame, pos, length = meta
    r0 = frame + 3 * (pos + aa_start - 1)          # 0-based start in the (possibly reversed) sequence
    r1 = frame + 3 * (pos + aa_end)                # exclusive end
    if strand == "+":
        return r0 + 1, r1
    return length - r1 + 1, length - r0


def parse_gff3(path):
    """return (features, regions) where regions maps gff3 seqid -> length from ##sequence-region."""
    feats, regions = [], {}
    for line in open(path):
        if line.startswith("##sequence-region"):
            parts = line.split()
            if len(parts) >= 4:
                regions[parts[1]] = int(parts[3])
            continue
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith(">"):
            break
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] in ("region", "gene"):
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        feats.append(dict(contig=f[0], type=f[2], start=int(f[3]), end=int(f[4]), strand=f[6],
                          locus_tag=attrs.get("locus_tag", attrs.get("ID", "")),
                          gene=attrs.get("gene", attrs.get("Name", "")),
                          product=unquote(attrs.get("product", ""))))
    return feats, regions


def map_contigs(feats, regions, fasta_len):
    """Bakta renames contigs (contig_1, contig_2 ...) unless told otherwise, so gff3 seqids may not match
    the assembly's names. Map gff3 ids to assembly ids by name, else by sequence length, else by order."""
    gff_ids = list(dict.fromkeys(f["contig"] for f in feats))
    if all(g in fasta_len for g in gff_ids):
        return {g: g for g in gff_ids}, "gff3 contig ids match the assembly"
    by_len = collections.defaultdict(list)
    for name, n in fasta_len.items():
        by_len[n].append(name)
    mapping, how = {}, "by sequence length (##sequence-region)"
    for g in gff_ids:
        cands = by_len.get(regions.get(g, -1), [])
        if len(cands) == 1:
            mapping[g] = cands[0]
    if len(mapping) < len(gff_ids):
        names = list(fasta_len)
        mapping, how = {}, "by contig order (contig_N = N-th contig of the assembly)"
        for g in gff_ids:
            tail = g.rsplit("_", 1)[-1]
            if tail.isdigit() and 1 <= int(tail) <= len(names):
                mapping[g] = names[int(tail) - 1]
    if len(mapping) < len(gff_ids):
        sys.exit("cannot map the gff3 contig ids onto the assembly (names, lengths and order all fail); "
                 f"gff3 ids start {gff_ids[:3]}, assembly ids start {list(fasta_len)[:3]}")
    for f in feats:  # sanity: every feature must fit inside the contig it was mapped to
        if f["end"] > fasta_len[mapping[f["contig"]]]:
            sys.exit(f"mapped contig {mapping[f['contig']]} is shorter than feature {f['locus_tag']}; mapping is wrong")
    return mapping, how


def overlap_label(feats_by_contig, contig, start, end):
    hits = [f for f in feats_by_contig.get(contig, []) if f["start"] <= end and f["end"] >= start]
    if not hits:
        return "no annotated feature"
    return "; ".join(f"{f['locus_tag']}|{f['gene'] or '-'}|{f['product'] or '-'}" for f in hits)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="reference proteins (headers 'gene|accession')")
    ap.add_argument("--assembly", required=True, help="assembly fasta (contigs)")
    ap.add_argument("--gff3", help="Bakta gff3 for the same assembly")
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-neighbourhood")
    ap.add_argument("--window-genes", default="", help="comma list of query gene names that define the window")
    ap.add_argument("--flank", type=int, default=5000, help="nt added each side of the window (default 5000)")
    ap.add_argument("--db_dir", required=True, help="working dir for the segment fasta, db and raw hits")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-aa", type=int, default=30, help="shortest stop-free segment searched")
    ap.add_argument("--merge-gap", type=int, default=3000, help="hits closer than this (nt) form one locus")
    ap.add_argument("--top-n", type=int, default=3)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    for f in (args.query, args.assembly) + ((args.gff3,) if args.gff3 else ()):
        if not os.path.isfile(f):
            sys.exit(f"Input file not found: {f}")
    work = Path(args.db_dir)
    work.mkdir(parents=True, exist_ok=True)
    for p in (args.out_loci, args.out_neighbourhood):
        if p:
            Path(p).parent.mkdir(parents=True, exist_ok=True)

    # 1. six-frame segments
    meta, seg_fa = {}, work / "segments.faa"
    n_contigs, fasta_len = 0, {}
    with open(seg_fa, "w") as fo:
        for contig, dna in read_fasta(args.assembly):
            n_contigs += 1
            fasta_len[contig] = len(dna)
            for sid, prot, m in segments(contig, dna, args.min_aa):
                meta[sid] = m
                fo.write(f">{sid}\n{prot}\n")
    print(f"{n_contigs} contigs, {len(meta)} stop-free segments >= {args.min_aa} aa")

    qlen = {name: len(seq) for name, seq in read_fasta(args.query)}
    qgene = {name: name.split("|")[0] for name in qlen}

    # 2. diamond
    db, raw = str(work / "segments_db"), work / "locus_hits.raw.tsv"
    subprocess.run(["diamond", "makedb", "--in", str(seg_fa), "-d", db, "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "--query", args.query, "--db", db, "--out", str(raw),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qstart", "qend",
                    "sstart", "send", "evalue", "bitscore",
                    "--sensitive", "--evalue", str(args.evalue), "--max-target-seqs", "100",
                    "--max-hsps", "0", "--threads", str(args.threads), "--quiet"], check=True)

    # 3. group HSPs into loci
    hsps = collections.defaultdict(list)
    for line in open(raw):
        q, s, pid, length, qs, qe, ss, se, ev, bs = line.rstrip("\n").split("\t")
        nt0, nt1 = aa_to_nt(meta[s], int(ss), int(se))
        contig, strand, frame = meta[s][0], meta[s][1], meta[s][2]
        hsps[(q, contig, strand)].append(dict(nt0=nt0, nt1=nt1, qs=int(qs), qe=int(qe), pid=float(pid),
                                              length=int(length), ev=float(ev), bs=float(bs),
                                              frame=f"{strand}{frame}"))
    loci = collections.defaultdict(list)
    for (q, contig, strand), hs in hsps.items():
        hs.sort(key=lambda h: h["nt0"])
        cur = [hs[0]]
        for h in hs[1:]:
            covered = {i for x in cur for i in range(x["qs"], x["qe"] + 1)}
            shared = len(covered & set(range(h["qs"], h["qe"] + 1)))
            same_gene = shared <= 0.2 * min(len(covered), h["qe"] - h["qs"] + 1)
            # same locus: close in nt AND covering a different part of the query
            # (a frameshifted gene hits in pieces; a neighbouring paralog hits the same part)
            if h["nt0"] - max(x["nt1"] for x in cur) <= args.merge_gap and same_gene:
                cur.append(h)
            else:
                loci[q].append((contig, strand, cur)); cur = [h]
        loci[q].append((contig, strand, cur))

    feats_by_contig = collections.defaultdict(list)
    feats = []
    if args.gff3:
        feats, regions = parse_gff3(args.gff3)
        if not feats:
            sys.exit(f"no features read from {args.gff3}")
        mapping, how = map_contigs(feats, regions, fasta_len)
        print(f"gff3: {len(feats)} features; {how}")
        for f in feats:
            f["contig"] = mapping[f["contig"]]       # from here on, assembly contig names
            feats_by_contig[f["contig"]].append(f)

    rows, found = [], collections.defaultdict(list)
    for q in qlen:
        scored = []
        for contig, strand, hs in loci.get(q, []):
            covered = set()
            for h in hs:
                covered.update(range(h["qs"], h["qe"] + 1))
            tot = sum(h["length"] for h in hs)
            scored.append((sum(h["bs"] for h in hs), contig, strand, hs, len(covered) / qlen[q],
                           sum(h["pid"] * h["length"] for h in hs) / tot))
        scored.sort(key=lambda x: -x[0])
        if not scored:
            rows.append([q, "-", "-", "-", "-", "-", "0", "-", "-", "-", "-", "-", "-", "-", "no hit on the assembly"])
            continue
        for rank, (bs, contig, strand, hs, cov, pid) in enumerate(scored[:args.top_n], 1):
            n0, n1 = min(h["nt0"] for h in hs), max(h["nt1"] for h in hs)
            frames = ",".join(sorted({h["frame"] for h in hs}))
            ev = min(h["ev"] for h in hs)
            rows.append([q, str(rank), contig, strand, str(n0), str(n1), str(len(hs)), frames,
                         f"{pid:.1f}", f"{100 * cov:.1f}", f"{bs:.0f}", f"{ev:.2e}",
                         str(fasta_len[contig]), str(min(n0 - 1, fasta_len[contig] - n1)),
                         overlap_label(feats_by_contig, contig, n0, n1) if feats else "(no gff3 given)"])
            found[qgene[q]].append((contig, n0, n1))

    with open(args.out_loci, "w") as out:
        out.write("\t".join(["query", "rank", "contig", "strand", "nt_start", "nt_end", "n_hsps", "frames",
                             "pident", "qcov", "bitscore", "evalue", "contig_len", "nt_to_contig_end",
                             "overlapping_feature"]) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"wrote {args.out_loci} ({len(rows)} rows)")

    # 5. neighbourhood
    if args.out_neighbourhood:
        if not feats:
            sys.exit("--out-neighbourhood needs --gff3")
        wg = [g for g in args.window_genes.split(",") if g]
        spans = [s for g in wg for s in found.get(g, [])[:1]]   # best locus of each window gene
        by_contig = collections.defaultdict(list)
        for contig, a, b in spans:
            by_contig[contig].append((a, b))
        if not by_contig:
            sys.exit("none of --window-genes was found on the assembly; no window to list")
        marks = collections.defaultdict(list)
        for g, lst in found.items():
            for contig, a, b in lst[:1]:
                marks[contig].append((a, b, g))
        n_listed = 0
        with open(args.out_neighbourhood, "w") as out:
            out.write("\t".join(["contig", "start", "end", "strand", "type", "locus_tag", "gene", "product",
                                 "query_hit", "mobile_element_word"]) + "\n")
            for contig, ab in by_contig.items():
                lo, hi = min(a for a, _ in ab) - args.flank, max(b for _, b in ab) + args.flank
                for f in sorted(feats_by_contig.get(contig, []), key=lambda x: x["start"]):
                    if f["end"] < lo or f["start"] > hi:
                        continue
                    hit = ",".join(g for a, b, g in marks[contig] if a <= f["end"] and b >= f["start"]) or "-"
                    mob = "yes" if MOBILE.search(f"{f['gene']} {f['product']}") else "-"
                    out.write("\t".join([contig, str(f["start"]), str(f["end"]), f["strand"], f["type"],
                                         f["locus_tag"], f["gene"] or "-", f["product"] or "-", hit, mob]) + "\n")
                    n_listed += 1
        if n_listed == 0:
            sys.exit("window is empty: no annotated feature lies near the hits (contig naming or coordinates differ from the gff3)")
        print(f"wrote {args.out_neighbourhood}")


if __name__ == "__main__":
    main()
