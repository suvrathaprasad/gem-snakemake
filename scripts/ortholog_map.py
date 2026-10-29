#!/usr/bin/env python3
"""
ortholog_map.py — map every reference protein (the iML1515 metabolic genes)
to its reciprocal best hit in one strain's proteome.

Adapted from get_rbh.py (dN-dS-Snakemake-pipeline). Same DIAMOND forward /
reverse search and the same reciprocal check, with these changes for the
genome scan:
  * ALL reference proteins are written out, not just the reciprocal pairs,
    with an `rbh` flag (yes / no / - for no hit), so the scan can tell
    "absent" from "best hit belongs to another gene".
  * No coverage filter at parse time. The scan needs truncated hits
    (down to its 50% coverage floor), so the floors are applied later by
    scan_metabolic_genes.py. Coverage is coverage of the reference protein.
  * The reverse search runs against the WHOLE K-12 proteome, not only the
    metabolic subset, so a strain protein whose real counterpart is a
    non-metabolic K-12 gene is not wrongly assigned to a metabolic one.
  * The reverse best hit's gene name and a score ratio are reported, so
    near-ties between close paralogs (e.g. gadA/gadB) are visible.
  * DIAMOND only (BLAST+ path dropped), --sensitive as in the scan.

Steps:
  1. forward: reference proteins -> strain proteins
  2. reverse: strain proteins -> whole K-12 proteome
  3. best hit per query by bitscore; a reference gene is a reciprocal hit
     if the reverse best hit of its forward target is the reference itself
  4. write one row per reference gene

Output columns: gene, ref_acc, best_hit_target, pident, qcovhsp, evalue,
bitscore, reverse_best_acc, reverse_best_gene, reverse_score_ratio, rbh
(reverse_score_ratio = reverse bitscore of the reference / best reverse
bitscore; 1.0 for a strict reciprocal hit). Two K-12 genes that are
near-identical (e.g. gadA/gadB, 98%) cannot be told apart by a strain
protein, so a reference whose reverse score is within --rbh-min-ratio
(default 0.95) of the best counts as reciprocal, otherwise a present gene
would be called absent. Genuine paralog losses scored well below that in
the test (argF 0.89, acrF 0.78, glcA 0.65). 0.95 is a starting heuristic,
not a validated cutoff; use 1.0 for the strict version.)
"""
import argparse
import collections
import os
import subprocess
import sys
from pathlib import Path

OUTFMT = "6 qseqid sseqid evalue bitscore pident qstart qend qlen sstart send slen"


def make_diamond_db(faa, db_path, threads):
    cmd = ["diamond", "makedb", "--in", faa, "--db", db_path,
           "--threads", str(threads), "--quiet"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"diamond makedb failed for {faa}:\n{result.stderr}")


def run_diamond(query, db, out, evalue, threads):
    cmd = ["diamond", "blastp", "--query", query, "--db", db, "--out", out,
           "--evalue", str(evalue), "--threads", str(threads),
           "--outfmt", *OUTFMT.split(),
           "--sensitive", "--quiet"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"diamond blastp failed:\n{result.stderr}")


def parse_hits(path):
    """best[query] = (subject, evalue, bitscore, pident, qcov); scores[query][subject] = bitscore."""
    best, scores = {}, collections.defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            qid, sid, evalue, bitscore, pident, qstart, qend, qlen = parts[:8]
            bitscore, qstart, qend, qlen = float(bitscore), int(qstart), int(qend), int(qlen)
            qcov = (qend - qstart + 1) / qlen * 100 if qlen > 0 else 0.0
            if bitscore > scores[qid].get(sid, 0.0):
                scores[qid][sid] = bitscore
            if qid not in best or bitscore > best[qid][2]:
                best[qid] = (sid, evalue, bitscore, float(pident), qcov)
    return best, scores


def load_reference_ids(faa):
    """Headers look like b0002|thrA|NP_414543.1 -> (gene, accession)."""
    ids = []
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].split()[0]
                parts = header.split("|")
                ids.append((header, parts[0], parts[2] if len(parts) > 2 else ""))
    return ids


def load_acc_to_gene(feature_table):
    """protein accession -> 'symbol' (or locus tag) from the RefSeq feature table."""
    acc2gene = {}
    if not feature_table:
        return acc2gene
    with open(feature_table) as fh:
        header = fh.readline().lstrip("#").strip().split("\t")
        pi, si, li = (header.index(c) for c in ("product_accession", "symbol", "locus_tag"))
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if len(x) > max(pi, si, li) and x[pi]:
                acc2gene[x[pi]] = x[si] or x[li]
    return acc2gene


def main():
    ap = argparse.ArgumentParser(description="Reciprocal-best-hit map of reference proteins in one strain.")
    ap.add_argument("--query", required=True, help="reference proteins (files/iML1515_metabolic_genes.faa)")
    ap.add_argument("--reference-proteome", required=True, help="whole K-12 proteome (files/GCF_000005845.2_ASM584v2_protein.faa)")
    ap.add_argument("--target", required=True, help="strain proteins (Bakta .faa)")
    ap.add_argument("--feature-table", help="K-12 feature table, to name the reverse best hit")
    ap.add_argument("--out", required=True)
    ap.add_argument("--db_dir", required=True)
    ap.add_argument("--fwd_hits", required=True, help="forward DIAMOND output (kept)")
    ap.add_argument("--rev_hits", required=True, help="reverse DIAMOND output (kept)")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--rbh-min-ratio", type=float, default=0.95)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    for f in (args.query, args.reference_proteome, args.target):
        if not os.path.isfile(f):
            sys.exit(f"Input file not found: {f}")
    for p in (args.out, args.fwd_hits, args.rev_hits):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
    db_dir = Path(args.db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    strain_db, k12_db = str(db_dir / "strain_db"), str(db_dir / "k12_db")
    make_diamond_db(args.target, strain_db, args.threads)
    make_diamond_db(args.reference_proteome, k12_db, args.threads)
    run_diamond(args.query, strain_db, args.fwd_hits, args.evalue, args.threads)
    run_diamond(args.target, k12_db, args.rev_hits, args.evalue, args.threads)

    fwd, _ = parse_hits(args.fwd_hits)
    rev, rev_scores = parse_hits(args.rev_hits)
    acc2gene = load_acc_to_gene(args.feature_table)

    header = ["gene", "ref_acc", "best_hit_target", "pident", "qcovhsp", "evalue", "bitscore",
              "reverse_best_acc", "reverse_best_gene", "reverse_score_ratio", "rbh"]
    counts = collections.Counter()
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for qid, gene, acc in load_reference_ids(args.query):
            hit = fwd.get(qid)
            if hit is None:
                out.write("\t".join([gene, acc] + ["-"] * 8 + ["-"]) + "\n")
                counts["no hit"] += 1
                continue
            target, evalue, bitscore, pident, qcov = hit
            rb = rev.get(target)
            if rb is None:
                rb_acc, rb_gene, ratio, rbh = "-", "none", "-", "no"
            else:
                rb_acc = rb[0]
                rb_gene = acc2gene.get(rb_acc, rb_acc)
                ratio = rev_scores[target].get(acc, 0.0) / rb[2] if rb[2] > 0 else 0.0
                rbh = "yes" if (rb_acc == acc or ratio >= args.rbh_min_ratio) else "no"
                ratio = f"{ratio:.3f}"
            counts[rbh] += 1
            out.write("\t".join([gene, acc, target, f"{pident:.1f}", f"{qcov:.1f}", evalue,
                                 f"{bitscore:.1f}", rb_acc, rb_gene, str(ratio), rbh]) + "\n")
    print(f"ortholog map: {sum(counts.values())} reference genes | reciprocal: {counts['yes']} | "
          f"best hit is another gene: {counts['no']} | no hit: {counts['no hit']}")


if __name__ == "__main__":
    main()
