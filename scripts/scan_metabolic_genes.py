#!/usr/bin/env python3
"""
Genome-wide scan of E. coli K-12 iML1515 metabolic genes against one strain's
own proteome. Takes the reciprocal-best-hit map from ortholog_map.py (each
iML1515 reference protein's best hit in the strain's Bakta proteins, and
whether that hit is reciprocal) and calls each reference gene intact,
disrupted or absent, then attaches the reactions and subsystems each gene
supports in iML1515. A best hit that is not reciprocal (the strain protein
really belongs to another K-12 gene, e.g. a paralog) is called absent.

Generalises interpret_operon_hits.py from 14 genes to ~1,500. The intact
cutoffs are the same ones used there (>=80% query coverage and >=70% identity).
Two things are new at genome scale, and both are arguments so they can be
changed without editing the script:

  --min-identity / --min-coverage (default 40 / 50): a best hit below EITHER
      floor is treated as absent, not disrupted. Without a floor, a distant
      paralog (many metabolic genes have isozymes, and paralogs often align
      full-length at 30-40% identity) would be reported as a "disrupted" copy
      of a gene the strain has simply lost. 40/50 are the gene-presence
      thresholds used by Barth et al. 2020 (Section 4.3). The hit's identity
      and coverage are still printed, and `detail` says which floor failed, so
      a heavily truncated but otherwise identical gene is not hidden.
  n_refs_sharing_target: how many reference genes picked the same strain
      protein as their best hit (forward). With the reciprocal check, a
      shared protein is only an ortholog of one of them; the others are
      called absent. ortholog_target is the strain protein only when the hit
      is reciprocal.

Thresholds are a starting heuristic, not validated cutoffs.
"""
import argparse
import collections
import csv

COVERAGE_INTACT_THRESHOLD = 80.0  # % of reference gene length aligned
IDENTITY_INTACT_THRESHOLD = 70.0  # % identity over the aligned region


def load_gene_table(path):
    genes = collections.OrderedDict()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes[row["gene"]] = row
    return genes


def load_query_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])  # e.g. b0001|thrL|NP_414542.1
    return ids


def load_hits(map_path):
    """Read ortholog_map.tsv -> {gene: hit}; genes without a forward hit are left out."""
    hits = {}
    with open(map_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["best_hit_target"] == "-":
                continue
            hits[row["gene"]] = {
                "target": row["best_hit_target"],
                "pident": float(row["pident"]),
                "qcovhsp": float(row["qcovhsp"]),
                "evalue": row["evalue"],
                "bitscore": row["bitscore"],
                "rbh": row["rbh"],
                "reverse_gene": row["reverse_best_gene"],
            }
    return hits


def classify(hit, min_identity, min_coverage):
    """Return (call, detail)."""
    if hit is None:
        return "absent", "no hit"
    if hit["pident"] < min_identity:
        return "absent", f"weak hit only (identity < {min_identity:g}%)"
    if hit["qcovhsp"] < min_coverage:
        return "absent", f"weak hit only (coverage < {min_coverage:g}%)"
    if hit["rbh"] != "yes":
        return "absent", f"best hit is not an ortholog (reciprocal best hit: {hit['reverse_gene']})"
    cov_ok = hit["qcovhsp"] >= COVERAGE_INTACT_THRESHOLD
    id_ok = hit["pident"] >= IDENTITY_INTACT_THRESHOLD
    if cov_ok and id_ok:
        return "intact", "-"
    if id_ok:
        return "disrupted", "truncated (high identity, low coverage)"
    if cov_ok:
        return "disrupted", "divergent (full length, low identity)"
    return "disrupted", "truncated and divergent"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True, help="files/iML1515_metabolic_genes.faa")
    parser.add_argument("--gene-table", required=True, help="files/iML1515_gene_table.tsv")
    parser.add_argument("--ortholog-map", required=True, help="{sample}.ortholog_map.tsv from ortholog_map.py")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--min-identity", type=float, default=40.0)
    parser.add_argument("--min-coverage", type=float, default=50.0)
    args = parser.parse_args()

    table = load_gene_table(args.gene_table)
    hits = load_hits(args.ortholog_map)

    target_counts = collections.Counter(h["target"] for h in hits.values())

    header = ["gene", "symbol", "call", "detail", "best_hit_target", "pident", "qcovhsp",
              "evalue", "bitscore", "n_refs_sharing_target", "ortholog_target", "rbh",
              "reverse_best_gene", "n_reactions", "subsystems", "reactions"]
    rows = []
    for qid in load_query_ids(args.query_fasta):
        gene = qid.split("|")[0]
        info = table.get(gene, {})
        hit = hits.get(gene)
        call, detail = classify(hit, args.min_identity, args.min_coverage)
        rows.append([
            gene, info.get("symbol", "-"), call, detail,
            hit["target"] if hit else "-",
            f"{hit['pident']:.1f}" if hit else "-",
            f"{hit['qcovhsp']:.1f}" if hit else "-",
            hit["evalue"] if hit else "-",
            hit["bitscore"] if hit else "-",
            str(target_counts[hit["target"]]) if hit else "-",
            hit["target"] if hit and hit["rbh"] == "yes" else "-",
            hit["rbh"] if hit else "-",
            hit["reverse_gene"] if hit else "-",
            info.get("n_reactions", "-"), info.get("subsystems", "-"), info.get("reactions", "-"),
        ])

    # disrupted/absent first so the genes worth looking at sit at the top
    order = {"disrupted": 0, "absent": 1, "intact": 2}
    rows.sort(key=lambda r: (order[r[2]], r[0]))

    with open(args.output, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    counts = collections.Counter(r[2] for r in rows)
    print(f"scanned {len(rows)} genes: " + ", ".join(f"{counts[c]} {c}" for c in ("intact", "disrupted", "absent")))


if __name__ == "__main__":
    main()
