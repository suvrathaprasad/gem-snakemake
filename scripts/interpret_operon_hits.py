#!/usr/bin/env python3
"""
Turn raw diamond blastp hits (reference gene vs. this strain's own
proteome) into a three-way call per gene: intact / disrupted / absent.

Thresholds below are a starting heuristic, not a validated cutoffs.
"""
import argparse
import csv

COVERAGE_INTACT_THRESHOLD = 80.0  # % of reference gene length aligned
IDENTITY_INTACT_THRESHOLD = 70.0  # % identity over the aligned region


def load_query_gene_names(fasta_path):
    genes = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                genes.append(line[1:].split()[0])
    return genes


def load_hits(raw_hits_path):
    hits = {}
    with open(raw_hits_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            qseqid, sseqid, pident, qcovhsp, evalue, bitscore = row
            hits[qseqid] = {
                "target": sseqid,
                "pident": float(pident),
                "qcovhsp": float(qcovhsp),
                "evalue": evalue,
                "bitscore": bitscore,
            }
    return hits


def classify(hit):
    if hit is None:
        return "absent"
    if hit["qcovhsp"] >= COVERAGE_INTACT_THRESHOLD and hit["pident"] >= IDENTITY_INTACT_THRESHOLD:
        return "intact"
    return "disrupted (partial/low-identity hit)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-fasta", required=True, help="reference/rha_glc_genes.faa")
    parser.add_argument("--raw-hits", required=True, help="diamond blastp outfmt 6 output")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    genes = load_query_gene_names(args.query_fasta)
    hits = load_hits(args.raw_hits)

    header = ["gene", "call", "best_hit_target", "pident", "qcovhsp", "evalue", "bitscore"]
    rows = []
    for gene in genes:
        hit = hits.get(gene)
        call = classify(hit)
        if hit:
            rows.append((gene, call, hit["target"], f"{hit['pident']:.1f}",
                         f"{hit['qcovhsp']:.1f}", hit["evalue"], hit["bitscore"]))
        else:
            rows.append((gene, call, "-", "-", "-", "-", "-"))

    with open(args.output, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
            print("\t".join(row))


if __name__ == "__main__":
    main()
