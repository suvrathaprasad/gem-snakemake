#!/usr/bin/env python3
"""
Pull the rha/glc operon reference genes out of the full E. coli K-12
MG1655 RefSeq proteome (GCF_000005845.2), using feature_table.txt to
map gene symbol -> protein accession, then extracting those sequences
from protein.faa.
"""

import argparse

GENES = [
    "rhaA", "rhaB", "rhaD", "rhaT", "rhaR", "rhaS", "rhaM",
    "glcA", "glcB", "glcC", "glcD", "glcE", "glcF", "glcG",
]


def load_gene_to_protein(feature_table_path):
    gene_to_protein = {}
    with open(feature_table_path) as f:
        header = f.readline().lstrip("#").strip().split("\t")
        symbol_idx = header.index("symbol")
        protein_idx = header.index("product_accession")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            symbol = fields[symbol_idx]
            protein_acc = fields[protein_idx]
            if symbol in GENES and protein_acc:
                gene_to_protein[symbol] = protein_acc
    return gene_to_protein


def extract_sequences(protein_faa_path, wanted_accessions):
    sequences = {}
    current_acc = None
    current_seq = []
    with open(protein_faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_acc in wanted_accessions:
                    sequences[current_acc] = "".join(current_seq)
                current_acc = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_acc in wanted_accessions:
            sequences[current_acc] = "".join(current_seq)
    return sequences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature-table", required=True)
    parser.add_argument("--protein-faa", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    gene_to_protein = load_gene_to_protein(args.feature_table)

    missing_genes = set(GENES) - set(gene_to_protein)
    if missing_genes:
        print(f"[warn] gene symbols not found in feature table: {sorted(missing_genes)}")

    accessions = set(gene_to_protein.values())
    sequences = extract_sequences(args.protein_faa, accessions)

    with open(args.output, "w") as out:
        for gene, acc in gene_to_protein.items():
            seq = sequences.get(acc)
            if seq is None:
                print(f"[warn] {gene} ({acc}) had no matching sequence in protein.faa")
                continue
            out.write(f">{gene}|{acc}\n{seq}\n")

    print(f"Wrote {len(sequences)} sequences to {args.output}")


if __name__ == "__main__":
    main()
