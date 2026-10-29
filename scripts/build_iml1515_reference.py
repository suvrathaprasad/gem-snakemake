#!/usr/bin/env python3
"""
One-off: build the reference inputs for the genome-wide metabolic gene scan
from E. coli K-12 iML1515 (BiGG JSON) and the MG1655 RefSeq proteome
(GCF_000005845.2).

iML1515 gene IDs are MG1655 locus tags (b-numbers). The RefSeq feature table
maps each locus tag to its protein accession, so every iML1515 gene with a
protein in the feature table gets its sequence pulled from protein.faa.

Writes two files:
  - FASTA of the reference proteins, header ">b0001|thrL|NP_414542.1"
  - TSV with one row per iML1515 gene (scanned or excluded): symbol, protein,
    reaction IDs, subsystems, and a reference_status column that records why a
    gene was excluded (for example the spontaneous pseudo-gene s0001, or a
    locus tag with no protein in the feature table).

Uses only the standard library (no cobra needed).
"""
import argparse
import csv
import gzip
import json
import re


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_iml1515(json_path):
    with open_text(json_path) as f:
        model = json.load(f)
    genes = {g["id"]: g.get("name", "") for g in model["genes"]}

    gene_reactions = {gid: [] for gid in genes}
    gene_subsystems = {gid: set() for gid in genes}
    for rxn in model["reactions"]:
        # the gene_reaction_rule is e.g. "b0001 and (b0002 or b0003)"
        tokens = set(re.findall(r"[A-Za-z0-9_]+", rxn.get("gene_reaction_rule", "")))
        for gid in tokens & genes.keys():
            gene_reactions[gid].append(rxn["id"])
            if rxn.get("subsystem"):
                gene_subsystems[gid].add(rxn["subsystem"])
    return genes, gene_reactions, gene_subsystems


def load_locus_to_protein(feature_table_path):
    locus_to_protein = {}
    with open(feature_table_path) as f:
        header = f.readline().lstrip("#").strip().split("\t")
        locus_idx = header.index("locus_tag")
        protein_idx = header.index("product_accession")
        symbol_idx = header.index("symbol")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(locus_idx, protein_idx, symbol_idx):
                continue
            if fields[protein_idx]:  # CDS rows with a protein product
                locus = fields[locus_idx]
                """
                a locus tag can list several products (e.g. b0484 copA: the 834 aa
                NP_ protein and a 70 aa ribosomal-slippage YP_ product): keep the
                canonical NP_ accession, otherwise the first one seen
                """
                if locus in locus_to_protein and (
                    locus_to_protein[locus][0].startswith("NP_")
                    or not fields[protein_idx].startswith("NP_")
                ):
                    continue
                locus_to_protein[fields[locus_idx]] = (fields[protein_idx], fields[symbol_idx])
    return locus_to_protein


def load_sequences(protein_faa_path, wanted):
    sequences, acc, seq = {}, None, []
    with open(protein_faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if acc in wanted:
                    sequences[acc] = "".join(seq)
                acc, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if acc in wanted:
        sequences[acc] = "".join(seq)
    return sequences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iml1515-json", required=True, help="iML1515.json (or .json.gz) from BiGG")
    parser.add_argument("--feature-table", required=True)
    parser.add_argument("--protein-faa", required=True)
    parser.add_argument("--out-faa", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    genes, gene_reactions, gene_subsystems = load_iml1515(args.iml1515_json)
    locus_to_protein = load_locus_to_protein(args.feature_table)
    sequences = load_sequences(
        args.protein_faa, {acc for acc, _ in locus_to_protein.values()}
    )

    n_written, excluded = 0, []
    with open(args.out_faa, "w") as faa, open(args.out_table, "w") as tsv:
        writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "symbol", "protein_acc", "n_reactions", "reactions",
                         "subsystems", "reference_status"])
        for gid in sorted(genes):
            acc, symbol = locus_to_protein.get(gid, (None, genes[gid]))
            if gid.startswith("s"):
                status = "excluded: spontaneous-reaction pseudo-gene"
            elif acc is None:
                status = "excluded: no protein for this locus tag in the feature table"
            elif acc not in sequences:
                status = "excluded: protein accession not found in protein.faa"
            else:
                status = "scanned"
                faa.write(f">{gid}|{symbol or genes[gid]}|{acc}\n{sequences[acc]}\n")
                n_written += 1
            if status != "scanned":
                excluded.append((gid, status))
            writer.writerow([
                gid, symbol or genes[gid], acc or "-",
                len(gene_reactions[gid]),
                ";".join(sorted(gene_reactions[gid])) or "-",
                ";".join(sorted(gene_subsystems[gid])) or "-",
                status,
            ])

    print(f"iML1515 genes: {len(genes)} | reference proteins written: {n_written} | excluded: {len(excluded)}")
    for gid, status in excluded:
        print(f"  [excluded] {gid}: {status}")


if __name__ == "__main__":
    main()
