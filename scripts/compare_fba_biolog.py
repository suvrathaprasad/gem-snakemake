#!/usr/bin/env python3
"""
Compare fba_media_test.py's growth/no-growth predictions against the
real Biolog phenotype microarray data (Supplement_TableS1_Biolog_AUC.xlsx,
Barth et al. 2020) for the same strain.

Only glyoxylate and rhamnose have a direct Biolog well to compare
against (F10/C06 on the Carbon1 PM1 plate) — glucose and no_carbon
are FBA-side sanity checks only (see fba_media_test.py) and are
reported but not scored against Biolog, since neither is one of the
paper's discriminatory substrates.

Biolog "positive"/"negative" call here is a simple AUC-above-negative-
-control threshold, NOT the paper's own statistical test (which used
replicate-level variance we don't have in this summary table). Treat
this as an approximate check, not a definitive one.
"""
import argparse
import csv

import openpyxl

BIOLOG_COLUMNS = {
    "glyoxylate": "F10 (Glyoxylic Acid)",
    "rhamnose": "C06 (L-Rhamnose)",
}
NEGATIVE_CONTROL_COLUMN = "A01 (Negative Control)"


def load_biolog_row(xlsx_path, labor_id, sheet="Carbon1"):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet]
    header = [c.value for c in ws[1]]
    id_col = header.index("LaborID")

    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[id_col] == labor_id:
            return dict(zip(header, row))

    raise ValueError(f"LaborID {labor_id!r} not found in sheet {sheet!r}")


def load_fba_results(tsv_path):
    with open(tsv_path) as f:
        return {row["carbon_source"]: row for row in csv.DictReader(f, delimiter="\t")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fba", required=True, help="fba_media_test.py output TSV")
    parser.add_argument("--biolog", required=True, help="Supplement_TableS1_Biolog_AUC.xlsx")
    parser.add_argument("--labor-id", required=True, help="Strain ID as it appears in the Biolog table, e.g. 13E0634")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    fba = load_fba_results(args.fba)
    biolog_row = load_biolog_row(args.biolog, args.labor_id)
    neg_ctrl_auc = biolog_row[NEGATIVE_CONTROL_COLUMN]

    rows = []
    for name, fba_row in fba.items():
        if name in BIOLOG_COLUMNS:
            biolog_auc = biolog_row[BIOLOG_COLUMNS[name]]
            biolog_call = "positive" if (biolog_auc - neg_ctrl_auc) >= 500 else "negative"
            fba_call = "growth" if fba_row["prediction"] == "growth" else "no growth"
            match = (fba_call == "growth") == (biolog_call == "positive")
            rows.append(
                (name, fba_row["growth_rate"], fba_call, f"{biolog_auc:.1f}",
                 f"{neg_ctrl_auc:.1f}", biolog_call, str(match))
            )
        else:
            # glucose / no_carbon: FBA-side sanity checks, no Biolog well to compare
            rows.append((name, fba_row["growth_rate"], fba_row["prediction"],
                         "NA", "NA", "NA (not in paper)", "NA"))

    header = ["condition", "fba_growth_rate", "fba_call", "biolog_auc",
              "biolog_negative_control_auc", "biolog_call", "match"]
    with open(args.output, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
            print("\t".join(row))


if __name__ == "__main__":
    main()

