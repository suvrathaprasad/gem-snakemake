#!/usr/bin/env python3
"""
Cross-check of knockout_flagged_genes.py: for every flagged gene whose
reciprocal-best-hit ortholog (ortholog_target in the scan table) is a gene
in the strain's own curated model, knock that locus tag out
with COBRApy on the strain model itself and compare the reactions
removed with the iML1515-GPR result in the knockout_single table.

Genes that cannot be checked this way are kept in the output with the
reason: no hit (no locus tag), best hit is not a reciprocal ortholog (the
strain protein belongs to another gene, e.g. glcA's best hit is lldP's
protein), or tag not in the strain model.

The two approaches use different gene rules (iML1515's vs the ones
CarveMe assigned to the strain model), so differences show where the
rule sets disagree. FBA is rerun after each knockout; there is no
gap-filling, so growth after a knockout is not forced back.
"""
import argparse
import csv
import io
from contextlib import redirect_stderr

import cobra

from fba_media_test import CARBON_SOURCES, run_fba_on_medium
from knockout_flagged_genes import SUBSTRATES, change, confidence


def growth(model):
    with redirect_stderr(io.StringIO()):
        return {s: run_fba_on_medium(model, CARBON_SOURCES[s])[0] for s in SUBSTRATES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", required=True, help="*.metabolic_gene_scan.tsv")
    ap.add_argument("--model", required=True, help="curated strain model (.xml)")
    ap.add_argument("--knockout-single", required=True, help="*.knockout_single.tsv")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    with open(args.scan) as f:
        scan = [r for r in csv.DictReader(f, delimiter="\t") if r["call"] != "intact"]
    with open(args.knockout_single) as f:
        ours = {r["gene"]: set(filter(None, r["rxns_removed_in_model"].split(";")))
                for r in csv.DictReader(f, delimiter="\t")}

    model = cobra.io.read_sbml_model(args.model)
    tags = {g.id for g in model.genes}
    orig = {r.id: r.bounds for r in model.reactions}
    base = growth(model)

    cols = (["gene", "symbol", "call", "confidence", "strain_locus_tag", "status",
             "rxns_removed_strain_model", "rxns_removed_iML1515", "agreement",
             "only_in_strain_model", "only_in_iML1515"]
            + [f"{s}_growth" for s in SUBSTRATES] + [f"{s}_change" for s in SUBSTRATES])
    n_checked = n_same = 0
    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in scan:
            tag = r["ortholog_target"] if r["ortholog_target"] not in ("", "-") else ""
            if r["best_hit_target"] in ("", "-"):
                status = "not checked: no hit, no locus tag"
            elif not tag:
                status = "not checked: best hit is not a reciprocal ortholog"
            elif tag not in tags:
                status = "not checked: locus tag not in strain model"
            else:
                status = "checked"
            head = [r["gene"], r["symbol"], r["call"], confidence(r), tag, status]
            if status != "checked":
                w.writerow(head + [""] * (len(cols) - len(head)))
                continue
            with model:
                model.genes.get_by_id(tag).knock_out()
                blocked = {x.id for x in model.reactions
                           if x.bounds == (0, 0) and orig[x.id] != (0, 0)}
                g = growth(model)
            mine = ours[r["gene"]]
            same = blocked == mine
            n_checked += 1
            n_same += same
            w.writerow(head + [";".join(sorted(blocked)), ";".join(sorted(mine)),
                               "same" if same else "different",
                               ";".join(sorted(blocked - mine)), ";".join(sorted(mine - blocked))]
                       + [f"{g[s]:.6f}" for s in SUBSTRATES]
                       + [change(base[s], g[s]) for s in SUBSTRATES])
    print(f"flagged genes: {len(scan)} | checked by locus tag: {n_checked} | "
          f"same reactions removed: {n_same} | different: {n_checked - n_same}")


if __name__ == "__main__":
    main()
