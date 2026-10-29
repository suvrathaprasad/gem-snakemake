#!/usr/bin/env python3
"""
Knock out the metabolic genes flagged by scan_metabolic_genes.py (call =
disrupted or absent) in a strain's curated model, and rerun FBA on
M9 + glucose / glyoxylate / rhamnose.

Reaction-level knockout via the iML1515 GPR rules: flagged genes are set
to lost, every iML1515 reaction's rule is re-evaluated (`and` = complex,
all subunits needed; `or` = isozymes, any one is enough), and reactions
whose rule comes out False are removed from the strain model by BiGG ID
(where present). Absent genes have no protein in the strain model, so
the strain model's own gene rules cannot be used for this.

Outputs
  --out-single    one row per flagged gene, knocked out on its own
  --out-combined  one row per carbon source, all flagged genes together
  --out-model     (optional) the combined-knockout model as SBML
"""
import argparse
import ast
import csv
import gzip
import io
import json
import sys
from contextlib import redirect_stderr

import cobra

from fba_media_test import CARBON_SOURCES, run_fba_on_medium

SUBSTRATES = ["glucose", "glyoxylate", "rhamnose"]  # glucose = positive control
TOL = 1e-6


def load_rules(path):
    """iML1515 JSON (gz ok) -> {reaction_id: rule string}, rules non-empty only."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        model = json.load(f)
    return {r["id"]: r["gene_reaction_rule"].strip()
            for r in model["reactions"] if r["gene_reaction_rule"].strip()}


def rule_holds(rule, lost):
    """ 
    If the genes in `lost` are non-functional, is this reaction's 
    GPR still satisfied? 
    True if the rule is still satisfied when the genes in `lost` are gone.
    """
    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.BoolOp):
            vals = [ev(v) for v in node.values]
            return all(vals) if isinstance(node.op, ast.And) else any(vals)
        if isinstance(node, ast.Name):
            return node.id not in lost
        raise ValueError(f"unexpected element in GPR rule: {rule!r}")
    return ev(ast.parse(rule, mode="eval"))


def reactions_lost(rules, gene_to_rxns, lost):
    """Reactions whose GPR fails once `lost` genes are removed. and category."""
    candidates = set()
    for g in lost:
        candidates |= gene_to_rxns.get(g, set())
    return {r for r in candidates if not rule_holds(rules[r], lost)}


def confidence(row):
    """low = weak hit, or best-hit protein shared with other reference genes."""
    shared = row["n_refs_sharing_target"]
    if row["detail"].startswith("weak"):
        return "low"
    if shared not in ("", "-") and int(shared) > 1:
        return "low"
    return "high"


def growth_all(model, remove_ids, quiet=True):
    """Growth per substrate with `remove_ids` removed (reverted afterwards)."""
    out = {}
    with model:
        model.remove_reactions([model.reactions.get_by_id(r) for r in remove_ids],
                               remove_orphans=False)
        for s in SUBSTRATES:
            if quiet:
                with redirect_stderr(io.StringIO()):
                    out[s] = run_fba_on_medium(model, CARBON_SOURCES[s])[0]
            else:
                out[s] = run_fba_on_medium(model, CARBON_SOURCES[s])[0]
    return out


def change(base, new):
    if abs(new - base) <= TOL * max(1.0, abs(base)):
        return "none"
    return "no growth" if new <= TOL else ("reduced" if new < base else "increased")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iml1515-json", required=True)
    ap.add_argument("--scan", required=True, help="*.metabolic_gene_scan.tsv")
    ap.add_argument("--model", required=True, help="curated strain model (.xml)")
    ap.add_argument("--out-single", required=True)
    ap.add_argument("--out-combined", required=True)
    ap.add_argument("--out-model")
    args = ap.parse_args()

    """Build reverse lookup - for every gene in iML1515, which reactions mention it?"""
    rules = load_rules(args.iml1515_json)
    gene_to_rxns = {}
    for rid, rule in rules.items():
        for node in ast.walk(ast.parse(rule, mode="eval")):
            if isinstance(node, ast.Name):
                gene_to_rxns.setdefault(node.id, set()).add(rid)

    with open(args.scan) as f:
        scan = list(csv.DictReader(f, delimiter="\t"))
    flagged = [r for r in scan if r["call"] in ("disrupted", "absent")]
    flagged_ids = {r["gene"] for r in flagged}

    model = cobra.io.read_sbml_model(args.model)
    in_model = {r.id for r in model.reactions}

    base = growth_all(model, [], quiet=False)
    print("baseline growth: " + ", ".join(f"{s}={base[s]:.4f}" for s in SUBSTRATES))
    if base["glucose"] <= TOL:
        sys.exit("baseline model does not grow on glucose - check the model before interpreting knockouts")

    # --- single knockouts -------------------------------------------------
    single_lost = {}
    cols = (["gene", "symbol", "call", "detail", "confidence", "subsystems",
             "rxns_lost_iML1515", "rxns_removed_in_model", "effect"]
            + [f"{s}_growth" for s in SUBSTRATES] + [f"{s}_change" for s in SUBSTRATES])
    with open(args.out_single, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in flagged:
            lost = reactions_lost(rules, gene_to_rxns, {r["gene"]}) # reactions predicted to fail based on GPR
            single_lost[r["gene"]] = lost
            removed = sorted(lost & in_model) # those predicted lost reactions that exist in the strain
            g = growth_all(model, removed) if removed else dict(base)
            if not lost:
                effect = "no reaction lost (isozyme intact or no reactions)"
            elif not removed:
                effect = "reactions lost but not in strain model"
            elif all(change(base[s], g[s]) == "none" for s in SUBSTRATES):
                effect = "reactions removed, growth unchanged"
            else:
                effect = "growth changed"
            w.writerow([r["gene"], r["symbol"], r["call"], r["detail"], confidence(r),
                        r["subsystems"], ";".join(sorted(lost)), ";".join(removed), effect]
                       + [f"{g[s]:.6f}" for s in SUBSTRATES]
                       + [change(base[s], g[s]) for s in SUBSTRATES])

    # --- combined knockout ------------------------------------------------
    lost_all = reactions_lost(rules, gene_to_rxns, flagged_ids)
    removed_all = sorted(lost_all & in_model) # both from and category
    only_combo = sorted((lost_all & in_model) - set().union(*single_lost.values())) # reactions lost when both in and category are lost
    g_all = growth_all(model, removed_all)
    with open(args.out_combined, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["carbon_source", "baseline_growth", "combined_growth", "change",
                    "n_flagged_genes", "n_rxns_lost_iML1515", "n_rxns_removed_in_model",
                    "rxns_removed_in_model", "rxns_lost_only_in_combination"])
        for s in SUBSTRATES:
            w.writerow([s, f"{base[s]:.6f}", f"{g_all[s]:.6f}", change(base[s], g_all[s]),
                        len(flagged_ids), len(lost_all), len(removed_all),
                        ";".join(removed_all), ";".join(only_combo)])

    if args.out_model:
        ko = model.copy()
        ko.remove_reactions([ko.reactions.get_by_id(r) for r in removed_all], remove_orphans=False)
        cobra.io.write_sbml_model(ko, args.out_model)

    print(f"flagged genes: {len(flagged_ids)} | iML1515 reactions lost: {len(lost_all)} | "
          f"removed from model: {len(removed_all)} | only lost in combination: {len(only_combo)}")
    print("combined: " + ", ".join(f"{s}={g_all[s]:.4f} ({change(base[s], g_all[s])})" for s in SUBSTRATES))


if __name__ == "__main__":
    main()
