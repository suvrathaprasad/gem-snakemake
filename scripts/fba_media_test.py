#!/usr/bin/env python3
"""
FBA growth prediction for a gap-filled CarveMe model under M9 + a
single defined carbon source, matching the media used in Barth et al.
2020 (Toxins 12(6):414, PRJNA559322) to independently verify Biolog
phenotype-microarray calls for glyoxylic acid and L-rhamnose.
"""
import argparse
import sys

import cobra

# Core M9 minimal-medium ions/inorganic nutrients, aerobic. Open
# (effectively unconstrained) bounds for everything except the single
# carbon source, which is varied per condition below. IDs follow BiGG
# universal naming (matches CarveMe's default/gramneg universe).
M9_BASE = {
    "EX_ca2_e": 1000.0,
    "EX_cl_e": 1000.0,
    "EX_co2_e": 1000.0,  # produced/exchanged, left open
    "EX_cobalt2_e": 1000.0,
    "EX_cu2_e": 1000.0,
    "EX_fe2_e": 1000.0,
    "EX_fe3_e": 1000.0,
    "EX_h_e": 1000.0,
    "EX_h2o_e": 1000.0,
    "EX_k_e": 1000.0,
    "EX_mg2_e": 1000.0,
    "EX_mn2_e": 1000.0,
    "EX_mobd_e": 1000.0,
    "EX_na1_e": 1000.0,
    "EX_nh4_e": 1000.0,  # nitrogen source, matches paper's M9 recipe
    "EX_ni2_e": 1000.0,
    "EX_o2_e": 1000.0,  # aerobic, matches paper's assay
    "EX_pi_e": 1000.0,
    "EX_so4_e": 1000.0,  # sulfate, matches paper's M9 recipe
    "EX_zn2_e": 1000.0,
}

# Carbon source variants. Uptake bound of 10 mmol/gDW/h is the standard
# modeling convention (matches e.g. iML1515's own default glucose
# bound) — a qualitative growth/no-growth check, not a literal
# translation of the paper's 0.4% w/v carbon source concentration.
# "no_carbon" is the FBA equivalent of Biolog's negative-control well —
# no carbon exchange open at all, ions/minerals only. Should always
# come back as "no growth"; if it doesn't, something's wrong with the
# model or the ion list below, independent of glyoxylate/rhamnose.
CARBON_SOURCES = {
    "no_carbon": {},                   # negative control
    "glucose": {"EX_glc__D_e": 10.0},  # positive control / sanity check
    "glyoxylate": {"EX_glx_e": 10.0},  # paper: STECper glc pathway often impaired
    "rhamnose": {"EX_rmn_e": 10.0},  # paper: STECper rha pathway often impaired
}


def run_fba_on_medium(model, carbon_exchanges):
    medium = dict(M9_BASE)
    medium.update(carbon_exchanges)

    # Only keep exchanges that actually exist in this model — CarveMe's
    # universe may not include every BiGG ID, and a missing key would
    # otherwise raise inside cobrapy's medium setter.
    missing = [rxn_id for rxn_id in medium if not model.reactions.has_id(rxn_id)]
    for rxn_id in missing:
        print(f"  [warn] {rxn_id} not in model — skipping", file=sys.stderr)
        medium = {k: v for k, v in medium.items() if k not in missing}

    with model:  # context manager: bounds revert after this block
        model.medium = medium
        solution = model.optimize()
        growth = solution.objective_value if solution.status == "optimal" else 0.0

    return growth, missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", help="Gap-filled SBML model (.xml)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(args.model)

    rows = []
    for name, exchanges in CARBON_SOURCES.items():
        growth, missing = run_fba_on_medium(model, exchanges)
        predicted = "growth" if growth > 1e-6 else "no growth"
        rows.append((name, f"{growth:.6f}", predicted, ";".join(missing)))
        print(f"{name}: growth_rate={growth:.6f} -> {predicted}")

    with open(args.output, "w") as f:
        f.write("carbon_source\tgrowth_rate\tprediction\tmissing_exchanges\n")
        for row in rows:
            f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
