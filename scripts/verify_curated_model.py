#!/usr/bin/env python3
"""
Verify the curated model's glyoxylate growth prediction.
Checks:
1. Mass balance on GLXt2r — the hand-written reaction.
   Confirms it does not silently create or destroy mass.
2. Flux trace under glyoxylate-only growth conditions.
   Shows which reactions carry nonzero flux and specifically checks
   the expected glyoxylate pathway reactions.
The verification report is automatically written next to the input
model.
Example:
    python scripts/verify_curated_model.py \
        PRJNA559322/models/SRR9972681/SRR9972681.curated.xml
This creates:
    PRJNA559322/models/SRR9972681/SRR9972681.verify_curated_model.txt
"""

import argparse
from pathlib import Path

import cobra


M9_BASE = {
    "EX_ca2_e": 1000.0,
    "EX_cl_e": 1000.0,
    "EX_co2_e": 1000.0,
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
    "EX_nh4_e": 1000.0,
    "EX_ni2_e": 1000.0,
    "EX_o2_e": 1000.0,
    "EX_pi_e": 1000.0,
    "EX_so4_e": 1000.0,
    "EX_zn2_e": 1000.0,
}


def run(model_path, output):
    """Run all verification checks and write the report to output."""

    model = cobra.io.read_sbml_model(model_path)

    # Check 1: mass balance on the hand-written reaction
    print("=" * 60, file=output)
    print("MASS BALANCE CHECK: GLXt2r", file=output)
    print("=" * 60, file=output)

    rxn = model.reactions.get_by_id("GLXt2r")
    imbalance = rxn.check_mass_balance()
    if imbalance:
        print(f"  [WARNING] Imbalanced: {imbalance}", file=output)
    else:
        print("  Balanced - no mass/charge created or destroyed.", file=output,)

    # Check 2: flux trace under glyoxylate-only condition
    print(file=output)
    print("=" * 60, file=output)
    print("FLUX TRACE: glyoxylate as sole carbon source", file=output)
    print("=" * 60, file=output)

    medium = dict(M9_BASE)
    medium["EX_glx_e"] = 10.0

    with model:
        model.medium = medium
        solution = model.optimize()

        print(f"  Growth rate: {solution.objective_value:.6f} 1/h", file=output,)
        print(f"  Status: {solution.status}", file=output)
        print(file=output)

        nonzero = solution.fluxes[solution.fluxes.abs() > 1e-6].sort_values(key=abs, ascending=False)

        print(f"  {len(nonzero)} reactions carry nonzero flux. Top 20 by magnitude:", file=output,)

        for rxn_id, flux in nonzero.head(20).items():
            rxn_obj = model.reactions.get_by_id(rxn_id)

            print(f"    {rxn_id:20s} {flux:12.4f}   {rxn_obj.name}", file=output,)

        # Check specifically expected glyoxylate pathway reactions
        print(file=output)
        print("  Specifically checking known candidates:", file=output)

        for candidate in [ "GLXt2r", "EX_glx_e", "GLXCL", "ACEA", "ACEB", "MALS", "ICL",]:
            if model.reactions.has_id(candidate):
                flux = solution.fluxes[candidate]
                print(f"    {candidate:20s} flux={flux:.6f}", file=output,)
            else:
                print(f"    {candidate:20s} not in model", file=output,)


def main():
    parser = argparse.ArgumentParser(description="Verify a curated COBRA model under glyoxylate-only growth.")
    parser.add_argument("model", help="Curated SBML model (.xml)",)
    args = parser.parse_args()

    model_path = Path(args.model)

    output_path = model_path.with_name(f"{model_path.stem.replace('.curated', '')}.verify_curated_model.txt")

    with output_path.open("w") as output:
        run(model_path, output)

    print(f"Verification written to: {output_path}")


if __name__ == "__main__":
    main()
