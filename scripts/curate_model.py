#!/usr/bin/env python3
"""
Add glyoxylate transport + exchange to a CarveMe-built model.

CarveMe's gram-negative universe has no glx_e metabolite or EX_glx_e
exchange reaction at all -- confirmed directly from the universe
build log, and confirmed NOT to be gene-driven (the same gap shows up
in strains where glcA is fully intact, e.g. SRR9972681, SRR9972673).
This is a limitation of the underlying BiGG-derived universe, not
something CarveMe's homology search could ever resolve by finding a
better gene hit for any strain.

Added uniformly to every strain's model based on genomic evidence
(glcA -- annotated glycolate/glyoxylate transporter -- is present, at
some identity, in every strain checked so far), NOT selectively to
make any one strain's FBA prediction match its Biolog result. See
build plan / chat log for the full reasoning.

Adds:
  - glx_e (extracellular glyoxylate) -- same formula/charge as glx_c
  - GLXt2r: glx_e + h_e <-> glx_c + h_c  (proton symport, matches
    GlcA's annotated symport mechanism; reversible)
  - EX_glx_e: glx_e <-> (boundary exchange, open bounds)
"""
import argparse

import cobra


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_model")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(args.input_model)

    if model.metabolites.has_id("glx_e"):
        print("glx_e already present -- nothing to add")
        cobra.io.write_sbml_model(model, args.output)
        return

    glx_c = model.metabolites.get_by_id("glx_c")
    h_e = model.metabolites.get_by_id("h_e")
    h_c = model.metabolites.get_by_id("h_c")

    # Derive the real external/cytoplasmic compartment IDs from an
    # existing metabolite rather than assuming "e"/"c" literally — this
    # model's SBML uses prefixed IDs (e.g. "C_e") instead.
    glx_e = cobra.Metabolite(
        "glx_e",
        formula=glx_c.formula,
        charge=glx_c.charge,
        name="Glyoxylate",
        compartment=h_e.compartment,
    )
    model.add_metabolites([glx_e])

    transport = cobra.Reaction(
        "GLXt2r", name="Glyoxylate transport, reversible proton symport"
    )
    transport.add_metabolites({glx_e: -1, h_e: -1, glx_c: 1, h_c: 1})
    transport.lower_bound = -1000
    transport.upper_bound = 1000
    model.add_reactions([transport])

    model.add_boundary(glx_e, type="exchange", reaction_id="EX_glx_e", lb=-1000, ub=1000)

    cobra.io.write_sbml_model(model, args.output)
    print(f"Added glx_e, GLXt2r, EX_glx_e -- wrote {args.output}")


if __name__ == "__main__":
    main()
