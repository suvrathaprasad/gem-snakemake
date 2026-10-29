"""Merge the per-strain scan, knockout and cross-check tables into one
cross-strain candidate table (one row per gene flagged in at least one strain).

Usage:
  candidate_table.py --out candidates.tsv \
      --strain NAME1 scan1.tsv knockout_single1.tsv knockout_crosscheck1.tsv \
      --strain NAME2 scan2.tsv knockout_single2.tsv knockout_crosscheck2.tsv ... \
      [--locus-check NAME1 glc_locus_check1.tsv --locus-check NAME2 ...]

Inputs per strain (all produced by the pipeline):
  {sample}.metabolic_gene_scan.tsv, {sample}.knockout_single.tsv,
  {sample}.knockout_crosscheck.tsv

Optional locus check: {sample}.glc_locus_check.tsv (the assembly-level search of
the rha/glc reference genes). For a gene whose symbol is one of its queries, the
table gets a `<strain>_locus_check` column with the top locus: "intact locus"
(>= 70% identity and >= 80% query coverage, one HSP), "split locus" (intact
identity/coverage over several HSPs, e.g. a frameshift), or "no intact locus"
(with the best hit and the feature it overlaps, e.g. a paralog). It does not
change the evidence tier.

Evidence tier, per strain (flagged genes only):
  weak   confidence = low in knockout_single (weak hit, or a hit shared with
         other reference genes)
  check  truncated call (could be a real pseudogene or a gene-calling split),
         or the strain-model knockout disagrees with the iML1515 rules
  strong everything else
There is no overall tier. Each gene carries its tier per strain, plus a
`strong_in` column ("strong in N of M strains", M = all strains in the run,
not only those in which the gene is flagged). Rows are sorted by the number of
strains with a strong call (most first), then by the number of strains flagged
(fewest first), then pattern and gene. These rules are a draft, not validated.

Growth effect: relative drop of each substrate's growth against the highest
growth in that strain's single-knockout table (the unaffected baseline),
reported in percent. A gene "affects a modelled substrate" when the drop is
>= --growth-threshold (default 0.01 = 1%) in any strain. Raw percentages are
always reported.
"""
import argparse
import csv
import sys

SUBSTRATES = ["glucose", "glyoxylate", "rhamnose"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def locus_note(row):
    pid, cov = fnum(row["pident"]), fnum(row["qcov"])
    feat = (row.get("overlapping_feature") or "").split(";")[0].split("|")
    feat = feat[-1].strip() if feat and feat[0] else "no annotated feature"
    desc = f"{pid:.1f}% identity, {cov:.1f}% coverage"
    if pid >= 70 and cov >= 80:
        if int(row["n_hsps"]) > 1:
            return f"split locus ({desc}, {row['n_hsps']} HSPs, possible frameshift)"
        return f"intact locus ({desc})"
    return f"no intact locus (best hit {desc}, overlaps: {feat})"


def tier_for(single, cross):
    if single.get("confidence") == "low":
        return "weak"
    if (single.get("detail") or "").startswith("truncated"):
        return "check"
    if cross and cross.get("agreement") == "different":
        return "check"
    return "strong"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strain", nargs=4, action="append", required=True,
                    metavar=("NAME", "SCAN", "SINGLE", "CROSSCHECK"))
    ap.add_argument("--locus-check", nargs=2, action="append", default=[],
                    metavar=("NAME", "LOCUS_CHECK"),
                    help="optional glc_locus_check.tsv for a strain named in --strain")
    ap.add_argument("--out", required=True)
    ap.add_argument("--growth-threshold", type=float, default=0.01)
    a = ap.parse_args()

    names = [s[0] for s in a.strain]
    if len(set(names)) != len(names):
        sys.exit("error: strain names must be unique")

    locus = {}
    for lname, lpath in a.locus_check:
        if lname not in names:
            sys.exit(f"error: --locus-check names an unknown strain: {lname}")
        top = {}
        for r in read_tsv(lpath):
            if r["rank"] == "1":
                top[r["query"].split("|")[0]] = r
        locus[lname] = top

    data = {}
    for name, scan_p, single_p, cross_p in a.strain:
        scan = {r["gene"]: r for r in read_tsv(scan_p)}
        single = {r["gene"]: r for r in read_tsv(single_p)}
        cross = {r["gene"]: r for r in read_tsv(cross_p)}
        flagged = {g for g, r in scan.items() if r["call"] in ("absent", "disrupted")}
        if flagged != set(single):
            sys.exit(f"error: {name}: flagged genes in the scan ({len(flagged)}) "
                     f"differ from knockout_single ({len(single)})")
        base = {}
        for s in SUBSTRATES:
            vals = [fnum(r.get(f"{s}_growth")) for r in single.values()]
            vals = [v for v in vals if v is not None]
            base[s] = max(vals) if vals else None
        data[name] = dict(scan=scan, single=single, cross=cross, flagged=flagged, base=base)

    all_flagged = sorted(set().union(*(d["flagged"] for d in data.values())))
    rows = []
    for gene in all_flagged:
        ref = next(d["scan"][gene] for d in data.values() if gene in d["scan"])
        in_strains = [n for n in names if gene in data[n]["flagged"]]
        row = {
            "gene": gene,
            "symbol": ref["symbol"],
            "subsystems": ref["subsystems"],
            "reactions_iML1515": ref["reactions"],
            "n_strains_flagged": len(in_strains),
        }
        if len(in_strains) == len(names):
            row["pattern"] = "all strains"
        elif len(in_strains) == 1:
            row["pattern"] = f"{in_strains[0]} only"
        else:
            row["pattern"] = " + ".join(in_strains)

        tiers, max_drop, max_where = [], 0.0, ""
        percol = {}
        for n in names:
            d = data[n]
            sc = d["scan"].get(gene, {})
            percol[f"{n}_call"] = sc.get("call", "")
            percol[f"{n}_identity"] = sc.get("pident", "")
            percol[f"{n}_coverage"] = sc.get("qcovhsp", "")
            if n in locus:
                lr = locus[n].get(ref["symbol"])
                percol[f"{n}_locus_check"] = locus_note(lr) if lr else ""
            if gene in d["flagged"]:
                sg = d["single"][gene]
                cr = d["cross"].get(gene)
                t = tier_for(sg, cr)
                tiers.append(t)
                percol[f"{n}_detail"] = sg.get("detail", "")
                percol[f"{n}_confidence"] = sg.get("confidence", "")
                percol[f"{n}_tier"] = t
                percol[f"{n}_effect"] = sg.get("effect", "")
                percol[f"{n}_strain_model_check"] = ((cr or {}).get("agreement")
                                                     or (cr or {}).get("status", ""))
                for s in SUBSTRATES:
                    g = fnum(sg.get(f"{s}_growth"))
                    b = d["base"][s]
                    if g is not None and b:
                        drop = (b - g) / b
                        if drop > max_drop:
                            max_drop, max_where = drop, f"{s} ({n})"
            else:
                for k in ("detail", "confidence", "tier", "effect", "strain_model_check"):
                    percol[f"{n}_{k}"] = ""
        row["n_strong"] = tiers.count("strong")
        row["strong_in"] = f"{row['n_strong']}/{len(names)}"
        row["max_growth_drop_pct"] = f"{max_drop * 100:.4f}"
        row["max_growth_drop_where"] = max_where
        row["affects_modelled_substrate"] = "yes" if max_drop >= a.growth_threshold else "no"
        row.update(percol)
        rows.append(row)

    rows.sort(key=lambda r: (-r["n_strong"], r["n_strains_flagged"],
                             r["pattern"], r["gene"]))

    head = ["gene", "symbol", "pattern", "n_strains_flagged", "strong_in",
            "affects_modelled_substrate", "max_growth_drop_pct", "max_growth_drop_where",
            "subsystems", "reactions_iML1515"]
    for n in names:
        head += [f"{n}_{k}" for k in ("call", "detail", "identity", "coverage", "confidence",
                                      "tier", "effect", "strain_model_check")]
        if any(n == l[0] for l in a.locus_check):
            head.append(f"{n}_locus_check")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=head, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} genes written to {a.out}", file=sys.stderr)
    for k in range(1, len(names) + 1):
        print(f"  flagged in {k} strain(s): {sum(r['n_strains_flagged'] == k for r in rows)}",
              file=sys.stderr)
    for k in range(len(names), -1, -1):
        print(f"  strong in {k} of {len(names)} strain(s): "
              f"{sum(r['n_strong'] == k for r in rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()

