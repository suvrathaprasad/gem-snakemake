#!/usr/bin/env python3
"""
Per-strain summary report: pulls together BUSCO completeness, Bakta
annotation stats, CarveMe draft/gapfilled model size, the FBA media
test, the FBA-vs-Biolog comparison, and the operon gene check into
one PDF with plots.

Everything shown is read directly from each stage's own output file —
nothing here is hardcoded per-strain, so the same script produces a
meaningful report for any sample with these files present (e.g. the
next, sporadic-strain run). The one piece of encoded domain knowledge
is GENE_RELEVANCE below (which genes matter for which substrate) —
that's pathway biology, not a per-strain conclusion.
"""
import argparse
import csv
import glob
import os
import re
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openpyxl
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

CALL_COLORS = {
    "intact": "#4caf50",
    "disrupted (partial/low-identity hit)": "#ff9800",
    "absent": "#e53935",
}

# Which reference genes are relevant to which substrate, for the
# integrated-findings page. Pathway biology, not a per-strain fact.
GENE_RELEVANCE = {
    "glyoxylate": ["glcA", "glcB", "glcC", "glcD", "glcE", "glcF", "glcG"],
    "rhamnose": ["rhaA", "rhaB", "rhaD", "rhaT", "rhaR", "rhaS", "rhaM"],
}


def find_busco_summary(busco_dir):
    matches = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
    if not matches:
        return None
    with open(matches[0]) as f:
        text = f.read()
    # Real BUSCO format: only S and D sit inside the brackets;
    # F and M come after, outside them.
    m = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    if not m:
        return None
    complete, single, dup, frag, missing, n = m.groups()
    return {
        "complete": float(complete),
        "single": float(single),
        "duplicated": float(dup),
        "fragmented": float(frag),
        "missing": float(missing),
        "n": int(n),
    }


def count_model_stats(xml_path):
    if not os.path.exists(xml_path):
        return None
    with open(xml_path) as f:
        text = f.read()
    return {
        "reactions": text.count("<reaction "),
        "species": text.count("<species "),
    }


def load_bakta_stats(bakta_txt_path):
    stats = {}
    if not os.path.exists(bakta_txt_path):
        return stats
    with open(bakta_txt_path) as f:
        for line in f:
            if ":" in line:
                key, _, value = line.partition(":")
                stats[key.strip()] = value.strip()
    return stats


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _find_column(header, keyword):
    """Case-insensitive grep"""
    for i, name in enumerate(header):
        if name and keyword.lower() in str(name).lower():
            return i
    return None


def load_biolog_strain_meta(xlsx_path, labor_id, sheet="Carbon1"):
    if not os.path.exists(xlsx_path):
        return {}
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet]
    header = [c.value for c in ws[1]]
    id_col = header.index("LaborID")
    geno_col = _find_column(header, "serotype")
    colonization_col = _find_column(header, "coloniz")
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[id_col] == labor_id:
            return {
                "geno_serotype": row[geno_col] if geno_col is not None else "NA",
                "colonization": row[colonization_col] if colonization_col is not None else "NA",
            }
    return{}

def _label_vertical_bars(ax, bars, fmt="{:.2f}", fontsize=8):
    """Label vertical bars using their actual height/value."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for bar in bars:
        height = bar.get_height()

        if height == 0:
            continue

        label = fmt.format(height)

        # Measure rendered bar height.
        x_center = bar.get_x() + bar.get_width() / 2
        bottom_display = ax.transData.transform((x_center, bar.get_y()))
        top_display = ax.transData.transform((x_center, bar.get_y() + height))

        bar_pixels = abs(top_display[1] - bottom_display[1])

        # Measure label height.
        text = ax.text(
            0,
            0,
            label,
            fontsize=fontsize,
            fontweight="bold",
            color="white",
        )
        bbox = text.get_window_extent(renderer=renderer)
        text_height = bbox.height
        text.remove()

        if bar_pixels < text_height + 4:
            label_size = max(5, fontsize - 2)
        elif bar_pixels < text_height + 10:
            label_size = max(6, fontsize - 1)
        else:
            label_size = fontsize

        if bar_pixels < 12:
            offset = 1
        elif bar_pixels < 20:
            offset = 3
        else:
            offset = 6

        ax.annotate(
            label,
            xy=(x_center, bar.get_y() + height),
            xytext=(0, -offset),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=label_size,
            color="white",
            fontweight="bold",
            clip_on=True,
        )


def _label_horizontal_bars(ax, bars, fmt="{:.0f}", fontsize=8):
    """Label horizontal bars using their actual width/value."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for bar in bars:
        width = bar.get_width()

        if width == 0:
            continue

        label = fmt.format(width)

        y_center = bar.get_y() + bar.get_height() / 2

        left_display = ax.transData.transform((bar.get_x(), y_center))
        right_display = ax.transData.transform((bar.get_x() + width, y_center))

        bar_pixels = abs(right_display[0] - left_display[0])

        text = ax.text(
            0,
            0,
            label,
            fontsize=fontsize,
            fontweight="bold",
            color="white",
        )
        bbox = text.get_window_extent(renderer=renderer)
        text_width = bbox.width
        text.remove()

        if bar_pixels < text_width + 4:
            label_size = max(5, fontsize - 2)
        elif bar_pixels < text_width + 10:
            label_size = max(6, fontsize - 1)
        else:
            label_size = fontsize

        if bar_pixels < 18:
            offset = 1
        elif bar_pixels < 30:
            offset = 3
        else:
            offset = 6

        ax.annotate(
            label,
            xy=(bar.get_x() + width, y_center),
            xytext=(-offset, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=label_size,
            color="white",
            fontweight="bold",
            clip_on=True,
        )


def plot_busco(busco, out_path):
    labels = ["Single", "Duplicated", "Fragmented", "Missing"]
    values = [busco["single"], busco["duplicated"], busco["fragmented"], busco["missing"]]
    bar_colors = ["#4caf50", "#2196f3", "#ff9800", "#e53935"]

    fig, ax = plt.subplots(figsize=(6, 2.8))
    left = 0
    for label, value, color in zip(labels, values, bar_colors):
        ax.barh(["BUSCO"], [value], left=left, color=color, label=f"{label} ({value:.1f}%)")
        left += value
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of BUSCO genes")
    ax.set_title(f"BUSCO completeness (n={busco['n']} genes)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.45), ncol=2, fontsize=8, frameon=True,)
    fig.subplots_adjust(bottom=0.42, top=0.82)
    fig.savefig(out_path, dpi=150, bbox_inches="tight",)
    plt.close(fig)


def plot_fba_growth(fba_rows, out_path):
    conditions = [r["carbon_source"] for r in fba_rows]
    rates = [float(r["growth_rate"]) for r in fba_rows]
    colors_list = ["#4caf50" if r["prediction"] == "growth" else "#e53935" for r in fba_rows]

    fig, ax = plt.subplots(figsize=(6, 3.5))
    bars = ax.bar(conditions, rates, color=colors_list)
    _label_vertical_bars(ax, bars)
    ax.set_ylabel("FBA growth rate (1/h)")
    ax.set_ylim(0, max(rates + [1.0]) * 1.25)
    ax.set_title("FBA growth prediction by condition")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_operon_genes(operon_rows, out_path):
    genes = [r["gene"].split("|")[0] for r in operon_rows]
    coverage = [float(r["qcovhsp"]) if r["qcovhsp"] != "-" else 0.0 for r in operon_rows]
    bar_colors = [CALL_COLORS.get(r["call"], "#9e9e9e") for r in operon_rows]

    fig, ax = plt.subplots(figsize=(6, 4.2))
    bars = ax.barh(genes, coverage, color=bar_colors)
    _label_horizontal_bars(ax, bars, fmt="{:.0f}%")
    ax.set_xlabel("Query coverage (%)")
    ax.set_xlim(0, 115)
    ax.set_title("Reference gene coverage vs. this strain's proteome")
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CALL_COLORS.values()]
    ax.legend(legend_handles, CALL_COLORS.keys(), loc="upper center",
              bbox_to_anchor=(0.5, -0.15), fontsize=7, ncol=1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_operon_scatter(operon_rows, out_path):
    """Identity vs. coverage — separates the two different failure
    modes seen: divergence (high coverage, low identity) vs.
    truncation (high identity, low coverage)."""
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for r in operon_rows:
        if r["pident"] == "-":
            continue
        gene = r["gene"].split("|")[0]
        x, y = float(r["pident"]), float(r["qcovhsp"])
        ax.scatter(x, y, color=CALL_COLORS.get(r["call"], "#9e9e9e"), s=60, zorder=3)
        ax.annotate(gene, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Identity (%)")
    ax.set_ylabel("Coverage (%)")
    ax.set_xlim(40, 105)
    ax.set_ylim(40, 105)
    ax.set_title("Divergence vs. truncation across the 14 reference genes")
    ax.grid(alpha=0.3)
    legend_handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8)
                       for c in CALL_COLORS.values()]
    ax.legend(legend_handles, CALL_COLORS.keys(), loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_fba_vs_biolog_grouped(compare_rows, out_path):
    """Compare FBA growth rate with Biolog AUC using separate y-axes.

    FBA growth rate is shown on the left axis and Biolog AUC on the
    right axis because the two measurements have different units and
    numerical scales.
    """
    rows = [r for r in compare_rows if r["biolog_auc"] != "NA"]
    if not rows:
        return None

    conditions = [r["condition"] for r in rows]
    fba_rates = [float(r["fba_growth_rate"]) for r in rows]
    biolog_auc = [float(r["biolog_auc"]) for r in rows]
    neg_ctrl = float(rows[0]["biolog_negative_control_auc"])
    x = range(len(conditions))
    width = 0.35
    fig, ax1 = plt.subplots(figsize=(5.8, 4.0))
    ax2 = ax1.twinx()
    # FBA growth rate — left axis
    bars_fba = ax1.bar(
        [i - width / 2 for i in x],
        fba_rates, width, color="#2196f3", label="FBA growth rate", )
    # Biolog AUC — right axis
    bars_biolog = ax2.bar(
        [i + width / 2 for i in x], 
        biolog_auc, width, color="#9c27b0", alpha=0.85, label="Biolog AUC", )
    # Negative-control reference on Biolog scale
    ax2.axhline(neg_ctrl, color="black", linestyle="--", linewidth=1, 
                label=f"Biolog negative control ({neg_ctrl:.0f})", )
    # Value labels
    _label_vertical_bars(ax1, bars_fba, fmt="{:.2f}", )
    _label_vertical_bars(ax2, bars_biolog, fmt="{:.0f}", )
    # X axis
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(conditions)
    # Separate scales
    ax1.set_ylabel("FBA growth rate (1/h)", color="#1976d2", )
    ax2.set_ylabel("Biolog AUC", color="#7b1fa2", )
    ax1.set_title("FBA prediction vs. Biolog result")
    # Both axes start at zero, but each has its own scale.
    max_fba = max(fba_rates + [0.01])
    ax1.set_ylim(0, max_fba * 1.35, )
    max_biolog = max(biolog_auc + [neg_ctrl, 0.01])
    ax2.set_ylim(0, max_biolog * 1.20, )
    # Combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, fontsize=7, loc="upper left", )
    # Subtle horizontal grid using the FBA axis.
    ax1.grid(axis="y", alpha=0.25, linestyle="-",)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150,)
    plt.close(fig)

    return out_path


def plot_curated_model_vs_biolog(compare_rows, curated_rows, out_path):
    """Show the effect of model curation alongside the experimental Biolog result.

    Original and curated FBA growth rates use the left axis.
    Biolog AUC uses the right axis because it is on a different scale.
    """
    curated_by_condition = {
        r["condition"]: r for r in curated_rows
    }

    rows = []
    for r in compare_rows:
        curated = curated_by_condition.get(r["condition"])
        if curated is None or r["biolog_auc"] == "NA":
            continue
        rows.append((r, curated))

    if not rows:
        return None

    conditions = [r["condition"] for r, _ in rows]
    original_rates = [float(r["fba_growth_rate"]) for r, _ in rows]
    curated_rates = [float(c["fba_growth_rate"]) for _, c in rows]
    biolog_auc = [float(r["biolog_auc"]) for r, _ in rows]
    neg_ctrl = float(rows[0][0]["biolog_negative_control_auc"])
    x = range(len(conditions))
    width = 0.25
    fig, ax1 = plt.subplots(figsize=(6.5, 4.2))
    ax2 = ax1.twinx()
    # Original model
    bars_original = ax1.bar(
        [i - width for i in x],
        original_rates, width, color="#b0bec5", label="Original FBA", )
    # Curated model
    bars_curated = ax1.bar(list(x), curated_rates, width, color="#4caf50", label="Curated FBA", )
    # Experimental Biolog result
    bars_biolog = ax2.bar(
        [i + width for i in x],
        biolog_auc, width, color="#9c27b0", alpha=0.85, label="Biolog AUC", )
    # Negative control reference
    ax2.axhline(neg_ctrl, color="black", linestyle="--", linewidth=1, 
                label=f"Biolog negative control ({neg_ctrl:.0f})", )
    # Value labels
    _label_vertical_bars(ax1, bars_original, fmt="{:.2f}")
    _label_vertical_bars(ax1, bars_curated, fmt="{:.2f}")
    _label_vertical_bars(ax2, bars_biolog, fmt="{:.0f}")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(conditions)
    ax1.set_ylabel("FBA growth rate (1/h)", color="#455a64", )
    ax2.set_ylabel("Biolog AUC", color="#7b1fa2", )
    ax1.set_title("Effect of model curation on FBA prediction")
    # Keep the FBA axis anchored at zero.
    max_fba = max(original_rates + curated_rates + [0.01])
    ax1.set_ylim(0, max_fba * 1.35)
    max_biolog = max(biolog_auc + [neg_ctrl, 0.01])
    ax2.set_ylim(0, max_biolog * 1.20)
    # Make the two legends one combined legend.
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, fontsize=7, loc="upper left", )
    ax1.grid(axis="y", alpha=0.25, linestyle="-", )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    return out_path


def plot_model_comparison(draft_stats, gapfilled_stats, out_path):
    categories = ["reactions", "species"]
    draft_vals = [draft_stats[c] for c in categories]
    gapfilled_vals = [gapfilled_stats[c] for c in categories]

    x = range(len(categories))
    width = 0.35
    fig, ax = plt.subplots(figsize=(5, 3.2))
    bars1 = ax.bar([i - width / 2 for i in x], draft_vals, width, color="#607d8b", label="draft")
    bars2 = ax.bar([i + width / 2 for i in x], gapfilled_vals, width, color="#4caf50", label="gap-filled")
    _label_vertical_bars(ax, bars1, fmt="{:.0f}")
    _label_vertical_bars(ax, bars2, fmt="{:.0f}")
    max_value = max(draft_vals + gapfilled_vals)
    ax.set_ylim(0, max_value * 1.20)
    ax.set_xticks(list(x))
    ax.set_xticklabels(["reactions", "metabolites"])
    ax.set_title("Draft vs. gap-filled model size")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_bakta_features(bakta_stats, out_path):
    keys = ["CDSs", "tRNAs", "rRNAs", "ncRNAs", "pseudogenes", "hypotheticals", "sORFs", "tmRNAs"]
    labels, values = [], []
    for k in keys:
        if k in bakta_stats:
            try:
                values.append(int(bakta_stats[k]))
                labels.append(k)
            except ValueError:
                continue
    if not values:
        return None

    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    bars = ax.bar(labels, values, color="#3f51b5")
    _label_vertical_bars(ax, bars, fmt="{:.0f}")
    ax.set_yscale("symlog")
    ax.set_ylim(0, max(values) * 1.35)
    ax.set_ylabel("Count")
    ax.set_title("Bakta annotation feature counts")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_biolog_compare(compare_rows, out_path):
    rows = [r for r in compare_rows if r["biolog_auc"] != "NA"]
    if not rows:
        return None
    conditions = [r["condition"] for r in rows]
    auc = [float(r["biolog_auc"]) for r in rows]
    neg_ctrl = float(rows[0]["biolog_negative_control_auc"])
    bar_colors = ["#4caf50" if r["biolog_call"] == "positive" else "#e53935" for r in rows]

    fig, ax = plt.subplots(figsize=(5, 3.5))
    bars = ax.bar(conditions, auc, color=bar_colors)
    _label_vertical_bars(ax, bars, fmt="{:.0f}")
    ax.axhline(neg_ctrl, color="black", linestyle="--", linewidth=1,
               label=f"negative control ({neg_ctrl:.0f})")
    ax.set_ylabel("Biolog AUC")
    ax.set_title("Biolog AUC vs. negative control")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def make_table(data, col_widths=None):
    t = Table(data, colWidths=col_widths)
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c3e50")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f5f5")]),
            ]
        )
    )
    return t


def build_integrated_findings(fba_rows, compare_rows, operon_rows):
    """For each mismatched substrate, list the status of its relevant
    genes — the connective narrative the rest of the report doesn't
    otherwise spell out.

    IMPORTANT distinction: a substrate whose FBA condition has a
    non-empty `missing_exchanges` means the model has NO representation
    of that pathway at all — the exchange reaction isn't in CarveMe's
    universe, full stop, regardless of any gene's status. Attributing
    that kind of gap to a specific gene's disruption is misleading
    (confirmed by comparing multiple strains: the same exchange was
    missing even where the relevant gene was fully intact). Only a
    substrate the model could actually evaluate (empty
    missing_exchanges) gets a gene-based candidate explanation.
    """
    fba_by_condition = {r["carbon_source"]: r for r in fba_rows}
    operon_by_gene = {r["gene"].split("|")[0]: r for r in operon_rows}
    findings = []
    for row in compare_rows:
        if row.get("match") != "False":
            continue
        condition = row["condition"]
        fba_row = fba_by_condition.get(condition, {})
        missing = fba_row.get("missing_exchanges", "").strip()

        if missing:
            note = (f"FBA predicted '{row['fba_call']}' but Biolog shows "
                    f"'{row['biolog_call']}'. This is NOT a gene-level "
                    f"finding — the model has no representation of this "
                    f"pathway at all (missing exchange reaction: {missing}), "
                    f"regardless of any gene's status. The operon gene check "
                    f"below cannot explain this kind of gap; treat it as a "
                    f"limitation of the model's reaction universe instead.")
        else:
            genes = GENE_RELEVANCE.get(condition, [])
            flagged = [(g, operon_by_gene[g]["call"]) for g in genes
                       if g in operon_by_gene and operon_by_gene[g]["call"] != "intact"]
            if flagged:
                gene_desc = "; ".join(f"{g}: {call}" for g, call in flagged)
                note = (f"FBA predicted '{row['fba_call']}' but Biolog shows "
                        f"'{row['biolog_call']}'. Candidate explanation — "
                        f"{gene_desc}.")
            else:
                note = (f"FBA predicted '{row['fba_call']}' but Biolog shows "
                        f"'{row['biolog_call']}'. No disruption found among the "
                        f"checked genes for this pathway — mismatch unexplained "
                        f"by this analysis.")
        findings.append((condition, note))
    return findings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="SRR sample ID")
    parser.add_argument("--labor-id", required=True, help="Strain ID as in the Biolog table")
    parser.add_argument("--busco-dir", required=True)
    parser.add_argument("--bakta-txt", required=True)
    parser.add_argument("--draft-model", required=True)
    parser.add_argument("--gapfilled-model", required=True)
    parser.add_argument("--fba-results", required=True)
    parser.add_argument("--biolog-compare", required=True)
    parser.add_argument("--operon-check", required=True)
    parser.add_argument("--biolog-xlsx", required=True)
    parser.add_argument("--curated-compare", help="fba_biolog_compare_curated.tsv, if the glyoxylate-curation step has run")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    busco = find_busco_summary(args.busco_dir)
    bakta_stats = load_bakta_stats(args.bakta_txt)
    draft_stats = count_model_stats(args.draft_model)
    gapfilled_stats = count_model_stats(args.gapfilled_model)
    fba_rows = load_tsv(args.fba_results)
    compare_rows = load_tsv(args.biolog_compare)
    operon_rows = load_tsv(args.operon_check)
    curated_rows = load_tsv(args.curated_compare) if args.curated_compare else []
    strain_meta = load_biolog_strain_meta(args.biolog_xlsx, args.labor_id)

    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph(f"Genome scale metabolic modelling of strain: {args.labor_id} ({args.sample})", styles["Title"]))
    story.append(Spacer(1, 12))

    meta_lines = [
        f"Geno-serotype: {strain_meta.get('geno_serotype', 'NA')}",
        f"Colonization type: {strain_meta.get('colonization', 'NA')}",
    ]
    for line in meta_lines:
        story.append(Paragraph(line, styles["Normal"]))
    story.append(Spacer(1, 12))

    story.append(Paragraph("Pipeline stage summary", styles["Heading2"]))
    stage_rows = [["Stage", "Result"]]
    if busco:
        stage_rows.append(["BUSCO completeness",
                            f"{busco['complete']:.1f}% (n={busco['n']} genes, "
                            f"{busco['duplicated']:.1f}% duplicated, {busco['missing']:.1f}% missing)"])
    if bakta_stats:
        cds = bakta_stats.get("CDS") or bakta_stats.get("CDSs") or "NA"
        stage_rows.append(["Bakta CDS count", cds])
    if draft_stats:
        stage_rows.append(["CarveMe draft model", f"{draft_stats['reactions']} reactions, {draft_stats['species']} metabolites"])
    if gapfilled_stats:
        stage_rows.append(["CarveMe gap-filled model", f"{gapfilled_stats['reactions']} reactions, {gapfilled_stats['species']} metabolites"])
    story.append(make_table(stage_rows, col_widths=[180, 300]))
    story.append(Spacer(1, 12))

    if bakta_stats:
        png = plot_bakta_features(bakta_stats, tempfile.NamedTemporaryFile(suffix=".png", delete=False).name)
        if png:
            story.append(Image(png, width=380, height=240))
            story.append(Spacer(1, 12))

    if busco:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_busco(busco, tmp.name)
            story.append(Image(tmp.name, width=380, height=235))

    if draft_stats and gapfilled_stats:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_model_comparison(draft_stats, gapfilled_stats, tmp.name)
            story.append(Image(tmp.name, width=320, height=210))
    story.append(PageBreak())

    # Integrated findings — the connective narrative
    findings = build_integrated_findings(fba_rows, compare_rows, operon_rows)
    if findings:
        story.append(Paragraph("Integrated findings", styles["Heading2"]))
        story.append(Paragraph(
            "Substrates where the FBA prediction disagreed with the real Biolog "
            "result, cross-referenced against the operon gene check to see "
            "whether a disrupted gene explains the mismatch.",
            styles["Normal"],
        ))
        story.append(Spacer(1, 8))
        for condition, note in findings:
            story.append(Paragraph(f"<b>{condition}:</b> {note}", styles["Normal"]))
            story.append(Spacer(1, 6))
        story.append(PageBreak())

    if fba_rows:
        story.append(Paragraph("FBA growth predictions", styles["Heading2"]))
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_fba_growth(fba_rows, tmp.name)
            story.append(Image(tmp.name, width=380, height=220))
        story.append(Spacer(1, 8))
        fba_table = [["Condition", "Growth rate", "Prediction", "Missing exchanges"]]
        for r in fba_rows:
            fba_table.append([r["carbon_source"], r["growth_rate"], r["prediction"], r.get("missing_exchanges", "")])
        story.append(make_table(fba_table, col_widths=[100, 90, 90, 160]))
        story.append(PageBreak())

    if compare_rows:
        story.append(Paragraph("FBA vs. Biolog comparison", styles["Heading2"]))

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            png = plot_biolog_compare(compare_rows, tmp.name)
            if png:
                story.append(Image(png, width=320, height=220))

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            png = plot_fba_vs_biolog_grouped(compare_rows, tmp.name)
            if png:
                story.append(Image(png, width=320, height=220))

        story.append(Spacer(1, 8))

        compare_table = [
            [
                "Condition",
                "FBA call",
                "Biolog AUC",
                "Negative ctrl AUC",
                "Biolog call",
                "Match",
            ]
        ]

        for r in compare_rows:
            compare_table.append([
                r["condition"],
                r["fba_call"],
                r["biolog_auc"],
                r["biolog_negative_control_auc"],
                r["biolog_call"],
                r["match"],
            ])

        story.append(make_table(compare_table, col_widths=[70, 70, 70, 90, 90, 60], ))
        story.append(PageBreak())


        if curated_rows:
            story.append(Paragraph(
                "Curated model re-test (glyoxylate transport added)",
                styles["Heading2"], ))

            story.append(
                Paragraph(
                    "CarveMe's universe has no glyoxylate exchange reaction, "
                    "which was confirmed directly from the universe build log. "
                    "This was therefore treated as a limitation of the model's "
                    "reaction universe rather than as evidence of a disrupted gene. "
                    "A transport + exchange reaction was subsequently added based "
                    "on genomic evidence and applied uniformly. The graph below "
                    "shows the original FBA prediction, the curated-model "
                    "prediction, and the corresponding Biolog result.",
                    styles["Normal"],
                )
            )

            story.append(Spacer(1, 8))

            # Curated-model comparison plot
            with tempfile.NamedTemporaryFile(
                suffix=".png",
                delete=False,
            ) as tmp:
                png = plot_curated_model_vs_biolog(
                    compare_rows,
                    curated_rows,
                    tmp.name,
                )

                if png:
                    story.append(Image(png, width=430, height=280, ))

            story.append(Spacer(1, 10))

            # Before/after table
            curated_by_condition = {
                r["condition"]: r
                for r in curated_rows
            }

            before_after = [[
                "Condition",
                "Original FBA",
                "Curated FBA",
                "Biolog",
                "Match?",
            ]]

            for r in compare_rows:
                cond = r["condition"]
                after = curated_by_condition.get(cond)

                if after is None:
                    continue

                before_after.append([
                    cond,
                    r["fba_call"],
                    after["fba_call"],
                    r["biolog_call"],
                    after["match"],
                ])

            story.append(make_table(before_after, col_widths=[80, 95, 95, 80, 70], ))
            story.append(PageBreak())
   
    if operon_rows:
        story.append(Paragraph("Operon gene check (rha/glc)", styles["Heading2"]))
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_operon_genes(operon_rows, tmp.name)
            story.append(Image(tmp.name, width=380, height=280))
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_operon_scatter(operon_rows, tmp.name)
            story.append(Image(tmp.name, width=340, height=280))
        story.append(PageBreak())
        operon_table = [["Gene", "Call", "Identity %", "Coverage %"]]
        for r in operon_rows:
            operon_table.append([r["gene"].split("|")[0], r["call"], r["pident"], r["qcovhsp"]])
        story.append(make_table(operon_table, col_widths=[70, 200, 80, 80]))

    doc = SimpleDocTemplate(args.output, pagesize=letter)
    doc.build(story)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
