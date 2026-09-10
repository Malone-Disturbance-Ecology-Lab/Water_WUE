#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure_Sx_5pct_threshold_evaluation.py

Creates a three‑panel supplementary figure and a cleaned supplementary table
from the saved 5‑percent threshold diagnostic outputs and the published Table 5 medians.

Reads:
- effect_size_threshold_sensitivity_5_to_75pct.csv  (for panel c)
- 5pct_threshold_month_averaged_summary.csv         (for panel b counts)

No analysis is rerun; this is pure downstream plotting and table formatting.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ============================================================================
# PATHS AND OUTPUT
# ============================================================================
BASE_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\5_percent_threshold_justification"

INPUT_SENS = os.path.join(BASE_DIR, "effect_size_threshold_sensitivity_5_to_75pct.csv")
INPUT_5PCT = os.path.join(BASE_DIR, "5pct_threshold_month_averaged_summary.csv")
INPUT_SUPP = os.path.join(BASE_DIR, "Table_Sx_5pct_threshold_justification.csv")

OUTPUT_FIG_PNG = os.path.join(BASE_DIR, "Figure_Sx_5pct_threshold_evaluation.png")
OUTPUT_TABLE = os.path.join(BASE_DIR, "Table_Sx_5pct_threshold_evaluation.csv")

# ============================================================================
# FORMATTING VARIABLES
# ============================================================================
FONT_FAMILY = "Arial"
BASE_FONT_SIZE = 12
AXIS_LABEL_SIZE = 16
TICK_SIZE = 13
PANEL_LABEL_SIZE = 18
FIG_WIDTH = 15
FIG_HEIGHT = 5.6
LINE_WIDTH = 1.3
MARKER_SIZE = 6
DPI = 600

# Colors (colorblind‑friendly)
COLOR_LT = "#4B8BBE"   # light blue
COLOR_GE = "#D55E00"   # orange
COLOR_BAR = "#4B8BBE"  # single color for panel (b)
COLOR_LINE = "#4B8BBE" # line colour for panel c

# ============================================================================
# READ DATA
# ============================================================================
print("Reading saved diagnostic CSVs...")
df_sens = pd.read_csv(INPUT_SENS)
df_5pct = pd.read_csv(INPUT_5PCT)
df_supp = pd.read_csv(INPUT_SUPP)
print(f"  sensitivity: {len(df_sens)} rows")
print(f"  5% summary:  {len(df_5pct)} rows")
print(f"  supp table:  {len(df_supp)} rows")

# ============================================================================
# PANEL A DATA – Table 5 medians
# ============================================================================
panel_a_data = {
    "Metric": ["Stability", "Plasticity range"],
    "<5%": [6.03, 2.11],
    "≥5%": [1.10, 3.13],
    "p": [0.03, 0.03]
}
df_a = pd.DataFrame(panel_a_data)

# ============================================================================
# PANEL B DATA – statistical evidence categories
# ============================================================================
status_counts = df_5pct["status"].value_counts()
supported_at_5 = status_counts.get("supported_at_5pct", 0)
support_before = status_counts.get("support_before_5pct", 0)
no_stat = status_counts.get("no_stat_support", 0)

# Clearer, concise labels
panel_b_categories = [
    "Significant at\n5% response",
    "Significant before\n5% response",
    "Not significant"
]
panel_b_counts = [supported_at_5, support_before, no_stat]

# ============================================================================
# PANEL C DATA – response magnitudes evaluated
# ============================================================================
threshold_levels = [5, 10, 15, 20, 25, 30, 35, 40, 50, 75]
counts = df_sens.groupby("threshold")["n_crossings"].sum().reindex(threshold_levels, fill_value=0)
x_c = counts.index
y_c = counts.values

# ============================================================================
# CREATE FIGURE
# ============================================================================
plt.rcParams["font.family"] = FONT_FAMILY

fig, axes = plt.subplots(1, 3, figsize=(FIG_WIDTH, FIG_HEIGHT))
ax_a, ax_b, ax_c = axes

# ---- Panel (a) ----
x = np.arange(len(df_a["Metric"]))
width = 0.35

ax_a.bar(x - width/2, df_a["<5%"], width, label="<5%", color=COLOR_LT, edgecolor="black", linewidth=0.8)
ax_a.bar(x + width/2, df_a["≥5%"], width, label="≥5%", color=COLOR_GE, edgecolor="black", linewidth=0.8)

ax_a.set_ylabel("Site-level metric", fontsize=AXIS_LABEL_SIZE)
ax_a.set_xticks(x)
ax_a.set_xticklabels(df_a["Metric"], fontsize=TICK_SIZE)
ax_a.tick_params(axis="both", labelsize=TICK_SIZE)
ax_a.legend(loc="upper right", fontsize=BASE_FONT_SIZE, framealpha=0.9)
ax_a.grid(axis="y", linestyle="--", alpha=0.3)
ax_a.set_ylim(0, 7.5)

for i, (pval, y1, y2) in enumerate(zip(df_a["p"], df_a["<5%"], df_a["≥5%"])):
    ymax = max(y1, y2) + 0.3
    ax_a.text(i, ymax, f"p = {pval:.2f}", ha="center", va="bottom",
              fontsize=BASE_FONT_SIZE, fontweight="bold")

ax_a.text(-0.12, 1.02, "(a)", transform=ax_a.transAxes, fontsize=PANEL_LABEL_SIZE,
          fontweight="bold", va="bottom", ha="left")

# ---- Panel (b) ----
bars_b = ax_b.bar(panel_b_categories, panel_b_counts, color=COLOR_BAR,
                  edgecolor="black", linewidth=0.8)
ax_b.set_ylabel("Number of modelled dry-side responses", fontsize=AXIS_LABEL_SIZE)
ax_b.tick_params(axis="both", labelsize=TICK_SIZE)
ax_b.yaxis.set_major_locator(MaxNLocator(integer=True))
ax_b.grid(axis="y", linestyle="--", alpha=0.3)
ax_b.set_ylim(0, 25)

for bar, cnt in zip(bars_b, panel_b_counts):
    if cnt > 0:
        ax_b.text(bar.get_x() + bar.get_width()/2, cnt + 0.4, str(cnt),
                  ha="center", va="bottom", fontsize=BASE_FONT_SIZE, fontweight="bold")

ax_b.text(-0.12, 1.02, "(b)", transform=ax_b.transAxes, fontsize=PANEL_LABEL_SIZE,
          fontweight="bold", va="bottom", ha="left")

# ---- Panel (c) ----
ax_c.plot(x_c, y_c, marker="o", markersize=MARKER_SIZE, color=COLOR_LINE,
          linewidth=LINE_WIDTH, zorder=2)
ax_c.set_xlabel("WUE$_T$ response magnitude (%)", fontsize=AXIS_LABEL_SIZE)
ax_c.set_ylabel("Number of modelled dry-side responses", fontsize=AXIS_LABEL_SIZE)
ax_c.set_xticks(x_c)
ax_c.set_xticklabels([str(v) for v in x_c], fontsize=TICK_SIZE, rotation=45, ha="right")
ax_c.tick_params(axis="both", labelsize=TICK_SIZE)
ax_c.yaxis.set_major_locator(MaxNLocator(integer=True))
ax_c.grid(axis="y", linestyle="--", alpha=0.3, zorder=1)
ax_c.set_ylim(-1, 46)

for xi, yi in zip(x_c, y_c):
    if yi > 0:
        ax_c.text(xi + 0.8, yi + 0.5, str(int(yi)), ha="left", va="bottom",
                  fontsize=BASE_FONT_SIZE, fontweight="bold")

ax_c.text(-0.12, 1.02, "(c)", transform=ax_c.transAxes, fontsize=PANEL_LABEL_SIZE,
          fontweight="bold", va="bottom", ha="left")

# ---- Adjust layout and save ----
plt.tight_layout(pad=2.0)

print("\nSaving figure...")
plt.savefig(OUTPUT_FIG_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")

print("\nDisplaying figure...")
plt.show()

# ============================================================================
# SUPPLEMENTARY TABLE – cleaned and formatted
# ============================================================================
print("\nCreating supplementary table...")

df_table = df_supp.copy()

df_table = df_table.rename(columns={
    "Coast": "Coastal region",
    "Ecosystem": "Ecosystem",
    "SPEI timescale": "SPEI timescale",
    "Response direction": "Response direction",
    "SPEI at 5% response": "SPEI at 5% WUE_T response",
    "SPEI at statistical-support threshold": "SPEI where response became statistically distinguishable from near-normal",
    "Difference in SPEI": "Difference in SPEI",
    "95% CI supports departure from near-normal at 5% crossing": "Statistically distinguishable from near-normal at 5% response"
})

df_table["Statistically distinguishable from near-normal at 5% response"] = df_table[
    "Statistically distinguishable from near-normal at 5% response"
].replace({True: "Yes", False: "No"})

for col in ["SPEI at 5% WUE_T response",
            "SPEI where response became statistically distinguishable from near-normal",
            "Difference in SPEI"]:
    df_table[col] = df_table[col].round(2)

coast_mapping = {
    "AK Coast": "Alaska",
    "Atlantic Coast": "Atlantic",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf"
}
df_table["Coastal region"] = df_table["Coastal region"].replace(coast_mapping)

df_table["SPEI timescale"] = df_table["SPEI timescale"].str.replace("_", "-")

df_table["SPEI where response became statistically distinguishable from near-normal"] = df_table[
    "SPEI where response became statistically distinguishable from near-normal"
].fillna("—")

df_table = df_table.sort_values(["Coastal region", "Ecosystem", "SPEI timescale", "Response direction"])

df_table.to_csv(OUTPUT_TABLE, index=False,encoding="utf-8-sig")
print(f"  Table saved: {OUTPUT_TABLE}")

# ============================================================================
# CONSOLE SUMMARY
# ============================================================================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Panel A – ecological relevance (Table 5 medians):")
print("  Stability: <5% = 6.03, ≥5% = 1.10, p = 0.03")
print("  Plasticity range: <5% = 2.11, ≥5% = 3.13, p = 0.03")
print(f"\nPanel B – statistical evidence categories:")
print(f"  Significant at 5%: {supported_at_5}")
print(f"  Significant before 5%: {support_before}")
print(f"  Not significant: {no_stat}")
print(f"\nPanel C – response counts:")
for thr, cnt in zip(threshold_levels, counts.values):
    print(f"  {thr}%: {int(cnt)} responses")
print(f"\nGenerated files:")
print(f"  Figure PNG: {OUTPUT_FIG_PNG}")
print(f"  Table CSV:  {OUTPUT_TABLE}")
print("="*60)