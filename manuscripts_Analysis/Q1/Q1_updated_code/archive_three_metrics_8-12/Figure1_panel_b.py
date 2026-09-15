import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)

smooth_file = os.path.join(output_dir, "Q1_final_smooth_gam_predictions_TET.csv")
gam_comp_file = os.path.join(output_dir, "Q1_final_gam_model_comparison.csv")
edf_file = os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv")

fig_output = os.path.join(figure_dir, "Q1_panel_B_smooth_GAM.png")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
df = pd.read_csv(smooth_file)
gam_comp = pd.read_csv(gam_comp_file)
edf_df = pd.read_csv(edf_file)

required = ["difference_label", "ecosystem_label", "Trans_ratio",
            "prediction", "lower", "upper"]
if not all(col in df.columns for col in required):
    raise ValueError("Missing required columns in smooth predictions CSV")

# ------------------------------------------------------------------------------
# Settings (identical to Panel A)
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500",
}

diff_label_mapping = {
    "WUE_ET - WUE_T": r"WUE$_{ET}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_ET - WUE_E": r"WUE$_{ET}$ – WUE$_{E}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_E - WUE_T": r"WUE$_{E}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
}
diff_labels = list(diff_label_mapping.keys())

df = df[df["difference_label"].isin(diff_labels)]
df["difference_label"] = pd.Categorical(df["difference_label"],
                                        categories=diff_labels, ordered=True)
df["ecosystem_label"] = pd.Categorical(df["ecosystem_label"],
                                       categories=ecosystem_order, ordered=True)

# ------------------------------------------------------------------------------
# Extract adjusted R² for smooth GAM (no "Nonlinear GAM" text)
# ------------------------------------------------------------------------------
smooth_row = gam_comp[gam_comp["model"] == "q1_smooth_gam"].iloc[0]
smooth_r2 = smooth_row["adjusted_r_squared"]

# ------------------------------------------------------------------------------
# Parse edf values and get p‑value for T:ET smooths
# ------------------------------------------------------------------------------
edf_df = edf_df[edf_df["term"].str.contains(r"s\(Trans_ratio\)", regex=True)].copy()

# Get p-value across the nine T:ET smooth terms
p_col = [c for c in edf_df.columns if "p" in c.lower()][0]
max_tet_p = edf_df[p_col].max()

if max_tet_p < 0.001:
    smooth_p_text = "p < 0.001"
else:
    smooth_p_text = f"p = {max_tet_p:.3f}"

model_text = f"Adjusted R² = {smooth_r2:.2f}; {smooth_p_text}"

# ------------------------------------------------------------------------------
# Parse edf values for plotting (ecosystem‑specific)
# ------------------------------------------------------------------------------
def parse_smooth_term(term):
    clean = term.replace("s(Trans_ratio):diff_ecosystem", "")
    diff_type, ecosystem = clean.split("__")
    diff_map = {
        "WUE_ET_minus_WUE_T": "WUE_ET - WUE_T",
        "WUE_ET_minus_WUE_E": "WUE_ET - WUE_E",
        "WUE_E_minus_WUE_T": "WUE_E - WUE_T",
    }
    return pd.Series({
        "difference_label": diff_map[diff_type],
        "ecosystem_label": ecosystem
    })

edf_parsed = edf_df["term"].apply(parse_smooth_term)
edf_df = pd.concat([edf_df, edf_parsed], axis=1)

edf_lookup = {
    (row["difference_label"], row["ecosystem_label"]): row["edf"]
    for _, row in edf_df.iterrows()
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(24, 7), sharey=False)

plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

for ax in axes:
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)

# Loop over panels
for i, diff_label in enumerate(diff_labels):
    ax = axes[i]
    subset = df[df["difference_label"] == diff_label]

    # Plot curves and confidence bands
    for eco in ecosystem_order:
        eco_sub = subset[subset["ecosystem_label"] == eco]
        if eco_sub.empty:
            continue
        eco_sub = eco_sub.sort_values("Trans_ratio")
        x = eco_sub["Trans_ratio"].values
        y = eco_sub["prediction"].values
        lower = eco_sub["lower"].values
        upper = eco_sub["upper"].values
        color = ecosystem_colours[eco]

        ax.fill_between(x, lower, upper, color=color, alpha=0.15)
        ax.plot(x, y, color=color, linewidth=3)

    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)

    ax.set_xlabel("T:ET ratio", fontsize=28)
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=28, labelpad=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Panel label (d), (e), (f) – top‑left corner
    ax.text(-0.15, 1.02, f"{chr(100 + i)})", transform=ax.transAxes,
            fontsize=34, fontweight='bold', va='bottom', ha='left')

    # --------------------------------------------------------------
    # edf text box – bottom‑right, widened (rect_width = 0.70)
    # --------------------------------------------------------------
    lines = []
    for eco in ecosystem_order:
        edf_val = edf_lookup.get((diff_label, eco), np.nan)
        if not np.isnan(edf_val):
            lines.append((f"{eco}: edf = {edf_val:.2f}", ecosystem_colours[eco]))
        else:
            lines.append((f"{eco}: edf = NA", ecosystem_colours[eco]))

    # Position: shift right and increase width
    x0 = 0.90               # moved right
    y0 = 0.05
    line_height = 0.06
    rect_width = 0.70       # wider box
    rect_height = len(lines) * line_height + 0.04

    rect = patches.Rectangle(
        (x0 - rect_width, y0),
        rect_width,
        rect_height,
        transform=ax.transAxes,
        facecolor='white',
        edgecolor='black',
        linewidth=1.2,
        alpha=0.85,
        clip_on=False
    )
    ax.add_patch(rect)

    # Right‑aligned coloured text (fontsize 23)
    for j, (line, color) in enumerate(lines):
        ax.text(
            x0 - 0.02,
            y0 + rect_height - j * line_height - 0.02,
            line,
            transform=ax.transAxes,
            ha='right',
            va='top',
            fontsize=23,
            color=color,
            clip_on=False
        )

# --------------------------------------------------------------
# Shared R² & p‑value label – moved high (y=0.985)
# --------------------------------------------------------------
fig.text(
    0.53, 0.985,
    model_text,
    ha="center",
    va="top",
    fontsize=24,
    bbox=dict(
        boxstyle="round,pad=0.35",
        facecolor="white",
        edgecolor="black",
        linewidth=1.5,
        alpha=0.9
    )
)

# Adjust margins: top=0.92 gives space for the label
plt.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.12, wspace=0.35)

# Save
fig.savefig(fig_output, dpi=600, bbox_inches='tight', pad_inches=0.5, facecolor='white')
print(f"Panel B (Smooth GAM) figure saved to: {fig_output}")