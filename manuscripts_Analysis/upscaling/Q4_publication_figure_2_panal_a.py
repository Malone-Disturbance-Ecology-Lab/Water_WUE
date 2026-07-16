# -*- coding: utf-8 -*-
"""
Chunk 1: Panel a – Event-frequency maps (with 5 x‑ticks)
- Two maps: negative events (WUE_T <= -5%) and positive events (WUE_T >= +5%)
- Panel label "a)" at (-0.10, 1.02)
- High resolution (900 dpi), wspace=0
- Saves as 'EDI_panel_A_maps.png' in talib_publication
"""

import csv
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import netCDF4 as nc
import numpy as np

# ---------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------
SPATIAL_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
OUT_DIR = SPATIAL_DIR / "talib_publication"
OUT_DIR.mkdir(parents=True, exist_ok=True)

neg_path = SPATIAL_DIR / "EDI_response_v2_freq_strong_decrease_upland.nc"
pos_path = SPATIAL_DIR / "EDI_response_v2_freq_strong_increase_upland.nc"
OUT_FIG = OUT_DIR / "EDI_panel_A_maps.png"

# ---------------------------------------------------------------------
# COAST REGION LABEL POSITIONS
# ---------------------------------------------------------------------
REGION_LABEL_POS = {
    "AK Coast": (-167.0, 73.0, "Alaska"),
    "Pacific Coast": (-137.0, 41.0, "Pacific"),
    "Gulf Coast": (-101.0, 22.0, "Gulf"),
    "Atlantic Coast": (-82.0, 50.0, "Atlantic"),
}

# ---------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------
def read_raster(path):
    ds = nc.Dataset(path)
    lat = np.array(ds.variables["lat"][:])
    lon = np.array(ds.variables["lon"][:])
    var_name = [v for v in ds.variables if v not in ("lat", "lon")][0]
    arr = np.array(ds.variables[var_name][:], dtype=float)
    fill = getattr(ds.variables[var_name], "_FillValue", None)
    ds.close()
    if fill is not None:
        arr[arr == fill] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    return lon, lat, arr

def style_map_axis(ax, label_fontsize=10):
    ax.set_xlim(-180, -60)
    ax.set_ylim(17, 77)
    ax.set_xlabel("Longitude", fontsize=label_fontsize)
    ax.set_ylabel("Latitude", fontsize=label_fontsize)
    ax.set_facecolor("white")
    ax.grid(False)
    for region, (x, y, label) in REGION_LABEL_POS.items():
        ax.text(
            x, y, label,
            ha="center", va="center",
            fontsize=7, fontweight="bold", color="black",
        )

# ---------------------------------------------------------------------
# READ DATA
# ---------------------------------------------------------------------
lon, lat, neg = read_raster(neg_path)
_, _, pos = read_raster(pos_path)

# ---------------------------------------------------------------------
# CREATE PANEL A
# ---------------------------------------------------------------------
fig, (ax_neg, ax_pos) = plt.subplots(1, 2, figsize=(6.5, 3), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)

# Colormaps
neg_cmap = mcolors.LinearSegmentedColormap.from_list(
    "negative_events", ["#FFFFFF", "#D7ECF7", "#73B3D8", "#2879B9", "#08306B"]
)
pos_cmap = mcolors.LinearSegmentedColormap.from_list(
    "positive_events", ["#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#7F0000"]
)
neg_cmap.set_bad((1, 1, 1, 0))
pos_cmap.set_bad((1, 1, 1, 0))

neg_vmin, neg_vmax = np.nanpercentile(neg, [2, 98])
pos_vmin, pos_vmax = np.nanpercentile(pos, [2, 98])

im_neg = ax_neg.pcolormesh(lon, lat, neg, cmap=neg_cmap, vmin=neg_vmin, vmax=neg_vmax, shading="auto")
im_pos = ax_pos.pcolormesh(lon, lat, pos, cmap=pos_cmap, vmin=pos_vmin, vmax=pos_vmax, shading="auto")

style_map_axis(ax_neg, label_fontsize=10)
style_map_axis(ax_pos, label_fontsize=10)
ax_neg.tick_params(labelsize=10)
ax_pos.tick_params(labelsize=10)
ax_pos.set_ylabel("")

# --- SET X‑AXIS TICKS: at least 5 ticks ---
ax_neg.set_xticks(np.arange(-180, -30, 30))  # gives -180, -150, -120, -90, -60
# ------------------------------------------

ax_neg.set_title("Negative events: WUE$_T$ ≤ -5%", fontweight="bold", fontsize=8)
ax_pos.set_title("Positive events: WUE$_T$ ≥ +5%", fontweight="bold", fontsize=8)

# Color bars – keep original tick formatting
cb_neg = fig.colorbar(im_neg, ax=ax_neg, fraction=0.045, pad=0.02)
cb_pos = fig.colorbar(im_pos, ax=ax_pos, fraction=0.045, pad=0.02)
cb_neg.set_label("% months", fontsize=10, labelpad=2)
cb_pos.set_label("% months", fontsize=10, labelpad=2)
cb_neg.ax.tick_params(labelsize=8)
cb_pos.ax.tick_params(labelsize=8)

# Panel label – moved to (-0.10, 1.02)
ax_neg.text(-0.10, 1.02, "a)", transform=ax_neg.transAxes, fontsize=10, fontweight="bold")

fig.tight_layout()

# ---------------------------------------------------------------------
# SAVE AND DISPLAY
# ---------------------------------------------------------------------
fig.savefig(OUT_FIG, dpi=900, bbox_inches="tight", facecolor="white")
print(f"Panel A saved to: {OUT_FIG}")

plt.show(block=True)
print("Panel A displayed.")