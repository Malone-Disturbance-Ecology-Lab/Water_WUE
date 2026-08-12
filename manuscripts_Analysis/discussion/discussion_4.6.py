# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 13:24:07 2026

@author: ammar
"""

import os
import pandas as pd
import re

OUT = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_outputs"

summary_file = os.path.join(
    OUT,
    "EDI_logistic_model_summary.txt"
)

smooth_file = os.path.join(
    OUT,
    "EDI_response_gam_smooth_terms.csv"
)

print("=" * 90)
print("CHECKING POSSIBLE SOURCE OF THE 4.2% RESULT")
print("=" * 90)

# ============================================================
# 1. RESPONSE GAM SUMMARY
# ============================================================

print("\n1. RESPONSE-BASED EDI GAM SUMMARY")
print("-" * 90)

if not os.path.exists(summary_file):
    print("Summary file NOT FOUND:")
    print(summary_file)
else:
    with open(summary_file, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    # Find the response-GAM section only
    start = text.find("=== Response GAM summary ===")

    if start >= 0:
        section = text[start:]
    else:
        section = text

    lines = section.splitlines()

    wanted = [
        "R-sq.(adj)",
        "Deviance explained",
        "Scale est.",
        "n =",
        "s(SPEI_3",
        "s(site_name",
        "Approximate significance of smooth terms"
    ]

    for i, line in enumerate(lines):
        if any(x.lower() in line.lower() for x in wanted):
            lo = max(0, i - 1)
            hi = min(len(lines), i + 2)

            for j in range(lo, hi):
                print(lines[j])

            print("-" * 60)

# ============================================================
# 2. SAVED SMOOTH TERMS
# ============================================================

print("\n2. RESPONSE GAM SMOOTH TERMS")
print("-" * 90)

if not os.path.exists(smooth_file):
    print("Smooth-term file NOT FOUND:")
    print(smooth_file)
else:
    smooth = pd.read_csv(smooth_file)

    print("Columns:")
    print(list(smooth.columns))

    print("\nAll smooth terms:")
    print(smooth.to_string(index=False))

    print("\nSite random-effect term:")
    site = smooth[
        smooth.iloc[:, 0].astype(str).str.contains(
            "site_name",
            case=False,
            na=False
        )
    ]

    if len(site):
        print(site.to_string(index=False))
    else:
        print("No site term found.")

print("\n" + "=" * 90)
print("DONE. NOTHING WAS SAVED.")
print("=" * 90)