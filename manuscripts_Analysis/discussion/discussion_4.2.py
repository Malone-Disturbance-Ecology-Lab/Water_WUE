# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 15:55:49 2026

@author: ammar
"""

import os
import pandas as pd
import numpy as np

d = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2_anomalies\Q2_anomalies_outputs"

summary = pd.read_csv(os.path.join(d, "Q2_anomalies_summary_by_class_metric.csv"))
signed = pd.read_csv(os.path.join(d, "Q2_anomalies_signed_response_fixed_effects.csv"))
magnitude = pd.read_csv(os.path.join(d, "Q2_anomalies_response_magnitude_fixed_effects.csv"))
wue = pd.read_csv(os.path.join(d, "Q2_anomalies_wue_value_fixed_effects.csv"))

print("\n=== SPEI-48 DRY MEANS ===")
print(
    summary[
        (summary["timescale"] == "SPEI_48") &
        (summary["anomaly_class"] == "Dry")
    ].to_string(index=False)
)

print("\n=== SIGNED RESPONSE: strongest / relevant coefficients ===")
print(
    signed[
        signed["term"].str.contains(
            "metric|timescale|SPEI_48|Saline",
            case=False, regex=True, na=False
        )
    ].sort_values("t_value").to_string(index=False)
)

print("\n=== RESPONSE MAGNITUDE: relevant coefficients ===")
print(
    magnitude[
        magnitude["term"].str.contains(
            "metric|timescale|SPEI_48|Saline",
            case=False, regex=True, na=False
        )
    ].sort_values("t_value").to_string(index=False)
)

print("\n=== WUE VALUE: relevant coefficients ===")
print(
    wue[
        wue["term"].str.contains(
            "metric|timescale|SPEI_48|Saline",
            case=False, regex=True, na=False
        )
    ].sort_values("t_value").to_string(index=False)
)




import os
import pandas as pd

d = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2_anomalies\Q2_anomalies_outputs"

files = {
    "SIGNED RESPONSE": "Q2_anomalies_signed_response_fixed_effects.csv",
    "RESPONSE MAGNITUDE": "Q2_anomalies_response_magnitude_fixed_effects.csv",
    "WUE VALUE": "Q2_anomalies_wue_value_fixed_effects.csv",
}

for label, fname in files.items():

    print("\n" + "=" * 100)
    print(label)
    print("=" * 100)

    df = pd.read_csv(os.path.join(d, fname))

    print("\nColumns:")
    print(list(df.columns))

    # Find term/name column
    term_col = None
    for c in ["term", "Term", "variable", "effect"]:
        if c in df.columns:
            term_col = c
            break

    if term_col is None:
        raise ValueError(f"Could not find term column in {fname}")

    # Show anything involving metric, SPEI48, timescale, or saline
    mask = (
        df[term_col].astype(str).str.contains(
            "metric|timescale|SPEI_48|SPEI48|Saline|saline",
            case=False,
            regex=True,
            na=False
        )
    )

    relevant = df.loc[mask].copy()

    print("\nRelevant coefficients:")
    if relevant.empty:
        print("No matching rows.")
    else:
        print(relevant.to_string(index=False))

print("\nDONE. Nothing was saved.")