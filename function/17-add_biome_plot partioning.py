# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 14:26:43 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
Merge site metadata (IGBP, climate, water_class) into ET partitioning files
No plotting - just add columns to existing CSV files
"""

import os
import glob
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

PRINT = "[ET_Merge]"

def load_site_metadata(site_info_csv):
    """Load site metadata from CSV file"""
    # Try different encodings
    encodings = ['utf-8', 'latin1', 'cp1252', 'iso-8859-1', 'windows-1252']
    
    info = None
    for encoding in encodings:
        try:
            info = pd.read_csv(site_info_csv, encoding=encoding)
            print(PRINT, f"Successfully read metadata with encoding: {encoding}")
            break
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(PRINT, f"Error with encoding {encoding}: {e}")
            continue
    
    if info is None:
        raise ValueError(f"Could not read {site_info_csv} with any common encoding")
    
    # Select only needed columns
    metadata_cols = ['site_name', 'IGBP', 'climate', 'water_class']
    available_cols = [col for col in metadata_cols if col in info.columns]
    
    if 'site_name' not in available_cols:
        raise ValueError("site_name column not found in metadata file")
    
    metadata = info[available_cols].copy()
    print(PRINT, f"Loaded metadata for {len(metadata)} sites")
    print(PRINT, f"Columns to add: {[col for col in available_cols if col != 'site_name']}")
    
    return metadata

def merge_metadata_into_files(input_dir, metadata):
    """Merge IGBP, climate, water_class into each ET partitioning CSV file"""
    files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    print(PRINT, f"Found {len(files)} files to process")
    
    updated_count = 0
    skipped_count = 0
    
    for file_path in files:
        site_name = Path(file_path).stem
        print(f"\n{PRINT} Processing {site_name}...")
        
        # Find matching metadata
        site_metadata = metadata[metadata['site_name'] == site_name]
        
        if site_metadata.empty:
            print(f"  ⚠ No metadata found for {site_name}, skipping")
            skipped_count += 1
            continue
        
        try:
            # Read the ET partitioning file
            df = pd.read_csv(file_path, engine='python')
            original_cols = df.columns.tolist()
            print(f"  Original shape: {df.shape}")
            print(f"  Original columns: {len(original_cols)}")
            
            # Add metadata columns (as constants)
            added_cols = []
            for col in ['IGBP', 'climate', 'water_class']:
                if col in site_metadata.columns:
                    value = site_metadata[col].iloc[0]
                    if col not in df.columns:
                        df[col] = value
                        added_cols.append(col)
                        print(f"  + Added column '{col}': {value}")
                    else:
                        print(f"  Column '{col}' already exists, skipping")
            
            # Save back to same file (overwrite)
            df.to_csv(file_path, index=False)
            print(f"  ✓ Updated {site_name}.csv (added {len(added_cols)} columns)")
            updated_count += 1
            
        except Exception as e:
            print(f"  ✗ Error processing {site_name}: {e}")
            skipped_count += 1
    
    print(f"\n{PRINT} {'='*50}")
    print(f"Merge completed!")
    print(f"  Updated: {updated_count} files")
    print(f"  Skipped: {skipped_count} files")
    print(f"{'='*50}")
    
    return updated_count

def verify_update(input_dir, sample_site=None):
    """Verify that columns were added correctly"""
    if sample_site is None:
        # Get first CSV file as sample
        files = glob.glob(os.path.join(input_dir, "*.csv"))
        if not files:
            print(PRINT, "No files found to verify")
            return
        sample_file = files[0]
        site_name = Path(sample_file).stem
    else:
        sample_file = os.path.join(input_dir, f"{sample_site}.csv")
        site_name = sample_site
    
    try:
        df = pd.read_csv(sample_file, engine='python')
        print(f"\n{PRINT} Verification for {site_name}:")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {df.columns.tolist()}")
        
        # Check if new columns exist
        for col in ['IGBP', 'climate', 'water_class']:
            if col in df.columns:
                unique_values = df[col].dropna().unique()
                print(f"  ✓ '{col}' present: {unique_values[0] if len(unique_values) > 0 else 'All NaN'}")
            else:
                print(f"  ✗ '{col}' NOT found")
    except Exception as e:
        print(f"  Error verifying {site_name}: {e}")

def main():
    """Main function to run the metadata merging"""
    
    # Define paths
    site_info_csv = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
    partitioning_dir = r"M:\Research\WUE_CUE\ameri_data\ET_partitioning"
    
    print(PRINT, "="*60)
    print(PRINT, "Merging metadata into ET partitioning files")
    print(PRINT, "="*60)
    print(PRINT, f"Metadata file: {site_info_csv}")
    print(PRINT, f"Target directory: {partitioning_dir}")
    print(PRINT, "="*60)
    
    # Load site metadata
    metadata = load_site_metadata(site_info_csv)
    
    # Merge metadata into existing files
    merge_metadata_into_files(partitioning_dir, metadata)
    
    # Verify with a sample site (US-A03 as example)
    verify_update(partitioning_dir, "US-A03")
    
    print(PRINT, "\n" + "="*60)
    print(PRINT, "COMPLETED! All CSV files have been updated with IGBP, climate, and water_class columns")
    print(PRINT, "="*60)

if __name__ == "__main__":
    main()
    
#########################################################################################################

# -*- coding: utf-8 -*-
"""
Boxplots of yearly evaporation/transpiration fractions from partitioned ET.
Data is already in growing season (no need for additional filtering)

Inputs per CSV (one per site, in the ET_partitioning folder):
  - DateTime, ET, Evap_pen, Trans_pen, IGBP, climate, water_class (half-hourly)

Processing:
  - Aggregate to annual sums per site-year for ET, Evap_pen, Trans_pen
  - Compute yearly fractions:
        evap_frac  = Evap_pen / ET
        trans_frac = Trans_pen / ET
    (handle div-by-zero / inf → NaN, clip to [0, 1])

Outputs (saved to ET_partitioning/plots folder):
  - evap_trans_by_site_box.png
  - evap_trans_by_biome_box.png
  - evap_trans_by_climate_box.png
  - evap_trans_by_water_class_box.png
  - summary_statistics.csv
"""

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

PRINT = "[ET_Partition_Plots]"

# Update these paths to your new data location
ET_PARTITIONING_DIR = r"M:\Research\WUE_CUE\ameri_data\ET_partitioning"
PLOTS_OUTPUT_DIR = r"M:\Research\WUE_CUE\ameri_data\ET_partitioning\plots"

# -------------------------------
# Data cleaning helpers
# -------------------------------
def clean_metadata_value(value):
    """Clean and standardize metadata values (remove extra spaces, standardize case)"""
    if pd.isna(value) or value == "Unknown":
        return "Unknown"
    
    # Convert to string and strip
    value = str(value).strip()
    
    # Standardize water classes
    water_class_mapping = {
        'freshwater': 'Freshwater',
        'fresh water': 'Freshwater',
        'Freshwater': 'Freshwater',
        'upland': 'Upland',
        'Upland': 'Upland',
        'UPLAND': 'Upland',
        'saline': 'Saline',
        'Saline': 'Saline',
        'SALINE': 'Saline'
    }
    
    # Apply mapping if it's a water class
    value_lower = value.lower()
    if value_lower in water_class_mapping:
        return water_class_mapping[value_lower]
    
    # Standardize IGBP classes
    igbp_mapping = {
        'bsv': 'BSV',
        'wet': 'WET',
        'osh': 'OSH',
        'gra': 'GRA',
        'csh': 'CSH',
        'mf': 'MF',
        'enf': 'ENF',
        'cro': 'CRO',
        'dbf': 'DBF'
    }
    
    if value_lower in igbp_mapping:
        return igbp_mapping[value_lower]
    
    return value

def clean_dataframe_metadata(df):
    """Clean metadata columns in the dataframe"""
    if 'IGBP' in df.columns:
        df['IGBP'] = df['IGBP'].apply(clean_metadata_value)
    
    if 'climate' in df.columns:
        df['climate'] = df['climate'].apply(clean_metadata_value)
    
    if 'water_class' in df.columns:
        df['water_class'] = df['water_class'].apply(clean_metadata_value)
    
    return df

# -------------------------------
# I/O and aggregation helpers
# -------------------------------
def _read_partition_long(et_dir):
    """Read all CSV files and return long df with site, DateTime, ET, Evap_pen, Trans_pen, and metadata."""
    et_dir = Path(et_dir)
    files = sorted(glob.glob(str(et_dir / "*.csv")))
    
    if not files:
        print(PRINT, f"No CSV files found in {et_dir}")
        return pd.DataFrame()
    
    print(PRINT, f"Found {len(files)} files to process")
    
    rows = []
    skipped_sites = []
    
    for f in files:
        site = Path(f).stem
        try:
            df = pd.read_csv(f)
            
            # Check required columns
            need = ["DateTime", "ET", "Evap_pen", "Trans_pen"]
            if any(col not in df.columns for col in need):
                print(PRINT, f"  Skipping {site} (missing one of {need})")
                skipped_sites.append(site)
                continue
            
            # Convert to numeric and datetime
            dt = pd.to_datetime(df["DateTime"], errors="coerce")
            et = pd.to_numeric(df["ET"], errors="coerce")
            e_p = pd.to_numeric(df["Evap_pen"], errors="coerce")
            t_p = pd.to_numeric(df["Trans_pen"], errors="coerce")
            
            # Get metadata if available - clean them
            igbp = df["IGBP"].iloc[0] if "IGBP" in df.columns else "Unknown"
            climate = df["climate"].iloc[0] if "climate" in df.columns else "Unknown"
            water_class = df["water_class"].iloc[0] if "water_class" in df.columns else "Unknown"
            
            # Clean metadata values
            igbp = clean_metadata_value(igbp)
            climate = clean_metadata_value(climate)
            water_class = clean_metadata_value(water_class)
            
            # Only keep rows with valid data
            m = dt.notna() & et.notna() & e_p.notna() & t_p.notna()
            if not m.any():
                print(PRINT, f"  Skipping {site} (no valid rows after cleaning)")
                skipped_sites.append(site)
                continue
            
            rows.append(pd.DataFrame({
                "site": site,
                "IGBP": igbp,
                "climate": climate,
                "water_class": water_class,
                "DateTime": dt[m],
                "ET": et[m],
                "Evap_pen": e_p[m],
                "Trans_pen": t_p[m]
            }))
            print(PRINT, f"  Loaded {site}: {m.sum()} valid rows")
            
        except Exception as e:
            print(PRINT, f"  [ERROR] {site}: {e}")
            skipped_sites.append(site)
    
    if skipped_sites:
        print(PRINT, f"Skipped {len(skipped_sites)} sites: {skipped_sites[:5]}..." if len(skipped_sites) > 5 else f"Skipped sites: {skipped_sites}")
    
    if rows:
        result = pd.concat(rows, ignore_index=True)
        print(PRINT, f"Total loaded: {len(result)} rows from {len(rows)} sites")
        return result
    
    return pd.DataFrame(columns=["site", "IGBP", "climate", "water_class", "DateTime", "ET", "Evap_pen", "Trans_pen"])

def _annual_sums(df_long):
    """Annual sums per site for ET, Evap_pen, Trans_pen."""
    if df_long.empty:
        return pd.DataFrame()
    
    d = df_long.copy()
    d["year"] = pd.to_datetime(d["DateTime"]).dt.year
    
    # Aggregate by site and year, keeping metadata
    agg = (d.groupby(["site", "IGBP", "climate", "water_class", "year"], as_index=False)
             .agg(ET=("ET", "sum"),
                  Evap_pen=("Evap_pen", "sum"),
                  Trans_pen=("Trans_pen", "sum"),
                  n_obs=("ET", "count")))
    
    # Fractions with safe division
    with np.errstate(divide="ignore", invalid="ignore"):
        agg["evap_frac"] = agg["Evap_pen"] / agg["ET"]
        agg["trans_frac"] = agg["Trans_pen"] / agg["ET"]
    
    # Clean inf/-inf and clip to [0,1]
    for c in ("evap_frac", "trans_frac"):
        v = agg[c].replace([np.inf, -np.inf], np.nan)
        agg[c] = v.clip(lower=0.0, upper=1.0)
    
    return agg

# -------------------------------
# Plotting helpers
# -------------------------------
def _wrap_labels(labels, max_chars=14, max_lines=3):
    """Wrap long labels onto <= max_lines with ~max_chars per line."""
    wrapped = []
    for lab in labels:
        parts = str(lab).split()
        if not parts:
            wrapped.append(lab)
            continue
        lines, cur = [], ""
        for w in parts:
            add = (w if cur == "" else " " + w)
            if len(cur + add) <= max_chars:
                cur += add
            else:
                lines.append(cur.strip())
                cur = w
                if len(lines) >= max_lines - 1:  # leave last line for rest
                    break
        if cur and len(lines) < max_lines:
            lines.append(cur.strip())
        # If overflowed, join rest of words into last line with an ellipsis
        if len(lines) == max_lines and len(parts) > 0 and (" ".join(parts) != " ".join(" ".join(lines).split())):
            lines[-1] = (lines[-1] + " …")
        wrapped.append("\n".join(lines))
    return wrapped

def _two_boxplots(df, group_col, out_path, title_prefix, x_wrap=False):
    """
    Make a single figure with two boxplots (evap_frac, trans_frac) grouped by group_col.
    """
    if df.empty:
        print(PRINT, f"No data for {group_col} plot")
        return
    
    # Order groups by median evap_frac (descending) for readability
    med = (df.groupby(group_col)["evap_frac"]
             .median()
             .sort_values(ascending=False))
    groups = med.index.tolist()
    
    # Filter out groups with too few observations
    min_obs = 3
    valid_groups = []
    evap_data_list = []
    trans_data_list = []
    
    for g in groups:
        evap_vals = df.loc[df[group_col] == g, "evap_frac"].dropna()
        trans_vals = df.loc[df[group_col] == g, "trans_frac"].dropna()
        if len(evap_vals) >= min_obs:
            valid_groups.append(g)
            evap_data_list.append(evap_vals.values)
            trans_data_list.append(trans_vals.values)
        else:
            print(PRINT, f"  Excluding {g} (only {len(evap_vals)} observations)")
    
    if not valid_groups:
        print(PRINT, f"No valid groups for {group_col} plot (need ≥{min_obs} observations)")
        return
    
    fig = plt.figure(figsize=(max(10, 0.5*len(valid_groups)+6), 8))
    
    # Evaporation fraction
    ax1 = fig.add_subplot(2,1,1)
    b1 = ax1.boxplot(evap_data_list, patch_artist=True, widths=0.6, showfliers=False)
    for box in b1['boxes']:
        box.set_facecolor('lightblue')
        box.set_alpha(0.7)
    ax1.set_title(f"{title_prefix}: Evaporation Fraction (Evap_pen / ET)", fontsize=12)
    ax1.set_ylabel("Fraction of ET", fontsize=10)
    ax1.set_xticklabels([])
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Transpiration fraction
    ax2 = fig.add_subplot(2,1,2)
    b2 = ax2.boxplot(trans_data_list, patch_artist=True, widths=0.6, showfliers=False)
    for box in b2['boxes']:
        box.set_facecolor('lightgreen')
        box.set_alpha(0.7)
    ax2.set_title(f"{title_prefix}: Transpiration Fraction (Trans_pen / ET)", fontsize=12)
    ax2.set_ylabel("Fraction of ET", fontsize=10)
    
    # X tick labels (shared order)
    labels = valid_groups
    if x_wrap:
        labels = _wrap_labels(valid_groups, max_chars=16, max_lines=3)
    ax2.set_xticklabels(labels, rotation=0, ha="center")
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=220, bbox_inches='tight')
    plt.close()
    print(PRINT, f"  Saved: {out_path.name}")

def _round_to_significant_digits(value, min_digits=2, max_digits=4):
    """Round a number to between min_digits and max_digits significant digits"""
    if pd.isna(value):
        return value
    
    # Handle zero
    if value == 0:
        return 0.0
    
    # Calculate number of significant digits to use (between min and max)
    # Use more digits for smaller numbers to maintain precision
    abs_val = abs(value)
    if abs_val < 0.01:
        sig_digits = max_digits
    elif abs_val < 0.1:
        sig_digits = max_digits - 1
    else:
        sig_digits = min(max_digits, max(min_digits, max_digits - int(np.log10(abs_val))))
    
    # Round to specified significant digits
    return round(value, sig_digits - int(np.floor(np.log10(abs_val))) - 1)

def _save_summary_statistics(df_annual, output_dir):
    """Save summary statistics to CSV - ROBUST version with proper length handling"""
    if df_annual.empty:
        return
    
    # Calculate statistics by different groupings
    stats_list = []
    
    # Helper function to create stats for a subset
    def create_stats_rows(group_name, subset):
        rows = []
        
        # Add Evap_Fraction row if data exists
        if subset['evap_frac'].notna().any():
            rows.append({
                'Group': group_name,
                'Variable': 'Evap_Fraction',
                'Mean': _round_to_significant_digits(subset['evap_frac'].mean()),
                'Median': _round_to_significant_digits(subset['evap_frac'].median()),
                'Std': _round_to_significant_digits(subset['evap_frac'].std()),
                'Min': _round_to_significant_digits(subset['evap_frac'].min()),
                'Max': _round_to_significant_digits(subset['evap_frac'].max()),
                'N': subset['evap_frac'].count()
            })
        
        # Add Trans_Fraction row if data exists
        if subset['trans_frac'].notna().any():
            rows.append({
                'Group': group_name,
                'Variable': 'Trans_Fraction',
                'Mean': _round_to_significant_digits(subset['trans_frac'].mean()),
                'Median': _round_to_significant_digits(subset['trans_frac'].median()),
                'Std': _round_to_significant_digits(subset['trans_frac'].std()),
                'Min': _round_to_significant_digits(subset['trans_frac'].min()),
                'Max': _round_to_significant_digits(subset['trans_frac'].max()),
                'N': subset['trans_frac'].count()
            })
        
        return rows
    
    # Overall statistics
    overall_rows = create_stats_rows('Overall', df_annual)
    stats_list.extend(overall_rows)
    
    # Statistics by IGBP
    for igbp in df_annual['IGBP'].dropna().unique():
        subset = df_annual[df_annual['IGBP'] == igbp]
        if len(subset) > 0:
            rows = create_stats_rows(igbp, subset)
            stats_list.extend(rows)
    
    # Statistics by Climate
    for climate in df_annual['climate'].dropna().unique():
        subset = df_annual[df_annual['climate'] == climate]
        if len(subset) > 0:
            rows = create_stats_rows(climate, subset)
            stats_list.extend(rows)
    
    # Statistics by Water Class
    for wc in df_annual['water_class'].dropna().unique():
        subset = df_annual[df_annual['water_class'] == wc]
        if len(subset) > 0:
            rows = create_stats_rows(wc, subset)
            stats_list.extend(rows)
    
    # Convert to DataFrame and save
    if stats_list:
        all_stats = pd.DataFrame(stats_list)
        output_file = output_dir / "summary_statistics.csv"
        all_stats.to_csv(output_file, index=False)
        print(PRINT, f"  Saved summary statistics: {output_file.name}")
        print(PRINT, f"  Total rows in summary: {len(all_stats)}")
    else:
        print(PRINT, "  No statistics to save")

# -------------------------------
# Main function
# -------------------------------
def make_partition_boxplots(et_dir=ET_PARTITIONING_DIR, plots_dir=PLOTS_OUTPUT_DIR):
    """Main function to create all plots from ET partitioning data"""
    
    et_dir = Path(et_dir)
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    print(PRINT, "="*60)
    print(PRINT, "Creating ET Partitioning Plots")
    print(PRINT, f"Input directory: {et_dir}")
    print(PRINT, f"Output directory: {plots_dir}")
    print(PRINT, "="*60)
    
    # Read all data
    print(PRINT, "\nStep 1: Reading partition files...")
    long_df = _read_partition_long(et_dir)
    
    if long_df.empty:
        print(PRINT, "No data found. Exiting.")
        return
    
    # Clean metadata in the combined dataframe
    long_df = clean_dataframe_metadata(long_df)
    
    # Aggregate to annual sums
    print(PRINT, "\nStep 2: Aggregating to annual sums...")
    ann = _annual_sums(long_df)
    
    if ann.empty:
        print(PRINT, "No annual data could be aggregated. Exiting.")
        return
    
    # Clean metadata in annual dataframe
    ann = clean_dataframe_metadata(ann)
    
    print(PRINT, f"  Annual data: {len(ann)} site-years from {ann['site'].nunique()} sites")
    print(PRINT, f"  Years range: {ann['year'].min()} to {ann['year'].max()}")
    
    # Print summary of data availability by group (cleaned)
    print(PRINT, "\nData availability by IGBP:")
    for igbp in sorted(ann['IGBP'].dropna().unique()):
        n = len(ann[ann['IGBP'] == igbp])
        n_sites = ann[ann['IGBP'] == igbp]['site'].nunique()
        print(PRINT, f"  {igbp}: {n} site-years from {n_sites} sites")
    
    print(PRINT, "\nData availability by Water Class:")
    for wc in sorted(ann['water_class'].dropna().unique()):
        n = len(ann[ann['water_class'] == wc])
        n_sites = ann[ann['water_class'] == wc]['site'].nunique()
        print(PRINT, f"  {wc}: {n} site-years from {n_sites} sites")
    
    # Create plots
    print(PRINT, "\nStep 3: Creating plots...")
    
    # 1) Boxplots by SITE (each site across years)
    print(PRINT, "  Creating site boxplots...")
    out1 = plots_dir / "evap_trans_by_site_box.png"
    _two_boxplots(ann, group_col="site", out_path=out1,
                  title_prefix="Yearly Fractions by Site", x_wrap=False)
    
    # 2) Boxplots by IGBP (biome)
    print(PRINT, "  Creating IGBP biome boxplots...")
    ann_igbp = ann.dropna(subset=["IGBP"])
    ann_igbp = ann_igbp[ann_igbp['IGBP'] != "Unknown"]
    if not ann_igbp.empty:
        out2 = plots_dir / "evap_trans_by_igbp_box.png"
        _two_boxplots(ann_igbp, group_col="IGBP", out_path=out2,
                      title_prefix="Yearly Fractions by IGBP Biome", x_wrap=False)
    else:
        print(PRINT, "  No IGBP data available for plotting")
    
    # 3) Boxplots by CLIMATE
    print(PRINT, "  Creating climate boxplots...")
    ann_clim = ann.dropna(subset=["climate"])
    ann_clim = ann_clim[ann_clim['climate'] != "Unknown"]
    if not ann_clim.empty:
        out3 = plots_dir / "evap_trans_by_climate_box.png"
        _two_boxplots(ann_clim, group_col="climate", out_path=out3,
                      title_prefix="Yearly Fractions by Climate", x_wrap=True)
    else:
        print(PRINT, "  No climate data available for plotting")
    
    # 4) Boxplots by WATER CLASS (cleaned)
    print(PRINT, "  Creating water class boxplots...")
    ann_wc = ann.dropna(subset=["water_class"])
    ann_wc = ann_wc[ann_wc['water_class'] != "Unknown"]
    if not ann_wc.empty:
        out4 = plots_dir / "evap_trans_by_water_class_box.png"
        _two_boxplots(ann_wc, group_col="water_class", out_path=out4,
                      title_prefix="Yearly Fractions by Water Class", x_wrap=False)
    else:
        print(PRINT, "  No water_class data available for plotting")
    
    # Save summary statistics
    print(PRINT, "\nStep 4: Saving summary statistics...")
    _save_summary_statistics(ann, plots_dir)
    
    # Print final summary
    print(PRINT, "\n" + "="*60)
    print(PRINT, "COMPLETED!")
    print(PRINT, f"Plots saved to: {plots_dir}")
    print(PRINT, "="*60)
    print(PRINT, "\nGenerated plots:")
    for p in plots_dir.glob("*.png"):
        print(PRINT, f"  - {p.name}")
    print(PRINT, f"\nSummary statistics: {plots_dir / 'summary_statistics.csv'}")

# Run directly
if __name__ == "__main__":
    make_partition_boxplots()