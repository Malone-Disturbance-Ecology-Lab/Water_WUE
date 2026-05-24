# -*- coding: utf-8 -*-
"""
ET Partitioning Data Merger - PRODUCTION VERSION (ENCODING FIXED)
Merges half-hourly ET partitioning data (already in growing season) with metadata and drought indices

FIXES APPLIED:
1. Fixed metadata encoding issue (tries multiple encodings)
2. Handles problematic CSV files with error catching
3. Drought merge uses site_name, Year, month for monthly; site_name, Year for yearly
4. safe_divide rejects missing numerator AND denominator
5. WUE metrics require GPP > 0 before calculation
6. Missing column logging per site
7. Metadata deduplication after cleaning
8. water_class preserved in final outputs
9. Year/month consistency check with logging
10. WUE_tra and WUE_eva require Trans_ratio >= 0.10 and Evap_ratio >= 0.10
"""

import os
import glob
import pandas as pd
import numpy as np
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input paths
ET_PARTITIONING_DIR = r"M:\Research\WUE_CUE\ameri_data\ET_partitioning"
METADATA_FILE = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
DROUGHT_DIR = r"M:\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\drought"

# Output paths
OUTPUT_DIR = r"M:\Research\WUE_CUE\data_products"
MONTHLY_OUTPUT = os.path.join(OUTPUT_DIR, "WUE_CUE_monthly_merged_indices.csv")
YEARLY_OUTPUT = os.path.join(OUTPUT_DIR, "WUE_CUE_yearly_merged_indices.csv")

# Expected columns in ET partitioning files
EXPECTED_COLS = [
    'DateTime', 'NEE', 'GPP', 'Reco', 'NEP', 'ET', 'Evap_pen', 'Trans_pen',
    'Tair_f', 'VPD_f', 'Rg_f', 'NETRAD_f', 'WS', 'WD', 'H_f', 'PAR_f', 'RH', 'PA',
    'lai', 'avg_sos', 'avg_eos'
]

# Positive-only variables (negative values set to NaN)
POSITIVE_VARS = ['GPP', 'Reco', 'ET', 'precip_mm', 'Evap_pen', 'Trans_pen', 'lai']

# Site name aliases for known inconsistencies
SITE_ALIASES = {
    'US-STS': 'US-StS',
    'US-StS_FILL': 'US-StS',
    'US-STS_FILL': 'US-StS',
    'US-Cms': 'US-Cms',  # Keep as is
}

# Small threshold for denominator checks
EPS = 0.001

# ============================================================================
# LOGGING CLASS
# ============================================================================

class Logger:
    """Simple logger to collect and print messages"""
    def __init__(self):
        self.messages = []
        self.indent = 0
    
    def add(self, msg, level='INFO'):
        indent_str = "  " * self.indent
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_msg = f"{indent_str}[{timestamp}] {msg}"
        self.messages.append(log_msg)
        print(log_msg)
    
    def indent_increase(self):
        self.indent += 1
    
    def indent_decrease(self):
        self.indent -= 1
    
    def save(self, filepath):
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write("\n".join(self.messages))
    
    def get_summary_stats(self, df, metrics, level='monthly'):
        """Generate summary statistics for key metrics"""
        self.add(f"\n{'='*60}")
        self.add(f"SUMMARY STATISTICS ({level.upper()})")
        self.add(f"{'='*60}")
        
        for metric in metrics:
            if metric in df.columns:
                values = df[metric].dropna()
                if len(values) > 0:
                    self.add(f"\n{metric}:")
                    self.add(f"  Mean: {values.mean():.3f}")
                    self.add(f"  Median: {values.median():.3f}")
                    self.add(f"  Std: {values.std():.3f}")
                    self.add(f"  Min: {values.min():.3f}")
                    self.add(f"  Max: {values.max():.3f}")
                    self.add(f"  N: {len(values)}")
                else:
                    self.add(f"\n{metric}: No valid data")
        
        # Top 5 highest values for QC
        self.add(f"\n{'='*60}")
        self.add(f"TOP 5 HIGHEST VALUES ({level.upper()})")
        self.add(f"{'='*60}")
        
        for metric in ['WUE', 'WUE_tra', 'WUE_eva']:
            if metric in df.columns:
                top5 = df.nlargest(5, metric)[['site_name', 'Year'] + (['month'] if level == 'monthly' else []) + [metric]]
                self.add(f"\n{metric} top 5:")
                for _, row in top5.iterrows():
                    if level == 'monthly':
                        self.add(f"  {row['site_name']} {row['Year']}-{row['month']:02d}: {row[metric]:.2f}")
                    else:
                        self.add(f"  {row['site_name']} {row['Year']}: {row[metric]:.2f}")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def normalize_site_name(name: str) -> str:
    """Apply site name aliases if needed"""
    if pd.isna(name):
        return name
    name = str(name).strip()
    return SITE_ALIASES.get(name, name)


def safe_divide(num, denom, require_positive_gpp=False, gpp=None):
    """
    Safe division with NaN for invalid denominators or missing numerators
    
    Parameters:
    -----------
    num : numerator
    denom : denominator
    require_positive_gpp : bool
    gpp : GPP value (required if require_positive_gpp=True)
    """
    # Check for missing values
    if num is None or pd.isna(num):
        return np.nan
    if denom is None or pd.isna(denom):
        return np.nan
    
    # Check GPP > 0 if required
    if require_positive_gpp:
        if gpp is None or pd.isna(gpp) or gpp <= 0:
            return np.nan
    
    # Check denominator threshold
    if denom <= EPS:
        return np.nan
    
    result = num / denom
    if np.isinf(result):
        return np.nan
    return result


def circular_mean_deg(angles, weights=None):
    """Calculate circular mean for wind direction (degrees)"""
    if len(angles) == 0 or angles.isna().all():
        return np.nan
    
    angles_rad = np.deg2rad(angles)
    if weights is not None and len(weights) == len(angles):
        weights = np.array(weights)
        cos_mean = np.average(np.cos(angles_rad), weights=weights)
        sin_mean = np.average(np.sin(angles_rad), weights=weights)
    else:
        cos_mean = np.mean(np.cos(angles_rad))
        sin_mean = np.mean(np.sin(angles_rad))
    
    mean_rad = np.arctan2(sin_mean, cos_mean)
    mean_deg = np.rad2deg(mean_rad) % 360
    return mean_deg


def clean_positive_vars(df, var_list, logger=None, site_name=None):
    """Set negative values to NaN for positive-only variables and log counts"""
    for var in var_list:
        if var in df.columns:
            n_neg = (df[var] < 0).sum()
            if n_neg > 0 and logger:
                logger.add(f"  {var}: {n_neg} negative values set to NaN")
            df.loc[df[var] < 0, var] = np.nan
    return df


def check_year_month_consistency(df, logger=None, site_name=None):
    """Check consistency between original Year/month columns and parsed DateTime"""
    if 'Year' in df.columns and 'month' in df.columns:
        # Get parsed values
        parsed_year = df['DateTime'].dt.year
        parsed_month = df['DateTime'].dt.month
        
        year_mismatch = (df['Year'] != parsed_year).sum()
        month_mismatch = (df['month'] != parsed_month).sum()
        
        if year_mismatch > 0 and logger:
            logger.add(f"  ⚠️  {year_mismatch} rows have Year mismatch (using parsed values)")
        if month_mismatch > 0 and logger:
            logger.add(f"  ⚠️  {month_mismatch} rows have month mismatch (using parsed values)")
    
    return df


def compute_metrics(df, logger=None, level='monthly'):
    """
    Compute WUE metrics, ratios, and apply bounds
    
    CRITICAL FIX: 
        - WUE requires GPP > 0
        - WUE_tra requires GPP > 0 AND Trans_ratio >= 0.10 (avoids tiny denominator inflation)
        - WUE_eva requires GPP > 0 AND Evap_ratio >= 0.10 (avoids tiny denominator inflation)
        - Existing bounds: WUE_tra > 80 -> NaN, WUE_eva > 80 -> NaN
    """
    
    # ================================================================
    # FIRST: Compute Evap_ratio and Trans_ratio if they don't exist
    # ================================================================
    if 'Evap_ratio' not in df.columns and 'Evap_pen' in df.columns and 'ET' in df.columns:
        df['Evap_ratio'] = df.apply(
            lambda row: safe_divide(row['Evap_pen'], row['ET'], require_positive_gpp=False),
            axis=1
        )
        # Bound between 0 and 1.2
        df.loc[(df['Evap_ratio'] < 0) | (df['Evap_ratio'] > 1.2), 'Evap_ratio'] = np.nan
    
    if 'Trans_ratio' not in df.columns and 'Trans_pen' in df.columns and 'ET' in df.columns:
        df['Trans_ratio'] = df.apply(
            lambda row: safe_divide(row['Trans_pen'], row['ET'], require_positive_gpp=False),
            axis=1
        )
        # Bound between 0 and 1.2
        df.loc[(df['Trans_ratio'] < 0) | (df['Trans_ratio'] > 1.2), 'Trans_ratio'] = np.nan
    
    # ================================================================
    # WUE (standard) - unchanged
    # ================================================================
    if 'GPP' in df.columns and 'ET' in df.columns:
        df['WUE'] = df.apply(
            lambda row: safe_divide(row['GPP'], row['ET'], require_positive_gpp=True, gpp=row['GPP']), 
            axis=1
        )
        # WUE only needs to be > 0 (no upper limit)
        df.loc[df['WUE'] <= 0, 'WUE'] = np.nan
    
    # ================================================================
    # WUE_tra (transpiration-based) - requires Trans_ratio >= 0.10
    # ================================================================
    if 'GPP' in df.columns and 'Trans_pen' in df.columns and 'Trans_ratio' in df.columns:
        # Count excluded records for logging
        if logger is not None:
            low_ratio_mask = (df['Trans_ratio'] < 0.10) & (df['Trans_ratio'].notna())
            n_low_ratio = low_ratio_mask.sum()
            if n_low_ratio > 0:
                logger.add(f"  {level}: {n_low_ratio} records excluded from WUE_tra (Trans_ratio < 0.10)")
        
        df['WUE_tra'] = df.apply(
            lambda row: safe_divide(row['GPP'], row['Trans_pen'], require_positive_gpp=True, gpp=row['GPP'])
            if pd.notna(row.get('Trans_ratio')) and row['Trans_ratio'] >= 0.10
            else np.nan,
            axis=1
        )
        # Bound: >0 and <= 80
        df.loc[(df['WUE_tra'] <= 0) | (df['WUE_tra'] > 80), 'WUE_tra'] = np.nan
    
    # ================================================================
    # WUE_eva (evaporation-based) - requires Evap_ratio >= 0.10
    # ================================================================
    if 'GPP' in df.columns and 'Evap_pen' in df.columns and 'Evap_ratio' in df.columns:
        # Count excluded records for logging
        if logger is not None:
            low_ratio_mask = (df['Evap_ratio'] < 0.10) & (df['Evap_ratio'].notna())
            n_low_ratio = low_ratio_mask.sum()
            if n_low_ratio > 0:
                logger.add(f"  {level}: {n_low_ratio} records excluded from WUE_eva (Evap_ratio < 0.10)")
        
        df['WUE_eva'] = df.apply(
            lambda row: safe_divide(row['GPP'], row['Evap_pen'], require_positive_gpp=True, gpp=row['GPP'])
            if pd.notna(row.get('Evap_ratio')) and row['Evap_ratio'] >= 0.10
            else np.nan,
            axis=1
        )
        # Bound: >0 and <= 80
        df.loc[(df['WUE_eva'] <= 0) | (df['WUE_eva'] > 80), 'WUE_eva'] = np.nan
    
    # ================================================================
    # CUE (Carbon Use Efficiency) - allow negative NEP (unchanged)
    # ================================================================
    if 'NEP' in df.columns and 'GPP' in df.columns:
        df['CUE'] = df.apply(
            lambda row: safe_divide(row['NEP'], row['GPP'], require_positive_gpp=False),
            axis=1
        )
        # Bound: between -10 and 1.5
        df.loc[(df['CUE'] < -10) | (df['CUE'] > 1.5), 'CUE'] = np.nan
    
    # ================================================================
    # Partition consistency check (warning only) - unchanged
    # ================================================================
    if all(col in df.columns for col in ['Evap_pen', 'Trans_pen', 'ET']):
        partition_sum = (df['Evap_pen'] + df['Trans_pen']) / df['ET']
        n_inconsistent = ((partition_sum < 0.5) | (partition_sum > 1.5)).sum()
        if n_inconsistent > 0 and logger:
            logger.add(f"  ⚠️  {n_inconsistent} {level} records have Evap_pen+Trans_pen inconsistent with ET")
    
    return df


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def load_metadata(metadata_file: str, logger: Logger) -> Optional[pd.DataFrame]:
    """Load and clean site metadata, deduplicate by site_name with encoding handling"""
    
    # Try multiple encodings
    encodings = ['utf-8', 'latin1', 'cp1252', 'ISO-8859-1', 'utf-16']
    
    meta = None
    for encoding in encodings:
        try:
            logger.add(f"  Trying encoding: {encoding}")
            meta = pd.read_csv(metadata_file, encoding=encoding)
            logger.add(f"  ✅ Successfully read metadata with {encoding} encoding")
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
        except Exception as e:
            logger.add(f"  Error with {encoding}: {e}")
            continue
    
    if meta is None:
        try:
            # Last resort: let pandas guess
            logger.add(f"  Trying automatic encoding detection...")
            meta = pd.read_csv(metadata_file, encoding_errors='ignore')
            logger.add(f"  ✅ Read with automatic encoding")
        except Exception as e:
            logger.add(f"❌ Error loading metadata: {e}")
            return None
    
    meta.columns = meta.columns.str.strip()
    
    # Keep all columns but we'll use specific ones
    keep_cols = ['site_name', 'lat', 'long', 'IGBP', 'climate', 'water_class']
    available_cols = [col for col in keep_cols if col in meta.columns]
    meta = meta[available_cols].copy()
    
    # Clean string columns
    for col in ['site_name', 'IGBP', 'climate', 'water_class']:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).str.strip()
    
    # Apply site name normalization
    meta['site_name'] = meta['site_name'].apply(normalize_site_name)
    
    # DEDUPLICATE: Keep first occurrence of each site_name
    before_dedup = len(meta)
    meta = meta.drop_duplicates(subset=['site_name'], keep='first')
    after_dedup = len(meta)
    
    if before_dedup > after_dedup:
        logger.add(f"  Deduplicated metadata: {before_dedup} → {after_dedup} rows")
    
    logger.add(f"✅ Loaded metadata: {len(meta)} unique sites")
    logger.add(f"   Sites: {sorted(meta['site_name'].unique())}")
    return meta


def load_drought_site(site_name: str, drought_dir: str, logger: Logger) -> Optional[pd.DataFrame]:
    """Load and process drought data for a single site to monthly values"""
    drought_file = os.path.join(drought_dir, f"{site_name}.csv")
    
    if not os.path.exists(drought_file):
        logger.add(f"  ⚠️  No drought file found")
        return None
    
    try:
        drought_df = pd.read_csv(drought_file)
        drought_df.columns = drought_df.columns.str.strip()
        
        # Parse DateTime
        if 'DateTime' not in drought_df.columns:
            logger.add(f"  ⚠️  No DateTime column")
            return None
        
        drought_df['DateTime'] = pd.to_datetime(drought_df['DateTime'], errors='coerce')
        drought_df = drought_df.dropna(subset=['DateTime'])
        
        if len(drought_df) == 0:
            logger.add(f"  ⚠️  No valid dates")
            return None
        
        drought_df['Year'] = drought_df['DateTime'].dt.year
        drought_df['month'] = drought_df['DateTime'].dt.month
        
        # SPEI columns
        spei_cols = [f'SPEI_{s}' for s in [1, 3, 6, 12, 24, 36, 48]]
        existing_spei = [col for col in spei_cols if col in drought_df.columns]
        
        if not existing_spei:
            logger.add(f"  ⚠️  No SPEI columns found")
            return None
        
        # Aggregate to monthly (take first non-null value per month)
        monthly_drought = drought_df.groupby(['Year', 'month'])[existing_spei].first().reset_index()
        monthly_drought['site_name'] = site_name
        
        logger.add(f"  ✅ Loaded drought: {len(monthly_drought)} months")
        return monthly_drought
        
    except Exception as e:
        logger.add(f"  ❌ Error loading drought: {e}")
        return None


def load_et_partitioning_file(filepath: str, logger: Logger) -> Optional[Tuple[pd.DataFrame, str]]:
    """Load and validate a single ET partitioning file with error handling"""
    site_name = Path(filepath).stem
    original_name = site_name
    site_name = normalize_site_name(site_name)
    
    if original_name != site_name:
        logger.add(f"  🔄 Site name normalized: {original_name} → {site_name}")
    
    try:
        # Try different engines for problematic files
        try:
            df = pd.read_csv(filepath)
        except Exception as e:
            logger.add(f"  ⚠️  Standard read failed: {e}, trying python engine")
            df = pd.read_csv(filepath, engine='python', on_bad_lines='skip')
        
        df.columns = df.columns.str.strip()
        
        # Check for missing expected columns
        missing_cols = [col for col in EXPECTED_COLS if col not in df.columns]
        if missing_cols:
            logger.add(f"  ⚠️  Missing columns: {missing_cols}")
            # Create missing columns as NaN
            for col in missing_cols:
                df[col] = np.nan
        
        # Parse DateTime
        if 'DateTime' not in df.columns:
            logger.add(f"  ❌ No DateTime column")
            return None
        
        df['DateTime'] = pd.to_datetime(df['DateTime'], errors='coerce')
        df = df.dropna(subset=['DateTime'])
        
        if len(df) == 0:
            logger.add(f"  ❌ No valid dates")
            return None
        
        # Extract time components
        df['Year'] = df['DateTime'].dt.year
        df['month'] = df['DateTime'].dt.month
        df['day'] = df['DateTime'].dt.day
        df['DoY'] = df['DateTime'].dt.dayofyear
        df['Hour'] = df['DateTime'].dt.hour + df['DateTime'].dt.minute / 60
        
        # Rename precip column if needed
        if 'precip_mm_per_30min' in df.columns and 'precip_mm' not in df.columns:
            df = df.rename(columns={'precip_mm_per_30min': 'precip_mm'})
            logger.add(f"  Renamed precip_mm_per_30min → precip_mm")
        
        # Check Year/month consistency
        df = check_year_month_consistency(df, logger, site_name)
        
        # Clean positive-only variables and log counts
        df = clean_positive_vars(df, POSITIVE_VARS, logger, site_name)
        
        logger.add(f"  ✅ Loaded: {len(df)} half-hourly records")
        return df, site_name
        
    except Exception as e:
        logger.add(f"  ❌ Error loading: {e}")
        import traceback
        logger.add(f"  Traceback: {traceback.format_exc()}")
        return None


# ============================================================================
# AGGREGATION FUNCTIONS
# ============================================================================

def aggregate_monthly(site_df: pd.DataFrame, site_name: str, metadata: pd.DataFrame, 
                      drought_monthly: pd.DataFrame, logger: Logger) -> pd.DataFrame:
    """Aggregate half-hourly data to monthly"""
    
    # Define aggregation rules
    sum_cols = ['NEE', 'GPP', 'Reco', 'NEP', 'ET', 'precip_mm', 'Evap_pen', 'Trans_pen']
    mean_cols = ['Tair_f', 'VPD_f', 'Rg_f', 'NETRAD_f', 'WS', 'H_f', 'PAR_f', 'RH', 'PA', 'lai']
    
    # Check which columns exist
    existing_sum = [col for col in sum_cols if col in site_df.columns]
    existing_mean = [col for col in mean_cols if col in site_df.columns]
    
    if not existing_sum and not existing_mean:
        logger.add(f"  ⚠️  No columns to aggregate")
        return pd.DataFrame()
    
    # Create aggregation dictionary
    agg_dict = {col: 'sum' for col in existing_sum}
    agg_dict.update({col: 'mean' for col in existing_mean})
    
    # Add avg_sos, avg_eos as constants (take first value)
    if 'avg_sos' in site_df.columns:
        agg_dict['avg_sos'] = 'first'
    if 'avg_eos' in site_df.columns:
        agg_dict['avg_eos'] = 'first'
    
    # Aggregate by year and month
    monthly = site_df.groupby(['Year', 'month']).agg(agg_dict).reset_index()
    
    if len(monthly) == 0:
        logger.add(f"  ⚠️  No monthly aggregates created")
        return pd.DataFrame()
    
    # Add site_name
    monthly['site_name'] = site_name
    
    # Calculate avg_length
    if 'avg_sos' in monthly.columns and 'avg_eos' in monthly.columns:
        monthly['avg_length'] = monthly['avg_eos'] - monthly['avg_sos'] + 1
    else:
        monthly['avg_length'] = np.nan
    
    # Add metadata
    if metadata is not None:
        site_meta = metadata[metadata['site_name'] == site_name]
        if not site_meta.empty:
            for col in ['lat', 'long', 'IGBP', 'climate', 'water_class']:
                if col in site_meta.columns:
                    monthly[col] = site_meta.iloc[0][col]
        else:
            logger.add(f"  ⚠️  No metadata found for {site_name}")
    
    # Add drought data - merge on site_name, Year, month
    if drought_monthly is not None and len(drought_monthly) > 0:
        monthly = monthly.merge(
            drought_monthly[['site_name', 'Year', 'month'] + [f'SPEI_{s}' for s in [1, 3, 6, 12, 24, 36, 48]]],
            on=['site_name', 'Year', 'month'],
            how='left'
        )
    else:
        # Add empty SPEI columns
        for s in [1, 3, 6, 12, 24, 36, 48]:
            monthly[f'SPEI_{s}'] = np.nan
    
    # Handle wind direction (circular mean)
    if 'WD' in site_df.columns:
        wd_monthly = site_df.groupby(['Year', 'month']).apply(
            lambda x: circular_mean_deg(x['WD'])
        ).reset_index(name='WD')
        wd_monthly.columns = ['Year', 'month', 'WD']
        monthly = monthly.merge(wd_monthly, on=['Year', 'month'], how='left')
    else:
        monthly['WD'] = np.nan
    
    # Compute metrics
    monthly = compute_metrics(monthly, logger, level='monthly')
    
    logger.add(f"  ✅ Monthly: {len(monthly)} records")
    return monthly


def aggregate_yearly(site_df: pd.DataFrame, site_name: str, metadata: pd.DataFrame, 
                     drought_monthly: pd.DataFrame, logger: Logger) -> pd.DataFrame:
    """Aggregate half-hourly data to yearly"""
    
    # Define aggregation rules
    sum_cols = ['NEE', 'GPP', 'Reco', 'NEP', 'ET', 'precip_mm', 'Evap_pen', 'Trans_pen']
    mean_cols = ['Tair_f', 'VPD_f', 'Rg_f', 'NETRAD_f', 'WS', 'H_f', 'PAR_f', 'RH', 'PA', 'lai']
    
    # Check which columns exist
    existing_sum = [col for col in sum_cols if col in site_df.columns]
    existing_mean = [col for col in mean_cols if col in site_df.columns]
    
    if not existing_sum and not existing_mean:
        logger.add(f"  ⚠️  No columns to aggregate")
        return pd.DataFrame()
    
    # Create aggregation dictionary
    agg_dict = {col: 'sum' for col in existing_sum}
    agg_dict.update({col: 'mean' for col in existing_mean})
    
    # Add avg_sos, avg_eos as constants (take first value)
    if 'avg_sos' in site_df.columns:
        agg_dict['avg_sos'] = 'first'
    if 'avg_eos' in site_df.columns:
        agg_dict['avg_eos'] = 'first'
    
    # Aggregate by year
    yearly = site_df.groupby(['Year']).agg(agg_dict).reset_index()
    
    if len(yearly) == 0:
        logger.add(f"  ⚠️  No yearly aggregates created")
        return pd.DataFrame()
    
    # Add site_name
    yearly['site_name'] = site_name
    
    # Calculate avg_length
    if 'avg_sos' in yearly.columns and 'avg_eos' in yearly.columns:
        yearly['avg_length'] = yearly['avg_eos'] - yearly['avg_sos'] + 1
    else:
        yearly['avg_length'] = np.nan
    
    # Add metadata
    if metadata is not None:
        site_meta = metadata[metadata['site_name'] == site_name]
        if not site_meta.empty:
            for col in ['lat', 'long', 'IGBP', 'climate', 'water_class']:
                if col in site_meta.columns:
                    yearly[col] = site_meta.iloc[0][col]
        else:
            logger.add(f"  ⚠️  No metadata found for {site_name}")
    
    # Aggregate drought data to yearly - merge on site_name and Year
    if drought_monthly is not None and len(drought_monthly) > 0:
        spei_cols = [f'SPEI_{s}' for s in [1, 3, 6, 12, 24, 36, 48]]
        existing_spei = [col for col in spei_cols if col in drought_monthly.columns]
        
        if existing_spei:
            yearly_drought = drought_monthly.groupby(['site_name', 'Year'])[existing_spei].mean().reset_index()
            yearly = yearly.merge(yearly_drought, on=['site_name', 'Year'], how='left')
        else:
            for s in [1, 3, 6, 12, 24, 36, 48]:
                yearly[f'SPEI_{s}'] = np.nan
    else:
        for s in [1, 3, 6, 12, 24, 36, 48]:
            yearly[f'SPEI_{s}'] = np.nan
    
    # Handle wind direction (circular mean)
    if 'WD' in site_df.columns:
        wd_yearly = site_df.groupby(['Year']).apply(
            lambda x: circular_mean_deg(x['WD'])
        ).reset_index(name='WD')
        wd_yearly.columns = ['Year', 'WD']
        yearly = yearly.merge(wd_yearly, on=['Year'], how='left')
    else:
        yearly['WD'] = np.nan
    
    # Compute metrics
    yearly = compute_metrics(yearly, logger, level='yearly')
    
    logger.add(f"  ✅ Yearly: {len(yearly)} records")
    return yearly


# ============================================================================
# PROCESSING FUNCTION
# ============================================================================

def process_site(site_file: str, metadata: pd.DataFrame, logger: Logger) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Process a single ET partitioning site file"""
    
    # Load ET partitioning data
    result = load_et_partitioning_file(site_file, logger)
    if result is None:
        return None, None
    
    site_df, site_name = result
    
    # Load drought data
    drought_monthly = load_drought_site(site_name, DROUGHT_DIR, logger)
    
    # Aggregate
    monthly = aggregate_monthly(site_df, site_name, metadata, drought_monthly, logger)
    yearly = aggregate_yearly(site_df, site_name, metadata, drought_monthly, logger)
    
    return monthly, yearly


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main processing function"""
    
    print("=" * 80)
    print("ET PARTITIONING DATA MERGER - PRODUCTION VERSION")
    print("=" * 80)
    
    # Initialize logger
    logger = Logger()
    logger.add(f"Started: {datetime.now()}")
    logger.add(f"ET Partitioning Dir: {ET_PARTITIONING_DIR}")
    logger.add(f"Output Dir: {OUTPUT_DIR}")
    
    # Find all ET partitioning files
    site_files = glob.glob(os.path.join(ET_PARTITIONING_DIR, "*.csv"))
    logger.add(f"\n📁 Found {len(site_files)} site files")
    
    if not site_files:
        logger.add("❌ No site files found!")
        return
    
    # Load metadata
    logger.add("\n📋 Loading metadata...")
    metadata = load_metadata(METADATA_FILE, logger)
    
    # Process each site
    all_monthly = []
    all_yearly = []
    successful_sites = 0
    failed_sites = []
    
    for i, site_file in enumerate(site_files, 1):
        site_name = Path(site_file).stem
        logger.add(f"\n[{i}/{len(site_files)}] Processing {site_name}...")
        logger.indent_increase()
        
        monthly, yearly = process_site(site_file, metadata, logger)
        
        if monthly is not None and len(monthly) > 0:
            all_monthly.append(monthly)
            successful_sites += 1
        else:
            failed_sites.append(site_name)
        
        if yearly is not None and len(yearly) > 0:
            all_yearly.append(yearly)
        
        logger.indent_decrease()
    
    # Combine and save results
    logger.add(f"\n{'='*60}")
    logger.add("SAVING RESULTS")
    logger.add(f"{'='*60}")
    
    # Define final column order (including water_class)
    final_cols = [
        'site_name', 'Year', 'month', 'NEE', 'GPP', 'Reco', 'ET', 'NEP', 'IGBP', 'water_class',
        'WUE', 'WUE_tra', 'WUE_eva', 'CUE', 'Tair_f', 'VPD_f', 'Rg_f',
        'avg_length', 'avg_sos', 'avg_eos', 'lat', 'long', 'climate',
        'NETRAD_f', 'WS', 'WD', 'H_f', 'PAR_f', 'RH', 'PA', 'precip_mm',
        'Evap_pen', 'Trans_pen', 'Evap_ratio', 'Trans_ratio', 'lai',
        'SPEI_1', 'SPEI_3', 'SPEI_6', 'SPEI_12', 'SPEI_24', 'SPEI_36', 'SPEI_48'
    ]
    
    # Save monthly
    if all_monthly:
        final_monthly = pd.concat(all_monthly, ignore_index=True)
        
        # Reorder columns
        existing_cols = [col for col in final_cols if col in final_monthly.columns]
        # Add 'month' only for monthly output
        if 'month' in final_monthly.columns and 'month' not in existing_cols:
            existing_cols.insert(2, 'month')
        final_monthly = final_monthly[existing_cols]
        
        # Save
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        final_monthly.to_csv(MONTHLY_OUTPUT, index=False)
        logger.add(f"\n💾 Saved monthly: {MONTHLY_OUTPUT}")
        logger.add(f"   Rows: {len(final_monthly)}")
        logger.add(f"   Columns: {len(final_monthly.columns)}")
        
        # Summary statistics
        monthly_metrics = ['WUE', 'WUE_tra', 'WUE_eva', 'CUE', 'Evap_ratio', 'Trans_ratio']
        logger.get_summary_stats(final_monthly, monthly_metrics, level='monthly')
        
    else:
        logger.add("\n❌ No monthly data generated!")
    
    # Save yearly (without month column)
    if all_yearly:
        final_yearly = pd.concat(all_yearly, ignore_index=True)
        
        # Remove month column for yearly output
        yearly_cols = [col for col in final_cols if col != 'month' and col in final_yearly.columns]
        final_yearly = final_yearly[yearly_cols]
        
        # Save
        final_yearly.to_csv(YEARLY_OUTPUT, index=False)
        logger.add(f"\n💾 Saved yearly: {YEARLY_OUTPUT}")
        logger.add(f"   Rows: {len(final_yearly)}")
        logger.add(f"   Columns: {len(final_yearly.columns)}")
        
        # Summary statistics
        yearly_metrics = ['WUE', 'WUE_tra', 'WUE_eva', 'CUE', 'Evap_ratio', 'Trans_ratio']
        logger.get_summary_stats(final_yearly, yearly_metrics, level='yearly')
        
    else:
        logger.add("\n❌ No yearly data generated!")
    
    # Final summary
    logger.add(f"\n{'='*60}")
    logger.add("PROCESSING SUMMARY")
    logger.add(f"{'='*60}")
    logger.add(f"Total sites found: {len(site_files)}")
    logger.add(f"Successfully processed: {successful_sites}")
    if failed_sites:
        logger.add(f"Failed sites: {failed_sites}")
    
    # Save log file
    log_file = os.path.join(OUTPUT_DIR, "merge_log.txt")
    logger.save(log_file)
    logger.add(f"\n✅ Log saved to: {log_file}")
    
    return


# ============================================================================
# RUN
# ============================================================================

if __name__ == "__main__":
    main()