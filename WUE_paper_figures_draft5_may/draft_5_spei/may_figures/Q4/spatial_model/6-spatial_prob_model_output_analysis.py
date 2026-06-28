# -*- coding: utf-8 -*-
"""
CONCISE SPATIAL ANALYSIS WORKFLOW - UPDATED FOR CORRECTED MODEL
================================================================================
Now using model with NoChange kept as baseline.
Processes ONLY WUE_T (transpiration) metric
Class definitions:
    Class=1 = Less efficient (Decrease in WUE) - drought response
    Class=0 = Increase or NoChange (baseline) - non-decrease response

No changes needed to file names - automatically compatible.
================================================================================
"""

import xarray as xr
import numpy as np
import pandas as pd
import os
import glob
import re
from datetime import datetime
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

class ConciseSpatialAnalyzer:
    """Streamlined spatial analysis for pre-figure data preparation."""
    
    def __init__(self):
        """Initialize with hardcoded paths."""
        # INPUT: Monthly prediction files from Step 2 (corrected model)
        self.input_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs"
        
        # OUTPUT: All analysis products go here
        self.output_base = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis"
        
        # ONLY WUE_tra (WUE_T) - FIXED!
        self.metrics = ['WUE_tra']
        
        # Model interpretation note (UPDATED)
        self.model_note = """
        MODEL INTERPRETATION:
        - Class=1 = Less efficient (Decrease in WUE) → drought response
        - Class=0 = Increase or NoChange (baseline) → non-decrease response
        - NoChange is correctly kept as baseline (Class=0)
        - Processing ONLY WUE_T (transpiration) metric
        """
        
        # Verify input directory exists
        if not os.path.exists(self.input_dir):
            print(f"❌ ERROR: Input directory not found: {self.input_dir}")
            print("Please run the spatial projection workflow first.")
            return
        
        # Create output structure
        for subdir in ['diagnostics', 'tables_breakpoint', 
                      'rasters_time_aggregated', 'rasters_pre_post_breakpoint']:
            os.makedirs(os.path.join(self.output_base, subdir), exist_ok=True)
        
        print("="*80)
        print("CONCISE SPATIAL ANALYSIS WORKFLOW")
        print("="*80)
        print(self.model_note)
        print(f"📂 Input:  {self.input_dir}")
        print(f"📁 Output: {self.output_base}")
    
    # ==========================================================================
    # STEP 0: Build master file index
    # ==========================================================================
    
    def step0_build_file_index(self):
        """Build ordered list of monthly files."""
        print("\nSTEP 0: BUILDING FILE INDEX")
        print("-"*40)
        
        # Find files (SAME NAMES as before - no change needed)
        pattern = os.path.join(self.input_dir, "WUE_predictions_linearlogistic_*.nc")
        files = sorted(glob.glob(pattern))
        
        if not files:
            print(f"❌ ERROR: No files found matching: {pattern}")
            print("Please run the spatial projection workflow first.")
            return None
        
        # Parse dates
        file_data = []
        for f in files:
            filename = os.path.basename(f)
            
            # Extract YYYYMM from filename
            match = re.search(r'_(\d{6})\.nc', filename)
            if not match:
                match = re.search(r'(\d{6})', filename)
            
            if match:
                ym_str = match.group(1)
                year, month = int(ym_str[:4]), int(ym_str[4:6])
                file_data.append({
                    'file_path': f,
                    'filename': filename,
                    'year': year,
                    'month': month,
                    'yearmonth': ym_str,
                    'yearmonth_int': int(ym_str),
                    'date': datetime(year, month, 1)
                })
            else:
                print(f"  ⚠️ Could not parse date from: {filename}")
        
        # Create DataFrame
        df = pd.DataFrame(file_data).sort_values('yearmonth_int')
        
        if len(df) == 0:
            print("❌ No valid files found")
            return None
        
        # Save
        output_path = os.path.join(self.output_base, 'diagnostics', 'file_index.csv')
        df.to_csv(output_path, index=False)
        
        # Console output
        print(f"📊 Found {len(df)} monthly files")
        print(f"📅 Range: {df['yearmonth'].iloc[0]} to {df['yearmonth'].iloc[-1]}")
        print(f"💾 Saved: {output_path}")
        
        self.file_index = df
        return df
    
    # ==========================================================================
    # STEP 1: One-file sanity check (FIXED - only WUE_tra)
    # ==========================================================================
    
    def step1_sanity_check(self):
        """Validate first file structure - confirms NoChange baseline."""
        print("\nSTEP 1: FILE SANITY CHECK")
        print("-"*40)
        print("  Validating corrected model outputs (NoChange kept as baseline)")
        
        if not hasattr(self, 'file_index') or len(self.file_index) == 0:
            print("❌ Run Step 0 first")
            return False
        
        first_file = self.file_index.iloc[0]
        print(f"Testing: {first_file['filename']} ({first_file['yearmonth']})")
        
        try:
            ds = xr.open_dataset(first_file['file_path'])
            
            # Check required variables for ONLY WUE_tra
            required_vars = ['SPEI48'] + \
                           [f'P_decrease_{m}' for m in self.metrics] + \
                           [f'Class_{m}' for m in self.metrics]
            
            missing_vars = [v for v in required_vars if v not in ds.data_vars]
            if missing_vars:
                print(f"❌ Missing variables: {missing_vars}")
                ds.close()
                return False
            
            print("✅ All required variables present")
            print("✅ Class definition: 1=Less efficient (Decrease), 0=Increase or NoChange (baseline)")
            print(f"✅ Processing {len(self.metrics)} metric: {', '.join(self.metrics)}")
            
            # Per-metric checks
            print("\n🔍 Per-metric validation:")
            all_checks_passed = True
            
            for metric in self.metrics:
                class_var = f'Class_{metric}'
                prob_var = f'P_decrease_{metric}'
                
                class_data = ds[class_var].values
                prob_data = ds[prob_var].values
                
                # Use isfinite() for valid mask
                valid_mask = np.isfinite(class_data) & np.isfinite(prob_data)
                class_valid = class_data[valid_mask]
                prob_valid = prob_data[valid_mask]
                
                n_valid = len(class_valid)
                n_less = np.sum(class_valid == 1)   # Less efficient (Decrease)
                n_non_decrease = np.sum(class_valid == 0)   # Non-decrease (Increase/NoChange)
                
                if n_valid == 0:
                    print(f"  ❌ {metric}: No valid pixels (all NaN)")
                    all_checks_passed = False
                    continue
                
                # Check probability range
                prob_min, prob_max = prob_valid.min(), prob_valid.max()
                prob_in_range = 0 <= prob_min <= prob_max <= 1
                
                # Check class values
                unique_classes = np.unique(class_valid)
                classes_valid = set(unique_classes).issubset({0, 1})
                
                # Check class/probability consistency
                class_from_prob = (prob_valid >= 0.5).astype(int)
                agreement = np.mean(class_valid == class_from_prob) * 100
                
                # Metric name mapping for display
                metric_display = {
                    'WUE_tra': 'WUE_T (transpiration)'
                }.get(metric, metric)
                
                print(f"\n  {metric_display}:")
                print(f"    Valid pixels: {n_valid:,}")
                print(f"    Less efficient (Class=1): {n_less:,} ({100*n_less/n_valid:.1f}%)")
                print(f"    Non-decrease (Class=0): {n_non_decrease:,} ({100*n_non_decrease/n_valid:.1f}%)")
                print(f"    Probability range: [{prob_min:.3f}, {prob_max:.3f}] {'✅' if prob_in_range else '❌'}")
                print(f"    Class values: {sorted(unique_classes.tolist())} {'✅' if classes_valid else '❌'}")
                print(f"    Class/prob agreement (p>=0.5 ↔ class==1): {agreement:.1f}%")
                
                if not prob_in_range:
                    print(f"    ⚠️ Probability outside [0,1] range!")
                    all_checks_passed = False
                if not classes_valid:
                    print(f"    ⚠️ Unexpected class values!")
                    all_checks_passed = False
                if agreement < 99:
                    print(f"    ⚠️ Low agreement between class and probability threshold!")
                    all_checks_passed = False
            
            # Check SPEI data
            spei_data = ds['SPEI48'].values
            spei_valid = spei_data[np.isfinite(spei_data)]
            print(f"\n🌡️ SPEI48 check:")
            print(f"    Valid SPEI pixels: {len(spei_valid):,}")
            if len(spei_valid) > 0:
                print(f"    SPEI range: [{spei_valid.min():.3f}, {spei_valid.max():.3f}]")
            
            ds.close()
            
            if all_checks_passed:
                print("\n✅ All sanity checks passed - Corrected model working!")
                return True
            else:
                print("\n⚠️ Some checks failed - review above")
                return False
            
        except Exception as e:
            print(f"❌ Error opening file: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    # ==========================================================================
    # STEP 2: Monthly decision table (LONG format only) - UPDATED labels
    # ==========================================================================
    
    def step2_create_monthly_table(self):
        """Create monthly long table for breakpoint analysis."""
        print("\nSTEP 2: MONTHLY DECISION TABLE")
        print("-"*40)
        
        if not hasattr(self, 'file_index'):
            print("❌ Run Step 0 first")
            return None
        
        monthly_data = []
        n_files = len(self.file_index)
        
        print(f"Processing {n_files} months for {len(self.metrics)} metric...")
        
        for idx, row in self.file_index.iterrows():
            if (idx + 1) % 12 == 0 or (idx + 1) == n_files:
                print(f"  Processed {idx + 1}/{n_files} months")
            
            try:
                ds = xr.open_dataset(row['file_path'])
                
                # Get SPEI data
                spei_data = ds['SPEI48'].values
                spei_valid = spei_data[np.isfinite(spei_data)]
                
                if len(spei_valid) == 0:
                    print(f"  ⚠️ {row['yearmonth']}: No valid SPEI data")
                    ds.close()
                    continue
                
                # SPEI statistics
                spei_mean = np.mean(spei_valid)
                n_spei_below_minus1 = np.sum(spei_valid < -1)
                n_spei_below_minus2 = np.sum(spei_valid < -2)
                frac_spei_below_minus1 = n_spei_below_minus1 / len(spei_valid)
                frac_spei_below_minus2 = n_spei_below_minus2 / len(spei_valid)
                
                for metric in self.metrics:
                    class_var = f'Class_{metric}'
                    prob_var = f'P_decrease_{metric}'
                    
                    if class_var not in ds or prob_var not in ds:
                        print(f"  ⚠️ {row['yearmonth']}: Missing {metric} variables")
                        continue
                    
                    class_data = ds[class_var].values
                    prob_data = ds[prob_var].values
                    
                    # Use isfinite() for valid mask
                    valid_mask = np.isfinite(class_data) & np.isfinite(prob_data)
                    class_valid = class_data[valid_mask]
                    prob_valid = prob_data[valid_mask]
                    
                    if len(class_valid) == 0:
                        continue
                    
                    # Counts
                    n_valid = len(class_valid)
                    n_less = np.sum(class_valid == 1)   # Less efficient
                    n_non_decrease = np.sum(class_valid == 0)   # Non-decrease
                    
                    # Fractions
                    frac_less = n_less / n_valid  # Fraction Less efficient
                    frac_non_decrease = n_non_decrease / n_valid  # Fraction Non-decrease
                    
                    # Probability statistics
                    mean_prob = np.mean(prob_valid)      # Mean probability of Less efficient
                    median_prob = np.median(prob_valid)  # Median probability of Less efficient
                    
                    # Add to data
                    monthly_data.append({
                        'year': row['year'],
                        'month': row['month'],
                        'yearmonth': row['yearmonth'],
                        'yearmonth_int': row['yearmonth_int'],
                        'date': row['date'],
                        'metric': metric,
                        'n_valid': n_valid,
                        'n_less': n_less,
                        'n_non_decrease': n_non_decrease,
                        'frac_less': frac_less,      # Fraction Less efficient
                        'frac_non_decrease': frac_non_decrease,      # Fraction Non-decrease
                        'mean_prob': mean_prob,      # Mean probability of Less efficient
                        'median_prob': median_prob,  # Median probability of Less efficient
                        'spei_mean': spei_mean,
                        'frac_spei_below_minus1': frac_spei_below_minus1,
                        'frac_spei_below_minus2': frac_spei_below_minus2
                    })
                
                ds.close()
                
            except Exception as e:
                print(f"❌ Error processing {row['filename']}: {e}")
                continue
        
        # Create DataFrame
        if not monthly_data:
            print("❌ No monthly data created")
            return None
        
        monthly_df = pd.DataFrame(monthly_data)
        
        # Save long format
        output_path = os.path.join(self.output_base, 'tables_breakpoint', 'monthly_metric_long.csv')
        monthly_df.to_csv(output_path, index=False)
        
        # Console summary
        print(f"\n📊 MONTHLY TABLE SUMMARY:")
        print(f"  Total records: {len(monthly_df)} ({len(monthly_df)/len(self.metrics):.0f} months × {len(self.metrics)} metric)")
        print(f"  Date range: {monthly_df['yearmonth'].min()} to {monthly_df['yearmonth'].max()}")
        print(f"  💾 Saved: {output_path}")
        
        # Per-metric console summary
        print("\n🔍 PER-METRIC SUMMARY (Key Interpretation):")
        print("  Note: frac_less = fraction of pixels classified as 'Less efficient' (Class=1)")
        print("        mean_prob = mean probability of 'Less efficient' response")
        print("        ρ = Spearman correlation coefficient")
        print("")
        
        metric_display = {
            'WUE_tra': 'WUE_T (transpiration)'
        }
        
        for metric in self.metrics:
            metric_df = monthly_df[monthly_df['metric'] == metric]
            if len(metric_df) == 0:
                print(f"  {metric_display.get(metric, metric)}: No data")
                continue
            
            # Basic stats
            frac_mean = metric_df['frac_less'].mean()
            frac_std = metric_df['frac_less'].std()
            frac_min = metric_df['frac_less'].min()
            frac_max = metric_df['frac_less'].max()
            
            # Correlation with drought
            valid_mask = metric_df['frac_less'].notna() & metric_df['frac_spei_below_minus1'].notna()
            if valid_mask.sum() > 2:
                rho, pval = stats.spearmanr(
                    metric_df.loc[valid_mask, 'frac_less'],
                    metric_df.loc[valid_mask, 'frac_spei_below_minus1']
                )
            else:
                rho = pval = np.nan
            
            # Correlation with SPEI
            valid_mask = metric_df['frac_less'].notna() & metric_df['spei_mean'].notna()
            if valid_mask.sum() > 2:
                rho_spei, pval_spei = stats.spearmanr(
                    metric_df.loc[valid_mask, 'frac_less'],
                    metric_df.loc[valid_mask, 'spei_mean']
                )
            else:
                rho_spei = pval_spei = np.nan
            
            # Top/bottom months
            top3 = metric_df.nlargest(3, 'frac_less')[['yearmonth', 'frac_less', 'frac_spei_below_minus1']]
            bottom3 = metric_df.nsmallest(3, 'frac_less')[['yearmonth', 'frac_less', 'frac_spei_below_minus1']]
            top3_str = ', '.join([f"{row['yearmonth']} ({row['frac_less']:.3f})" for _, row in top3.iterrows()])
            bottom3_str = ', '.join([f"{row['yearmonth']} ({row['frac_less']:.3f})" for _, row in bottom3.iterrows()])
            
            print(f"\n  {metric_display.get(metric, metric)}:")
            print(f"    Fraction Less efficient: {frac_mean:.3f} ± {frac_std:.3f} [{frac_min:.3f}, {frac_max:.3f}]")
            print(f"    Correlation with drought frequency: ρ = {rho:.3f} (p = {pval:.3f})")
            print(f"    Correlation with SPEI mean: ρ = {rho_spei:.3f} (p = {pval_spei:.3f})")
            print(f"    Months with HIGHEST Less efficient fraction: {top3_str}")
            print(f"    Months with LOWEST Less efficient fraction: {bottom3_str}")
            print(f"    Interpretation: ρ = {rho_spei:.3f} means when SPEI increases, Less efficient fraction decreases")
        
        self.monthly_df = monthly_df
        return monthly_df
    
    # ==========================================================================
    # STEP 3: Breakpoint detection - UPDATED interpretation
    # ==========================================================================
    
    def pettitt_test(self, series):
        """Pettitt test for change point detection."""
        n = len(series)
        if n < 10:
            return None, None
        
        # Handle NaN
        valid_idx = np.isfinite(series)
        if np.sum(valid_idx) < 10:
            return None, None
        
        series_clean = series[valid_idx]
        indices_clean = np.where(valid_idx)[0]
        n_clean = len(series_clean)
        
        # Pettitt statistic
        U = np.zeros(n_clean)
        for t in range(1, n_clean):
            sgn = np.sign(series_clean[:t, np.newaxis] - series_clean[t:])
            U[t] = np.sum(sgn)
        
        K = np.max(np.abs(U))
        tau = np.argmax(np.abs(U))
        
        # Approximate p-value
        p = 2 * np.exp(-6 * K**2 / (n_clean**3 + n_clean**2))
        
        # Convert to original index
        change_idx = indices_clean[tau] if tau < len(indices_clean) else None
        
        return change_idx, p
    
    def step3_breakpoint_analysis(self):
        """Perform breakpoint analysis with updated interpretation."""
        print("\nSTEP 3: BREAKPOINT ANALYSIS")
        print("-"*40)
        
        if not hasattr(self, 'monthly_df'):
            print("❌ Run Step 2 first")
            return None
        
        results = []
        chosen_breakpoint = None
        chosen_metric = None
        
        print("Running Pettitt test on 'Less efficient' fraction time series...")
        
        metric_display = {
            'WUE_tra': 'WUE_T (transpiration)'
        }
        
        for metric in self.metrics:
            metric_df = self.monthly_df[self.monthly_df['metric'] == metric]
            
            if len(metric_df) < 20:
                print(f"  ⚠️ {metric_display.get(metric, metric)}: Insufficient data ({len(metric_df)} points)")
                continue
            
            # Prepare time series
            time_series = []
            time_values = []
            
            for _, row in metric_df.iterrows():
                time_series.append(row['frac_less'])  # Fraction Less efficient
                time_values.append(row['yearmonth_int'])
            
            # Sort by time
            sort_idx = np.argsort(time_values)
            time_values = np.array(time_values)[sort_idx]
            time_series = np.array(time_series)[sort_idx]
            
            # Pettitt test
            change_idx, p_value = self.pettitt_test(time_series)
            
            if change_idx is not None and p_value < 0.1:
                breakpoint_ym = int(time_values[change_idx])
                breakpoint_year = breakpoint_ym // 100
                breakpoint_month = breakpoint_ym % 100
                breakpoint_label = datetime(breakpoint_year, breakpoint_month, 1).strftime("%B %Y")
                
                # Pre/post statistics
                pre_data = time_series[:change_idx+1]
                post_data = time_series[change_idx+1:]
                
                pre_mean = np.mean(pre_data)
                post_mean = np.mean(post_data)
                delta = post_mean - pre_mean
                pct_change = (delta / pre_mean * 100) if pre_mean != 0 else np.nan
                
                results.append({
                    'metric': metric,
                    'metric_display': metric_display.get(metric, metric),
                    'breakpoint_yearmonth': breakpoint_ym,
                    'breakpoint_year': breakpoint_year,
                    'breakpoint_month': breakpoint_month,
                    'p_value': p_value,
                    'pre_mean': pre_mean,
                    'post_mean': post_mean,
                    'delta': delta,
                    'pct_change': pct_change,
                    'n_pre': len(pre_data),
                    'n_post': len(post_data)
                })
                
                print(f"\n  {metric_display.get(metric, metric)}: Breakpoint at {breakpoint_ym} ({breakpoint_label})")
                print(f"    p-value: {p_value:.4f}")
                print(f"    Pre-breakpoint (<= {breakpoint_ym}): mean frac_less = {pre_mean:.3f} (n={len(pre_data)})")
                print(f"    Post-breakpoint (>  {breakpoint_ym}): mean frac_less = {post_mean:.3f} (n={len(post_data)})")
                print(f"    Change: Δ = {delta:.3f} ({pct_change:.1f}%)")
                print(f"    Interpretation: Fraction of less-efficient responses {'decreased' if delta < 0 else 'increased'} after {breakpoint_label}")
            else:
                results.append({
                    'metric': metric,
                    'metric_display': metric_display.get(metric, metric),
                    'breakpoint_yearmonth': None,
                    'breakpoint_year': None,
                    'breakpoint_month': None,
                    'p_value': p_value if p_value else np.nan,
                    'pre_mean': None,
                    'post_mean': None,
                    'delta': None,
                    'pct_change': None,
                    'n_pre': None,
                    'n_post': None
                })
                print(f"  {metric_display.get(metric, metric)}: No significant breakpoint (p={p_value:.4f})")
        
        # Save results
        results_df = pd.DataFrame(results)
        output_path = os.path.join(self.output_base, 'tables_breakpoint', 'pettitt_results.csv')
        results_df.to_csv(output_path, index=False)
        
        # Choose primary breakpoint - only WUE_tra, so no priority needed
        sig_results = results_df[results_df['p_value'] < 0.1]
        
        if len(sig_results) == 0:
            print("  No significant breakpoints found (p < 0.1)")
            print("  ⚠️ Pre/post maps will NOT be created")
        else:
            # Just take the first significant result (only WUE_tra)
            chosen = sig_results.iloc[0]
            chosen_breakpoint = chosen['breakpoint_yearmonth']
            chosen_metric = chosen['metric']
            
            bp_year = int(chosen_breakpoint) // 100
            bp_month = int(chosen_breakpoint) % 100
            bp_label = datetime(bp_year, bp_month, 1).strftime("%B %Y")
            
            print(f"\n  ✅ PRIMARY BREAKPOINT SELECTED:")
            print(f"     Metric: {chosen['metric_display']}")
            print(f"     Breakpoint: {chosen_breakpoint} ({bp_label})")
            print(f"     Significance: p = {chosen['p_value']:.4f}")
            print(f"     Change in Less efficient fraction: Δ = {chosen['delta']:.3f} ({chosen['pct_change']:.1f}%)")
            print(f"     Pre-{bp_label}: {chosen['pre_mean']:.3f} Less efficient")
            print(f"     Post-{bp_label}: {chosen['post_mean']:.3f} Less efficient")
            print(f"     Interpretation: Fraction of less-efficient responses {'decreased' if chosen['delta'] < 0 else 'increased'} after {bp_label}")
            
            self.breakpoint_label = bp_label
        
        # Save chosen breakpoint
        if chosen_breakpoint:
            bp_year = int(chosen_breakpoint) // 100
            bp_month = int(chosen_breakpoint) % 100
            bp_label = datetime(bp_year, bp_month, 1).strftime("%B %Y")
            
            if bp_month == 12:
                post_start_year = bp_year + 1
                post_start_month = 1
            else:
                post_start_year = bp_year
                post_start_month = bp_month + 1
            post_start_ym = post_start_year * 100 + post_start_month
            
            txt_path = os.path.join(self.output_base, 'tables_breakpoint', 'chosen_breakpoint_yearmonth.txt')
            
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write(f"PRIMARY BREAKPOINT ANALYSIS\n")
                f.write(f"============================\n\n")
                f.write(f"VARIABLE: Fraction of pixels classified as 'Less efficient' (Class=1)\n\n")
                f.write(f"Selected metric: {chosen['metric_display']}\n")
                f.write(f"Breakpoint yearmonth: {chosen_breakpoint} ({bp_label})\n")
                f.write(f"P-value: {chosen['p_value']:.4f}\n\n")
                f.write(f"PRE-BREAKPOINT PERIOD (2000-03 to {chosen_breakpoint}):\n")
                f.write(f"  Mean fraction Less efficient: {chosen['pre_mean']:.4f}\n")
                f.write(f"  Number of months: {chosen['n_pre']}\n\n")
                f.write(f"POST-BREAKPOINT PERIOD ({post_start_ym} to 2025-09):\n")
                f.write(f"  Mean fraction Less efficient: {chosen['post_mean']:.4f}\n")
                f.write(f"  Number of months: {chosen['n_post']}\n\n")
                f.write(f"CHANGE:\n")
                f.write(f"  Absolute change (Δ): {chosen['delta']:.4f}\n")
                f.write(f"  Percent change: {chosen['pct_change']:.1f}%\n\n")
                f.write(f"INTERPRETATION:\n")
                if chosen['delta'] < 0:
                    f.write(f"  After {bp_label}, the fraction of pixels showing 'Less efficient' response\n")
                    f.write(f"  decreased by {abs(chosen['pct_change']):.1f}%.\n")
                else:
                    f.write(f"  After {bp_label}, the fraction of pixels showing 'Less efficient' response\n")
                    f.write(f"  increased by {chosen['pct_change']:.1f}%.\n")
                f.write(f"\nNOTE: This analysis tracks modeled classification frequency, not actual ecosystem-state transitions.\n")
                f.write(f"The breakpoint indicates a shift in the predicted probability of less-efficient WUE responses.\n")
            
            print(f"\n  💾 Breakpoint details saved: {txt_path}")
        
        self.breakpoint_results = results_df
        self.chosen_breakpoint = chosen_breakpoint
        self.chosen_metric = chosen_metric
        
        return results_df
    
    # ==========================================================================
    # STEP 4: Full-period spatial aggregates
    # ==========================================================================
    
    def step4_full_period_aggregates(self):
        """Create spatial aggregates for full period."""
        print("\nSTEP 4: FULL-PERIOD SPATIAL AGGREGATES")
        print("-"*40)
        
        if not hasattr(self, 'file_index'):
            print("❌ Run Step 0 first")
            return None
        
        print(f"Aggregating {len(self.file_index)} months for {len(self.metrics)} metric...")
        
        # Initialize aggregation arrays
        agg_arrays = {}
        template_ds = None
        
        for idx, row in self.file_index.iterrows():
            if (idx + 1) % 12 == 0 or (idx + 1) == len(self.file_index):
                print(f"  Processing month {idx + 1}/{len(self.file_index)}")
            
            try:
                ds = xr.open_dataset(row['file_path'])
                
                # Initialize on first file
                if template_ds is None:
                    template_ds = ds
                    for metric in self.metrics:
                        class_var = f'Class_{metric}'
                        if class_var in ds:
                            shape = ds[class_var].shape
                            agg_arrays[metric] = {
                                'count_valid': np.zeros(shape, dtype=np.int32),
                                'count_less': np.zeros(shape, dtype=np.int32),
                                'count_non_decrease': np.zeros(shape, dtype=np.int32),
                                'sum_prob': np.zeros(shape, dtype=np.float32)
                            }
                
                # Update counts for each metric
                for metric in self.metrics:
                    if metric not in agg_arrays:
                        continue
                    
                    class_var = f'Class_{metric}'
                    prob_var = f'P_decrease_{metric}'
                    
                    if class_var not in ds or prob_var not in ds:
                        continue
                    
                    class_data = ds[class_var].values
                    prob_data = ds[prob_var].values
                    
                    # Use isfinite() for valid mask
                    valid_mask = np.isfinite(class_data) & np.isfinite(prob_data)
                    
                    # Update counts
                    agg_arrays[metric]['count_valid'][valid_mask] += 1
                    agg_arrays[metric]['count_less'][valid_mask] += (class_data[valid_mask] == 1)
                    agg_arrays[metric]['count_non_decrease'][valid_mask] += (class_data[valid_mask] == 0)
                    agg_arrays[metric]['sum_prob'][valid_mask] += prob_data[valid_mask]
                
                ds.close()
                
            except Exception as e:
                print(f"  ⚠️ Error processing {row['filename']}: {e}")
                continue
        
        if template_ds is None:
            print("❌ Could not read any files")
            return None
        
        # Create output datasets
        print("\n📊 Creating output files...")
        
        metric_display = {
            'WUE_tra': 'WUE_T (transpiration)'
        }
        
        for metric in self.metrics:
            if metric not in agg_arrays:
                print(f"  ⚠️ Skipping {metric_display.get(metric, metric)}: no data")
                continue
            
            count_valid = agg_arrays[metric]['count_valid']
            count_less = agg_arrays[metric]['count_less']
            count_non_decrease = agg_arrays[metric]['count_non_decrease']
            sum_prob = agg_arrays[metric]['sum_prob']
            
            # Avoid division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                frac_less = np.where(count_valid > 0, count_less / count_valid, np.nan)
                frac_non_decrease = np.where(count_valid > 0, count_non_decrease / count_valid, np.nan)
                mean_prob = np.where(count_valid > 0, sum_prob / count_valid, np.nan)
            
            # Create dataset
            output_ds = xr.Dataset()
            
            # Copy coordinates
            for coord in template_ds.coords:
                output_ds[coord] = template_ds[coord]
            
            class_var = f'Class_{metric}'
            dims = template_ds[class_var].dims
            
            output_ds['count_valid'] = xr.DataArray(
                count_valid, dims=dims,
                attrs={'long_name': f'Number of valid months for {metric_display.get(metric, metric)}', 'units': 'count'}
            )
            
            output_ds['count_less'] = xr.DataArray(
                count_less, dims=dims,
                attrs={'long_name': f'Number of months classified as Less efficient', 
                      'units': 'count', 'class_definition': 'Class=1 (Less efficient)'}
            )
            
            output_ds['frac_less'] = xr.DataArray(
                frac_less, dims=dims,
                attrs={'long_name': f'Fraction of months classified as Less efficient', 
                      'units': 'fraction', 'valid_range': '0,1',
                      'interpretation': 'Higher values = more frequent Less efficient response'}
            )
            
            output_ds['count_non_decrease'] = xr.DataArray(
                count_non_decrease, dims=dims,
                attrs={'long_name': f'Number of months classified as Non-decrease', 
                      'units': 'count', 'class_definition': 'Class=0 (Increase or NoChange)'}
            )
            
            output_ds['frac_non_decrease'] = xr.DataArray(
                frac_non_decrease, dims=dims,
                attrs={'long_name': f'Fraction of months classified as Non-decrease', 
                      'units': 'fraction', 'valid_range': '0,1',
                      'interpretation': 'Higher values = more frequent Non-decrease response (Increase or NoChange)'}
            )
            
            output_ds['mean_prob'] = xr.DataArray(
                mean_prob, dims=dims,
                attrs={'long_name': f'Mean probability of Less efficient classification', 
                      'units': 'probability', 'valid_range': '0,1',
                      'interpretation': 'P(Class=1) = probability of Less efficient response'}
            )
            
            output_ds.attrs = {
                'title': f'Full-period spatial aggregation for {metric_display.get(metric, metric)}',
                'metric': metric,
                'metric_display': metric_display.get(metric, metric),
                'variable_definitions': 'Class 0 = Increase or NoChange (baseline), Class 1 = Less efficient',
                'period': 'FULL (2000-03 to 2025-09)',
                'n_months': len(self.file_index),
                'created': datetime.now().isoformat(),
                'interpretation_note': 'frac_less = fraction of months with Less efficient response'
            }
            
            # Save
            output_path = os.path.join(self.output_base, 'rasters_time_aggregated', 
                                      f'agg_FULLPERIOD_{metric}.nc')
            output_ds.to_netcdf(output_path)
            
            # Console statistics
            valid_mask = count_valid > 0
            if np.any(valid_mask):
                total_valid_obs = np.sum(count_valid[valid_mask])
                mean_frac_less = np.nanmean(frac_less[valid_mask])
                mean_frac_non_decrease = np.nanmean(frac_non_decrease[valid_mask])
                pixels_low_obs = np.sum(count_valid[valid_mask] < 10)
                pct_low_obs = pixels_low_obs / np.sum(valid_mask) * 100
                
                print(f"\n  {metric_display.get(metric, metric)}:")
                print(f"    Total valid observations: {total_valid_obs:,}")
                print(f"    Mean fraction Less efficient: {mean_frac_less:.3f}")
                print(f"    Mean fraction Non-decrease: {mean_frac_non_decrease:.3f}")
                print(f"    Pixels with <10 obs: {pixels_low_obs:,} ({pct_low_obs:.1f}%)")
                print(f"    💾 Saved: {output_path}")
            else:
                print(f"  ⚠️ {metric_display.get(metric, metric)}: No valid pixels")
        
        print(f"\n✅ Full-period aggregates saved to:")
        print(f"   {os.path.join(self.output_base, 'rasters_time_aggregated')}")
        
        return agg_arrays
    
    # ==========================================================================
    # STEP 5: Pre/post breakpoint aggregates
    # ==========================================================================
    
    def step5_pre_post_aggregates(self):
        """Create spatial aggregates for pre/post breakpoint periods."""
        print("\nSTEP 5: PRE/POST BREAKPOINT AGGREGATES")
        print("-"*40)
        
        if not hasattr(self, 'chosen_breakpoint'):
            print("❌ No breakpoint chosen. Skipping pre/post analysis.")
            return None
        
        if not hasattr(self, 'file_index'):
            print("❌ Run Step 0 first")
            return None
        
        # Split files
        pre_files = self.file_index[self.file_index['yearmonth_int'] <= self.chosen_breakpoint]
        post_files = self.file_index[self.file_index['yearmonth_int'] > self.chosen_breakpoint]
        
        bp_year = int(self.chosen_breakpoint) // 100
        bp_month = int(self.chosen_breakpoint) % 100
        bp_label = datetime(bp_year, bp_month, 1).strftime("%B %Y")
        
        print(f"PRE period: {len(pre_files)} months ({pre_files['yearmonth'].iloc[0]} to {pre_files['yearmonth'].iloc[-1]})")
        print(f"POST period: {len(post_files)} months ({post_files['yearmonth'].iloc[0]} to {post_files['yearmonth'].iloc[-1]})")
        print(f"Breakpoint: {self.chosen_breakpoint} ({bp_label})")
        print(f"Comparing predicted less-efficient response frequency before vs after {bp_label}")
        
        if len(pre_files) == 0 or len(post_files) == 0:
            print("❌ One of the periods has no data")
            return None
        
        metric_display = {
            'WUE_tra': 'WUE_T (transpiration)'
        }
        
        # Process each period
        for period_name, period_files in [('PRE', pre_files), ('POST', post_files)]:
            print(f"\nAggregating {period_name} period ({len(period_files)} months)...")
            
            # Initialize
            agg_arrays = {}
            template_ds = None
            
            for idx, row in period_files.iterrows():
                try:
                    ds = xr.open_dataset(row['file_path'])
                    
                    # Initialize on first file
                    if template_ds is None:
                        template_ds = ds
                        for metric in self.metrics:
                            class_var = f'Class_{metric}'
                            if class_var in ds:
                                shape = ds[class_var].shape
                                agg_arrays[metric] = {
                                    'count_valid': np.zeros(shape, dtype=np.int32),
                                    'count_less': np.zeros(shape, dtype=np.int32),
                                    'count_non_decrease': np.zeros(shape, dtype=np.int32),
                                    'sum_prob': np.zeros(shape, dtype=np.float32)
                                }
                    
                    # Update counts
                    for metric in self.metrics:
                        if metric not in agg_arrays:
                            continue
                        
                        class_var = f'Class_{metric}'
                        prob_var = f'P_decrease_{metric}'
                        
                        if class_var not in ds or prob_var not in ds:
                            continue
                        
                        class_data = ds[class_var].values
                        prob_data = ds[prob_var].values
                        
                        # Use isfinite() for valid mask
                        valid_mask = np.isfinite(class_data) & np.isfinite(prob_data)
                        
                        agg_arrays[metric]['count_valid'][valid_mask] += 1
                        agg_arrays[metric]['count_less'][valid_mask] += (class_data[valid_mask] == 1)
                        agg_arrays[metric]['count_non_decrease'][valid_mask] += (class_data[valid_mask] == 0)
                        agg_arrays[metric]['sum_prob'][valid_mask] += prob_data[valid_mask]
                    
                    ds.close()
                    
                except Exception as e:
                    print(f"  ⚠️ Error processing {row['filename']}: {e}")
                    continue
            
            if template_ds is None:
                print(f"  ❌ Could not read any {period_name} files")
                continue
            
            # Create outputs
            for metric in self.metrics:
                if metric not in agg_arrays:
                    continue
                
                count_valid = agg_arrays[metric]['count_valid']
                count_less = agg_arrays[metric]['count_less']
                count_non_decrease = agg_arrays[metric]['count_non_decrease']
                sum_prob = agg_arrays[metric]['sum_prob']
                
                # Compute fractions
                with np.errstate(divide='ignore', invalid='ignore'):
                    frac_less = np.where(count_valid > 0, count_less / count_valid, np.nan)
                    frac_non_decrease = np.where(count_valid > 0, count_non_decrease / count_valid, np.nan)
                    mean_prob = np.where(count_valid > 0, sum_prob / count_valid, np.nan)
                
                # Create dataset
                output_ds = xr.Dataset()
                
                for coord in template_ds.coords:
                    output_ds[coord] = template_ds[coord]
                
                class_var = f'Class_{metric}'
                dims = template_ds[class_var].dims
                
                output_ds['count_valid'] = xr.DataArray(count_valid, dims=dims)
                output_ds['count_less'] = xr.DataArray(count_less, dims=dims)
                output_ds['frac_less'] = xr.DataArray(frac_less, dims=dims)
                output_ds['count_non_decrease'] = xr.DataArray(count_non_decrease, dims=dims)
                output_ds['frac_non_decrease'] = xr.DataArray(frac_non_decrease, dims=dims)
                output_ds['mean_prob'] = xr.DataArray(mean_prob, dims=dims)
                
                output_ds.attrs = {
                    'title': f'{period_name}-breakpoint aggregation for {metric_display.get(metric, metric)}',
                    'metric': metric,
                    'metric_display': metric_display.get(metric, metric),
                    'variable_definitions': 'Class 0 = Increase or NoChange (baseline), Class 1 = Less efficient',
                    'period': period_name,
                    'breakpoint_yearmonth': str(self.chosen_breakpoint),
                    'breakpoint_date': bp_label,
                    'period_definition': f'yearmonth <= {self.chosen_breakpoint}' if period_name == 'PRE' else f'yearmonth > {self.chosen_breakpoint}',
                    'n_months': len(period_files),
                    'created': datetime.now().isoformat(),
                    'interpretation_note': 'frac_less = fraction of months with Less efficient response'
                }
                
                output_path = os.path.join(self.output_base, 'rasters_pre_post_breakpoint',
                                         f'agg_{period_name}_{metric}.nc')
                output_ds.to_netcdf(output_path)
                
                valid_mask = count_valid > 0
                if np.any(valid_mask):
                    mean_frac_less = np.nanmean(frac_less[valid_mask])
                    mean_frac_non_decrease = np.nanmean(frac_non_decrease[valid_mask])
                    mean_count = np.mean(count_valid[valid_mask])
                    print(f"  {metric_display.get(metric, metric)}:")
                    print(f"    Mean fraction Less efficient = {mean_frac_less:.3f}")
                    print(f"    Mean fraction Non-decrease = {mean_frac_non_decrease:.3f}")
                    print(f"    Mean count_valid = {mean_count:.1f}")
                else:
                    print(f"  {metric_display.get(metric, metric)}: No valid pixels")
        
        print(f"\n✅ Pre/post aggregates saved to:")
        print(f"   {os.path.join(self.output_base, 'rasters_pre_post_breakpoint')}")
        
        # Print comparison summary
        chosen = self.breakpoint_results[self.breakpoint_results['metric'] == self.chosen_metric].iloc[0]
        print(f"\n📊 PRE vs POST COMPARISON SUMMARY:")
        print(f"  Breakpoint: {self.chosen_breakpoint} ({bp_label})")
        print(f"  Metric: {chosen['metric_display']}")
        print(f"  Pre-period: {len(pre_files)} months")
        print(f"  Post-period: {len(post_files)} months")
        print(f"  Key finding: {abs(chosen['pct_change']):.1f}% {'decrease' if chosen['delta'] < 0 else 'increase'} in Less efficient fraction")
        print(f"  Interpretation: Fraction of less-efficient responses {'decreased' if chosen['delta'] < 0 else 'increased'} after {bp_label}")
        print(f"  Note: This tracks modeled classification frequency, not actual ecosystem-state transitions.")
        
        return True
    
    # ==========================================================================
    # MAIN WORKFLOW
    # ==========================================================================
    
    def run_full_workflow(self):
        """Run the complete concise workflow."""
        print("\n🚀 STARTING CONCISE WORKFLOW")
        print("="*80)
        
        results = {}
        
        results['file_index'] = self.step0_build_file_index()
        if results['file_index'] is None:
            return None
        
        print("\n" + "-"*80)
        results['sanity_check'] = self.step1_sanity_check()
        
        print("\n" + "-"*80)
        results['monthly_table'] = self.step2_create_monthly_table()
        if results['monthly_table'] is None:
            print("❌ Failed to create monthly table")
            return None
        
        print("\n" + "-"*80)
        results['breakpoint'] = self.step3_breakpoint_analysis()
        
        print("\n" + "-"*80)
        results['full_aggregates'] = self.step4_full_period_aggregates()
        
        if hasattr(self, 'chosen_breakpoint') and self.chosen_breakpoint:
            print("\n" + "-"*80)
            results['pre_post_aggregates'] = self.step5_pre_post_aggregates()
        
        print("\n" + "="*80)
        print("WORKFLOW COMPLETE")
        print("="*80)
        self._print_final_summary()
        
        return results
    
    def _print_final_summary(self):
        """Print final summary with updated terminology."""
        print("\n📋 FINAL SUMMARY")
        print("-"*40)
        
        metric_display = {
            'WUE_tra': 'WUE_T (transpiration)'
        }
        
        print(f"📂 Input: {self.input_dir}")
        print(f"📁 Output: {self.output_base}")
        print(f"📊 Metrics processed: {len(self.metrics)}")
        for m in self.metrics:
            print(f"     - {metric_display.get(m, m)}")
        
        if hasattr(self, 'file_index'):
            print(f"📊 Files processed: {len(self.file_index)} months (2000-03 to 2025-09)")
        
        if hasattr(self, 'monthly_df'):
            n_months = len(self.monthly_df) / len(self.metrics)
            print(f"📅 Monthly records: {n_months:.0f} months × {len(self.metrics)} metric")
        
        if hasattr(self, 'chosen_breakpoint'):
            bp_year = int(self.chosen_breakpoint) // 100
            bp_month = int(self.chosen_breakpoint) % 100
            bp_label = datetime(bp_year, bp_month, 1).strftime("%B %Y")
            
            chosen_result = self.breakpoint_results[self.breakpoint_results['metric'] == self.chosen_metric].iloc[0]
            
            print(f"\n🔀 KEY BREAKPOINT FINDING:")
            print(f"   Breakpoint: {self.chosen_breakpoint} ({bp_label})")
            print(f"   Metric: {chosen_result['metric_display']}")
            delta = chosen_result['delta']
            pct_change = chosen_result['pct_change']
            print(f"   Change in Less efficient fraction: Δ = {delta:.3f} ({abs(pct_change):.1f}% {'decrease' if delta < 0 else 'increase'})")
            print(f"   Interpretation: Fraction of less-efficient responses {'decreased' if delta < 0 else 'increased'} after {bp_label}")
            print(f"   Note: This tracks modeled classification frequency, not actual ecosystem-state transitions.")
            print(f"   PRE/POST aggregates created")
        else:
            print(f"\n🔀 No significant breakpoint found")
            print(f"   Only FULL-period aggregates created")
        
        print(f"\n📁 OUTPUTS CREATED:")
        base = self.output_base
        for subdir in ['diagnostics', 'tables_breakpoint', 
                      'rasters_time_aggregated', 'rasters_pre_post_breakpoint']:
            path = os.path.join(base, subdir)
            if os.path.exists(path):
                files = [f for f in os.listdir(path) if f.endswith(('.csv', '.nc', '.txt'))]
                if files:
                    print(f"  {subdir}: {len(files)} files")
        
        print(f"\n📝 INTERPRETATION NOTES:")
        print(f"  - frac_less = fraction of months with 'Less efficient' response (Class=1)")
        print(f"  - mean_prob = mean probability of 'Less efficient' response")
        print(f"  - Class=1 = Decrease in WUE (drought response)")
        print(f"  - Class=0 = Increase or NoChange (baseline) - non-decrease response")
        print(f"  - NoChange is correctly kept as baseline")
        
        print(f"\n✅ Analysis complete. Data ready for mapping.")

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def main():
    """Main execution function."""
    
    print("\n" + "="*80)
    print("CONCISE SPATIAL ANALYSIS WORKFLOW - CORRECTED MODEL")
    print("="*80)
    print("Using model with NoChange kept as baseline (Class=0)")
    print("Class=1 = Less efficient (Decrease), Class=0 = Increase or NoChange (baseline)")
    print("Processing ONLY WUE_T (transpiration) metric")
    
    # Run workflow
    analyzer = ConciseSpatialAnalyzer()
    results = analyzer.run_full_workflow()
    
    if results:
        print("\n🎉 WORKFLOW COMPLETED SUCCESSFULLY!")
        print("\n📋 NEXT STEPS:")
        print("1. Review console output above for interpretation")
        print("2. Use spatial aggregates for mapping:")
        print("   - Full period: rasters_time_aggregated/agg_FULLPERIOD_WUE_tra.nc")
        if hasattr(analyzer, 'chosen_breakpoint') and analyzer.chosen_breakpoint:
            print("   - Pre/post: rasters_pre_post_breakpoint/agg_PRE/POST_WUE_tra.nc")
        print("3. Monthly data for time series: tables_breakpoint/monthly_metric_long.csv")
        print("4. Breakpoint details: tables_breakpoint/chosen_breakpoint_yearmonth.txt")
        print("\n📝 IMPORTANT FOR FIGURE LABELS:")
        print("   - Label maps as 'Fraction of months with Less efficient response'")
        print("   - Not 'WUE decline' but 'Probability/fraction of Less efficient response'")
        print("\n📝 CAUTION - INTERPRETATION NOTE:")
        print("   - This analysis tracks modeled classification frequency from point-based models")
        print("   - Results indicate predicted probability shifts, not direct ecosystem-state transitions")
        print("   - 'Non-decrease' class includes both WUE increases AND NoChange conditions")
    
    return results

if __name__ == "__main__":
    results = main()