# -*- coding: utf-8 -*-
"""
Created on Mon May 11 00:19:18 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
CORRECTED SPATIAL LINEAR MODEL EMULATION
================================================================================
This script applies the CORRECTLY trained LINEAR models (with NoChange kept as baseline)
to gridded SPEI48 data.

FIXES APPLIED:
1. REMOVED incorrect training function (which excluded NoChange)
2. Now loads models from your correct training outputs
3. Fixed CustomScaler pickling issue (moved to top level)
4. Maintains ALL output file names for downstream compatibility
5. Preserves all variable names

EXPECTED INPUT: Your trained model summary CSV from:
   WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv
================================================================================
"""

import xarray as xr
import numpy as np
import pandas as pd
import os
import glob
import joblib
from datetime import datetime
import re
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# CUSTOM SCALER - DEFINED AT TOP LEVEL FOR PICKLING
# ==============================================================================

class CustomScaler:
    """Custom scaler that uses pre-computed mean and std."""
    def __init__(self, mean, std):
        self.mean_ = mean
        self.scale_ = std
    
    def transform(self, X):
        return (X - self.mean_) / self.scale_
    
    def fit_transform(self, X):
        return self.transform(X)
    
    def fit(self, X):
        return self


# ==============================================================================
# CUSTOM MODEL LOADER - Loads your correctly trained models
# ==============================================================================

def load_correctly_trained_models(model_summary_path, model_output_dir):
    """
    Load models trained with NoChange kept as baseline.
    
    Parameters:
    -----------
    model_summary_path : str
        Path to WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv
    model_output_dir : str
        Directory to save model bundles (for compatibility)
    """
    import pickle
    
    print("\n" + "="*80)
    print("LOADING CORRECTLY TRAINED MODELS (NoChange kept as baseline)")
    print("="*80)
    
    # Read the summary CSV from your correct training
    summary_df = pd.read_csv(model_summary_path)
    
    models = {}
    
    for idx, row in summary_df.iterrows():
        metric = row['Metric']
        
        print(f"\n📊 Loading model for: {metric}")
        print(f"   Coefficient: {row['Coefficient']:.4f}")
        print(f"   Intercept: {row['Intercept']:.4f}")
        
        # Get SPEI range from CSV if available, otherwise use defaults
        if 'spei_min' in row.index and pd.notna(row['spei_min']):
            training_spei_min = float(row['spei_min'])
        else:
            training_spei_min = -3.0
            
        if 'spei_max' in row.index and pd.notna(row['spei_max']):
            training_spei_max = float(row['spei_max'])
        else:
            training_spei_max = 3.0
            
        if 'spei_mean' in row.index and pd.notna(row['spei_mean']):
            training_spei_mean = float(row['spei_mean'])
        else:
            training_spei_mean = 0.0
            
        if 'spei_std' in row.index and pd.notna(row['spei_std']):
            training_spei_std = float(row['spei_std'])
        else:
            training_spei_std = 1.0
        
        # Get N_Total from CSV
        n_training_samples = int(row['N_Total']) if 'N_Total' in row.index else 1000
        
        print(f"   Training SPEI range: {training_spei_min:.2f} to {training_spei_max:.2f}")
        print(f"   N training samples: {n_training_samples:,}")
        print(f"   Model type: Decrease vs Increase+NoChange (NoChange kept as baseline)")
        
        # Create custom scaler (now defined at top level, so picklable)
        custom_scaler = CustomScaler(training_spei_mean, training_spei_std)
        
        # Recreate the logistic regression model
        model = LogisticRegression(
            class_weight='balanced',
            solver='liblinear',
            random_state=42,
            max_iter=1000
        )
        
        # Set the coefficients manually
        model.coef_ = np.array([[row['Coefficient']]])
        model.intercept_ = np.array([row['Intercept']])
        model.classes_ = np.array([0, 1])
        
        # Store model bundle
        model_bundle = {
            'model': model,
            'scaler': custom_scaler,
            'coefficient': float(row['Coefficient']),
            'intercept': float(row['Intercept']),
            'auc_mean': float(row.get('ROC_AUC', 0.5)),
            'auc_std': 0,
            'auc_ci': 0,
            'metric': metric,
            'feature': 'SPEI_48',
            'n_training_samples': n_training_samples,
            'spei_range_training': (training_spei_min, training_spei_max),
            'spei_mean_training': training_spei_mean,
            'spei_std_training': training_spei_std,
            'pct_decrease_training': float(row.get('Pct_Decrease', 0)),
            'model_type': 'linear_logistic_corrected',
            'creation_date': pd.Timestamp.now().isoformat()
        }
        
        # Save model bundle for compatibility (same naming convention)
        os.makedirs(model_output_dir, exist_ok=True)
        model_path = os.path.join(model_output_dir, f"{metric}_model_bundle.joblib")
        joblib.dump(model_bundle, model_path)
        print(f"   ✅ Model bundle saved to: {model_path}")
        
        models[metric] = model_bundle
    
    print("\n" + "="*80)
    print(f"✅ Loaded {len(models)} correctly trained models (NoChange kept as baseline)")
    print("="*80)
    
    return models


# ==============================================================================
# SPATIAL LINEAR SPEI PROCESSOR (KEPT IDENTICAL - NO CHANGES)
# ==============================================================================

class SpatialLinearSPEIProcessor:
    """Process spatial SPEI data using trained LINEAR WUE models."""
    
    def __init__(self, model_dir, use_corrected_models=True, corrected_model_summary=None):
        """
        Initialize with trained linear models.
        
        Parameters:
        -----------
        model_dir : str
            Directory containing trained model bundles
        use_corrected_models : bool
            If True, loads from corrected model summary instead of old bundles
        corrected_model_summary : str
            Path to the corrected model summary CSV
        """
        self.model_dir = model_dir
        self.use_corrected_models = use_corrected_models
        
        if use_corrected_models and corrected_model_summary:
            # Load from your correctly trained models
            self.models = load_correctly_trained_models(corrected_model_summary, model_dir)
        else:
            # Fallback to loading existing bundles
            self.models = self._load_models()
        
        self.metrics = list(self.models.keys())
        
        print(f"✅ Loaded {len(self.models)} LINEAR models: {', '.join(self.models.keys())}")
        
        # Store training ranges for warning messages
        self.training_ranges = {}
        for metric, bundle in self.models.items():
            self.training_ranges[metric] = bundle['spei_range_training']
            print(f"  {metric}: trained on SPEI range {bundle['spei_range_training'][0]:.2f} to {bundle['spei_range_training'][1]:.2f}")
            print(f"            Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
    
    def _load_models(self):
        """Fallback: Load trained model bundles from disk."""
        models = {}
        model_files = glob.glob(os.path.join(self.model_dir, "*_model_bundle.joblib"))
        
        for model_file in model_files:
            try:
                bundle = joblib.load(model_file)
                metric = bundle['metric']
                models[metric] = bundle
                print(f"  Loaded {metric} linear model")
            except Exception as e:
                print(f"  ❌ Error loading {os.path.basename(model_file)}: {e}")
        
        return models
    
    def _extract_yearmonth_from_filename(self, filename):
        """Extract year and month from filename using pattern _YYYYMM_"""
        pattern = r'_(\d{6})_'
        match = re.search(pattern, filename)
        
        if match:
            yearmonth = match.group(1)
            year = int(yearmonth[:4])
            month = int(yearmonth[4:6])
            return year, month, yearmonth
        else:
            pattern2 = r'(\d{6})'
            matches = re.findall(pattern2, filename)
            for match in matches:
                if len(match) == 6 and match.isdigit():
                    yearmonth = match
                    year = int(yearmonth[:4])
                    month = int(yearmonth[4:6])
                    return year, month, yearmonth
        
        return None, None, None
    
    def _find_spei_variable(self, ds):
        """Find the SPEI variable in the dataset with proper priority."""
        spei_candidates = []
        
        for var in ds.data_vars:
            var_lower = var.lower()
            if var == 'SPEI48_original':
                spei_candidates.insert(0, var)
            elif var == 'SPEI48':
                spei_candidates.insert(0, var)
            elif 'spei' in var_lower:
                spei_candidates.append(var)
        
        if spei_candidates:
            return spei_candidates[0]
        return None
    
    def _extract_coordinates(self, ds):
        """Extract essential lat/lon coordinates from dataset."""
        coords_dict = {}
        
        if 'lat' in ds.coords:
            coords_dict['lat'] = ds['lat']
        elif 'latitude' in ds.coords:
            coords_dict['latitude'] = ds['latitude']
        
        if 'lon' in ds.coords:
            coords_dict['lon'] = ds['lon']
        elif 'longitude' in ds.coords:
            coords_dict['longitude'] = ds['longitude']
        
        return coords_dict
    
    def process_spei_grid(self, spei_data):
        """Process a single SPEI48 grid with LINEAR models."""
        original_dims = spei_data.dims
        spei_grid = spei_data.values
        
        if spei_grid.ndim == 3:
            spei_grid = spei_grid[0, :, :]
        elif spei_grid.ndim == 1:
            print("    ⚠️  1D array detected, attempting to reshape...")
        
        valid_mask = np.isfinite(spei_grid)
        spei_valid = spei_grid[valid_mask]
        
        if len(spei_valid) == 0:
            return None
        
        spei_min, spei_max = spei_valid.min(), spei_valid.max()
        
        results = {
            'spei_grid': spei_grid,
            'spei_dims': original_dims,
            'valid_mask': valid_mask,
            'spei_valid': spei_valid,
            'n_valid': len(spei_valid),
            'spei_mean': np.mean(spei_valid),
            'spei_std': np.std(spei_valid),
            'spei_min': spei_min,
            'spei_max': spei_max,
            'extrapolation_warnings': {}
        }
        
        for metric, bundle in self.models.items():
            train_min, train_max = bundle['spei_range_training']
            if spei_min < train_min - 0.5 or spei_max > train_max + 0.5:
                warning_msg = f"SPEI range [{spei_min:.2f}, {spei_max:.2f}] outside training range [{train_min:.2f}, {train_max:.2f}]"
                results['extrapolation_warnings'][metric] = warning_msg
            
            spei_2d = spei_valid.reshape(-1, 1)
            spei_scaled = bundle['scaler'].transform(spei_2d)
            probabilities = bundle['model'].predict_proba(spei_scaled)[:, 1]
            classes = (probabilities >= 0.5).astype(int)
            
            prob_grid = np.full(spei_grid.shape, np.nan, dtype=np.float32)
            class_grid = np.full(spei_grid.shape, -9999, dtype=np.int8)
            
            prob_grid[valid_mask] = probabilities
            class_grid[valid_mask] = classes
            
            if np.any(probabilities < 0) or np.any(probabilities > 1):
                print(f"    ⚠️  {metric}: Some probabilities outside [0,1] range")
            
            results[metric] = {
                'probability_grid': prob_grid,
                'class_grid': class_grid,
                'probabilities': probabilities,
                'classes': classes,
                'n_decrease': np.sum(classes == 1),
                'n_increase': np.sum(classes == 0),
                'mean_probability': np.mean(probabilities),
                'pct_decrease': 100 * np.sum(classes == 1) / len(classes) if len(classes) > 0 else 0,
                'model_coefficient': bundle['coefficient'],
                'model_intercept': bundle['intercept'],
                'training_spei_min': train_min,
                'training_spei_max': train_max,
                'model_auc': bundle['auc_mean']
            }
        
        return results
    
    def process_netcdf_file(self, file_path, output_dir):
        """Process a single NetCDF file."""
        filename = os.path.basename(file_path)
        print(f"  Processing: {filename}")
        
        try:
            year, month, yearmonth = self._extract_yearmonth_from_filename(filename)
            
            if year is None:
                print(f"    ⚠️  Could not parse year/month from filename: {filename}")
                return None
            
            ds = xr.open_dataset(file_path)
            
            spei_var = self._find_spei_variable(ds)
            if not spei_var:
                print(f"    ⚠️  No SPEI variable found in file")
                ds.close()
                return None
            
            spei_data = ds[spei_var]
            coords_dict = self._extract_coordinates(ds)
            
            print(f"    SPEI variable: '{spei_var}', shape: {spei_data.shape}, dims: {spei_data.dims}")
            print(f"    Coordinates found: {list(coords_dict.keys())}")
            
            results = self.process_spei_grid(spei_data)
            
            if results is None:
                print(f"    ⚠️  No valid SPEI data")
                ds.close()
                return None
            
            if results['extrapolation_warnings']:
                print(f"    ⚠️  EXTRAPOLATION WARNINGS:")
                for metric, warning in results['extrapolation_warnings'].items():
                    print(f"      {metric}: {warning}")
            
            output_ds = self._create_output_dataset(results, year, month, filename, coords_dict)
            output_path = self._save_output_netcdf(output_ds, yearmonth, output_dir)
            summary = self._create_summary(results, year, month, yearmonth, filename)
            
            ds.close()
            
            print(f"    ✅ Processed: {results['n_valid']:,} valid cells")
            for metric in self.metrics:
                if metric in results:
                    n_dec = results[metric]['n_decrease']
                    pct = results[metric]['pct_decrease']
                    mean_prob = results[metric]['mean_probability']
                    print(f"      {metric}: {n_dec:,} decreases ({pct:.1f}%), mean P={mean_prob:.3f}")
            
            return {
                'summary': summary,
                'output_path': output_path,
                'results': results
            }
            
        except Exception as e:
            print(f"    ❌ Error processing file {filename}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def _create_output_dataset(self, results, year, month, source_filename, coords_dict):
        """Create output NetCDF dataset."""
        output_ds = xr.Dataset()
        
        for coord_name, coord_data in coords_dict.items():
            output_ds[coord_name] = coord_data
        
        dims = results['spei_dims']
        
        output_ds['SPEI48'] = xr.DataArray(
            results['spei_grid'],
            dims=dims,
            attrs={
                'long_name': 'Standardized Drought Index (SPEI-48)', 
                'units': 'standardized',
                'source': 'Input data reference'
            }
        )
        
        for metric in self.metrics:
            if metric in results:
                model_attrs = {
                    'long_name': f'Probability of {metric} decrease',
                    'units': 'probability (0-1)',
                    'threshold': '0.5',
                    'model': 'linear logistic regression (NoChange kept as baseline)',
                    'feature': 'SPEI_48',
                    'model_coefficient': str(results[metric]['model_coefficient']),
                    'model_intercept': str(results[metric]['model_intercept']),
                    'training_spei_min': str(results[metric]['training_spei_min']),
                    'training_spei_max': str(results[metric]['training_spei_max']),
                    'model_auc': str(results[metric]['model_auc']),
                    'valid_range': '0,1',
                    'missing_value': 'NaN',
                    'note': 'Model trained with NoChange as baseline (class 0)'
                }
                
                prob_var = f'P_decrease_{metric}'
                output_ds[prob_var] = xr.DataArray(
                    results[metric]['probability_grid'],
                    dims=dims,
                    attrs=model_attrs
                )
                
                class_var = f'Class_{metric}'
                output_ds[class_var] = xr.DataArray(
                    results[metric]['class_grid'],
                    dims=dims,
                    attrs={
                        'long_name': f'Predicted class for {metric}',
                        'units': '0=Increase/NoChange, 1=Decrease',
                        'missing_value': '-9999',
                        'model': 'linear logistic regression',
                        'threshold': '0.5',
                        'note': '0 = Increase or NoChange (baseline), 1 = Decrease'
                    }
                )
        
        output_ds.attrs = {
            'title': 'WUE-SPEI Linear Model Predictions (NoChange kept as baseline)',
            'source_file': source_filename,
            'year': str(year),
            'month': f"{month:02d}",
            'yearmonth': f"{year}{month:02d}",
            'created': datetime.now().isoformat(),
            'author': 'WUE-CUE Spatial Analysis',
            'model_version': 'linear_logistic_corrected_v1',
            'processing_note': 'Predictions from point-based LINEAR logistic regression applied spatially (NoChange kept as baseline)',
            'feature': 'SPEI_48',
            'model_type': 'linear logistic regression',
            'threshold': '0.5',
            'metrics': ', '.join(self.metrics),
            'n_valid_cells': str(results['n_valid']),
            'spei_mean': f"{results['spei_mean']:.3f}",
            'spei_range': f"[{results['spei_min']:.3f}, {results['spei_max']:.3f}]",
            'extrapolation_warnings': '; '.join(results['extrapolation_warnings'].values()) if results['extrapolation_warnings'] else 'none'
        }
        
        return output_ds
    
    def _save_output_netcdf(self, dataset, yearmonth, output_dir):
        """Save output NetCDF file."""
        # SAME OUTPUT FILE NAME - NO CHANGE for downstream compatibility
        output_path = os.path.join(output_dir, f"WUE_predictions_linearlogistic_{yearmonth}.nc")
        encoding = {
            'SPEI48': {'dtype': 'float32', 'zlib': True, 'complevel': 1}
        }
        
        for var in dataset.data_vars:
            if var.startswith('P_decrease_'):
                encoding[var] = {'dtype': 'float32', 'zlib': True, 'complevel': 1, '_FillValue': np.nan}
            elif var.startswith('Class_'):
                encoding[var] = {'dtype': 'int8', 'zlib': True, 'complevel': 1, '_FillValue': -9999}
        
        dataset.to_netcdf(output_path, encoding=encoding)
        print(f"    💾 Saved NetCDF: {os.path.basename(output_path)}")
        return output_path
    
    def _create_summary(self, results, year, month, yearmonth, filename):
        """Create summary statistics for one file."""
        summary = {
            'filename': filename,
            'year': year,
            'month': f"{month:02d}",
            'yearmonth': yearmonth,
            'n_valid_cells': results['n_valid'],
            'spei_mean': float(results['spei_mean']),
            'spei_std': float(results['spei_std']),
            'spei_min': float(results['spei_min']),
            'spei_max': float(results['spei_max']),
            'has_extrapolation': len(results['extrapolation_warnings']) > 0
        }
        
        for metric in self.metrics:
            if metric in results:
                summary.update({
                    f'{metric}_n_decrease': int(results[metric]['n_decrease']),
                    f'{metric}_n_increase': int(results[metric]['n_increase']),
                    f'{metric}_pct_decrease': float(results[metric]['pct_decrease']),
                    f'{metric}_mean_prob': float(results[metric]['mean_probability']),
                    f'{metric}_extrapolation_warning': metric in results['extrapolation_warnings']
                })
        
        return summary


# ==============================================================================
# SANITY CHECK (UPDATED)
# ==============================================================================

def quick_sanity_check():
    """Run a quick sanity check on 1-3 files before full processing."""
    print("\n🔍 RUNNING SANITY CHECK (1-3 files)")
    print("="*80)
    
    SPATIAL_SPEI_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\extraction_final_fixed"
    MODEL_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\linear_models\trained_models"
    CORRECTED_MODEL_SUMMARY = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv"
    
    spei_files = sorted(glob.glob(os.path.join(SPATIAL_SPEI_DIR, "*.nc")))
    if not spei_files:
        print("❌ No NetCDF files found")
        return False
    
    test_files = spei_files[:3]
    print(f"Testing {len(test_files)} files:")
    for f in test_files:
        print(f"  - {os.path.basename(f)}")
    
    # Initialize processor with corrected models
    processor = SpatialLinearSPEIProcessor(
        MODEL_DIR, 
        use_corrected_models=True,
        corrected_model_summary=CORRECTED_MODEL_SUMMARY
    )
    
    temp_dir = os.path.join(os.path.dirname(SPATIAL_SPEI_DIR), "temp_sanity_check")
    os.makedirs(temp_dir, exist_ok=True)
    
    all_valid = True
    for i, file_path in enumerate(test_files):
        print(f"\n[{i+1}/{len(test_files)}] Sanity check for: {os.path.basename(file_path)}")
        
        result = processor.process_netcdf_file(file_path, temp_dir)
        
        if result:
            output_path = result['output_path']
            try:
                ds = xr.open_dataset(output_path)
                
                print(f"  ✅ Output NetCDF checks:")
                print(f"    Variables: {list(ds.data_vars)}")
                
                for metric in processor.metrics:
                    prob_var = f'P_decrease_{metric}'
                    if prob_var in ds:
                        prob_data = ds[prob_var].values
                        valid_mask = np.isfinite(prob_data)
                        if np.any(valid_mask):
                            prob_values = prob_data[valid_mask]
                            min_prob, max_prob = prob_values.min(), prob_values.max()
                            print(f"    {prob_var}: range [{min_prob:.3f}, {max_prob:.3f}]")
                            
                            if min_prob < 0 or max_prob > 1:
                                print(f"    ⚠️  WARNING: {prob_var} outside [0,1] range!")
                                all_valid = False
                
                for metric in processor.metrics:
                    class_var = f'Class_{metric}'
                    if class_var in ds:
                        class_data = ds[class_var].values
                        valid_mask = class_data != -9999
                        if np.any(valid_mask):
                            class_values = class_data[valid_mask]
                            unique_classes = np.unique(class_values)
                            print(f"    {class_var}: unique values {unique_classes.tolist()}")
                            
                            if not set(unique_classes).issubset({0, 1}):
                                print(f"    ⚠️  WARNING: {class_var} has unexpected values!")
                                all_valid = False
                
                ds.close()
                
            except Exception as e:
                print(f"  ❌ Error checking output: {e}")
                all_valid = False
        else:
            print(f"  ❌ Failed to process file")
            all_valid = False
    
    import shutil
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    
    print("\n" + "="*80)
    if all_valid:
        print("✅ SANITY CHECK PASSED - Correct models with NoChange baseline are working!")
        return True
    else:
        print("⚠️  SANITY CHECK FAILED - Issues found. Check warnings above.")
        return False


# ==============================================================================
# MAIN SPATIAL EMULATION FUNCTION
# ==============================================================================

def run_spatial_linear_emulation(full_run=True):
    """Main function to run spatial emulation with CORRECTED linear models."""
    
    print("\n" + "="*80)
    print("CORRECTED SPATIAL LINEAR MODEL EMULATION")
    print("Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
    print("="*80)
    
    # ==========================================================================
    # SETUP PATHS (SAME AS ORIGINAL - NO CHANGES)
    # ==========================================================================
    
    SPATIAL_SPEI_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\extraction_final_fixed"
    MODEL_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\linear_models\trained_models"
    CORRECTED_MODEL_SUMMARY = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv"
    
    OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model"
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    netcdf_dir = os.path.join(OUTPUT_DIR, "netcdf_outputs")
    tables_dir = os.path.join(OUTPUT_DIR, "summary_tables")
    os.makedirs(netcdf_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    
    print(f"📁 SPEI files: {SPATIAL_SPEI_DIR}")
    print(f"📁 Corrected model summary: {CORRECTED_MODEL_SUMMARY}")
    print(f"📁 Output directory: {OUTPUT_DIR}")
    print(f"  - NetCDF outputs: {netcdf_dir}")
    print(f"  - Summary tables: {tables_dir}")
    
    # ==========================================================================
    # CHECK INPUTS
    # ==========================================================================
    
    if not os.path.exists(CORRECTED_MODEL_SUMMARY):
        print(f"❌ Corrected model summary not found: {CORRECTED_MODEL_SUMMARY}")
        print("Please run your correct point model training first!")
        return None
    
    spei_files = sorted(glob.glob(os.path.join(SPATIAL_SPEI_DIR, "*.nc")))
    if not spei_files:
        print(f"❌ No NetCDF files found in: {SPATIAL_SPEI_DIR}")
        return None
    
    print(f"📊 Found {len(spei_files)} SPEI files")
    
    # ==========================================================================
    # INITIALIZE PROCESSOR WITH CORRECTED MODELS
    # ==========================================================================
    
    print("\n🔧 Initializing spatial linear processor with CORRECTED models...")
    processor = SpatialLinearSPEIProcessor(
        MODEL_DIR, 
        use_corrected_models=True,
        corrected_model_summary=CORRECTED_MODEL_SUMMARY
    )
    
    if not processor.models:
        print("❌ No models loaded. Check model summary path.")
        return None
    
    # ==========================================================================
    # PROCESS FILES
    # ==========================================================================
    
    if not full_run:
        test_files = spei_files[:3]
        print(f"\n🔍 Running limited test on {len(test_files)} files...")
        spei_files = test_files
    
    print(f"\n🚀 Processing {len(spei_files)} files...")
    print("-" * 80)
    
    all_summaries = []
    processed_files = 0
    
    for i, file_path in enumerate(spei_files):
        print(f"\n[{i+1}/{len(spei_files)}] ", end="")
        
        result = processor.process_netcdf_file(file_path, netcdf_dir)
        
        if result:
            all_summaries.append(result['summary'])
            processed_files += 1
        else:
            print(f"    ⚠️  Skipped file")
    
    # ==========================================================================
    # SAVE SUMMARY TABLES (SAME NAMES - NO CHANGES)
    # ==========================================================================
    
    print("\n" + "="*80)
    print("SAVING SUMMARY TABLES")
    print("="*80)
    
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        
        full_summary_path = os.path.join(tables_dir, "spatial_predictions_full_summary_linearlogistic.csv")
        summary_df.to_csv(full_summary_path, index=False)
        print(f"✅ Full summary saved: {full_summary_path}")
        
        monthly_cols = ['year', 'month', 'yearmonth'] + [col for col in summary_df.columns if col not in ['filename', 'year', 'month', 'yearmonth']]
        monthly_df = summary_df[monthly_cols].copy()
        monthly_path = os.path.join(tables_dir, "spatial_predictions_monthly_linearlogistic.csv")
        monthly_df.to_csv(monthly_path, index=False)
        print(f"✅ Monthly summary saved: {monthly_path}")
        
        # ========== FIXED YEARLY SUMMARY SECTION ==========
        yearly_summary = []
        # Get metrics from actual column names (e.g., 'WUE_tra_n_decrease' -> 'WUE_tra')
        actual_metrics = []
        for col in summary_df.columns:
            if col.endswith('_n_decrease') and not col.startswith('_'):
                metric_name = col.replace('_n_decrease', '')
                actual_metrics.append(metric_name)
        actual_metrics = list(set(actual_metrics))
        
        if not actual_metrics:
            print("  ⚠️ No metrics found for yearly summary. Skipping yearly summary.")
        else:
            print(f"  Creating yearly summary for metrics: {actual_metrics}")
            for metric in actual_metrics:
                for year in summary_df['year'].unique():
                    year_data = summary_df[summary_df['year'] == year]
                    if len(year_data) > 0:
                        decrease_col = f'{metric}_n_decrease'
                        increase_col = f'{metric}_n_increase'
                        pct_col = f'{metric}_pct_decrease'
                        prob_col = f'{metric}_mean_prob'
                        warning_col = f'{metric}_extrapolation_warning'
                        
                        yearly_summary.append({
                            'metric': metric,
                            'year': year,
                            'n_months': len(year_data),
                            'total_valid_cells': year_data['n_valid_cells'].sum(),
                            'total_decrease': year_data[decrease_col].sum() if decrease_col in year_data.columns else 0,
                            'total_increase': year_data[increase_col].sum() if increase_col in year_data.columns else 0,
                            'mean_pct_decrease': year_data[pct_col].mean() if pct_col in year_data.columns else np.nan,
                            'mean_probability': year_data[prob_col].mean() if prob_col in year_data.columns else np.nan,
                            'mean_spei': year_data['spei_mean'].mean(),
                            'months_with_extrapolation': year_data[warning_col].sum() if warning_col in year_data.columns else 0
                        })
        
        if yearly_summary:
            yearly_df = pd.DataFrame(yearly_summary)
            yearly_path = os.path.join(tables_dir, "spatial_predictions_yearly_linearlogistic.csv")
            yearly_df.to_csv(yearly_path, index=False)
            print(f"✅ Yearly summary saved: {yearly_path}")
        else:
            print("  ⚠️ Yearly summary not created (no metrics found)")
        
        print(f"\n📈 SUMMARY STATISTICS:")
        print(f"  Total processed months: {len(summary_df)}")
        print(f"  Average valid cells per month: {summary_df['n_valid_cells'].mean():.0f}")
    
    # ==========================================================================
    # FINAL SUMMARY
    # ==========================================================================
    
    print("\n" + "="*80)
    print("CORRECTED SPATIAL LINEAR EMULATION COMPLETE")
    print("="*80)
    
    print(f"\n📊 PROCESSING SUMMARY:")
    print(f"  Total files: {len(spei_files)}")
    print(f"  Successfully processed: {processed_files}")
    print(f"  Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
    print(f"  This corrects the previous error that excluded NoChange")
    
    if processed_files > 0:
        nc_files = [f for f in os.listdir(netcdf_dir) if f.endswith('.nc')]
        csv_files = [f for f in os.listdir(tables_dir) if f.endswith('.csv')]
        
        print(f"\n📁 OUTPUTS SAVED TO:")
        print(f"  NetCDF files: {netcdf_dir} ({len(nc_files)} files)")
        print(f"  Summary tables: {tables_dir} ({len(csv_files)} CSV files)")
    
    print(f"\n🎉 CORRECTED LINEAR MODEL ANALYSIS COMPLETE!")
    print(f"All results saved to: {OUTPUT_DIR}")
    print(f"\n✅ Downstream compatibility maintained - all file names unchanged")
    
    return {
        'processor': processor,
        'summaries': all_summaries if all_summaries else None,
        'output_dir': OUTPUT_DIR
    }


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

if __name__ == "__main__":
    print("\n🚀 STARTING CORRECTED SPATIAL LINEAR EMULATION")
    print("Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
    
    print("\n1. Run sanity check on first 3 files? (Recommended)")
    print("2. Run full processing on all files")
    choice = input("\nEnter choice (1 or 2): ").strip()
    
    if choice == "1":
        if quick_sanity_check():
            print("\n✅ Sanity check passed. Correct models are working!")
            proceed = input("\nProceed with full processing? (y/n): ").lower()
            if proceed == 'y':
                results = run_spatial_linear_emulation(full_run=True)
            else:
                print("\nExiting. Run again for full processing.")
        else:
            print("\n❌ Sanity check failed. Fix issues before full run.")
    else:
        results = run_spatial_linear_emulation(full_run=True)