#!/usr/bin/env python3
"""
Q4_Ecological_Impacts.py

Exact replication of Malone's Ecological_Impacts.R workflow using Ammara's Python Q2 and Q3 outputs.

IMPORTANT PATH RULES:
- Read Q3 inputs from: M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs
- Write Q4 outputs to: M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_outputs
- Write Q4 figures to: M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_figures
- Write Q4 temp files to: M:\Research\WUE_CUE\WUE_manuscript_version6\Q4\Q4_Ecological_Impacts_temp
- NEVER write to or modify the Python Q3 output folder

Primary result: Q3-aligned upland coast-specific EDI using predicted WUE_T change from near-normal SPEI.
"""

import os
import sys
import warnings
import subprocess
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

# Suppress warnings
warnings.filterwarnings('ignore')

# ============================================================================
# DIRECTORY SETUP - Q4 only
# ============================================================================

# Q4 Base directory - all outputs go here
BASE_OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q4"

# Q4 output directories - write only to these
OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "Q4_Ecological_Impacts_outputs")
FIGURE_DIR = os.path.join(BASE_OUTPUT_DIR, "Q4_Ecological_Impacts_figures")
TEMP_DIR = os.path.join(BASE_OUTPUT_DIR, "Q4_Ecological_Impacts_temp")

# ============================================================================
# INPUT DIRECTORIES - read only from these
# ============================================================================

# Python Q3 outputs - read from here, NEVER write to here
PYTHON_Q3_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"

# Python Q2 outputs - read from here
PYTHON_Q2_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs"

# Monthly data
MONTHLY_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# ============================================================================
# CONSTANTS - matching Malone exactly
# ============================================================================

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]

# EDI severity labels (logistic model)
EDI_LABELS = ["None", "Watch", "Stress", "Impact"]
EDI_BREAKS = [0, 0.30, 0.50, 0.70, 1.0]

# Color maps
ECOSYSTEM_COLORS = {"Upland": "#800080", "Freshwater": "#0000FF", "Saline": "#FFA500"}
COAST_COLORS = {
    "Atlantic Coast": "#008B8B",
    "Pacific Coast": "#C44E52",
    "Gulf Coast": "#8C564B",
    "AK Coast": "#4D4D4D"
}
Q3_CLASS_COLORS = {
    "No meaningful change (<5%)": "#F2F2F2",
    "Watch (5-10%)": "#BFD3E6",
    "Stress (10-20%)": "#F4A582",
    "Impact (>=20%)": "#B2182B"
}
Q3_PI_COLORS = {
    "PI supports decrease": "#2166AC",
    "PI overlaps no-change": "#D9D9D9",
    "PI supports increase": "#B2182B"
}
SEVERITY_COLORS = {"None": "#2166AC", "Watch": "#92C5DE", "Stress": "#F4A582", "Impact": "#D6604D"}

# ============================================================================
# R Configuration
# ============================================================================

R_PATHS = [
    r"C:\Program Files\R\R-4.4.2\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.4.0\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.3\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.1\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.0\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.2.3\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.2.2\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.2.1\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.2.0\bin\x64\Rscript.exe",
]

RSCRIPT_EXE = None
for path in R_PATHS:
    if os.path.exists(path):
        RSCRIPT_EXE = path
        break

if RSCRIPT_EXE is None:
    RSCRIPT_EXE = shutil.which("Rscript")

R_AVAILABLE = RSCRIPT_EXE is not None and os.path.exists(RSCRIPT_EXE)

if R_AVAILABLE:
    print(f"R found at: {RSCRIPT_EXE}")
else:
    print("ERROR: R not found. Please install R and add Rscript.exe to PATH.")
    sys.exit(1)


def r_path(path):
    """Convert Windows path to R-safe forward slashes."""
    return os.path.abspath(path).replace("\\", "/")


# ============================================================================
# Helper Functions
# ============================================================================

def create_directories():
    """Create Q4 output directories only."""
    for d in [OUTPUT_DIR, FIGURE_DIR, TEMP_DIR]:
        os.makedirs(d, exist_ok=True)
    print(f"Q4 Output directory: {OUTPUT_DIR}")
    print(f"Q4 Figure directory: {FIGURE_DIR}")
    print(f"Q4 Temp directory: {TEMP_DIR}")
    print(f"Q3 Input directory: {PYTHON_Q3_DIR}")
    print(f"Q2 Input directory: {PYTHON_Q2_DIR}")


def clean_output_directories():
    """Clean Q4 output directories only - never touch Q3."""
    for d in [OUTPUT_DIR, FIGURE_DIR, TEMP_DIR]:
        if os.path.exists(d):
            for f in os.listdir(d):
                fpath = os.path.join(d, f)
                if os.path.isfile(fpath):
                    try:
                        os.remove(fpath)
                    except:
                        pass
                elif os.path.isdir(fpath):
                    shutil.rmtree(fpath)
    print("Cleaned Q4 output directories.")


def write_table(df, filename):
    """Write to Q4 output directory only."""
    if df is None or df.empty:
        print(f"Warning: Empty dataframe for {filename}")
        return
    outpath = os.path.join(OUTPUT_DIR, filename)
    df.to_csv(outpath, index=False)
    print(f"  Wrote Q4 output: {filename} ({len(df)} rows)")


def save_plot(fig, filename, width=11, height=7, dpi=300):
    """Save figure to Q4 figure directory only."""
    outpath = os.path.join(FIGURE_DIR, filename)
    fig.savefig(outpath, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved Q4 figure: {filename}")


def normalize_malone_column_names(df):
    """Rename columns to match Malone's exact column names for Q3 tables only."""
    if df is None:
        return None
    rename_dict = {
        'p-value': 'p.value',
        'Std. Error': 'Std..Error',
        't value': 't.value',
        'Pr(>|t|)': 'Pr...t..'
    }
    # Only rename columns that exist
    rename_dict = {k: v for k, v in rename_dict.items() if k in df.columns}
    if rename_dict:
        return df.rename(columns=rename_dict)
    return df


def logistic(x, b0, b1):
    return 1 / (1 + np.exp(-(b0 + b1 * x)))


def spei_at_p(p, b0, b1):
    return (np.log(p / (1 - p)) - b0) / b1


def classify_q3_response_magnitude(x):
    if pd.isna(x) or not np.isfinite(x):
        return "No meaningful change (<5%)"
    if x < 5:
        return "No meaningful change (<5%)"
    elif x < 10:
        return "Watch (5-10%)"
    elif x < 20:
        return "Stress (10-20%)"
    else:
        return "Impact (>=20%)"


def classify_q3_response_direction(x):
    if pd.isna(x) or not np.isfinite(x):
        return "No meaningful change"
    if abs(x) < 5:
        return "No meaningful change"
    elif x < 0:
        return "WUE_T decrease"
    else:
        return "WUE_T increase"


def classify_response_magnitude(x):
    if pd.isna(x) or not np.isfinite(x):
        return "Stable (<0.5%)"
    if x < 0.5:
        return "Stable (<0.5%)"
    elif x < 1:
        return "Very mild (0.5-1%)"
    elif x < 2.5:
        return "Mild (1-2.5%)"
    elif x < 5:
        return "Moderate (2.5-5%)"
    else:
        return "Strong (>=5%)"


def classify_response_direction(x):
    if pd.isna(x) or not np.isfinite(x):
        return "Stable"
    if abs(x) < 0.5:
        return "Stable"
    elif x < 0:
        return "WUE_T decrease"
    else:
        return "WUE_T increase"


def classify_coast_region(lat, long):
    if lat > 50:
        return "AK Coast"
    elif long > -100:
        return "Atlantic Coast"
    elif long < -120:
        return "Pacific Coast"
    else:
        return "Gulf Coast"


def theme_wue(ax, base_size=13):
    ax.set_facecolor('white')
    ax.grid(True, linestyle='-', alpha=0.15, color='gray')
    ax.grid(True, which='minor', linestyle='-', alpha=0.05)
    ax.tick_params(labelsize=base_size-1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    return ax


def find_first_threshold(group_data, threshold, side, direction='any'):
    if direction == 'any':
        response_filter = abs(group_data['predicted_pct_change']) >= threshold
    elif direction == 'increase':
        response_filter = group_data['predicted_pct_change'] >= threshold
    else:
        response_filter = group_data['predicted_pct_change'] <= -threshold
    
    if side == 'dry':
        side_filter = group_data['SPEI_3'] < -1
    else:
        side_filter = group_data['SPEI_3'] > 1
    
    threshold_values = group_data.loc[side_filter & response_filter, 'SPEI_3'].values
    
    if len(threshold_values) == 0:
        return np.nan
    
    if side == 'dry':
        return np.max(threshold_values)
    else:
        return np.min(threshold_values)


def find_spei_threshold(group_data, pct_column, threshold_pct, anomaly_side, impact_direction):
    pct_values = group_data[pct_column].values
    
    if anomaly_side == 'dry':
        side_filter = group_data['SPEI_value'] < -1
    else:
        side_filter = group_data['SPEI_value'] > 1
    
    if impact_direction == 'decrease':
        direction_filter = pct_values <= -threshold_pct
    else:
        direction_filter = pct_values >= threshold_pct
    
    threshold_values = group_data.loc[side_filter & direction_filter, 'SPEI_value'].values
    
    if len(threshold_values) == 0:
        return np.nan
    
    if anomaly_side == 'dry':
        return np.max(threshold_values)
    else:
        return np.min(threshold_values)


# ============================================================================
# R GLM Wrapper
# ============================================================================

def fit_glm_in_r(data, formula, family="binomial(link = 'logit')"):
    """Fit GLM using R glm."""
    
    if not R_AVAILABLE:
        raise RuntimeError("R not available")
    
    temp_r_script = tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False)
    temp_r_script_path = temp_r_script.name
    
    temp_data_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    temp_coef_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    temp_summary_path = tempfile.NamedTemporaryFile(suffix='.txt', delete=False).name
    
    temp_data_path_r = r_path(temp_data_path)
    temp_coef_path_r = r_path(temp_coef_path)
    temp_summary_path_r = r_path(temp_summary_path)
    
    data.to_csv(temp_data_path, index=False)
    
    r_script = f'''
    suppressPackageStartupMessages({{
        library(dplyr)
    }})
    
    df <- read.csv("{temp_data_path_r}", stringsAsFactors = FALSE)
    glm_formula <- as.formula("{formula}")
    glm_model <- glm(glm_formula, data = df, family = {family})
    
    sink("{temp_summary_path_r}")
    print(summary(glm_model))
    sink()
    
    coef_df <- data.frame(
        term = names(coef(glm_model)),
        estimate = coef(glm_model)
    )
    write.csv(coef_df, "{temp_coef_path_r}", row.names = FALSE)
    
    cat(paste("AIC:", AIC(glm_model), "\\n"))
    '''
    
    with open(temp_r_script_path, 'w') as f:
        f.write(r_script)
    
    try:
        result = subprocess.run([RSCRIPT_EXE, temp_r_script_path], capture_output=True, text=True, check=False)
        
        if result.returncode != 0:
            raise RuntimeError(f"R GLM failed: {result.stderr}")
        
        coef_df = pd.read_csv(temp_coef_path)
        aic = None
        for line in result.stdout.split('\n'):
            if 'AIC:' in line:
                aic = float(line.split(':')[1].strip())
        
        with open(temp_summary_path, 'r') as f:
            summary_text = f.read()
        
        return {'coef_df': coef_df, 'aic': aic, 'summary': summary_text}
        
    finally:
        for f in [temp_r_script_path, temp_data_path, temp_coef_path, temp_summary_path]:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except:
                pass


# ============================================================================
# R GAM Wrapper - Single fit, multiple predictions
# ============================================================================

def fit_gam_in_r(data, formula, method="REML", select=True, 
                  exclude_terms=None, predict_newdata_list=None,
                  rds_output_name=None):
    """
    Fit GAM using R mgcv once, predict on multiple datasets.
    
    Args:
        data: pandas DataFrame
        formula: R formula string
        method: "REML" or "ML"
        select: boolean
        exclude_terms: list of terms to exclude
        predict_newdata_list: list of (name, DataFrame) pairs
        rds_output_name: if provided, save RDS to Q4 output dir with this name
    """
    
    if not R_AVAILABLE:
        raise RuntimeError("R not available")
    
    temp_r_script = tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False)
    temp_r_script_path = temp_r_script.name
    
    temp_data_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    temp_rds_path = tempfile.NamedTemporaryFile(suffix='.rds', delete=False).name
    temp_smooth_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    temp_param_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    temp_summary_path = tempfile.NamedTemporaryFile(suffix='.txt', delete=False).name
    
    temp_data_path_r = r_path(temp_data_path)
    temp_rds_path_r = r_path(temp_rds_path)
    temp_smooth_path_r = r_path(temp_smooth_path)
    temp_param_path_r = r_path(temp_param_path)
    temp_summary_path_r = r_path(temp_summary_path)
    
    data.to_csv(temp_data_path, index=False)
    
    exclude_str = ""
    if exclude_terms:
        terms_str = ", ".join([f"'{t}'" for t in exclude_terms])
        exclude_str = f", exclude = c({terms_str})"
    
    factor_code = """
    if("water_class" %in% colnames(df)) {
        df$water_class <- factor(df$water_class, levels = c("Upland", "Freshwater", "Saline"))
    }
    if("month_f" %in% colnames(df)) {
        df$month_f <- factor(df$month_f, levels = 1:12)
    }
    if("site_name" %in% colnames(df)) {
        df$site_name <- factor(df$site_name)
    }
    if("SPEI_timescale" %in% colnames(df)) {
        df$SPEI_timescale <- factor(df$SPEI_timescale, levels = c("SPEI_1","SPEI_3","SPEI_6","SPEI_12","SPEI_24","SPEI_36","SPEI_48"))
    }
    if("coast_region" %in% colnames(df)) {
        df$coast_region <- factor(df$coast_region, levels = c("Atlantic Coast","Pacific Coast","Gulf Coast","AK Coast"))
    }
    """
    
    spei_coast_code = """
    if("SPEI_timescale" %in% colnames(df) && "coast_region" %in% colnames(df)) {
        df$spei_coast <- interaction(df$SPEI_timescale, df$coast_region, sep = "__", drop = TRUE)
    }
    """
    
    pred_code = ""
    pred_paths = {}
    
    if predict_newdata_list:
        for name, newdata_df in predict_newdata_list:
            temp_newdata_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
            temp_pred_path = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
            
            temp_newdata_path_r = r_path(temp_newdata_path)
            temp_pred_path_r = r_path(temp_pred_path)
            
            newdata_df.to_csv(temp_newdata_path, index=False)
            pred_paths[name] = {'pred_path': temp_pred_path, 'newdata_path': temp_newdata_path}
            
            pred_code += f'''
            newdata_{name} <- read.csv("{temp_newdata_path_r}", stringsAsFactors = FALSE)
            if("water_class" %in% colnames(newdata_{name})) {{
                newdata_{name}$water_class <- factor(newdata_{name}$water_class, levels = levels(df$water_class))
            }}
            if("month_f" %in% colnames(newdata_{name})) {{
                newdata_{name}$month_f <- factor(newdata_{name}$month_f, levels = levels(df$month_f))
            }}
            if("site_name" %in% colnames(newdata_{name})) {{
                newdata_{name}$site_name <- factor(newdata_{name}$site_name, levels = levels(df$site_name))
            }}
            if("SPEI_timescale" %in% colnames(newdata_{name})) {{
                newdata_{name}$SPEI_timescale <- factor(newdata_{name}$SPEI_timescale, levels = levels(df$SPEI_timescale))
            }}
            if("coast_region" %in% colnames(newdata_{name})) {{
                newdata_{name}$coast_region <- factor(newdata_{name}$coast_region, levels = levels(df$coast_region))
            }}
            if("SPEI_timescale" %in% colnames(newdata_{name}) && "coast_region" %in% colnames(newdata_{name})) {{
                newdata_{name}$spei_coast <- interaction(newdata_{name}$SPEI_timescale, newdata_{name}$coast_region, sep = "__", drop = TRUE)
                newdata_{name}$spei_coast <- factor(newdata_{name}$spei_coast, levels = levels(df$spei_coast))
            }}
            
            pred_obj_{name} <- predict(gam_model, newdata = newdata_{name}, type = "link", se.fit = TRUE{exclude_str})
            newdata_{name}$predicted <- as.numeric(pred_obj_{name}$fit)
            newdata_{name}$se_fit <- as.numeric(pred_obj_{name}$se.fit)
            newdata_{name}$predicted_lower <- newdata_{name}$predicted - 1.96 * newdata_{name}$se_fit
            newdata_{name}$predicted_upper <- newdata_{name}$predicted + 1.96 * newdata_{name}$se_fit
            
            write.csv(newdata_{name}, "{temp_pred_path_r}", row.names = FALSE)
            '''
    
    r_script = f'''
    suppressPackageStartupMessages({{
        library(mgcv)
        library(dplyr)
    }})
    
    df <- read.csv("{temp_data_path_r}", stringsAsFactors = FALSE)
    {factor_code}
    {spei_coast_code}
    
    gam_formula <- as.formula("{formula}")
    gam_model <- gam(gam_formula, data = df, method = "{method}", select = {str(select).upper()})
    
    saveRDS(gam_model, "{temp_rds_path_r}")
    
    sink("{temp_summary_path_r}")
    print(summary(gam_model))
    sink()
    
    smooth_table <- as.data.frame(summary(gam_model)$s.table)
    smooth_table$smooth_term <- rownames(smooth_table)
    rownames(smooth_table) <- NULL
    smooth_table <- smooth_table[, c("smooth_term", setdiff(names(smooth_table), "smooth_term"))]
    write.csv(smooth_table, "{temp_smooth_path_r}", row.names = FALSE)
    
    param_table <- as.data.frame(summary(gam_model)$p.table)
    param_table$term <- rownames(param_table)
    rownames(param_table) <- NULL
    param_table <- param_table[, c("term", setdiff(names(param_table), "term"))]
    write.csv(param_table, "{temp_param_path_r}", row.names = FALSE)
    
    {pred_code}
    
    cat(paste("AIC:", AIC(gam_model), "\\n"))
    cat(paste("BIC:", BIC(gam_model), "\\n"))
    cat(paste("logLik:", as.numeric(logLik(gam_model)), "\\n"))
    cat(paste("dev_expl:", summary(gam_model)$dev.expl, "\\n"))
    cat(paste("r_sq:", summary(gam_model)$r.sq, "\\n"))
    '''
    
    with open(temp_r_script_path, 'w') as f:
        f.write(r_script)
    
    try:
        result = subprocess.run([RSCRIPT_EXE, temp_r_script_path], capture_output=True, text=True, check=False)
        
        if result.returncode != 0:
            raise RuntimeError(f"R GAM failed: {result.stderr}")
        
        stdout_lines = result.stdout.strip().split('\n')
        model_info = {}
        for line in stdout_lines:
            if ':' in line and not line.startswith(' '):
                key, value = line.split(':', 1)
                model_info[key.strip()] = value.strip()
        
        smooth_df = pd.read_csv(temp_smooth_path) if os.path.exists(temp_smooth_path) else None
        param_df = pd.read_csv(temp_param_path) if os.path.exists(temp_param_path) else None
        
        with open(temp_summary_path, 'r') as f:
            summary_text = f.read()
        
        predictions = {}
        for name, paths in pred_paths.items():
            if os.path.exists(paths['pred_path']):
                predictions[name] = pd.read_csv(paths['pred_path'])
        
        # Save RDS to Q4 output directory if requested
        if rds_output_name and os.path.exists(temp_rds_path):
            rds_output_path = os.path.join(OUTPUT_DIR, rds_output_name)
            shutil.copy2(temp_rds_path, rds_output_path)
            print(f"  Saved RDS model to Q4: {rds_output_name}")
        
        return {
            'predictions': predictions,
            'smooth_df': smooth_df,
            'param_df': param_df,
            'summary': summary_text,
            'model_info': model_info
        }
        
    finally:
        temp_files = [temp_r_script_path, temp_data_path, temp_rds_path, 
                      temp_smooth_path, temp_param_path, temp_summary_path]
        for f in temp_files:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except:
                pass
        for name, paths in pred_paths.items():
            try:
                if os.path.exists(paths['pred_path']):
                    os.remove(paths['pred_path'])
                if os.path.exists(paths['newdata_path']):
                    os.remove(paths['newdata_path'])
            except:
                pass


# ============================================================================
# Step 1: Load Data
# ============================================================================

def load_data():
    print("\n" + "="*60)
    print("Step 1: Loading Data")
    print("="*60)
    print(f"  Reading monthly data from: {MONTHLY_FILE}")
    print(f"  Reading Q2 from: {PYTHON_Q2_DIR}")
    print(f"  Reading Q3 from: {PYTHON_Q3_DIR}")
    
    monthly = pd.read_csv(MONTHLY_FILE)
    for col in monthly.select_dtypes(include=['object']).columns:
        monthly[col] = monthly[col].str.strip()
    print(f"  Monthly data: {len(monthly)} rows, {len(monthly['site_name'].unique())} sites")
    
    resistance = pd.read_csv(os.path.join(PYTHON_Q2_DIR, "Q2_site_resistance.csv"))
    recovery = pd.read_csv(os.path.join(PYTHON_Q2_DIR, "Q2_site_recovery.csv"))
    print(f"  Resistance: {len(resistance)} sites, Recovery: {len(recovery)} sites")
    
    # Read Q3 outputs from Python Q3 folder
    q3_pred = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv"))
    q3_threshold = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv"))
    print(f"  Q3 predicted: {len(q3_pred)} rows, Q3 threshold: {len(q3_threshold)} rows")
    
    q3_smooth = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")) if os.path.exists(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")) else None
    q3_param = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv")) if os.path.exists(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv")) else None
    q3_model_comp = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_model_comparison.csv")) if os.path.exists(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_model_comparison.csv")) else None
    q3_site_sens = pd.read_csv(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv")) if os.path.exists(os.path.join(PYTHON_Q3_DIR, "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv")) else None
    
    if q3_smooth is not None:
        print(f"  Q3 smooth terms: {len(q3_smooth)} rows")
    if q3_param is not None:
        print(f"  Q3 parametric terms: {len(q3_param)} rows")
    if q3_model_comp is not None:
        print(f"  Q3 model comparison: {len(q3_model_comp)} rows")
    if q3_site_sens is not None:
        print(f"  Q3 site sensitivity: {len(q3_site_sens)} rows")
    
    return {'monthly': monthly, 'resistance': resistance, 'recovery': recovery,
            'q3_pred': q3_pred, 'q3_threshold': q3_threshold,
            'q3_smooth': q3_smooth, 'q3_param': q3_param,
            'q3_model_comp': q3_model_comp, 'q3_site_sens': q3_site_sens}


# ============================================================================
# Step 2: Prepare Base Data
# ============================================================================

def prepare_base_data(monthly):
    print("\n" + "="*60)
    print("Step 2: Preparing Base Data")
    print("="*60)
    
    required_cols = ["site_name", "Year", "month", "water_class", "WUE_tra", "SPEI_3"]
    base_df = monthly[required_cols + ["lat", "long"]].copy()
    base_df = base_df.dropna(subset=required_cols)
    base_df = base_df[(base_df['site_name'] != "")]
    base_df = base_df[base_df['water_class'].isin(ECOSYSTEM_CLASSES)]
    
    base_df['WUE_T'] = base_df['WUE_tra']
    base_df['water_class'] = pd.Categorical(base_df['water_class'], categories=ECOSYSTEM_CLASSES)
    base_df['month_f'] = pd.Categorical(base_df['month'], categories=range(1, 13))
    base_df = base_df.sort_values(['site_name', 'Year', 'month'])
    
    print(f"  Base data: {len(base_df)} rows, {len(base_df['site_name'].unique())} sites")
    return base_df


# ============================================================================
# Step 3: Compute Baseline
# ============================================================================

def compute_baseline(base_df):
    print("\n" + "="*60)
    print("Step 3: Computing Near-Normal Baseline")
    print("="*60)
    
    near_normal = base_df[(base_df['SPEI_3'] >= -1) & (base_df['SPEI_3'] <= 1)].copy()
    baselines = near_normal.groupby(['site_name', 'month_f']).agg(
        baseline_WUE_T=('WUE_T', 'mean'),
        baseline_sd_WUE_T=('WUE_T', 'std'),
        n_baseline_months=('WUE_T', 'count')
    ).reset_index()
    baselines = baselines[baselines['n_baseline_months'] >= 3]
    
    print(f"  Baselines: {len(baselines)} site-month combos, {len(baselines['site_name'].unique())} sites")
    
    base_df = base_df.merge(baselines, on=['site_name', 'month_f'], how='left')
    base_df['monthly_resistance'] = base_df['WUE_T'] / base_df['baseline_WUE_T']
    base_df['WUE_T_pct_change'] = 100 * (base_df['WUE_T'] - base_df['baseline_WUE_T']) / base_df['baseline_WUE_T']
    
    model_ready = base_df[base_df['baseline_WUE_T'].notna()].copy()
    print(f"  Rows with valid baseline: {len(model_ready)}")
    
    return base_df, model_ready


# ============================================================================
# Step 4: Logistic EDI (R glm)
# ============================================================================

def fit_logistic_models_r(model_ready):
    print("\n" + "="*60)
    print("Step 4: Fitting Diagnostic Logistic EDI Models (R glm)")
    print("="*60)
    
    model_df = model_ready[
        (model_ready['monthly_resistance'].notna()) &
        (np.isfinite(model_ready['monthly_resistance'])) &
        (model_ready['SPEI_3'] < 0)
    ].copy()
    
    if len(model_df) == 0:
        raise RuntimeError(f"No logistic calibration rows. model_ready has {len(model_ready)} rows.")
    
    model_df['y'] = (model_df['monthly_resistance'] < 1).astype(int)
    print(f"  Logistic data: {len(model_df)} rows, impacted: {model_df['y'].sum()}")
    
    # Pooled model
    print("  Fitting pooled logistic model in R...")
    pooled_result = fit_glm_in_r(
        data=model_df[['SPEI_3', 'y']],
        formula="y ~ SPEI_3",
        family="binomial(link = 'logit')"
    )
    
    b0_pooled = pooled_result['coef_df'][pooled_result['coef_df']['term'] == '(Intercept)']['estimate'].values[0]
    b1_pooled = pooled_result['coef_df'][pooled_result['coef_df']['term'] == 'SPEI_3']['estimate'].values[0]
    aic_pooled = pooled_result['aic']
    
    print(f"  Pooled: intercept={b0_pooled:.4f}, slope={b1_pooled:.4f}, AIC={aic_pooled:.1f}")
    
    # Ecosystem-specific models
    eco_results = {}
    eco_coefs = {}
    
    for ec in ECOSYSTEM_CLASSES:
        sub_df = model_df[model_df['water_class'] == ec].copy()
        n_impacted = sub_df['y'].sum()
        n_non_impacted = len(sub_df) - n_impacted
        
        if (n_impacted < 5) or (n_non_impacted < 5):
            print(f"  {ec}: skipped (insufficient events)")
            eco_results[ec] = None
            eco_coefs[ec] = None
            continue
        
        try:
            print(f"  Fitting {ec} logistic model in R...")
            eco_result = fit_glm_in_r(
                data=sub_df[['SPEI_3', 'y']],
                formula="y ~ SPEI_3",
                family="binomial(link = 'logit')"
            )
            eco_results[ec] = eco_result
            b0_ec = eco_result['coef_df'][eco_result['coef_df']['term'] == '(Intercept)']['estimate'].values[0]
            b1_ec = eco_result['coef_df'][eco_result['coef_df']['term'] == 'SPEI_3']['estimate'].values[0]
            eco_coefs[ec] = (b0_ec, b1_ec)
            print(f"  {ec}: intercept={b0_ec:.4f}, slope={b1_ec:.4f}, AIC={eco_result['aic']:.1f}")
        except Exception as e:
            print(f"  {ec}: failed: {e}")
            eco_results[ec] = None
            eco_coefs[ec] = None
    
    return {'pooled': {'b0': b0_pooled, 'b1': b1_pooled, 'aic': aic_pooled, 'result': pooled_result},
            'eco_results': eco_results, 'eco_coefs': eco_coefs, 'model_df': model_df}


def apply_logistic_edi(base_df, logistic_results):
    print("\n" + "="*60)
    print("Applying Logistic EDI")
    print("="*60)
    
    b0, b1 = logistic_results['pooled']['b0'], logistic_results['pooled']['b1']
    base_df['EDI'] = logistic(base_df['SPEI_3'], b0, b1)
    
    for ec in ECOSYSTEM_CLASSES:
        if logistic_results['eco_coefs'].get(ec):
            b0_ec, b1_ec = logistic_results['eco_coefs'][ec]
            idx = base_df['water_class'] == ec
            base_df.loc[idx, 'EDI'] = logistic(base_df.loc[idx, 'SPEI_3'], b0_ec, b1_ec)
    
    base_df['EDI_class'] = pd.cut(base_df['EDI'], bins=EDI_BREAKS, labels=EDI_LABELS, include_lowest=True)
    
    print(f"  EDI applied to {len(base_df)} rows")
    for cls in EDI_LABELS:
        count = (base_df['EDI_class'] == cls).sum()
        print(f"    {cls}: {count} ({count/len(base_df)*100:.1f}%)")
    
    return base_df


# ============================================================================
# Step 5: Response EDI GAM
# ============================================================================

def fit_response_edi_gam(model_ready):
    print("\n" + "="*60)
    print("Step 5: Fitting Response-Based EDI GAM")
    print("="*60)
    
    # Ensure row_id is present (not __row_id - R mangles column names with underscores)
    if 'row_id' not in model_ready.columns:
        raise RuntimeError("model_ready is missing row_id. Add row ID before calling fit_response_edi_gam.")
    
    response_df = model_ready[
        (np.isfinite(model_ready['WUE_T_pct_change'])) &
        (np.isfinite(model_ready['SPEI_3'])) &
        (model_ready['baseline_WUE_T'].notna()) &
        (model_ready['n_baseline_months'] >= 3)
    ].copy()
    
    if len(response_df) == 0:
        raise RuntimeError(f"No response model rows. model_ready has {len(model_ready)} rows.")
    
    print(f"  Response data: {len(response_df)} rows, {len(response_df['site_name'].unique())} sites")
    
    response_df['month_f'] = response_df['month'].astype(int)
    response_df['site_name'] = response_df['site_name'].astype(str)
    response_df['water_class'] = response_df['water_class'].astype(str)
    
    formula = ("WUE_T_pct_change ~ water_class + month_f + "
               "s(SPEI_3, by = water_class, k = 6) + "
               "s(site_name, bs = 're')")
    
    # Prediction grid
    spei_seq = np.linspace(response_df['SPEI_3'].quantile(0.02), response_df['SPEI_3'].quantile(0.98), 250)
    grid_dfs = []
    for ec in ECOSYSTEM_CLASSES:
        grid_df = pd.DataFrame({
            'water_class': [ec] * len(spei_seq),
            'SPEI_3': spei_seq,
            'month_f': [7] * len(spei_seq),
            'site_name': [response_df['site_name'].iloc[0]] * len(spei_seq)
        })
        grid_dfs.append(grid_df)
    grid_data = pd.concat(grid_dfs, ignore_index=True)
    grid_data['month_f'] = grid_data['month_f'].astype(int)
    grid_data['site_name'] = grid_data['site_name'].astype(str)
    grid_data['water_class'] = grid_data['water_class'].astype(str)
    
    # All rows for prediction - include row_id for safe merge-back (not __row_id)
    all_pred_df = model_ready.copy()
    all_pred_df['month_f'] = all_pred_df['month'].astype(int)
    all_pred_df['site_name'] = all_pred_df['site_name'].astype(str)
    all_pred_df['water_class'] = all_pred_df['water_class'].astype(str)
    all_pred_subset = all_pred_df[
        ['row_id', 'site_name', 'Year', 'month', 'water_class', 'month_f', 'SPEI_3']
    ].copy()
    
    print("  Fitting response GAM in R...")
    result = fit_gam_in_r(
        data=response_df,
        formula=formula,
        method="REML",
        select=True,
        exclude_terms=["s(site_name)"],
        predict_newdata_list=[
            ('grid', grid_data),
            ('all_rows', all_pred_subset)
        ],
        rds_output_name=None  # Don't save RDS for response GAM
    )
    
    grid_pred = result['predictions'].get('grid')
    all_pred = result['predictions'].get('all_rows')
    
    if grid_pred is not None:
        # Remove se_fit column if present (Malone's output does not have it)
        if 'se_fit' in grid_pred.columns:
            grid_pred = grid_pred.drop(columns=['se_fit'])
        
        grid_pred = grid_pred.rename(columns={
            'predicted': 'predicted_pct_change',
            'predicted_lower': 'predicted_pct_change_lower',
            'predicted_upper': 'predicted_pct_change_upper'
        })
        grid_pred['predicted_magnitude'] = grid_pred['predicted_pct_change'].abs()
        grid_pred['response_class'] = grid_pred['predicted_magnitude'].apply(classify_response_magnitude)
        grid_pred['response_direction'] = grid_pred['predicted_pct_change'].apply(classify_response_direction)
    
    return {
        'response_df': response_df,
        'grid_pred': grid_pred,
        'all_pred': all_pred,
        'smooth_df': result['smooth_df'],
        'param_df': result['param_df'],
        'summary': result['summary']
    }


# ============================================================================
# Step 6: Q3-Aligned EDI
# ============================================================================

def prepare_q3_aligned_data(monthly, q3_pred, q3_threshold, q3_smooth, q3_param, q3_model_comp, q3_site_sens):
    print("\n" + "="*60)
    print("Step 6: Primary Q3-Aligned EDI")
    print("="*60)
    
    # --- 6a: Q3-aligned prediction scores ---
    # This comes from Python Q3 output, only add Q4 columns
    q3_pred_edi = q3_pred.copy()
    q3_pred_edi['EDI_Q3_magnitude'] = q3_pred_edi['predicted_pct_change'].abs()
    q3_pred_edi['EDI_Q3_class'] = q3_pred_edi['predicted_pct_change'].apply(
        lambda x: classify_q3_response_magnitude(abs(x))
    )
    q3_pred_edi['EDI_Q3_direction'] = q3_pred_edi['predicted_pct_change'].apply(
        lambda x: classify_q3_response_direction(x)
    )
    q3_pred_edi['prediction_interval_support'] = np.where(
        q3_pred_edi['predicted_upper_pct'] <= -5, "PI supports decrease",
        np.where(q3_pred_edi['predicted_lower_pct'] >= 5, "PI supports increase", "PI overlaps no-change")
    )
    print(f"  Q3-aligned prediction scores: {len(q3_pred_edi)} rows")
    
    # --- 6b: Q3 threshold summary ---
    q3_threshold_edi = q3_threshold.copy()
    q3_threshold_edi = q3_threshold_edi.sort_values(
        ['SPEI_timescale', 'coast_region', 'anomaly_side', 'impact_direction', 'threshold_pct']
    )
    print(f"  Q3 threshold summary: {len(q3_threshold_edi)} rows")
    
    # --- 6c: Q3 model info ---
    q3_model_info = q3_model_comp[q3_model_comp['model'] == 'coast_threshold_gam'].copy() if q3_model_comp is not None else None
    
    # --- 6d: Q3 monthly scores ---
    # This requires fitting a Q3-like GAM in R for diagnostic monthly scores
    print("  Fitting Q3 coast GAM in R for monthly scores...")
    
    q3_required_cols = ["site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra"] + SPEI_COLS
    q3_base = monthly[q3_required_cols].copy()
    q3_base = q3_base.dropna(subset=q3_required_cols)
    q3_base = q3_base[(q3_base['site_name'] != "")]
    q3_base = q3_base[(q3_base['water_class'] != "")]
    q3_base = q3_base[q3_base['water_class'].isin(ECOSYSTEM_CLASSES)]
    q3_base = q3_base[(np.isfinite(q3_base['lat'])) & (np.isfinite(q3_base['long']))]
    q3_base = q3_base[np.isfinite(q3_base['WUE_tra'])]
    
    q3_base['coast_region'] = q3_base.apply(lambda row: classify_coast_region(row['lat'], row['long']), axis=1)
    q3_base['month_f'] = q3_base['month'].astype(int)
    q3_base['site_name'] = q3_base['site_name'].astype(str)
    q3_base['WUE_T'] = q3_base['WUE_tra']
    
    # Reshape to long
    q3_long = pd.melt(
        q3_base,
        id_vars=['site_name', 'Year', 'month', 'water_class', 'coast_region', 'lat', 'long', 'WUE_T', 'month_f'],
        value_vars=SPEI_COLS,
        var_name='SPEI_timescale',
        value_name='SPEI_value'
    )
    q3_long = q3_long.dropna(subset=['SPEI_value', 'WUE_T'])
    print(f"  Q3 long data: {len(q3_long)} rows")
    
    formula = ("WUE_T ~ SPEI_timescale * coast_region + water_class + month_f + "
               "s(SPEI_value, by = spei_coast, k = 6) + "
               "s(site_name, bs = 're')")
    
    # Prediction grid for near-normal
    spei_seq = np.linspace(q3_long['SPEI_value'].quantile(0.02), q3_long['SPEI_value'].quantile(0.98), 120)
    grid_dfs = []
    for ts in SPEI_COLS:
        for coast in COAST_REGION_LEVELS:
            grid_df = pd.DataFrame({
                'SPEI_timescale': [ts] * len(spei_seq),
                'coast_region': [coast] * len(spei_seq),
                'SPEI_value': spei_seq,
                'water_class': ['Upland'] * len(spei_seq),
                'month_f': [7] * len(spei_seq),
                'site_name': [q3_long['site_name'].iloc[0]] * len(spei_seq)
            })
            grid_dfs.append(grid_df)
    grid_data = pd.concat(grid_dfs, ignore_index=True)
    grid_data['month_f'] = grid_data['month_f'].astype(int)
    grid_data['site_name'] = grid_data['site_name'].astype(str)
    grid_data['water_class'] = grid_data['water_class'].astype(str)
    
    # Monthly upland predictions
    q3_monthly_upland = q3_base[q3_base['water_class'] == "Upland"].copy()
    q3_monthly_upland['SPEI_value'] = q3_monthly_upland['SPEI_3']
    q3_monthly_upland['SPEI_timescale'] = 'SPEI_3'
    pred_data = q3_monthly_upland[['site_name', 'Year', 'month', 'water_class', 'coast_region', 
                                    'lat', 'long', 'WUE_T', 'month_f', 'SPEI_value', 
                                    'SPEI_timescale']].copy()
    
    # Fit GAM - save RDS for Q3 coast GAM
    print("  Fitting Q3 coast GAM in R...")
    q3_result = fit_gam_in_r(
        data=q3_long,
        formula=formula,
        method="ML",
        select=True,
        exclude_terms=["s(site_name)"],
        predict_newdata_list=[
            ('near_normal_grid', grid_data),
            ('monthly_upland', pred_data)
        ],
        rds_output_name="EDI_Q3_aligned_coast_threshold_gam_model.rds"  # Save RDS for Q3 coast GAM
    )
    
    q3_near_normal_pred = q3_result['predictions'].get('near_normal_grid')
    q3_monthly_scores = q3_result['predictions'].get('monthly_upland')
    
    q3_near_normal = None
    if q3_near_normal_pred is not None:
        q3_near_normal = q3_near_normal_pred[
            (q3_near_normal_pred['SPEI_value'] >= -1) & (q3_near_normal_pred['SPEI_value'] <= 1)
        ].groupby(['SPEI_timescale', 'coast_region']).agg(
            predicted_near_normal_WUE_T=('predicted', 'mean')
        ).reset_index()
        print(f"  Near-normal WUE_T computed")
    
    if q3_monthly_scores is not None and q3_near_normal is not None:
        q3_nn_sp3 = q3_near_normal[q3_near_normal['SPEI_timescale'] == 'SPEI_3']
        q3_monthly_scores = q3_monthly_scores.merge(
            q3_nn_sp3[['coast_region', 'predicted_near_normal_WUE_T']],
            on='coast_region', how='left'
        )
        
        q3_monthly_scores['EDI_Q3_pct_change'] = 100 * (
            q3_monthly_scores['predicted'] - q3_monthly_scores['predicted_near_normal_WUE_T']
        ) / q3_monthly_scores['predicted_near_normal_WUE_T']
        q3_monthly_scores['EDI_Q3_pct_change_lower'] = 100 * (
            q3_monthly_scores['predicted_lower'] - q3_monthly_scores['predicted_near_normal_WUE_T']
        ) / q3_monthly_scores['predicted_near_normal_WUE_T']
        q3_monthly_scores['EDI_Q3_pct_change_upper'] = 100 * (
            q3_monthly_scores['predicted_upper'] - q3_monthly_scores['predicted_near_normal_WUE_T']
        ) / q3_monthly_scores['predicted_near_normal_WUE_T']
        
        q3_monthly_scores['EDI_Q3_magnitude'] = q3_monthly_scores['EDI_Q3_pct_change'].abs()
        q3_monthly_scores['EDI_Q3_class'] = q3_monthly_scores['EDI_Q3_pct_change'].apply(
            lambda x: classify_q3_response_magnitude(abs(x))
        )
        q3_monthly_scores['EDI_Q3_direction'] = q3_monthly_scores['EDI_Q3_pct_change'].apply(
            lambda x: classify_q3_response_direction(x)
        )
        q3_monthly_scores['prediction_interval_support'] = np.where(
            q3_monthly_scores['EDI_Q3_pct_change_upper'] <= -5, "PI supports decrease",
            np.where(q3_monthly_scores['EDI_Q3_pct_change_lower'] >= 5, "PI supports increase", "PI overlaps no-change")
        )
        
        output_cols = ['site_name', 'Year', 'month', 'water_class', 'coast_region', 'lat', 'long',
                       'WUE_T', 'SPEI_value', 'predicted_near_normal_WUE_T', 'predicted',
                       'EDI_Q3_pct_change', 'EDI_Q3_pct_change_lower', 'EDI_Q3_pct_change_upper',
                       'EDI_Q3_magnitude', 'EDI_Q3_direction', 'EDI_Q3_class', 'prediction_interval_support']
        q3_monthly_out = q3_monthly_scores[output_cols].copy()
        q3_monthly_out = q3_monthly_out.rename(columns={'SPEI_value': 'SPEI_3', 'predicted': 'predicted_WUE_T'})
        print(f"  Q3 monthly scores: {len(q3_monthly_out)} rows")
    else:
        q3_monthly_out = pd.DataFrame()
        print("  WARNING: Q3 monthly scores could not be generated")
    
    # --- Regional summary ---
    if not q3_monthly_out.empty:
        q3_regional_summary = q3_monthly_out.groupby(
            ['coast_region', 'EDI_Q3_class', 'EDI_Q3_direction', 'prediction_interval_support']
        ).agg(
            n_site_months=('site_name', 'count'),
            n_sites=('site_name', 'nunique'),
            mean_SPEI_3=('SPEI_3', 'mean'),
            mean_pct_change=('EDI_Q3_pct_change', 'mean'),
            median_abs_pct_change=('EDI_Q3_magnitude', 'median')
        ).reset_index()
        q3_regional_summary['prop_site_months'] = q3_regional_summary.groupby('coast_region')['n_site_months'].transform(
            lambda x: x / x.sum() * 100
        )
        
        q3_regional_validation = q3_monthly_out.groupby('coast_region').agg(
            n_site_months=('site_name', 'count'),
            n_sites=('site_name', 'nunique'),
            mean_abs_Q3_EDI_pct_change=('EDI_Q3_magnitude', 'mean'),
            median_abs_Q3_EDI_pct_change=('EDI_Q3_magnitude', 'median'),
            pct_months_5pct_or_more=('EDI_Q3_magnitude', lambda x: (x >= 5).mean() * 100 if len(x) > 0 else 0),
            pct_months_10pct_or_more=('EDI_Q3_magnitude', lambda x: (x >= 10).mean() * 100 if len(x) > 0 else 0),
            pct_months_PI_supported=('prediction_interval_support', 
                                      lambda x: (x != "PI overlaps no-change").mean() * 100 if len(x) > 0 else 0)
        ).reset_index()
        
        if q3_site_sens is not None:
            q3_site_sens_subset = q3_site_sens[q3_site_sens['SPEI_timescale'] == 'SPEI_3']
            q3_regional_validation = q3_regional_validation.merge(
                q3_site_sens_subset[['coast_region', 'n_sites', 'median_abs_pct_change',
                                     'q25_abs_pct_change', 'q75_abs_pct_change', 'mean_abs_pct_change']],
                on='coast_region', how='left', suffixes=('_EDI_monthly', '_Q3_site')
            )
        print(f"  Q3 regional validation: {len(q3_regional_validation)} rows")
    else:
        q3_regional_summary = pd.DataFrame()
        q3_regional_validation = pd.DataFrame()
    
    return {
        'q3_pred_edi': q3_pred_edi,
        'q3_threshold_edi': q3_threshold_edi,
        'q3_model_info': q3_model_info,
        'q3_smooth': q3_smooth,  # From Python Q3 source
        'q3_param': q3_param,    # From Python Q3 source
        'q3_monthly_out': q3_monthly_out,
        'q3_regional_summary': q3_regional_summary,
        'q3_regional_validation': q3_regional_validation,
        'q3_model_summary': q3_result['summary'],
        'q3_model_info_dict': q3_result['model_info'],
        'q3_near_normal': q3_near_normal
    }


# ============================================================================
# Step 7: Compute Response Thresholds
# ============================================================================

def compute_response_thresholds(grid_pred):
    if grid_pred is None or grid_pred.empty:
        return None
    
    thresholds_list = []
    
    for ec in ECOSYSTEM_CLASSES:
        group_data = grid_pred[grid_pred['water_class'] == ec].copy()
        if group_data.empty:
            continue
        
        for threshold in [0.5, 1, 2.5, 5]:
            pred_mag = group_data['predicted_magnitude']
            stable_values = group_data.loc[pred_mag < threshold, 'SPEI_3'].values
            
            row = {
                'water_class': ec,
                'response_threshold_pct': threshold,
                'stable_SPEI_min': stable_values.min() if len(stable_values) > 0 else np.nan,
                'stable_SPEI_max': stable_values.max() if len(stable_values) > 0 else np.nan,
                'dry_any_impact_threshold': find_first_threshold(group_data, threshold, 'dry', 'any'),
                'dry_decrease_threshold': find_first_threshold(group_data, threshold, 'dry', 'decrease'),
                'dry_increase_threshold': find_first_threshold(group_data, threshold, 'dry', 'increase'),
                'wet_any_impact_threshold': find_first_threshold(group_data, threshold, 'wet', 'any'),
                'wet_decrease_threshold': find_first_threshold(group_data, threshold, 'wet', 'decrease'),
                'wet_increase_threshold': find_first_threshold(group_data, threshold, 'wet', 'increase'),
                'min_predicted_pct_change': group_data['predicted_pct_change'].min(),
                'max_predicted_pct_change': group_data['predicted_pct_change'].max(),
                'max_predicted_magnitude': group_data['predicted_magnitude'].max()
            }
            thresholds_list.append(row)
    
    return pd.DataFrame(thresholds_list) if thresholds_list else None


# ============================================================================
# Step 8: Write Outputs
# ============================================================================

def write_all_outputs(data):
    print("\n" + "="*60)
    print("Step 8: Writing Q4 Output Tables")
    print("="*60)
    
    # Q3-aligned outputs - use Q3 smooth/param from Python Q3 source
    write_table(data.get('q3_pred_edi'), "EDI_Q3_aligned_prediction_scores_upland.csv")
    write_table(data.get('q3_threshold_edi'), "EDI_Q3_aligned_coast_thresholds_5_10_20pct.csv")
    write_table(data.get('q3_model_info'), "EDI_Q3_aligned_model_info.csv")
    
    # Normalize Q3 smooth/param column names to match Malone's exact format
    write_table(
        normalize_malone_column_names(data.get('q3_smooth')),
        "EDI_Q3_aligned_coast_GAM_smooth_terms.csv"
    )
    write_table(
        normalize_malone_column_names(data.get('q3_param')),
        "EDI_Q3_aligned_coast_GAM_parametric_terms.csv"
    )
    
    write_table(data.get('q3_monthly_out'), "EDI_Q3_aligned_monthly_scores.csv")
    write_table(data.get('q3_regional_summary'), "EDI_Q3_aligned_monthly_regional_summary.csv")
    write_table(data.get('q3_regional_validation'), "EDI_Q3_aligned_regional_validation.csv")
    
    # Logistic diagnostic outputs
    write_table(data.get('logistic_coefs'), "EDI_logistic_coefficients.csv")
    write_table(data.get('thresholds'), "EDI_threshold_summary.csv")
    write_table(data.get('monthly_out'), "EDI_monthly_site_scores.csv")
    write_table(data.get('validation_df'), "EDI_validation_by_class.csv")
    
    # Response GAM outputs - DO NOT normalize column names
    # These already match Malone exactly from the R output
    write_table(data.get('response_grid'), "EDI_response_prediction_curves.csv")
    write_table(data.get('response_thresholds'), "EDI_response_threshold_summary.csv")
    write_table(data.get('response_smooth'), "EDI_response_gam_smooth_terms.csv")
    write_table(data.get('response_param'), "EDI_response_gam_parametric_terms.csv")


# ============================================================================
# Step 9: Write Summary Files
# ============================================================================

def write_q3_summary(data):
    print("\n" + "="*60)
    print("Step 9: Writing Q3-Aligned Summary")
    print("="*60)
    
    q3_pred = data.get('q3_pred_edi')
    q3_threshold = data.get('q3_threshold_edi')
    q3_monthly = data.get('q3_monthly_out')
    q3_model_summary = data.get('q3_model_summary', '')
    q3_model_info = data.get('q3_model_info_dict', {})
    q3_near_normal = data.get('q3_near_normal')
    
    with open(os.path.join(OUTPUT_DIR, "EDI_Q3_aligned_model_summary.txt"), 'w') as f:
        f.write("Ecological_Impacts — Q3-aligned EDI run summary\n")
        f.write("="*60 + "\n\n")
        
        f.write("=== Q3-aligned purpose ===\n")
        f.write("Primary ecological impact score is predicted upland WUE_T percent change from near-normal SPEI.\n")
        f.write("Near-normal SPEI is defined as -1 <= SPEI <= 1 within each coastline and SPEI timescale.\n")
        f.write("The model matches the final Q3 coast-threshold GAM logic.\n\n")
        
        f.write("=== Model formula ===\n")
        f.write("WUE_T ~ SPEI_timescale * coast_region + water_class + month_f +\n")
        f.write("  s(SPEI_value, by = spei_coast, k = 6) + s(site_name, bs = 're')\n\n")
        
        f.write("=== Model information ===\n")
        if q3_model_info:
            for k, v in q3_model_info.items():
                f.write(f"  {k}: {v}\n")
        if q3_model_summary:
            f.write("\n=== Model summary ===\n")
            f.write(q3_model_summary)
        
        f.write("\n=== Upland coast-specific 5/10/20% SPEI thresholds ===\n")
        if q3_threshold is not None:
            f.write(q3_threshold.to_string(index=False))
        
        f.write("\n\n=== SPEI-3 monthly upland EDI class counts by coastline ===\n")
        if q3_monthly is not None and not q3_monthly.empty:
            counts = q3_monthly.groupby(['coast_region', 'EDI_Q3_class']).size().unstack(fill_value=0)
            f.write(counts.to_string())
        
        f.write("\n\n=== SPEI-3 monthly upland EDI direction counts by coastline ===\n")
        if q3_monthly is not None and not q3_monthly.empty:
            dir_counts = q3_monthly.groupby(['coast_region', 'EDI_Q3_direction']).size().unstack(fill_value=0)
            f.write(dir_counts.to_string())
        
        f.write("\n\n=== SPEI-3 prediction-interval support by coastline ===\n")
        if q3_monthly is not None and not q3_monthly.empty:
            pi_counts = q3_monthly.groupby(['coast_region', 'prediction_interval_support']).size().unstack(fill_value=0)
            f.write(pi_counts.to_string())
        
        if q3_near_normal is not None:
            f.write("\n\n=== Near-normal WUE_T by coast and timescale ===\n")
            f.write(q3_near_normal.to_string(index=False))
    
    print("  Wrote Q4 output: EDI_Q3_aligned_model_summary.txt")


def write_logistic_summary(data, logistic_results):
    print("\n" + "="*60)
    print("Step 10: Writing Logistic Model Summary")
    print("="*60)
    
    model_df = logistic_results['model_df']
    response_df = data.get('response_df')
    response_thresholds = data.get('response_thresholds')
    response_summary = data.get('response_summary')
    validation_df = data.get('validation_df')
    monthly_out = data.get('monthly_out')
    
    with open(os.path.join(OUTPUT_DIR, "EDI_logistic_model_summary.txt"), 'w') as f:
        f.write("Ecological_Impacts — run summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"SPEI timescale used: SPEI_3\n")
        f.write(f"Output: {OUTPUT_DIR}\n\n")
        
        f.write("=== Calibration data ===\n")
        f.write(f"Total site-months with SPEI_3 < 0 and valid baseline: {len(model_df)}\n")
        f.write(f"Of which resistance < 1 (impacted): {model_df['y'].sum()}\n")
        f.write(f"Proportion impacted: {model_df['y'].mean() * 100:.1f}%\n\n")
        
        f.write("=== Response-based EDI model ===\n")
        f.write("Response: WUE_T percent change from site-month near-normal baseline\n")
        if response_df is not None:
            f.write(f"Rows with valid response baseline: {len(response_df)}\n")
            f.write(f"Sites: {len(response_df['site_name'].unique())}\n")
        
        if response_summary:
            f.write("\n=== Response GAM summary ===\n")
            f.write(str(response_summary))
            f.write("\n\n")
        
        if response_thresholds is not None:
            f.write("=== Response-based EDI thresholds ===\n")
            f.write(response_thresholds.to_string(index=False))
            f.write("\n\n")
        
        if monthly_out is not None:
            if 'response_EDI_class' in monthly_out.columns:
                f.write("=== Response class counts across all site-months ===\n")
                counts = monthly_out.groupby(['water_class', 'response_EDI_class']).size().unstack(fill_value=0)
                f.write(counts.to_string())
                f.write("\n\n")
                f.write("=== Response direction counts across all site-months ===\n")
                dir_counts = monthly_out.groupby(['water_class', 'response_EDI_direction']).size().unstack(fill_value=0)
                f.write(dir_counts.to_string())
                f.write("\n\n")
        
        f.write("=== Pooled logistic model ===\n")
        f.write(f"Intercept: {logistic_results['pooled']['b0']:.4f}\n")
        f.write(f"Slope: {logistic_results['pooled']['b1']:.4f}\n")
        f.write(f"AIC: {logistic_results['pooled']['aic']:.1f}\n\n")
        
        f.write("=== Ecosystem-stratified models ===\n")
        for ec in ECOSYSTEM_CLASSES:
            if logistic_results['eco_coefs'].get(ec):
                b0_ec, b1_ec = logistic_results['eco_coefs'][ec]
                aic_ec = logistic_results['eco_results'][ec]['aic']
                f.write(f"{ec}: intercept={b0_ec:.4f}, slope={b1_ec:.4f}, AIC={aic_ec:.1f}\n")
            else:
                f.write(f"{ec}: model not fitted (insufficient events)\n")
        
        f.write("\n=== SPEI-3 thresholds at EDI class boundaries ===\n")
        thresholds = data.get('thresholds')
        if thresholds is not None:
            f.write(thresholds.to_string(index=False))
        
        f.write("\n\n=== EDI validation against Q2 metrics (drought months, SPEI_3 < -1) ===\n")
        if validation_df is not None and not validation_df.empty:
            f.write(validation_df.to_string(index=False))
        
        f.write("\n\n=== Severity class counts across all site-months ===\n")
        if monthly_out is not None:
            severity_counts = monthly_out.groupby(['water_class', 'EDI_class']).size().unstack(fill_value=0)
            f.write(severity_counts.to_string())
    
    print("  Wrote Q4 output: EDI_logistic_model_summary.txt")


# ============================================================================
# Step 11: Figures
# ============================================================================

def create_all_figures(data, base_df, logistic_results):
    print("\n" + "="*60)
    print("Step 11: Creating Q4 Figures")
    print("="*60)
    
    create_figure_08_q3_composite(data)
    create_figure_00_response_curve(data)
    create_figure_00b_response_proportions(base_df)
    create_figure_01_logistic_curve(logistic_results)
    create_figure_02_resistance_calibration(base_df)
    create_figure_03_spei_edi_scatter(base_df)
    create_figure_04_edi_timeseries(base_df, data.get('resistance'))
    create_figure_05_severity_proportions(base_df)
    create_figure_06_threshold_validation(data.get('validation_df'))
    create_figure_07_composite(base_df, data)


def create_figure_08_q3_composite(data):
    print("  Creating Figure 08: Q3-Aligned Composite...")
    q3_pred = data.get('q3_pred_edi')
    q3_monthly = data.get('q3_monthly_out')
    q3_threshold = data.get('q3_threshold_edi')
    
    if q3_pred is None:
        print("    Skipping")
        return
    
    fig = plt.figure(figsize=(11, 15))
    
    gs_top = fig.add_gridspec(4, 3, top=0.68, bottom=0.35, left=0.08, right=0.92, hspace=0.15, wspace=0.15)
    selected_timescales = ['SPEI_1', 'SPEI_3', 'SPEI_48']
    q3_plot = q3_pred[q3_pred['SPEI_timescale'].isin(selected_timescales)].copy()
    
    panel_labels = {'Atlantic Coast': 'Atlantic', 'Pacific Coast': 'Pacific', 
                    'Gulf Coast': 'Gulf', 'AK Coast': 'AK'}
    
    threshold_markers = q3_threshold[q3_threshold['SPEI_timescale'].isin(selected_timescales)].copy() if q3_threshold is not None else None
    
    for i, coast in enumerate(COAST_REGION_LEVELS):
        for j, ts in enumerate(selected_timescales):
            ax = fig.add_subplot(gs_top[i, j])
            sub = q3_plot[(q3_plot['coast_region'] == coast) & (q3_plot['SPEI_timescale'] == ts)]
            
            if sub.empty:
                ax.set_xlim(-4, 4)
                ax.set_ylim(-30, 30)
                continue
            
            ax.axvline(x=-1, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=0, color='gray', alpha=0.5)
            for y in [-20, -10, -5, 5, 10, 20]:
                ax.axhline(y=y, color='gray', linestyle=':', alpha=0.3)
            
            ax.fill_between(sub['SPEI_value'], sub['predicted_lower_pct'], sub['predicted_upper_pct'],
                            alpha=0.15, color=COAST_COLORS.get(coast, 'gray'))
            ax.plot(sub['SPEI_value'], sub['predicted_pct_change'],
                    color=COAST_COLORS.get(coast, 'black'), linewidth=1.5)
            
            if threshold_markers is not None:
                markers = threshold_markers[
                    (threshold_markers['coast_region'] == coast) & 
                    (threshold_markers['SPEI_timescale'] == ts)
                ]
                if not markers.empty:
                    colors = {'5%': 'blue', '10%': 'orange', '20%': 'red'}
                    shapes = {'5%': 'o', '10%': 's', '20%': '^'}
                    for _, row in markers.iterrows():
                        if pd.notna(row['SPEI_threshold']):
                            ax.scatter(row['SPEI_threshold'], row['pct_change_threshold'],
                                      color=colors.get(str(row['threshold_pct']), 'gray'),
                                      s=80, marker=shapes.get(str(row['threshold_pct']), 'o'),
                                      edgecolor='black', linewidth=0.5)
                        if pd.notna(row['SPEI_threshold_lower']) and pd.notna(row['SPEI_threshold_upper']):
                            ax.plot([row['SPEI_threshold_lower'], row['SPEI_threshold_upper']],
                                   [row['pct_change_threshold'], row['pct_change_threshold']],
                                   color=colors.get(str(row['threshold_pct']), 'gray'),
                                   linewidth=2, alpha=0.6)
            
            if i == 3:
                ax.set_xlabel('SPEI')
            if j == 0:
                ax.set_ylabel('WUE$_T$ change (%)')
            if i == 0:
                ax.set_title(ts.replace('SPEI_', 'SPEI-'))
            if j == 0:
                ax.text(0.02, 0.95, panel_labels.get(coast, coast), transform=ax.transAxes, 
                        fontweight='bold', fontsize=11, verticalalignment='top')
            ax.set_xlim(-4, 4)
            ax.set_ylim(-30, 30)
            theme_wue(ax)
    
    gs_bottom = fig.add_gridspec(1, 2, top=0.28, bottom=0.08, left=0.08, right=0.92, hspace=0.15, wspace=0.15)
    
    ax_bl = fig.add_subplot(gs_bottom[0, 0])
    if q3_monthly is not None and not q3_monthly.empty:
        q3_class_props = q3_monthly.groupby(['coast_region', 'EDI_Q3_class']).size().reset_index(name='n')
        q3_class_props['prop'] = q3_class_props.groupby('coast_region')['n'].transform(lambda x: x / x.sum() * 100)
        pivot = q3_class_props.pivot(index='coast_region', columns='EDI_Q3_class', values='prop').fillna(0)
        for cls in ["No meaningful change (<5%)", "Watch (5-10%)", "Stress (10-20%)", "Impact (>=20%)"]:
            if cls not in pivot.columns:
                pivot[cls] = 0
        pivot = pivot[["No meaningful change (<5%)", "Watch (5-10%)", "Stress (10-20%)", "Impact (>=20%)"]]
        pivot.plot(kind='bar', stacked=True, ax=ax_bl, width=0.68,
                   color=[Q3_CLASS_COLORS.get(cls, '#D9D9D9') for cls in pivot.columns],
                   edgecolor='white', linewidth=0.25)
        ax_bl.set_xlabel(None)
        ax_bl.set_ylabel('% of upland site-months')
        ax_bl.set_title('Monthly SPEI-3 Upland EDI Classes')
        ax_bl.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
        ax_bl.legend(title='Q3 EDI class', loc='upper right', fontsize=8)
        theme_wue(ax_bl)
    
    ax_br = fig.add_subplot(gs_bottom[0, 1])
    if q3_monthly is not None and not q3_monthly.empty:
        q3_pi_props = q3_monthly.groupby(['coast_region', 'prediction_interval_support']).size().reset_index(name='n')
        q3_pi_props['prop'] = q3_pi_props.groupby('coast_region')['n'].transform(lambda x: x / x.sum() * 100)
        pivot_pi = q3_pi_props.pivot(index='coast_region', columns='prediction_interval_support', values='prop').fillna(0)
        pi_cats = ["PI supports decrease", "PI overlaps no-change", "PI supports increase"]
        for cat in pi_cats:
            if cat not in pivot_pi.columns:
                pivot_pi[cat] = 0
        pivot_pi = pivot_pi[pi_cats]
        pivot_pi.plot(kind='bar', stacked=True, ax=ax_br, width=0.68,
                      color=[Q3_PI_COLORS.get(cat, '#D9D9D9') for cat in pi_cats],
                      edgecolor='white', linewidth=0.25)
        ax_br.set_xlabel(None)
        ax_br.set_ylabel('% of upland site-months')
        ax_br.set_title('Prediction-Interval Support')
        ax_br.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
        ax_br.legend(title='95% interval', loc='upper right', fontsize=8)
        theme_wue(ax_br)
    
    fig.suptitle('Q3-Aligned Ecological Drought Index for Upland WUE$_T$', fontsize=18, fontweight='bold', y=0.98)
    fig.text(0.5, 0.32, 
             'Primary EDI is predicted WUE$_T$ change from near-normal SPEI; '
             'older logistic EDI outputs are retained as diagnostics.',
             ha='center', fontsize=11)
    plt.tight_layout()
    save_plot(fig, "EDI_plot_08_Q3_aligned_composite.png", width=11, height=15)


def create_figure_00_response_curve(data):
    print("  Creating Figure 00: Response Curve...")
    grid_df = data.get('response_grid')
    if grid_df is None or grid_df.empty:
        print("    Skipping")
        return
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    for idx, ec in enumerate(ECOSYSTEM_CLASSES):
        ax = axes[idx]
        sub = grid_df[grid_df['water_class'] == ec]
        if sub.empty:
            continue
        ax.axvspan(-1, 1, alpha=0.15, color='gray')
        ax.axhline(y=0, color='gray', linewidth=0.5)
        for y in [-5, -2.5, 2.5, 5]:
            ax.axhline(y=y, color='gray', linestyle=':', alpha=0.5)
        ax.fill_between(sub['SPEI_3'], sub['predicted_pct_change_lower'], sub['predicted_pct_change_upper'],
                        alpha=0.15, color=ECOSYSTEM_COLORS.get(ec, 'gray'))
        ax.plot(sub['SPEI_3'], sub['predicted_pct_change'],
                color=ECOSYSTEM_COLORS.get(ec, 'black'), linewidth=1.5)
        ax.set_xlabel('SPEI-3')
        ax.set_ylabel('Predicted WUE$_T$ change (%)')
        ax.set_title(f'{ec}')
        ax.set_xlim(-4, 4)
        theme_wue(ax)
    plt.suptitle('Response-Based EDI: Expected WUE$_T$ Change Across SPEI-3', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "EDI_plot_00_response_curve.png", width=10, height=12)


def create_figure_00b_response_proportions(base_df):
    print("  Creating Figure 00b: Response Class Proportions...")
    if 'response_EDI_class' not in base_df.columns:
        print("    Skipping")
        return
    
    response_df = base_df[base_df['response_EDI_class'].notna()].copy()
    if response_df.empty:
        print("    Skipping")
        return
    
    props = response_df.groupby(['water_class', 'response_EDI_class', 'response_EDI_direction']).size().reset_index(name='n')
    props['prop'] = props.groupby('water_class')['n'].transform(lambda x: x / x.sum() * 100)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    pivot = props.pivot_table(index='water_class', columns='response_EDI_class', values='prop', fill_value=0)
    pivot.plot(kind='bar', stacked=True, ax=ax, width=0.65)
    ax.set_ylabel('% of site-months')
    ax.set_xlabel('Ecosystem class')
    ax.set_title('Response-Based EDI Class Proportions')
    ax.legend(title='Response class', loc='upper right')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    theme_wue(ax)
    plt.tight_layout()
    save_plot(fig, "EDI_plot_00b_response_class_proportions.png", width=10, height=7)


def create_figure_01_logistic_curve(logistic_results):
    print("  Creating Figure 01: Logistic Curve...")
    spei_seq = np.linspace(-4, 2, 300)
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for ymin, ymax, label, color in [
        (0, 0.30, 'None', SEVERITY_COLORS['None']),
        (0.30, 0.50, 'Watch', SEVERITY_COLORS['Watch']),
        (0.50, 0.70, 'Stress', SEVERITY_COLORS['Stress']),
        (0.70, 1.0, 'Impact', SEVERITY_COLORS['Impact'])
    ]:
        ax.axhspan(ymin, ymax, alpha=0.10, color=color)
    for y in [0.30, 0.50, 0.70]:
        ax.axhline(y=y, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=-1, color='gray', linestyle=':', alpha=0.5)
    
    b0, b1 = logistic_results['pooled']['b0'], logistic_results['pooled']['b1']
    ax.plot(spei_seq, logistic(spei_seq, b0, b1), color='black', linewidth=1.5, label='All (pooled)')
    
    for ec in ECOSYSTEM_CLASSES:
        if logistic_results['eco_coefs'].get(ec):
            b0_ec, b1_ec = logistic_results['eco_coefs'][ec]
            ax.plot(spei_seq, logistic(spei_seq, b0_ec, b1_ec),
                    color=ECOSYSTEM_COLORS.get(ec, 'gray'), linewidth=1.5, linestyle='--', label=ec)
    
    for pos, label in zip([0.15, 0.40, 0.60, 0.85], ['None', 'Watch', 'Stress', 'Impact']):
        ax.text(-3.8, pos, label, fontsize=10, color='gray')
    
    ax.set_xlabel('SPEI-3')
    ax.set_ylabel('Ecological Drought Index (EDI)')
    ax.set_title('Logistic Calibration to SPEI-3')
    ax.legend(loc='upper left')
    ax.set_xlim(-4, 2)
    ax.set_ylim(0, 1)
    theme_wue(ax)
    plt.tight_layout()
    save_plot(fig, "EDI_plot_01_logistic_curve.png", width=10, height=7)


def create_figure_02_resistance_calibration(base_df):
    print("  Creating Figure 02: Resistance Calibration...")
    cal_df = base_df[
        (base_df['monthly_resistance'].notna()) &
        (np.isfinite(base_df['monthly_resistance'])) &
        (base_df['SPEI_3'] < 0) &
        (base_df['monthly_resistance'] < 5)
    ].copy()
    
    if cal_df.empty:
        print("    Skipping")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(13, 6))
    for idx, ec in enumerate(ECOSYSTEM_CLASSES):
        ax = axes[idx]
        sub = cal_df[cal_df['water_class'] == ec]
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
        ax.axvline(x=-1, color='gray', linestyle=':', alpha=0.5)
        ax.scatter(sub['SPEI_3'], sub['monthly_resistance'], alpha=0.2, s=10, 
                   color=ECOSYSTEM_COLORS.get(ec, 'gray'))
        if len(sub) > 10:
            sorted_sub = sub.sort_values('SPEI_3')
            window = max(10, len(sorted_sub) // 20)
            smoothed = sorted_sub['monthly_resistance'].rolling(window, center=True).mean()
            ax.plot(sorted_sub['SPEI_3'], smoothed, color=ECOSYSTEM_COLORS.get(ec, 'black'), linewidth=1.5)
        ax.set_title(ec)
        ax.set_xlabel('SPEI-3')
        if idx == 0:
            ax.set_ylabel('Monthly resistance')
        ax.set_xlim(-4, 0.5)
        theme_wue(ax)
    plt.suptitle('Monthly Resistance vs SPEI-3: Calibration Data', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "EDI_plot_02_monthly_resistance_calibration.png", width=13, height=6)


def create_figure_03_spei_edi_scatter(base_df):
    print("  Creating Figure 03: SPEI vs EDI Scatter...")
    scatter_df = base_df[np.isfinite(base_df['SPEI_3'])].copy()
    if scatter_df.empty:
        print("    Skipping")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(13, 6))
    for idx, ec in enumerate(ECOSYSTEM_CLASSES):
        ax = axes[idx]
        sub = scatter_df[scatter_df['water_class'] == ec]
        ax.scatter(sub['SPEI_3'], sub['EDI'], alpha=0.2, s=8, color=ECOSYSTEM_COLORS.get(ec, 'gray'))
        if len(sub) > 10:
            sorted_sub = sub.sort_values('SPEI_3')
            window = max(10, len(sorted_sub) // 20)
            smoothed = sorted_sub['EDI'].rolling(window, center=True).mean()
            ax.plot(sorted_sub['SPEI_3'], smoothed, color=ECOSYSTEM_COLORS.get(ec, 'black'), linewidth=1.5)
        for y in [0.30, 0.50, 0.70]:
            ax.axhline(y=y, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=-1, color='gray', linestyle=':', alpha=0.5)
        ax.set_title(ec)
        ax.set_xlabel('SPEI-3')
        if idx == 0:
            ax.set_ylabel('EDI')
        ax.set_xlim(-4, 2)
        ax.set_ylim(0, 1)
        theme_wue(ax)
    plt.suptitle('EDI vs SPEI-3: All Site-Months', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "EDI_plot_03_spei_edi_by_ecosystem.png", width=13, height=6)


def create_figure_04_edi_timeseries(base_df, resistance_df):
    print("  Creating Figure 04: EDI Time Series...")
    if resistance_df is None:
        print("    Skipping")
        return
    
    resistance_sorted = resistance_df.sort_values('resistance')
    featured = list(resistance_sorted['site_name'].head(4)) + list(resistance_sorted['site_name'].tail(4))
    featured = [s for s in featured if s in base_df['site_name'].unique()][:min(12, len(featured))]
    
    if not featured:
        print("    Skipping")
        return
    
    ts_df = base_df[base_df['site_name'].isin(featured)].copy()
    ts_df['date'] = pd.to_datetime(ts_df[['Year', 'month']].assign(day=15))
    
    n_cols = min(3, len(featured))
    n_rows = int(np.ceil(len(featured) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 10))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, site in enumerate(featured):
        if idx >= len(axes):
            break
        ax = axes[idx]
        site_df = ts_df[ts_df['site_name'] == site].sort_values('date')
        for cls in EDI_LABELS:
            sub = site_df[site_df['EDI_class'] == cls]
            if not sub.empty:
                ax.bar(sub['date'], sub['EDI'], width=25, color=SEVERITY_COLORS.get(cls, 'gray'), label=cls if idx == 0 else None)
        for y in [0.30, 0.50, 0.70]:
            ax.axhline(y=y, color='gray', linestyle='--', alpha=0.5, linewidth=0.5)
        ax.set_title(site)
        ax.set_ylim(0, 1)
        ax.set_xlim(site_df['date'].min(), site_df['date'].max())
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.tick_params(axis='x', rotation=30)
        if idx % n_cols == 0:
            ax.set_ylabel('EDI')
        theme_wue(ax)
    
    for idx in range(len(featured), len(axes)):
        axes[idx].set_visible(False)
    
    handles = [mpatches.Patch(color=SEVERITY_COLORS.get(cls, 'gray'), label=cls) for cls in EDI_LABELS]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=10)
    plt.suptitle('EDI Time Series: Sites with Lowest and Highest Resistance', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "EDI_plot_04_edi_timeseries_facet.png", width=14, height=10)


def create_figure_05_severity_proportions(base_df):
    print("  Creating Figure 05: Severity Class Proportions...")
    if 'EDI_class' not in base_df.columns:
        print("    Skipping")
        return
    
    props = base_df.groupby(['water_class', 'EDI_class']).size().reset_index(name='n')
    props['prop'] = props.groupby('water_class')['n'].transform(lambda x: x / x.sum() * 100)
    
    fig, ax = plt.subplots(figsize=(9, 7))
    pivot = props.pivot_table(index='water_class', columns='EDI_class', values='prop', fill_value=0)
    pivot = pivot[[c for c in EDI_LABELS if c in pivot.columns]]
    pivot.plot(kind='bar', stacked=True, ax=ax, width=0.6,
               color=[SEVERITY_COLORS.get(cls, 'gray') for cls in pivot.columns])
    ax.set_ylabel('% of site-months')
    ax.set_xlabel('Ecosystem class')
    ax.set_title('EDI Severity Class Proportions by Ecosystem')
    ax.legend(title='EDI class', loc='upper right')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    theme_wue(ax)
    plt.tight_layout()
    save_plot(fig, "EDI_plot_05_severity_class_proportions.png", width=9, height=7)


def create_figure_06_threshold_validation(validation_df):
    print("  Creating Figure 06: Threshold Validation...")
    if validation_df is None or validation_df.empty:
        print("    Skipping")
        return
    
    val_long = validation_df.melt(
        id_vars=['EDI_class', 'water_class'],
        value_vars=['mean_monthly_resist', 'mean_site_resistance', 'mean_recovery'],
        var_name='metric', value_name='value'
    )
    metric_labels = {
        'mean_monthly_resist': 'Monthly resistance\n(WUE_T / baseline, this month)',
        'mean_site_resistance': 'Site resistance\n(mean drought / near-normal, Q2)',
        'mean_recovery': 'Site recovery\n(post / pre-drought WUE_T, Q2)'
    }
    val_long['metric_label'] = val_long['metric'].map(metric_labels)
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 6))
    for idx, metric in enumerate(['mean_monthly_resist', 'mean_site_resistance', 'mean_recovery']):
        ax = axes[idx]
        sub = val_long[val_long['metric'] == metric]
        for ec in ECOSYSTEM_CLASSES:
            sub_ec = sub[sub['water_class'] == ec]
            if not sub_ec.empty:
                ax.plot(sub_ec['EDI_class'], sub_ec['value'], marker='o', linewidth=1.5,
                        color=ECOSYSTEM_COLORS.get(ec, 'gray'), label=ec if idx == 0 else None)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
        ax.set_title(metric_labels.get(metric, metric))
        ax.set_xlabel('EDI class')
        if idx == 0:
            ax.set_ylabel('Mean performance ratio')
        ax.set_ylim(0, 1.5)
        theme_wue(ax)
    
    handles = [mpatches.Patch(color=ECOSYSTEM_COLORS.get(ec, 'gray'), label=ec) for ec in ECOSYSTEM_CLASSES]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=10)
    plt.suptitle('EDI Validation Against Q2 Performance Metrics', fontsize=14, fontweight='bold')
    plt.tight_layout()
    save_plot(fig, "EDI_plot_06_threshold_validation.png", width=14, height=6)


def create_figure_07_composite(base_df, data):
    print("  Creating Figure 07: Composite...")
    
    fig = plt.figure(figsize=(16, 34))
    fig.suptitle('Ecological Drought Index (EDI): SPEI-3 and WUE$_T$ Response', fontsize=18, fontweight='bold')
    gs = fig.add_gridspec(4, 2, hspace=0.3, wspace=0.3)
    
    # Plot 1: Response curve
    ax1 = fig.add_subplot(gs[0, 0])
    grid_df = data.get('response_grid')
    if grid_df is not None and not grid_df.empty:
        for ec in ECOSYSTEM_CLASSES:
            sub = grid_df[grid_df['water_class'] == ec]
            if not sub.empty:
                ax1.plot(sub['SPEI_3'], sub['predicted_pct_change'], 
                         color=ECOSYSTEM_COLORS.get(ec, 'gray'), linewidth=1.5, label=ec)
                ax1.fill_between(sub['SPEI_3'], sub['predicted_pct_change_lower'], sub['predicted_pct_change_upper'],
                                 alpha=0.15, color=ECOSYSTEM_COLORS.get(ec, 'gray'))
        ax1.axvspan(-1, 1, alpha=0.15, color='gray')
        ax1.axhline(y=0, color='black', linewidth=0.5)
        ax1.set_xlabel('SPEI-3')
        ax1.set_ylabel('WUE$_T$ change (%)')
        ax1.set_title('Response-Based EDI Curves')
        ax1.legend(loc='best')
        theme_wue(ax1)
    else:
        ax1.text(0.5, 0.5, 'Response Curve\n(no data)', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_axis_off()
    
    # Plot 2: Response proportions
    ax2 = fig.add_subplot(gs[0, 1])
    if 'response_EDI_class' in base_df.columns:
        response_df = base_df[base_df['response_EDI_class'].notna()].copy()
        if not response_df.empty:
            props = response_df.groupby(['water_class', 'response_EDI_class']).size().reset_index(name='n')
            props['prop'] = props.groupby('water_class')['n'].transform(lambda x: x / x.sum() * 100)
            pivot = props.pivot(index='water_class', columns='response_EDI_class', values='prop').fillna(0)
            pivot.plot(kind='bar', stacked=True, ax=ax2, width=0.6)
            ax2.set_xlabel('Ecosystem class')
            ax2.set_ylabel('% of site-months')
            ax2.set_title('Response Class Proportions')
            ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
            theme_wue(ax2)
    else:
        ax2.text(0.5, 0.5, 'Response Proportions\n(no data)', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_axis_off()
    
    # Plot 3: Logistic curve
    ax3 = fig.add_subplot(gs[1, 0])
    if 'EDI' in base_df.columns and 'SPEI_3' in base_df.columns:
        ax3.scatter(base_df['SPEI_3'], base_df['EDI'], alpha=0.05, s=5)
        ax3.set_xlabel('SPEI-3')
        ax3.set_ylabel('EDI')
        ax3.set_title('Logistic EDI Distribution')
        theme_wue(ax3)
    else:
        ax3.text(0.5, 0.5, 'Logistic Curve\n(no data)', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_axis_off()
    
    # Plot 4: Severity proportions
    ax4 = fig.add_subplot(gs[1, 1])
    if 'EDI_class' in base_df.columns:
        props = base_df.groupby(['water_class', 'EDI_class']).size().reset_index(name='n')
        props['prop'] = props.groupby('water_class')['n'].transform(lambda x: x / x.sum() * 100)
        pivot = props.pivot(index='water_class', columns='EDI_class', values='prop').fillna(0)
        pivot = pivot[[c for c in EDI_LABELS if c in pivot.columns]]
        pivot.plot(kind='bar', stacked=True, ax=ax4, width=0.6,
                   color=[SEVERITY_COLORS.get(cls, 'gray') for cls in pivot.columns])
        ax4.set_xlabel('Ecosystem class')
        ax4.set_ylabel('% of site-months')
        ax4.set_title('EDI Severity Class Proportions')
        ax4.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
        theme_wue(ax4)
    else:
        ax4.text(0.5, 0.5, 'Severity Proportions\n(no data)', ha='center', va='center', transform=ax4.transAxes)
        ax4.set_axis_off()
    
    # Plot 5: Resistance calibration
    ax5 = fig.add_subplot(gs[2, :])
    cal_df = base_df[
        (base_df['monthly_resistance'].notna()) &
        (np.isfinite(base_df['monthly_resistance'])) &
        (base_df['SPEI_3'] < 0) &
        (base_df['monthly_resistance'] < 5)
    ].copy()
    if not cal_df.empty:
        for ec in ECOSYSTEM_CLASSES:
            sub = cal_df[cal_df['water_class'] == ec]
            if not sub.empty:
                ax5.scatter(sub['SPEI_3'], sub['monthly_resistance'], alpha=0.15, s=8,
                           color=ECOSYSTEM_COLORS.get(ec, 'gray'), label=ec)
        ax5.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
        ax5.axvline(x=-1, color='gray', linestyle=':', alpha=0.5)
        ax5.set_xlabel('SPEI-3')
        ax5.set_ylabel('Monthly resistance')
        ax5.set_title('Resistance Calibration')
        ax5.legend(loc='best')
        ax5.set_xlim(-4, 0.5)
        theme_wue(ax5)
    else:
        ax5.text(0.5, 0.5, 'Resistance Calibration\n(no data)', ha='center', va='center', transform=ax5.transAxes)
        ax5.set_axis_off()
    
    # Plot 6: Validation
    ax6 = fig.add_subplot(gs[3, :])
    val_df = data.get('validation_df')
    if val_df is not None and not val_df.empty:
        val_long = val_df.melt(
            id_vars=['EDI_class', 'water_class'],
            value_vars=['mean_monthly_resist', 'mean_site_resistance', 'mean_recovery'],
            var_name='metric', value_name='value'
        )
        for ec in ECOSYSTEM_CLASSES:
            sub = val_long[val_long['water_class'] == ec]
            if not sub.empty:
                ax6.plot(sub['EDI_class'], sub['value'], marker='o', linewidth=1.5,
                        color=ECOSYSTEM_COLORS.get(ec, 'gray'), label=ec)
        ax6.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
        ax6.set_xlabel('EDI class')
        ax6.set_ylabel('Mean performance ratio')
        ax6.set_title('EDI Validation')
        ax6.legend(loc='best')
        theme_wue(ax6)
    else:
        ax6.text(0.5, 0.5, 'Threshold Validation\n(no data)', ha='center', va='center', transform=ax6.transAxes)
        ax6.set_axis_off()
    
    plt.tight_layout()
    save_plot(fig, "EDI_plot_07_composite.png", width=16, height=34)


# ============================================================================
# Step 12: Verification
# ============================================================================

def verify_outputs():
    print("\n" + "="*60)
    print("Step 12: Verifying Q4 Outputs")
    print("="*60)
    
    required = [
        "EDI_Q3_aligned_prediction_scores_upland.csv",
        "EDI_Q3_aligned_coast_thresholds_5_10_20pct.csv",
        "EDI_Q3_aligned_model_info.csv",
        "EDI_Q3_aligned_coast_GAM_smooth_terms.csv",
        "EDI_Q3_aligned_coast_GAM_parametric_terms.csv",
        "EDI_Q3_aligned_monthly_scores.csv",
        "EDI_Q3_aligned_monthly_regional_summary.csv",
        "EDI_Q3_aligned_regional_validation.csv",
        "EDI_Q3_aligned_coast_threshold_gam_model.rds",
        "EDI_logistic_coefficients.csv",
        "EDI_threshold_summary.csv",
        "EDI_monthly_site_scores.csv",
        "EDI_validation_by_class.csv",
        "EDI_response_prediction_curves.csv",
        "EDI_response_threshold_summary.csv",
        "EDI_response_gam_smooth_terms.csv",
        "EDI_response_gam_parametric_terms.csv",
        "EDI_Q3_aligned_model_summary.txt",
        "EDI_logistic_model_summary.txt",
        "EDI_plot_08_Q3_aligned_composite.png",
        "EDI_plot_00_response_curve.png",
        "EDI_plot_00b_response_class_proportions.png",
        "EDI_plot_01_logistic_curve.png",
        "EDI_plot_02_monthly_resistance_calibration.png",
        "EDI_plot_03_spei_edi_by_ecosystem.png",
        "EDI_plot_04_edi_timeseries_facet.png",
        "EDI_plot_05_severity_class_proportions.png",
        "EDI_plot_06_threshold_validation.png",
        "EDI_plot_07_composite.png"
    ]
    
    missing = []
    for f in required:
        out_path = os.path.join(OUTPUT_DIR, f)
        fig_path = os.path.join(FIGURE_DIR, f)
        if not os.path.exists(out_path) and not os.path.exists(fig_path):
            missing.append(f)
    
    if missing:
        print(f"\nERROR: Missing {len(missing)} required Q4 outputs:")
        for m in missing:
            print(f"  - {m}")
        raise FileNotFoundError(f"Missing {len(missing)} required Q4 outputs.")
    
    print("  All required Q4 outputs verified!")


# ============================================================================
# Main
# ============================================================================

def main():
    print("\n" + "="*60)
    print("Q4_Ecological_Impacts.py")
    print("Exact Replication of Malone's Ecological_Impacts.R Workflow")
    print("="*60)
    print(f"Start time: {datetime.now()}")
    print(f"R available: {R_AVAILABLE}")
    print(f"R path: {RSCRIPT_EXE}")
    print("\n=== Q4 PATH RULES ===")
    print(f"Reading Q3 inputs from: {PYTHON_Q3_DIR}")
    print(f"Writing Q4 outputs to: {OUTPUT_DIR}")
    print(f"Writing Q4 figures to: {FIGURE_DIR}")
    print("NEVER write to or modify Python Q3 output folder")
    print("="*60)
    
    try:
        create_directories()
        clean_output_directories()
        
        loaded = load_data()
        monthly = loaded['monthly']
        resistance = loaded['resistance']
        recovery = loaded['recovery']
        q3_pred = loaded['q3_pred']
        q3_threshold = loaded['q3_threshold']
        q3_smooth = loaded['q3_smooth']
        q3_param = loaded['q3_param']
        q3_model_comp = loaded['q3_model_comp']
        q3_site_sens = loaded['q3_site_sens']
        
        # Prepare base data with row ID for safe merge-back
        base_df = prepare_base_data(monthly)
        base_df, model_ready = compute_baseline(base_df)
        
        # Add stable row ID for safe prediction merge-back
        # Use 'row_id' (not '__row_id') - R mangles column names with underscores
        base_df = base_df.reset_index(drop=True)
        base_df['row_id'] = np.arange(len(base_df))
        
        # Recreate model_ready with row ID
        model_ready = base_df[base_df['baseline_WUE_T'].notna()].copy()
        
        logistic_results = fit_logistic_models_r(model_ready)
        base_df = apply_logistic_edi(base_df, logistic_results)
        
        response_result = fit_response_edi_gam(model_ready)
        response_grid = response_result.get('grid_pred')
        response_thresholds = compute_response_thresholds(response_grid)
        
        # Merge response predictions using row_id (safe, avoids floating point issues)
        if response_result.get('all_pred') is not None:
            all_pred = response_result['all_pred']
            if 'row_id' not in all_pred.columns:
                raise RuntimeError("Response prediction output is missing row_id; cannot safely merge predictions back.")
            
            all_pred = all_pred.rename(columns={
                'predicted': 'response_EDI_pct_change',
                'predicted_lower': 'response_EDI_pct_change_lower',
                'predicted_upper': 'response_EDI_pct_change_upper'
            })
            
            base_df = base_df.merge(
                all_pred[
                    [
                        'row_id',
                        'response_EDI_pct_change',
                        'response_EDI_pct_change_lower',
                        'response_EDI_pct_change_upper'
                    ]
                ],
                on='row_id',
                how='left'
            )
            
            base_df['response_EDI_magnitude'] = base_df['response_EDI_pct_change'].abs()
            base_df['response_EDI_class'] = base_df['response_EDI_magnitude'].apply(
                lambda x: classify_response_magnitude(x) if pd.notna(x) else np.nan
            )
            base_df['response_EDI_direction'] = base_df['response_EDI_pct_change'].apply(
                lambda x: classify_response_direction(x) if pd.notna(x) else np.nan
            )
        
        q3_result = prepare_q3_aligned_data(
            monthly, q3_pred, q3_threshold, q3_smooth, q3_param, q3_model_comp, q3_site_sens
        )
        
        # Build outputs
        outputs = {
            'q3_pred_edi': q3_result['q3_pred_edi'],
            'q3_threshold_edi': q3_result['q3_threshold_edi'],
            'q3_model_info': q3_result['q3_model_info'],
            'q3_smooth': q3_smooth,  # From Python Q3 source
            'q3_param': q3_param,    # From Python Q3 source
            'q3_monthly_out': q3_result['q3_monthly_out'],
            'q3_regional_summary': q3_result['q3_regional_summary'],
            'q3_regional_validation': q3_result['q3_regional_validation'],
            'q3_model_summary': q3_result['q3_model_summary'],
            'q3_model_info_dict': q3_result['q3_model_info_dict'],
            'q3_near_normal': q3_result.get('q3_near_normal'),
            'response_grid': response_grid,
            'response_thresholds': response_thresholds,
            'response_smooth': response_result.get('smooth_df'),
            'response_param': response_result.get('param_df'),
            'response_df': response_result.get('response_df'),
            'response_summary': response_result.get('summary'),
            'resistance': resistance,
            'recovery': recovery
        }
        
        # Logistic coefficients
        coef_rows = [{'water_class': 'All (pooled)', 'intercept': logistic_results['pooled']['b0'],
                      'slope_SPEI3': logistic_results['pooled']['b1'], 'AIC': logistic_results['pooled']['aic'],
                      'n_obs': len(logistic_results['model_df']), 'n_impacted': logistic_results['model_df']['y'].sum()}]
        for ec in ECOSYSTEM_CLASSES:
            if logistic_results['eco_coefs'].get(ec):
                b0_ec, b1_ec = logistic_results['eco_coefs'][ec]
                sub = logistic_results['model_df'][logistic_results['model_df']['water_class'] == ec]
                coef_rows.append({'water_class': ec, 'intercept': b0_ec, 'slope_SPEI3': b1_ec,
                                  'AIC': logistic_results['eco_results'][ec]['aic'], 'n_obs': len(sub),
                                  'n_impacted': sub['y'].sum()})
        outputs['logistic_coefs'] = pd.DataFrame(coef_rows)
        
        # Thresholds
        thresh_rows = [{'model': 'Pooled', 'water_class': 'All',
                        'EDI_0.30_SPEI3': spei_at_p(0.30, logistic_results['pooled']['b0'], logistic_results['pooled']['b1']),
                        'EDI_0.50_SPEI3': spei_at_p(0.50, logistic_results['pooled']['b0'], logistic_results['pooled']['b1']),
                        'EDI_0.70_SPEI3': spei_at_p(0.70, logistic_results['pooled']['b0'], logistic_results['pooled']['b1'])}]
        for ec in ECOSYSTEM_CLASSES:
            if logistic_results['eco_coefs'].get(ec):
                b0_ec, b1_ec = logistic_results['eco_coefs'][ec]
                thresh_rows.append({'model': 'Ecosystem', 'water_class': ec,
                                    'EDI_0.30_SPEI3': spei_at_p(0.30, b0_ec, b1_ec),
                                    'EDI_0.50_SPEI3': spei_at_p(0.50, b0_ec, b1_ec),
                                    'EDI_0.70_SPEI3': spei_at_p(0.70, b0_ec, b1_ec)})
        outputs['thresholds'] = pd.DataFrame(thresh_rows)
        
        # Monthly site scores with response columns
        monthly_cols = ['site_name', 'Year', 'month', 'water_class', 'WUE_T', 'SPEI_3',
                        'baseline_WUE_T', 'monthly_resistance', 'WUE_T_pct_change', 'EDI', 'EDI_class']
        for col in ['response_EDI_pct_change', 'response_EDI_magnitude', 'response_EDI_direction', 'response_EDI_class']:
            if col in base_df.columns:
                monthly_cols.append(col)
        outputs['monthly_out'] = base_df[monthly_cols].copy()
        
        # Validation - use observed=True to only keep observed groups
        val_df = base_df[(base_df['SPEI_3'] < -1) & (base_df['EDI_class'].notna())].copy()
        if not val_df.empty and 'resistance' in resistance.columns:
            val_df = val_df.merge(resistance[['site_name', 'resistance']], on='site_name', how='left')
            if 'mean_recovery' in recovery.columns:
                val_df = val_df.merge(recovery[['site_name', 'mean_recovery']], on='site_name', how='left')
            val_df = val_df.dropna(subset=['resistance'])
            if not val_df.empty:
                outputs['validation_df'] = val_df.groupby(
                    ['EDI_class', 'water_class'],
                    observed=True
                ).agg(
                    n_months=('site_name', 'count'),
                    n_sites=('site_name', 'nunique'),
                    mean_SPEI_3=('SPEI_3', 'mean'),
                    mean_EDI=('EDI', 'mean'),
                    mean_WUE_T=('WUE_T', 'mean'),
                    mean_monthly_resist=('monthly_resistance', 'mean'),
                    mean_site_resistance=('resistance', 'mean'),
                    mean_recovery=('mean_recovery', 'mean') if 'mean_recovery' in val_df.columns else pd.NaT
                ).reset_index()
                # Filter to only rows with observations
                outputs['validation_df'] = outputs['validation_df'][
                    outputs['validation_df']['n_months'] > 0
                ].copy()
            else:
                outputs['validation_df'] = pd.DataFrame()
        else:
            outputs['validation_df'] = pd.DataFrame()
        
        write_all_outputs(outputs)
        write_q3_summary(outputs)
        write_logistic_summary(outputs, logistic_results)
        create_all_figures(outputs, base_df, logistic_results)
        verify_outputs()
        
        print("\n" + "="*60)
        print("Q4_Ecological_Impacts.py completed successfully!")
        print(f"End time: {datetime.now()}")
        print("="*60)
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()