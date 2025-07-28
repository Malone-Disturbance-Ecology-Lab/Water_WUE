# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 13:16:28 2025

@author: ammar
"""




### this code take one merger ameri monthly file that has all sites with monthly WUE and CUE 

####this code merge salinity data with ameri data 

#####################################################################################
import os
import pandas as pd
import numpy as np


def merge_con_salinity_ameri_gs(
    salinity_folder, ameri_monthly_csv, growing_season_csv, save_folder
):
    salinity_dfs = []

    # Step 1: Read and combine all salinity files
    for file in os.listdir(salinity_folder):
        if file.endswith(".csv"):
            site_name = file.replace(".csv", "")
            file_path = os.path.join(salinity_folder, file)
            try:
                df = pd.read_csv(file_path)
                df["DateTime"] = pd.to_datetime(df["DateTime"], errors="coerce")
                df["Year"] = df["DateTime"].dt.year
                df["month"] = df["DateTime"].dt.month
                df["DoY"] = df["DateTime"].dt.dayofyear
                df["salinity"] = pd.to_numeric(df["salinity"], errors="coerce")
                df.loc[df["salinity"] < 0, "salinity"] = np.nan
                df["site_name"] = site_name
                salinity_dfs.append(df[["site_name", "DateTime", "Year", "month", "DoY", "salinity"]])
                print(f"📄 Processed {file}, columns: {df.columns.tolist()}")
            except Exception as e:
                print(f"❌ Error reading {file}: {e}")

    if not salinity_dfs:
        raise ValueError("❌ No valid salinity CSV files found.")

    salinity_all = pd.concat(salinity_dfs, ignore_index=True)

    # Step 2: Clean site_name whitespace for all dataframes
    salinity_all["site_name"] = salinity_all["site_name"].str.strip()

    gs_df = pd.read_csv(growing_season_csv)
    gs_df["site_name"] = gs_df["site_name"].str.strip()

    ameri_df = pd.read_csv(ameri_monthly_csv)
    ameri_df.rename(columns={"year": "Year"}, inplace=True)
    ameri_df["site_name"] = ameri_df["site_name"].str.strip()

    # Step 3: Compute site-level avg_sos and avg_eos
    site_avg = gs_df.groupby("site_name")[["avg_sos", "avg_eos"]].mean().reset_index()

    # Step 4: Merge salinity with site-level SOS/EOS
    sal_merge = pd.merge(salinity_all, site_avg, on="site_name", how="inner")

    # Step 5: Filter by growing season DoY range
    sal_gs = sal_merge[
        (sal_merge["DoY"] >= sal_merge["avg_sos"]) & (sal_merge["DoY"] <= sal_merge["avg_eos"])
    ]

    # Step 6: Calculate monthly mean salinity within growing season
    monthly_salinity_gs = (
        sal_gs.groupby(["site_name", "Year", "month"])
        .agg(mean_salinity=("salinity", "mean"))
        .reset_index()
    )

    # Step 7: Merge with AmeriFlux
    merged_df = pd.merge(
        monthly_salinity_gs, ameri_df, on=["site_name", "Year", "month"], how="inner"
    )

    # Step 8: Report unmatched salinity sites
    unmatched_sites = set(monthly_salinity_gs["site_name"].unique()) - set(merged_df["site_name"].unique())
    if unmatched_sites:
        print(f"\n❓ Sites in salinity GS data but not in merged output: {unmatched_sites}")

    # Step 9: Save outputs
    merged_path = os.path.join(save_folder, "ameri_salinity_gs_monthly.csv")
    sal_path = os.path.join(save_folder, "monthly_salinity_gs.csv")

    merged_df.to_csv(merged_path, index=False)
    monthly_salinity_gs.to_csv(sal_path, index=False)

    print(f"\n✅ Merged GS salinity-AmeriFlux monthly data saved to: {merged_path}")
    print(f"📁 Monthly growing season salinity data saved to: {sal_path}")

    return merged_df, monthly_salinity_gs



# Import the function (if saved in a separate file, make sure to import it properly)
# from your_script import merge_con_salinity_ameri_gs

salinity_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\continous_salinity\salinity"
ameri_monthly_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
growing_season_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"

# Run the function
merged_df, monthly_gs_salinity = merge_con_salinity_ameri_gs(
    salinity_folder,
    ameri_monthly_csv,
    growing_season_csv,
    save_folder
)
###################################################################################################
###salinity versus ET, GPP, WUE, CUE graph, normalized data
## dont change this code. first see what happens with normalization by site 

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import linregress

def plot_salinity_relationships(input_path_csv, save_folder):
    df = pd.read_csv(input_path_csv)

    # Clip extreme low CUE values
    if 'CUE' in df.columns:
        df['CUE'] = df['CUE'].clip(lower=-0.5)

    target_vars = ['ET', 'GPP', 'WUE', 'CUE']

    ylabel_map = {
        'ET': 'ET (mm)',
        'GPP': 'GPP (g C m⁻² month⁻¹)',
        'WUE': 'WUE (g C per Kg H₂O)',
        'CUE': 'CUE (g C m⁻²)'
    }

    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'H', 'X']
    palette = sns.color_palette("tab10", n_colors=len(df['site_name'].unique()))
    site_names = sorted(df['site_name'].dropna().unique())
    site_styles = {site: (markers[i % len(markers)], palette[i % len(palette)]) for i, site in enumerate(site_names)}

    # Normalize ET and GPP by site to [-1, 1]
    df_norm = df.copy()
    for target in ['ET', 'GPP']:
        if target in df.columns:
            def min_max_norm(x):
                if x.max() == x.min():
                    return pd.Series([0] * len(x), index=x.index)
                return 2 * (x - x.min()) / (x.max() - x.min()) - 1
            df_norm[target] = df.groupby('site_name')[target].transform(min_max_norm)
    
    # Retain WUE and CUE without normalization
    for target in ['WUE', 'CUE']:
        if target in df.columns:
            df_norm[target] = df[target]

    for target in target_vars:
        if target not in df.columns:
            print(f"⚠️ Skipping {target} — not found in data.")
            continue

        df_clean = df_norm.dropna(subset=['mean_salinity', target])
        df_clean = df_clean[df_clean['mean_salinity'] >= 0]

        if df_clean.empty:
            print(f"❌ No valid data with non-negative salinity and {target} found.")
            continue

        for site in site_names:
            site_df = df_clean[df_clean['site_name'] == site]
            if len(site_df) < 2:
                continue

            x = site_df['mean_salinity'].values.reshape(-1, 1)
            y = site_df[target].values

            model = LinearRegression()
            model.fit(x, y)
            y_pred = model.predict(x)
            slope, intercept, r_value, p_value, std_err = linregress(site_df['mean_salinity'], site_df[target])

            color = 'royalblue' if slope > 0 else 'red'

            fig, ax = plt.subplots(figsize=(7, 5))
            marker, _ = site_styles[site]
            sns.scatterplot(x=site_df['mean_salinity'], y=site_df[target], ax=ax, s=120, color=color, marker=marker, label=site)
            ax.plot(site_df['mean_salinity'], y_pred, linestyle='--', linewidth=2.5, color=color)
            ax.set_xlabel("Salinity (ppt)", fontsize=18)
            ax.set_ylabel(ylabel_map[target], fontsize=18)
            title = f"{site} | {target} | Slope = {slope:.2f}, $R^2$ = {r2_score(y, y_pred):.2f}, p = {'< 0.01' if p_value < 0.01 else f'{p_value:.2f}'}"
            ax.set_title(title, fontdict={'fontsize': 16})
            ax.tick_params(axis='both', labelsize=16)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.legend(fontsize=22)
            fig.tight_layout()

            filename = f"{site}_{target}_salinity_plot.png"
            fig.savefig(os.path.join(save_folder, filename), dpi=300)
            plt.close()

        # Combined plot
        if len(df_clean) > 1:
            X_all = df_clean['mean_salinity'].values.reshape(-1, 1)
            y_all = df_clean[target].values
            model_all = LinearRegression()
            model_all.fit(X_all, y_all)
            y_all_pred = model_all.predict(X_all)
            r2_all = r2_score(y_all, y_all_pred)
            slope, intercept, r_value, p_value, std_err = linregress(df_clean['mean_salinity'], df_clean[target])

            fig, ax = plt.subplots(figsize=(7, 5))
            for site in site_names:
                site_df = df_clean[df_clean['site_name'] == site]
                marker, _ = site_styles[site]
                color = 'royalblue' if slope > 0 else 'red'
                ax.scatter(site_df['mean_salinity'], site_df[target], s=110, label=site, marker=marker, color=color)

            ax.plot(df_clean['mean_salinity'], y_all_pred, linestyle='--', linewidth=2.5, color='black')
            ax.set_xlabel("Salinity (ppt)", fontsize=18)
            ax.set_ylabel(ylabel_map[target], fontsize=18)
            ax.set_title(f"{target} | Slope = {slope:.2f}, $R^2$ = {r2_all:.2f}, p = {'< 0.01' if p_value < 0.01 else f'{p_value:.2f}'}", fontsize=18)
            ax.tick_params(axis='both', labelsize=16)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.legend(fontsize=14, title=None, loc='best', ncol=2)
            fig.tight_layout()

            filename_combined = f"combined_{target}_salinity_plot.png"
            fig.savefig(os.path.join(save_folder, filename_combined), dpi=300)
            plt.close()

# Example call
input_path_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\ameri_salinity_gs_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\continous_salinity"
plot_salinity_relationships(input_path_csv, save_folder)


#################################################################################################

## GPP plot to see how phenofit did 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_gpp_by_site(input_path_csv):
    df = pd.read_csv(input_path_csv)

    # Ensure necessary columns
    if not {'site_name', 'Year', 'month', 'GPP'}.issubset(df.columns):
        raise ValueError("The CSV must have 'site_name', 'Year', 'month', 'GPP'.")

    # Create datetime from Year and month
    df['date'] = pd.to_datetime(df[['Year', 'month']].assign(day=1))

    # Sort for safety
    df = df.sort_values(['site_name', 'date'])
    
    sns.set(style="whitegrid")
    site_names = sorted(df['site_name'].dropna().unique())

    for site in site_names:
        site_df = df[df['site_name'] == site][['date', 'GPP']].sort_values('date').reset_index(drop=True)

        # Identify gaps in time (where difference > 32 days)
        site_df['gap'] = site_df['date'].diff().dt.days > 32
        site_df['segment'] = site_df['gap'].cumsum()

        # Plot each continuous segment separately
        plt.figure(figsize=(10, 5))
        for _, segment_df in site_df.groupby('segment'):
            plt.plot(segment_df['date'], segment_df['GPP'], marker='o', linestyle='-', color='green')

        plt.title(f"GPP Over Time – {site}", fontsize=16)
        plt.xlabel("Date", fontsize=14)
        plt.ylabel("GPP (g C m⁻² month⁻¹)", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.show()



###############################################################################################
###salinity versus ET, GPP, WUE, CUE graph, not normalized data

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import linregress

def plot_salinity_relationships(input_path_csv, save_folder):
    df = pd.read_csv(input_path_csv)

    # Clip extreme low CUE values
    if 'CUE' in df.columns:
        df['CUE'] = df['CUE'].clip(lower=-0.5)

    target_vars = ['ET', 'GPP', 'WUE', 'CUE']

    ylabel_map = {
        'ET': 'ET (mm)',
        'GPP': 'GPP (g C m⁻² month⁻¹)',
        'WUE': 'WUE (g C per Kg H₂O)',
        'CUE': 'CUE (g C m⁻²)'
    }

    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'H', 'X']
    palette = sns.color_palette("tab10", n_colors=len(df['site_name'].unique()))
    site_names = sorted(df['site_name'].dropna().unique())
    site_styles = {site: (markers[i % len(markers)], palette[i % len(palette)]) for i, site in enumerate(site_names)}

    for target in target_vars:
        if target not in df.columns:
            print(f"⚠️ Skipping {target} — not found in data.")
            continue

        df_clean = df.dropna(subset=['salinity', target])
        df_clean = df_clean[df_clean['salinity'] >= 0]

        if df_clean.empty:
            print(f"❌ No valid data with non-negative salinity and {target} found.")
            continue

        for site in site_names:
            site_df = df_clean[df_clean['site_name'] == site]
            if len(site_df) < 2:
                continue

            x = site_df['salinity'].values.reshape(-1, 1)
            y = site_df[target].values

            model = LinearRegression()
            model.fit(x, y)
            y_pred = model.predict(x)
            slope, intercept, r_value, p_value, std_err = linregress(site_df['salinity'], site_df[target])

            # Color based on slope (royal blue for positive slope, red for negative slope)
            color = 'royalblue' if slope > 0 else 'red'

            fig, ax = plt.subplots(figsize=(7, 5))
            marker, _ = site_styles[site]
            sns.scatterplot(x=site_df['salinity'], y=site_df[target], ax=ax, s=120, color=color, marker=marker, label=site)
            ax.plot(site_df['salinity'], y_pred, linestyle='--', linewidth=2.5, color=color)
            ax.set_xlabel("Salinity (ppt)", fontsize=18)
            ax.set_ylabel(ylabel_map[target], fontsize=18)
            title = f"{site} | {target} | $R^2$ = {r2_score(y, y_pred):.2f}, p = {'< 0.01' if p_value < 0.01 else f'{p_value:.2f}'}"
            ax.set_title(title, fontsize=18)
            ax.tick_params(axis='both', labelsize=16)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.legend(fontsize=22)
            fig.tight_layout()

            filename = f"{site}_{target}_salinity_plot.png"
            fig.savefig(os.path.join(save_folder, filename), dpi=300)
            plt.close()

        # Combined plot
        if len(df_clean) > 1:
            X_all = df_clean['salinity'].values.reshape(-1, 1)
            y_all = df_clean[target].values
            model_all = LinearRegression()
            model_all.fit(X_all, y_all)
            y_all_pred = model_all.predict(X_all)
            r2_all = r2_score(y_all, y_all_pred)
            slope, intercept, r_value, p_value, std_err = linregress(df_clean['salinity'], df_clean[target])

            fig, ax = plt.subplots(figsize=(7, 5))
            for site in site_names:
                site_df = df_clean[df_clean['site_name'] == site]
                marker, _ = site_styles[site]
                # Use slope to decide the color (royal blue for positive slope, red for negative slope)
                color = 'royalblue' if slope > 0 else 'red'
                ax.scatter(site_df['salinity'], site_df[target], s=110, label=site, marker=marker, color=color)

            ax.plot(df_clean['salinity'], y_all_pred, linestyle='--', linewidth=2.5, color='black')
            ax.set_xlabel("Salinity (ppt)", fontsize=18)
            ax.set_ylabel(ylabel_map[target], fontsize=18)
            ax.set_title(f"{target} | $R^2$ = {r2_all:.2f}, p = {'< 0.01' if p_value < 0.01 else f'{p_value:.2f}'}", fontsize=22)
            ax.tick_params(axis='both', labelsize=16)
            ax.spines['top'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)
            ax.legend(fontsize=14, title=None, loc='best', ncol=2)
            fig.tight_layout()

            filename_combined = f"combined_{target}_salinity_plot.png"
            fig.savefig(os.path.join(save_folder, filename_combined), dpi=300)
            plt.close()

# Example call
input_path_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\ameri_salinity_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\continous_salinity"
plot_salinity_relationships(input_path_csv, save_folder)

#################################################################################################only chnage color baased on slope, no need to display slope

plot_normalized_salinity_relationships(input_path_csv, save_folder)

use slope to decide color. calculate linear fit slope. if positive slope use 
blue color, if negative slope use red color. dont change anything about labl
############################################################################################

## box plot 

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def boxplot_salinity(input_path_csv, save_folder):
    df = pd.read_csv(input_path_csv)
    os.makedirs(save_folder, exist_ok=True)

    var_labels = {
        'ET': 'ET (mm)',
        'GPP': 'GPP (g C m⁻² month⁻¹)',
        'WUE': 'WUE (g C per Kg H₂O)',
        'CUE': 'CUE (g C m⁻² month⁻¹)'
    }

    # Order sites by mean WUE
    wue_means = df.groupby('site_name')['WUE'].mean().sort_values()
    site_order = wue_means.index.tolist()

    # Color-blind friendly color palette
    base_palette = sns.color_palette("colorblind", len(site_order))
    site_colors = dict(zip(site_order, base_palette))

    # Marker styles for data points
    marker_styles = ['o', 's', 'D', '^', 'v', '>', '<', 'P', 'X', '*', 'H', '8', 'p']
    marker_dict = dict(zip(site_order, marker_styles * ((len(site_order) // len(marker_styles)) + 1)))

    for var in ['WUE', 'ET', 'GPP', 'CUE']:
        df_clean = df[['site_name', var, 'salinity_ppt']].dropna()
        if var == 'CUE':
            df_clean.loc[df_clean['CUE'] < -0.5, 'CUE'] = -0.5

        df_clean = df_clean[df_clean['site_name'].isin(site_order)]

        fig, ax = plt.subplots(figsize=(12, 6))

        # Draw boxplot
        box_colors = [sns.set_hls_values(site_colors[site], l=0.6) for site in site_order]
        sns.boxplot(
            data=df_clean,
            x='site_name',
            y=var,
            order=site_order,
            linewidth=2,
            fliersize=0,
            palette=box_colors,
            boxprops=dict(edgecolor='black'),
            medianprops=dict(color='black'),
            whiskerprops=dict(color='black'),
            capprops=dict(color='black'),
            ax=ax
        )

        # Overlay data points
        for i, site in enumerate(site_order):
            site_df = df_clean[df_clean['site_name'] == site]
            ax.scatter(
                x=[i] * len(site_df),
                y=site_df[var],
                color=site_colors[site],
                marker=marker_dict[site],
                s=60,
                edgecolor='black',
                alpha=0.8
            )

        # Add mean salinity text (colored)
        y_max = df_clean[var].max()
        y_range = y_max - df_clean[var].min()
        text_y = y_max + 0.08 * y_range

        for i, site in enumerate(site_order):
            mean_sal = df[df['site_name'] == site]['salinity_ppt'].mean()
            ax.text(
                i,
                text_y,
                f"Salinity = {mean_sal:.0f} ppt",
                ha='center',
                va='bottom',
                fontsize=14,
                color=site_colors[site]
            )

        # Color x-tick labels to match boxplot color
        ax.set_xticks(range(len(site_order)))
        ax.set_xticklabels(site_order, fontsize=16, rotation=45)
        for label, site in zip(ax.get_xticklabels(), site_order):
            label.set_color(site_colors[site])

        ax.set_xlabel("")  # Remove "Site"
        ax.set_ylabel(var_labels[var], fontsize=20)
        ax.tick_params(axis='y', labelsize=16)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

        plt.tight_layout()
        filename = os.path.join(save_folder, f"boxplot_{var}.png")
        fig.savefig(filename, dpi=300)
        plt.close()



input_path_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\ameri_salinity_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\box_plot_salinity"

boxplot_salinity(input_path_csv, save_folder)

###################################################################################################
### previous analysis
###############################################################################################


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pymannkendall import original_test
from sklearn.linear_model import LinearRegression

def plot_salinity_trends(csv_folder, save_folder):
    os.makedirs(save_folder, exist_ok=True)
    csv_files = [f for f in os.listdir(csv_folder) if f.endswith('.csv')]

    # Updated color-blind-friendly palette (CUD)
    color_palette = ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D']
    marker_styles = ['o', 's', '^', 'D', 'v', 'P', 'X']  # Different marker shapes

    plt.figure(figsize=(12, 8))

    for idx, file in enumerate(csv_files):
        file_path = os.path.join(csv_folder, file)
        site_df = pd.read_csv(file_path)
        site_name = os.path.splitext(file)[0]

        if 'DateTime' not in site_df.columns or 'salinity' not in site_df.columns:
            continue

        site_df['DateTime'] = pd.to_datetime(site_df['DateTime'], errors='coerce')
        site_df = site_df.dropna(subset=['DateTime', 'salinity'])
        site_df = site_df[site_df['salinity'] >= 0]
        site_df['Year'] = site_df['DateTime'].dt.year
        site_df['Month'] = site_df['DateTime'].dt.month
        seasonal_df = site_df[site_df['Month'].isin([4, 5, 6, 7, 8, 9])]

        if seasonal_df.empty:
            continue

        monthly_salinity = seasonal_df.groupby('Year')['salinity'].mean().reset_index()
        if len(monthly_salinity) < 2:
            continue

        X = monthly_salinity['Year'].values.reshape(-1, 1)
        y = monthly_salinity['salinity'].values
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)
        slope = model.coef_[0]

        mk_result = original_test(y)
        p_val = mk_result.p

        color = color_palette[idx % len(color_palette)]
        marker = marker_styles[idx % len(marker_styles)]

        bars = plt.bar(monthly_salinity['Year'], y, color=color, width=0.6, alpha=0.6,
                       label=f"{site_name} (p={p_val:.3f}, slope={slope:.2f})")

        for bar in bars:
            bar_x = bar.get_x() + bar.get_width() / 2
            bar_y = bar.get_height()

            # Value label in same color
            plt.text(bar_x, bar_y, f'{int(bar_y)}', ha='center', va='bottom',
                     color=color, fontsize=12)

            # Marker on top of the bar
            plt.plot(bar_x, bar_y, marker=marker, color=color, markersize=8)

        plt.plot(monthly_salinity['Year'], y_pred, linestyle='--', color=color, linewidth=2)

    plt.xlabel("Year", fontsize=18)
    plt.ylabel("Mean Salinity (ppt)", fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=14, loc='best', frameon=True)
    plt.grid(True)
    plt.title("Seasonal (Apr–Sep) Mean Salinity Trend by Site", fontsize=20)
    plt.tight_layout(pad=3.0)

    plot_path = os.path.join(save_folder, "salinity_trend_all_sites.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()


csv_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\continous_salinity\salinity\significant_salinity"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\salinity_trends"

plot_salinity_trends(csv_folder, save_folder)



#################################################################################################

### max wUE, CUE only for year with >=4 years for data with months 5,6,7,8



import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from pymannkendall import original_test

def plot_max_wue_cue_trends(csv_file, save_folder):
    df = pd.read_csv(csv_file)
    df['CUE'] = df['CUE'].clip(lower=-0.5)

    required_cols = ['site_name', 'Year', 'month', 'WUE', 'CUE', 'salinity_ppt']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns: {required_cols}")

    os.makedirs(save_folder, exist_ok=True)

    for site, site_df in df.groupby("site_name"):
        seasonal = site_df[site_df['month'].isin([5, 6, 7, 8])]
        valid_years = (
            seasonal.groupby('Year')['month']
            .nunique()
            .loc[lambda x: x == 4]
            .index
        )

        if len(valid_years) < 4:
            continue

        filtered = seasonal[seasonal['Year'].isin(valid_years)]
        yearly_max = (
            filtered.groupby('Year')[['WUE', 'CUE']]
            .max()
            .reset_index()
            .sort_values('Year')
        )

        years = yearly_max['Year'].values
        x = years.reshape(-1, 1)

        wue_y = yearly_max['WUE'].values
        wue_model = LinearRegression().fit(x, wue_y)
        wue_pred = wue_model.predict(x)
        wue_slope = wue_model.coef_[0]
        wue_p = original_test(wue_y).p

        cue_y = yearly_max['CUE'].values
        cue_model = LinearRegression().fit(x, cue_y)
        cue_pred = cue_model.predict(x)
        cue_slope = cue_model.coef_[0]
        cue_p = original_test(cue_y).p

        if wue_p >= 0.05 and cue_p >= 0.05:
            continue

        max_salinity = filtered['salinity_ppt'].max()

        fig, ax1 = plt.subplots(figsize=(10, 6))

        color_wue = '#CC79A7'  # Pink/Purple
        color_cue = '#56B4E9'  # Sky Blue
        marker_wue = 'D'       # Diamond
        marker_cue = '^'       # Triangle

        ax1.set_xlabel("Year", fontsize=15, fontweight='bold')
        ax1.set_ylabel("Max WUE (g C kg⁻¹ H₂O)", color=color_wue, fontsize=15, fontweight='bold')

        wue_line, = ax1.plot(
            years, wue_y,
            marker=marker_wue,
            linestyle='-',
            color=color_wue,
            linewidth=3,
            markersize=9,
            markeredgewidth=2,
            label='Max WUE'
        )
        ax1.plot(years, wue_pred, linestyle=':', color=color_wue, linewidth=3)
        ax1.tick_params(axis='y', labelcolor=color_wue)
        ax1.tick_params(axis='both', labelsize=13)
        ax1.grid(True, linestyle='--', alpha=0.4)
        ax1.yaxis.set_major_locator(MaxNLocator(nbins='auto', prune='lower'))

        ax2 = ax1.twinx()
        ax2.set_ylabel("Max CUE (g C m⁻²)", color=color_cue, fontsize=15, fontweight='bold')
        cue_line, = ax2.plot(
            years, cue_y,
            marker=marker_cue,
            linestyle='-',
            color=color_cue,
            linewidth=3,
            markersize=9,
            markeredgewidth=2,
            label='Max CUE'
        )
        ax2.plot(years, cue_pred, linestyle=':', color=color_cue, linewidth=3)
        ax2.tick_params(axis='y', labelcolor=color_cue)
        ax2.tick_params(axis='both', labelsize=13)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins='auto', prune='lower'))

        plt.title(f"{site} - Max WUE and CUE Trend\nSalinity = {max_salinity:.2f} ppt",
                  fontsize=15, fontweight='bold')

        # Custom color-coded legend at the bottom right
        legend_handles = [
            mpatches.Patch(color=color_wue, label='Max WUE'),
            mpatches.Patch(color=color_cue, label='Max CUE')
        ]
        ax1.legend(handles=legend_handles, loc='lower right', fontsize=12, frameon=True)

        annotation_text = ""
        if wue_p < 0.05:
            annotation_text += f"WUE: slope = {wue_slope:.2g}, p = {wue_p:.3f}\n"
        if cue_p < 0.05:
            annotation_text += f"CUE: slope = {cue_slope:.2g}, p = {cue_p:.3f}"
        ax1.text(0.01, 0.97, annotation_text.strip(),
                 transform=ax1.transAxes, fontsize=11, va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        plot_path = os.path.join(save_folder, f"{site}_max_WUE_CUE_trend.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()


# Example usage
csv_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\long_term_trend"

plot_max_wue_cue_trends(csv_file, save_folder)




#################################################################################################

## same function now instead of mankendall only see if linear fit trend, include p value 


x label should not be more than two significan digits, remove years if one of the following months  5,6,7,8 is missing in a given year




csv_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
save_folder= \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\long_term_trend

###########################################################################################


import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import matplotlib.colors as mcolors

def linear_long_term(csv_file, save_folder, significance_level=0.05):
    df = pd.read_csv(csv_file)

    # Ensure required columns
    required_cols = ['site_name', 'Year', 'month', 'ET', 'GPP', 'NEP', 'salini_coarse']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns: {required_cols}")

    os.makedirs(save_folder, exist_ok=True)

    # Y-axis variable definitions
    variables = {
        'WUE': {'label': 'WUE (g C per Kg H₂O)'},
        'CUE': {'label': 'CUE (g C m⁻²)'},
        'ET': {'label': 'ET (mm)'},
        'GPP': {'label': 'GPP (g C m⁻²)'}
    }

    colors = list(mcolors.TABLEAU_COLORS.values())
    markers = ['o', 'D', '^', 's', 'p', '*', 'H', 'X']
    color_map, marker_map = {}, {}
    c_idx, m_idx = 0, 0

    plot_data = {var: [] for var in variables}

    # Filter sites with 4+ full years (May–Aug)
    for site, site_df in df.groupby("site_name"):
        seasonal = site_df[site_df['month'].isin([5, 6, 7, 8])]
        year_counts = seasonal.groupby('Year')['month'].nunique()
        full_years = year_counts[year_counts == 4].index

        if len(full_years) < 4:
            continue  # skip site

        if site not in color_map:
            color_map[site] = colors[c_idx % len(colors)]
            marker_map[site] = markers[m_idx % len(markers)]
            c_idx += 1
            m_idx += 1

        # Aggregate by year
        yearly = (
            seasonal[seasonal['Year'].isin(full_years)]
            .groupby('Year')[['GPP', 'ET', 'NEP']]
            .sum()
            .reset_index()
            .sort_values('Year')
        )

        yearly['WUE'] = yearly['GPP'] / yearly['ET']
        yearly['CUE'] = yearly['NEP'] / yearly['GPP']
        salinity_category = site_df['salini_coarse'].iloc[0]

        for var in variables:
            # Exclude missing years from regression
            x = yearly['Year']
            y = yearly[var]

            # Exclude NaN years (years with missing data for the variable)
            valid_data = ~y.isna()
            x_valid = x[valid_data]
            y_valid = y[valid_data]

            if len(x_valid) >= 4:  # Need at least 4 valid years
                slope, intercept, r, p, std = stats.linregress(x_valid, y_valid)

                if p < significance_level:
                    plot_data[var].append({
                        'x': x_valid, 'y': y_valid,
                        'slope': slope, 'intercept': intercept, 'p': p,
                        'site': site,
                        'salinity': salinity_category,
                        'color': color_map[site],
                        'marker': marker_map[site]
                    })

    # Plotting
    for var, entries in plot_data.items():
        if not entries:
            continue

        plt.figure(figsize=(8, 5))
        for item in entries:
            label = f"{item['site']} (p={item['p']:.3f}, {item['salinity']})"
            plt.scatter(
                item['x'], item['y'],
                color=item['color'], marker=item['marker'],
                edgecolors='black', s=120, linewidths=1.5,
                label=label
            )
            x_vals = np.array(item['x'])
            y_fit = item['slope'] * x_vals + item['intercept']
            plt.plot(x_vals, y_fit, color=item['color'], linewidth=3, linestyle='--')

        plt.xlabel("Year", fontsize=16, fontweight='bold')
        plt.ylabel(variables[var]['label'], fontsize=16, fontweight='bold')
        plt.title(f"{var} Dynamics (May–Aug)", fontsize=18, fontweight='bold')

        # Reduce x-axis ticks
        all_years = df['Year'].unique()
        tick_spacing = 2 if len(all_years) <= 15 else 5
        xticks = np.arange(df['Year'].min(), df['Year'].max() + 1, tick_spacing)
        plt.xticks(xticks, fontsize=12)
        plt.yticks(fontsize=12)

        # Legend
        plt.legend(loc='best', fontsize=11, title_fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, f"{var}_dynamics.png"), dpi=400)
        plt.close()

# Example usage
csv_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\long_term_dynamics"
linear_long_term(csv_file, save_folder)


################################################################################################

## site map



import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

def plot_sites_from_info(info_csv, save_folder):
    # Load site info CSV
    info_df = pd.read_csv(info_csv)
    
    # Filter and clean
    df = info_df[['site_name', 'lat', 'long', 'salinity_con']].dropna()
    
    # Colors
    circle_color = '#2ca02c'  # Green for circles (salinity_con == "yes")
    triangle_color = '#9467bd'  # Purple for triangles (salinity_con == "no")
    
    # Marker size (same for both circle and triangle)
    marker_size = 12

    # Compute map extent based on site lat/lon
    min_lat = df['lat'].min() - 2
    max_lat = df['lat'].max() + 2
    min_lon = df['long'].min() - 2
    max_lon = df['long'].max() + 2

    # Plot setup
    fig = plt.figure(figsize=(16, 12))
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#f0e4d7')  # Light beige land color
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='#cce5ff')  # Softer blue water color
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor='#aad3df')  # Softer lake blue
    ax.add_feature(cfeature.BORDERS, linewidth=1.2, edgecolor='gray')  # Gray borders for countries
    ax.add_feature(cfeature.COASTLINE, linewidth=1)  # Coastlines in black
    ax.add_feature(cfeature.STATES, linewidth=0.8, edgecolor='lightgray')  # Light gray state boundaries

    # Add state boundaries and label once per state
    state_boundaries = shpreader.Reader(
        shpreader.natural_earth(category='cultural', name='admin_1_states_provinces')
    )
    shown_states = set()
    for _, row in df.iterrows():
        point = Point(row['long'], row['lat'])
        for state in state_boundaries.records():
            if state.geometry.contains(point):
                ax.add_geometries([state.geometry], crs=ccrs.PlateCarree(), edgecolor='lightgray', facecolor='none')
                state_name = state.attributes['name']
                if state_name not in shown_states:
                    ax.text(row['long'] + 0.5, row['lat'] + 0.5, state_name,
                            transform=ccrs.PlateCarree(), fontsize=10, fontweight='bold', color='black')
                    shown_states.add(state_name)
                break

    # Plot site points
    for _, row in df.iterrows():
        lat, lon = row['lat'], row['long']
        is_salinity = str(row['salinity_con']).strip().lower() == "yes"
        
        if is_salinity:
            # Plot circle marker for salinity_con == "yes"
            ax.plot(lon, lat, marker='o', color=circle_color, markersize=marker_size, markeredgewidth=3, transform=ccrs.PlateCarree())
        else:
            # Plot triangle marker for salinity_con == "no"
            ax.plot(lon, lat, marker='^', color=triangle_color, markersize=marker_size, transform=ccrs.PlateCarree())

    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    circle_handle = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=circle_color, markersize=marker_size, label="Salinity = Yes")
    triangle_handle = plt.Line2D([0], [0], marker='^', color='w', markerfacecolor=triangle_color, markersize=marker_size, label="Salinity = No")
    ax.legend(handles=[circle_handle, triangle_handle], loc='lower right', fontsize=14, title="Salinity Status", title_fontsize=16)

    plt.title("Site Locations with Salinity Info", fontsize=20, fontweight='bold')
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(save_folder, 'site_map_with_salinity_info.png')
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Map saved to: {output_path}")

# Example usage
info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_salinity.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\maps"

plot_sites_from_info(info_csv, save_folder)







#######################################################################################33

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from scipy.stats import linregress
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from scipy.stats import linregress
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

def plot_slope_maps(csv_file, save_folder):
    # Load and filter data
    df = pd.read_csv(csv_file)
    keep_cols = ['site_name', 'ET', 'GPP', 'WUE', 'lat', 'long']
    df = df[keep_cols].dropna()

    # Calculate slopes
    results = []
    for site in df['site_name'].unique():
        sub = df[df['site_name'] == site]
        if len(sub) < 2:
            continue
        try:
            et_slope, _, _, _, _ = linregress(sub['ET'], sub['WUE'])
            gpp_slope, _, _, _, _ = linregress(sub['GPP'], sub['WUE'])
            lat = sub['lat'].iloc[0]
            lon = sub['long'].iloc[0]
            results.append({
                'site_name': site,
                'ET_slope': et_slope,
                'GPP_slope': gpp_slope,
                'lat': lat,
                'lon': lon
            })
        except:
            continue

    slope_df = pd.DataFrame(results)

    # Extend map to include Florida more clearly by increasing the southern buffer
    buffer = 4  # Increased buffer to ensure Florida's bottom coast is fully visible
    min_lat = slope_df['lat'].min() - buffer
    max_lat = slope_df['lat'].max() + buffer
    min_lon = slope_df['lon'].min() - buffer
    max_lon = slope_df['lon'].max() + buffer

    # Map setup function
    def setup_map(title):
        fig = plt.figure(figsize=(14, 10))
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.Geodetic())

        # Updated map colors
        ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#f4e1c1')  # Lighter beige land color
        ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='#a3d3f5')  # Lighter water color (light blue)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor='#80deea')  # Light lake blue
        ax.add_feature(cfeature.BORDERS, linewidth=1.2, edgecolor='gray')  # Neutral gray borders
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        ax.add_feature(cfeature.STATES, linewidth=0.8, edgecolor='#e0e0e0', facecolor='#f5f5f5')  # Lighter state boundaries

        # State boundaries and labels for only involved states
        state_boundaries = shpreader.Reader(
            shpreader.natural_earth(category='cultural', name='admin_1_states_provinces')
        )
        shown_states = set()
        for _, row in slope_df.iterrows():
            for state in state_boundaries.records():
                if state.geometry.contains(Point(row['lon'], row['lat'])):
                    ax.add_geometries([state.geometry], crs=ccrs.PlateCarree(), edgecolor='#e0e0e0', facecolor='#f5f5f5')
                    state_name = state.attributes['name']
                    if state_name not in shown_states:
                        ax.text(row['lon'] + 0.3, row['lat'] + 0.3, state_name,
                                transform=ccrs.Geodetic(), fontsize=10, fontweight='bold', color='black')
                        shown_states.add(state_name)

        gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.4)
        gl.top_labels = gl.right_labels = False
        plt.title(title, fontsize=20, fontweight='bold', color='darkslategray', pad=20)
        return fig, ax

    # ET vs WUE slope map
    fig_et, ax_et = setup_map("ET vs WUE Slope (+/-) at Sites (USA Coasts)")
    for _, row in slope_df.iterrows():
        color = '#1976D2' if row['ET_slope'] > 0 else '#e53935'  # Darker blue for positive, red for negative
        ax_et.text(row['lon'], row['lat'], "+" if row['ET_slope'] > 0 else "-",
                   transform=ccrs.Geodetic(), fontsize=40, fontweight='bold', color=color, ha='center', va='center')

    et_path = os.path.join(save_folder, 'ET_WUE_Slope_Signs_Map.png')
    plt.tight_layout(pad=0)  # Remove extra white space around the map
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Remove the borders at the four corners
    plt.savefig(et_path, dpi=300)
    plt.close()

    # GPP vs WUE slope map
    fig_gpp, ax_gpp = setup_map("GPP vs WUE Slope (+/-) at Sites (USA Coasts)")
    for _, row in slope_df.iterrows():
        color = '#1976D2' if row['GPP_slope'] > 0 else '#e53935'  # Darker blue for positive, red for negative
        ax_gpp.text(row['lon'], row['lat'], "+" if row['GPP_slope'] > 0 else "-",
                    transform=ccrs.Geodetic(), fontsize=40, fontweight='bold', color=color, ha='center', va='center')

    gpp_path = os.path.join(save_folder, 'GPP_WUE_Slope_Signs_Map.png')
    plt.tight_layout(pad=0)  # Remove extra white space around the map
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Remove the borders at the four corners
    plt.savefig(gpp_path, dpi=300)
    plt.close()

    print(f"Maps saved:\n - {et_path}\n - {gpp_path}")

# Example usage
csv_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly.csv"
save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\plots\maps"

plot_slope_maps(csv_file, save_folder)

##########################################################################3 binned salinity

## axiallary code

###################################################################

import pandas as pd
import os

# Path to the file
folder_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\continous_salinity\orig_salinity"
filename = "list_of_hydro.csv"
full_path = os.path.join(folder_path, filename)


def convert_dms_csv_to_decimal(folder_path, filename):
    """
    Converts LAT/LON in ddmmss.sss format to decimal degrees for a CSV file.

    Parameters:
    - folder_path: str, path to the folder containing the CSV
    - filename: str, name of the CSV file (e.g., 'list_of_hydro.csv')

    Returns:
    - df: pandas DataFrame with added 'LAT_dd' and 'LON_dd' columns
    """

    # Construct full file path
    full_path = os.path.join(folder_path, filename)

    # Helper function to convert ddmmss.sss to decimal degrees
    def dms_to_dd(value):
        try:
            value = float(value)
            degrees = int(value // 10000)
            minutes = int((value % 10000) // 100)
            seconds = value % 100
            return degrees + (minutes / 60) + (seconds / 3600)
        except:
            return None

    # Load CSV
    df = pd.read_csv(full_path)

    # Strip whitespace from column names
    df.columns = [col.strip() for col in df.columns]

    # Convert LAT and LON columns
    if 'LAT (ddmmss.sss)' in df.columns and 'LON (ddmmss.sss)' in df.columns:
        df['LAT_dd'] = df['LAT (ddmmss.sss)'].apply(dms_to_dd)
        df['LON_dd'] = df['LON (ddmmss.sss)'].apply(dms_to_dd)
    else:
        raise KeyError("Expected columns 'LAT (ddmmss.sss)' and 'LON (ddmmss.sss)' not found.")

    # Save to a new file
    output_filename = filename.replace('.csv', '_converted.csv')
    output_path = os.path.join(folder_path, output_filename)
    df.to_csv(output_path, index=False)

    print(f"Converted file saved to: {output_path}")
    return df

df=convert_dms_csv_to_decimal(folder_path, filename)

##################################################################################################
## plot florida ameri sites and salinity sites

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import os

folder_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\continous_salinity\orig_salinity"


def plot_fl_salinity_ameriflux(folder_path, salinity_csv, ameri_csv):
    """
    Plots salinity and AmeriFlux sites on a Florida map using their lat/lon coordinates.

    Parameters:
    - folder_path: str, shared directory where the files are located
    - salinity_csv: str, filename of the salinity data
    - ameri_csv: str, filename of the AmeriFlux data
    """

    # Read CSVs
    salinity_df = pd.read_csv(os.path.join(folder_path, salinity_csv))
    ameri_df = pd.read_csv(os.path.join(folder_path, ameri_csv))

    # Clean column names
    salinity_df.columns = salinity_df.columns.str.strip()
    ameri_df.columns = ameri_df.columns.str.strip()

    # Check for required columns
    for col in ['LAT_dd', 'LON_dd']:
        if col not in salinity_df.columns or col not in ameri_df.columns:
            raise ValueError(f"Column '{col}' must be present in both CSVs with decimal degree coordinates.")

    # Create GeoDataFrames
    salinity_gdf = gpd.GeoDataFrame(salinity_df, geometry=gpd.points_from_xy(salinity_df['LON_dd'], salinity_df['LAT_dd']), crs="EPSG:4326")
    ameri_gdf = gpd.GeoDataFrame(ameri_df, geometry=gpd.points_from_xy(ameri_df['LON_dd'], ameri_df['LAT_dd']), crs="EPSG:4326")

    # Load USA state map and filter Florida
    usa = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    florida = usa[(usa['name'] == 'United States')]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 10))
    florida.to_crs(epsg=4326).plot(ax=ax, color='white', edgecolor='black')

    salinity_gdf.plot(ax=ax, color='blue', markersize=40, label='Salinity Sites')
    ameri_gdf.plot(ax=ax, color='red', markersize=40, label='AmeriFlux Sites')

    plt.title('Florida: Salinity vs AmeriFlux Sites')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    
plot_fl_salinity_ameriflux(folder_path, "salinity_FL_lat_long.csv", "ameri_FL_lat_long.csv")
    
#################################################################################################    


import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point
from adjustText import adjust_text
import numpy as np
import os

def plot_closest_salinity_to_ameriflux(folder_path, salinity_csv, ameri_csv):
    # Load data
    salinity_path = os.path.join(folder_path, salinity_csv)
    ameri_path = os.path.join(folder_path, ameri_csv)

    sal_df = pd.read_csv(salinity_path)
    ameri_df = pd.read_csv(ameri_path)

    print("Sample salinity coords:\n", sal_df[['LAT_dd', 'LON_dd']].head())
    print("Sample AmeriFlux coords:\n", ameri_df[['lat', 'long']].head())

    # Find closest salinity site for each AmeriFlux site
    closest_points = []
    for idx, ameri_row in ameri_df.iterrows():
        # Calculate distance (Euclidean approx) between current ameri site and all salinity sites
        dists = np.sqrt((sal_df['LAT_dd'] - ameri_row['lat'])**2 + (sal_df['LON_dd'] - ameri_row['long'])**2)
        min_idx = dists.idxmin()
        closest_points.append(min_idx)

    closest_salinity_df = sal_df.loc[closest_points].drop_duplicates()

    print("Closest salinity sites to AmeriFlux sites:")
    print(closest_salinity_df[['Station_ID', 'LAT_dd', 'LON_dd']])

    # Plotting
    fig = plt.figure(figsize=(12, 14))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-88, -79, 24, 32], crs=ccrs.PlateCarree())

    # Light color for Florida state background
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#f2f2f2')
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='aliceblue')
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black')
    ax.gridlines(draw_labels=True)

    # Plot closest salinity sites
    ax.scatter(
        closest_salinity_df['LON_dd'], closest_salinity_df['LAT_dd'],
        color='deepskyblue', s=160, marker='o',
        edgecolor='black', linewidth=1.2,
        label='Closest Salinity Sites', zorder=5
    )

    # Plot AmeriFlux sites
    ax.scatter(
        ameri_df['long'], ameri_df['lat'],
        color='orange', s=200, marker='^',
        edgecolor='black', linewidth=1.2,
        label='AmeriFlux Sites', zorder=6
    )

    # Add labels and use adjustText for better placement
    texts = []
    for _, row in closest_salinity_df.iterrows():
        texts.append(ax.text(row['LON_dd'], row['LAT_dd'], str(row.get('Station_ID', 'N/A')),
                             fontsize=9, color='blue', weight='bold', zorder=10))

    for _, row in ameri_df.iterrows():
        texts.append(ax.text(row['long'], row['lat'], row.get('site_name', 'N/A'),
                             fontsize=9, color='darkred', weight='bold', zorder=10))

    adjust_text(texts,
                expand_points=(1.2, 1.2),
                expand_text=(1.2, 1.2),
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.5),
                ax=ax)

    ax.set_title("Closest Salinity Sites to AmeriFlux Sites in Florida", fontsize=18, weight='bold')
    ax.legend(loc='lower left', fontsize=13)

    plt.tight_layout()
    plt.show()

   
 
plot_closest_salinity_to_ameriflux(
    folder_path,
    "salinity_FL_lat_long.csv",
    "ameri_FL_lat_long.csv"
)


