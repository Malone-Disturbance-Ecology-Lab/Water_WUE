# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 10:54:13 2025

@author: ammar
"""

################################################################################################################
import pandas as pd
import os

### not standardized code

def calculate_driven_components():
    """
    Calculate GPP-driven and ET-driven components, derivatives, and classification columns.
    Saves cleaned data to:
    \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_yearly.csv
    """
    # Define paths
    input_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_yearly.csv"
    output_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_yearly.csv"
    
    try:
        # Load and sort data
        df = pd.read_csv(input_path)
        print(f"Loaded data with {len(df)} rows")
        df = df.sort_values(by=["site_name", "year"])
        
        # Calculate yearly differences
        df["dt_GPP"] = df.groupby("site_name")["GPP"].diff()
        df["dt_ET"] = df.groupby("site_name")["ET"].diff()
        
        # Handle zero ET values
        et_nonzero = df["ET"].replace(0, float('nan'))
        
        # Calculate driven components
        df["GPP_driven"] = (1 / et_nonzero) * df["dt_GPP"]
        df["ET_driven"] = (-df["GPP"] / (et_nonzero**2)) * df["dt_ET"]
        
        # Calculate derivatives
        df["GPP_derivative"] = 1 / et_nonzero
        df["ET_derivative"] = -df["GPP"] / (et_nonzero**2)
        
        # Create cleaned version (remove rows with NaN in dt_GPP/dt_ET)
        cleaned_df = df.dropna(subset=['dt_GPP', 'dt_ET'])
        
        # Add 'driven' column based on magnitude comparison
        cleaned_df['driven'] = np.select(
            [
                cleaned_df['GPP_driven'].abs() > cleaned_df['ET_driven'].abs(),
                cleaned_df['GPP_driven'].abs() < cleaned_df['ET_driven'].abs()
            ],
            ['GPP_driven', 'ET_driven'],
            default='Equal'
        )
        
        # Add 'sensitivity' column based on magnitude comparison
        cleaned_df['sensitivity'] = np.select(
            [
                cleaned_df['GPP_derivative'].abs() > cleaned_df['ET_derivative'].abs(),
                cleaned_df['GPP_derivative'].abs() < cleaned_df['ET_derivative'].abs()
            ],
            ['GPP', 'ET'],
            default='Equal'
        )
        
        # Save cleaned data
        cleaned_df.to_csv(output_path, index=False)
        print("Processing completed successfully!")
        print(f"Cleaned data saved to: {output_path}")
        print(f"Original rows: {len(df)} | Cleaned rows: {len(cleaned_df)}")
        print(f"Rows removed: {len(df) - len(cleaned_df)} (first year per site)")
        
        # Show distribution of new classification columns
        print("\nDriven Component Distribution:")
        print(cleaned_df['driven'].value_counts(dropna=False))
        
        print("\nSensitivity Distribution:")
        print(cleaned_df['sensitivity'].value_counts(dropna=False))
        
        # Return cleaned DataFrame for Spyder console
        return cleaned_df
        
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return None
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None

# Execute and display in Spyder
if __name__ == "__main__":
    cleaned_df = calculate_driven_components()
    
    if cleaned_df is not None:
        print("\nCleaned DataFrame Preview:")
        print(cleaned_df.head())
        
        # Optional: Show summary stats
        print("\nSummary Statistics for New Columns:")
        print(cleaned_df[['dt_GPP', 'dt_ET', 'GPP_driven', 'ET_driven', 
                          'GPP_derivative', 'ET_derivative']].describe())
        
        print("\nThe cleaned DataFrame is available as 'cleaned_df' in Variable Explorer")


#############################################################################################################

## standardized code

def calculate_driven_components():
    """
    Calculate GPP-driven and ET-driven components, derivatives, and classification columns.
    Saves cleaned data to:
    \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_yearly.csv
    """
    # Define paths
    input_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_yearly.csv"
    output_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_yearly.csv"
    
    try:
        # Load and sort data
        df = pd.read_csv(input_path)
        print(f"Loaded data with {len(df)} rows")
        df = df.sort_values(by=["site_name", "year"])
        
        # Calculate yearly differences
        df["dt_GPP"] = df.groupby("site_name")["GPP"].diff()
        df["dt_ET"] = df.groupby("site_name")["ET"].diff()
        
        # Handle zero ET values
        et_nonzero = df["ET"].replace(0, float('nan'))
        
        # Calculate driven components
        df["GPP_driven"] = (1 / et_nonzero) * df["dt_GPP"]
        df["ET_driven"] = (-df["GPP"] / (et_nonzero**2)) * df["dt_ET"]
        
        # Calculate derivatives
        df["GPP_derivative"] = 1 / et_nonzero
        df["ET_derivative"] = -df["GPP"] / (et_nonzero**2)
        
        # Define min-max scaling function for [-1, 1] range
        def min_max_scale(series):
            """Scale a series to [-1, 1] range"""
            non_nan = series.dropna()
            if non_nan.empty:
                return series
            min_val = non_nan.min()
            max_val = non_nan.max()
            
            if min_val == max_val:
                # Handle constant series by setting to 0
                return series.mask(series.notna(), 0)
            else:
                # Scale to [-1, 1]
                return 2 * ((series - min_val) / (max_val - min_val)) - 1
        
        # Columns to standardize per site
        cols_to_scale = [
            'GPP_driven', 'ET_driven', 
            'GPP_derivative', 'ET_derivative'
        ]
        
        # Apply scaling per site for each column
        for col in cols_to_scale:
            df[col] = df.groupby('site_name')[col].transform(min_max_scale)
        
        # Create cleaned version (remove rows with NaN in dt_GPP/dt_ET)
        cleaned_df = df.dropna(subset=['dt_GPP', 'dt_ET'])
        
        # Define tolerance for floating point comparison
        tol = 1e-10  # Very small tolerance for floating point equality
        
        # Add 'driven' column based on magnitude comparison with tolerance
        cleaned_df['driven'] = np.where(
            cleaned_df['GPP_driven'].abs() - cleaned_df['ET_driven'].abs() > tol,
            'GPP_driven',
            np.where(
                cleaned_df['ET_driven'].abs() - cleaned_df['GPP_driven'].abs() > tol,
                'ET_driven',
                # If difference is within tolerance, classify based on raw values
                np.where(
                    cleaned_df['GPP_driven'].abs() > cleaned_df['ET_driven'].abs(),
                    'GPP_driven',
                    'ET_driven'
                )
            )
        )
        
        # Add 'sensitivity' column based on magnitude comparison with tolerance
        cleaned_df['sensitivity'] = np.where(
            cleaned_df['GPP_derivative'].abs() - cleaned_df['ET_derivative'].abs() > tol,
            'GPP',
            np.where(
                cleaned_df['ET_derivative'].abs() - cleaned_df['GPP_derivative'].abs() > tol,
                'ET',
                # If difference is within tolerance, classify based on raw values
                np.where(
                    cleaned_df['GPP_derivative'].abs() > cleaned_df['ET_derivative'].abs(),
                    'GPP',
                    'ET'
                )
            )
        )
        
        # Save cleaned data
        cleaned_df.to_csv(output_path, index=False)
        print("Processing completed successfully!")
        print(f"Cleaned data saved to: {output_path}")
        print(f"Original rows: {len(df)} | Cleaned rows: {len(cleaned_df)}")
        print(f"Rows removed: {len(df) - len(cleaned_df)} (first year per site)")
        
        # Show distribution of new classification columns
        print("\nDriven Component Distribution:")
        print(cleaned_df['driven'].value_counts(dropna=False))
        
        print("\nSensitivity Distribution:")
        print(cleaned_df['sensitivity'].value_counts(dropna=False))
        
        # Return cleaned DataFrame for Spyder console
        return cleaned_df
        
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return None
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None

# Execute and display in Spyder
if __name__ == "__main__":
    cleaned_df = calculate_driven_components()
    
    if cleaned_df is not None:
        print("\nCleaned DataFrame Preview:")
        print(cleaned_df.head())
        
        # Optional: Show summary stats
        print("\nSummary Statistics for New Columns:")
        print(cleaned_df[['dt_GPP', 'dt_ET', 'GPP_driven', 'ET_driven', 
                          'GPP_derivative', 'ET_derivative']].describe())
        
        print("\nThe cleaned DataFrame is available as 'cleaned_df' in Variable Explorer")


###############################################################################################################
import pandas as pd
import numpy as np
import os

def calculate_site_level_drivers():
    """
    Calculate site-level mean drivers and sensitivities from yearly data.
    Saves results to:
    \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_sensitivity.csv
    """
    # Define paths
    input_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_yearly.csv"
    output_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_sensitivity.csv"
    
    try:
        # Load yearly data
        df = pd.read_csv(input_path)
        print(f"Loaded yearly data with {len(df)} rows from {input_path}")
        
        # Calculate WUE if not present (GPP/ET)
        if 'WUE' not in df.columns:
            print("Note: Calculating WUE as GPP/ET")
            df['WUE'] = df['GPP'] / df['ET']
        
        # Columns to average per site
        metrics = ['GPP_derivative', 'ET_derivative', 'WUE', 'GPP_driven', 'ET_driven']
        
        # Calculate mean values per site
        site_df = df.groupby('site_name')[metrics].mean().reset_index()
        
        # Calculate magnitude of mean driven components
        site_df['GPP_driven_mag'] = site_df['GPP_driven'].abs()
        site_df['ET_driven_mag'] = site_df['ET_driven'].abs()
        
        # Classify driven component
        site_df['driven'] = np.select(
            [
                site_df['GPP_driven_mag'] > site_df['ET_driven_mag'],
                site_df['GPP_driven_mag'] < site_df['ET_driven_mag']
            ],
            ['GPP_driven', 'ET_driven'],
            default='Equal'
        )
        
        # Calculate magnitude of derivatives
        site_df['GPP_derivative_mag'] = site_df['GPP_derivative'].abs()
        site_df['ET_derivative_mag'] = site_df['ET_derivative'].abs()
        
        # Classify sensitivity
        site_df['sensitivity'] = np.select(
            [
                site_df['GPP_derivative_mag'] > site_df['ET_derivative_mag'],
                site_df['GPP_derivative_mag'] < site_df['ET_derivative_mag']
            ],
            ['GPP', 'ET'],
            default='Equal'
        )
        
        # Save results
        site_df.to_csv(output_path, index=False)
        print("Processing completed successfully!")
        print(f"Site-level data saved to: {output_path}")
        print(f"Processed {len(site_df)} sites")
        
        # Show classification distribution
        print("\nDriven Component Distribution:")
        print(site_df['driven'].value_counts(dropna=False))
        
        print("\nSensitivity Distribution:")
        print(site_df['sensitivity'].value_counts(dropna=False))
        
        return site_df
        
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return None
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None

# Execute and display in Spyder
if __name__ == "__main__":
    site_level_df = calculate_site_level_drivers()
    
    if site_level_df is not None:
        print("\nSite-level DataFrame Preview:")
        print(site_level_df.head())
        
        print("\nSummary Statistics:")
        print(site_level_df.describe())
        
        # Show sites with different classifications
        print("\nSites with GPP-driven WUE:")
        gpp_driven = site_level_df[site_level_df['driven'] == 'GPP_driven']
        print(gpp_driven[['site_name', 'GPP_driven_mag', 'ET_driven_mag']])
        
        print("\nThe site-level DataFrame is available as 'site_level_df' in Variable Explorer")


#################################################################################################

## function to show graph for sentivity to ET and GPP


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

def create_wue_plots():
    """
    Create and save bar plots of WUE values with annotations for driven and sensitivity.
    Saves plots to:
    \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\
    """
    # Define paths
    input_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\gpp_et_driven_sensitivity.csv"
    save_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products"
    
    try:
        # Load site-level data
        df = pd.read_csv(input_path)
        print(f"Loaded data with {len(df)} sites")

        # Global plot settings
        rcParams['font.size'] = 12
        plt.rc('axes', titlesize=16)
        plt.rc('axes', labelsize=14)
        plt.rc('xtick', labelsize=11)
        plt.rc('ytick', labelsize=11)

        sensitivity_colors = {'GPP': '#1f77b4', 'ET': '#ff7f0e', 'Equal': '#2ca02c'}
        driven_colors = {'GPP_driven': '#1f77b4', 'ET_driven': '#ff7f0e', 'Equal': '#2ca02c'}

        def create_plot(df, sort_column, title, filename, color_map=None):
            fig, ax = plt.subplots(figsize=(16, 10))

            if sort_column == 'WUE':
                df = df.sort_values(sort_column)
            else:
                df = df.sort_values([sort_column, 'WUE'])

            x = np.arange(len(df)) * 1.2  # Increase space between bars
            bar_width = 0.5  # Narrower bars

            # Bar colors
            if color_map:
                colors = df[sort_column].map(color_map)
                bars = ax.bar(x, df['WUE'], width=bar_width, color=colors)
            else:
                bars = ax.bar(x, df['WUE'], width=bar_width, color='skyblue')

            # Annotate bars
            for i, bar in enumerate(bars):
                height = bar.get_height()
                driven = df['driven'].iloc[i].replace("_driven", "")
                sensitivity = df['sensitivity'].iloc[i]

                ax.text(
                    bar.get_x() + bar.get_width() / 2.,
                    height + 0.05,
                    f"D: {driven}",
                    ha='center',
                    va='bottom',
                    fontsize=10,
                    fontweight='bold'
                )
                ax.text(
                    bar.get_x() + bar.get_width() / 2.,
                    height + 0.25,
                    f"S: {sensitivity}",
                    ha='center',
                    va='bottom',
                    fontsize=10
                )

            ax.set_xlabel('Site Name', fontsize=14, weight='bold')
            ax.set_ylabel('WUE (GPP/ET)', fontsize=14, weight='bold')
            ax.set_title(title, fontsize=16, weight='bold')

            ax.set_xticks(x)
            ax.set_xticklabels(df['site_name'], rotation=45, ha='right')

            ax.grid(True, linestyle='--', alpha=0.7)

            if color_map:
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor=color, label=label.replace("_driven", "")) 
                                   for label, color in color_map.items()]
                ax.legend(handles=legend_elements, loc='upper left', title=sort_column)

            plt.tight_layout()
            save_path = os.path.join(save_folder, filename)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved plot: {save_path}")
            plt.close()

        # Create all plots
        create_plot(df.copy(), 'WUE', 'Site WUE Values (Sorted Low to High)', 'wue_sorted.png')
        create_plot(df.copy(), 'sensitivity', 'Site WUE Values (Sorted by Sensitivity)', 'wue_sensitivity_sorted.png', sensitivity_colors)
        create_plot(df.copy(), 'driven', 'Site WUE Values (Sorted by Driven Component)', 'wue_driven_sorted.png', driven_colors)

        print("All plots created successfully!")
        return df

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return None
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None

# Run
if __name__ == "__main__":
    result_df = create_wue_plots()

    if result_df is not None:
        print("\nPreview of the data used for plotting:")
        print(result_df[['site_name', 'WUE', 'driven', 'sensitivity']].head())

        print("\nWUE Statistics:")
        print(result_df['WUE'].describe())

        print("\nSensitivity Distribution:")
        print(result_df['sensitivity'].value_counts())

        print("\nDriven Distribution:")
        print(result_df['driven'].value_counts())

        print("\nThe DataFrame is available as 'result_df' in Variable Explorer.")

######################################################################################################




# problem in GPP driven and ET driven when both zero 






#############################################################################################################










