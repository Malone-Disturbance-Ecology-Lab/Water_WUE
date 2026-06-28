# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:56:51 2026

@author: ammar
"""
# sites missing day time GPP, US-Elm, US-So2

## this code will check after gap filling what variables do we have
###################################################################################################
import os
import pandas as pd
from pathlib import Path

def check_csv_variables(main_folder, required_variables):
    """
    Check CSV files in the main folder for required variables.
    
    Parameters:
    main_folder: str, path to the main folder
    required_variables: list, list of variable names to check
    """
    
    # Convert to Path object for easier handling
    main_path = Path(main_folder)
    
    # Check if the folder exists
    if not main_path.exists():
        print(f"Error: Folder '{main_folder}' does not exist.")
        return
    
    # Get all CSV files in the main folder (not subfolders)
    csv_files = [f for f in main_path.iterdir() if f.is_file() and f.suffix.lower() == '.csv']
    
    if not csv_files:
        print(f"No CSV files found in '{main_folder}'.")
        return
    
    # Results storage
    results = []
    
    print(f"Checking {len(csv_files)} CSV file(s) in '{main_folder}'...\n")
    print("-" * 80)
    
    # Check each CSV file
    for csv_file in csv_files:
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Get columns in the file
            file_columns = set(df.columns)
            
            # Check which required variables are present
            present_vars = []
            missing_vars = []
            
            for var in required_variables:
                if var in file_columns:
                    present_vars.append(var)
                else:
                    missing_vars.append(var)
            
            # Store results
            results.append({
                'filename': csv_file.name,
                'present': present_vars,
                'missing': missing_vars,
                'all_present': len(missing_vars) == 0
            })
            
            # Print results for this file
            print(f"\nFile: {csv_file.name}")
            if len(missing_vars) == 0:
                print(f"  ✓ All required variables found: {', '.join(present_vars)}")
            else:
                print(f"  ✗ Missing variables: {', '.join(missing_vars)}")
                if present_vars:
                    print(f"  ✓ Present variables: {', '.join(present_vars)}")
            
        except Exception as e:
            print(f"\nFile: {csv_file.name}")
            print(f"  ✗ Error reading file: {str(e)}")
            results.append({
                'filename': csv_file.name,
                'present': [],
                'missing': required_variables,
                'all_present': False,
                'error': str(e)
            })
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print("=" * 80)
    
    files_with_all_vars = [r for r in results if r.get('all_present', False)]
    files_with_errors = [r for r in results if 'error' in r]
    
    print(f"\nTotal files checked: {len(results)}")
    print(f"Files with all required variables: {len(files_with_all_vars)}")
    print(f"Files missing some variables: {len(results) - len(files_with_all_vars) - len(files_with_errors)}")
    print(f"Files with errors: {len(files_with_errors)}")
    
    if files_with_all_vars:
        print("\nFiles that have all required variables:")
        for r in files_with_all_vars:
            print(f"  - {r['filename']}")
    
    # Optionally, save results to a text file
    save_summary = input("\n\nDo you want to save the summary to a text file? (y/n): ").lower()
    if save_summary == 'y':
        output_file = main_path / 'variable_check_summary.txt'
        with open(output_file, 'w') as f:
            f.write("CSV Variable Check Summary\n")
            f.write("=" * 50 + "\n\n")
            for r in results:
                f.write(f"File: {r['filename']}\n")
                if 'error' in r:
                    f.write(f"  ERROR: {r['error']}\n")
                else:
                    f.write(f"  All required variables present: {r['all_present']}\n")
                    if r['present']:
                        f.write(f"  Present: {', '.join(r['present'])}\n")
                    if r['missing']:
                        f.write(f"  Missing: {', '.join(r['missing'])}\n")
                f.write("\n")
        
        print(f"Summary saved to: {output_file}")

# Main execution
if __name__ == "__main__":
    # Define the main folder path
    main_folder = r"M:\Research\WUE_CUE\ameri_data\ameri_fill"
    
    # Define the required variables
    required_variables = ['GPP_DT', 'GPP_nt', 'Reco_DT', 'Reco_nt', 'NEE_f', 'LE_f']
    
    # Run the check
    check_csv_variables(main_folder, required_variables)
    
    
      
##############################################################################################
## this code will check start and end date for sites

##############################################################################################
import os
import pandas as pd
import glob

# Define the folder path
folder_path = r"M:\Research\WUE_CUE\ameri_data\ameri_fill"

# Find all CSV files in the main folder
csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

# Store results
results = []

# Process each CSV file
for file_path in sorted(csv_files):
    try:
        # Extract site name from filename
        filename = os.path.basename(file_path)
        site_name = filename.replace('.csv', '').replace('_fill', '')
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Convert DateTime column to datetime format
        df['DateTime_parsed'] = pd.to_datetime(df['DateTime'])
        
        # Get start and end dates
        start_date = df['DateTime_parsed'].min()
        end_date = df['DateTime_parsed'].max()
        
        # Store results
        results.append({
            'Site': site_name,
            'Start Date': start_date.strftime('%Y-%m-%d'),
            'End Date': end_date.strftime('%Y-%m-%d')
        })
        
    except Exception as e:
        print(f"Error processing {filename}: {e}")

# Create and display table
if results:
    summary_df = pd.DataFrame(results)
    print("\n" + summary_df.to_string(index=False))
else:
    print("No files processed successfully.")
    
##############################################################################################

### this code will help how many csv files in an input folder

import os
import glob

folder_path = r"M:\Research\WUE_CUE\drivers\ameri_drivers\PET_drought"

csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
print(len(csv_files))