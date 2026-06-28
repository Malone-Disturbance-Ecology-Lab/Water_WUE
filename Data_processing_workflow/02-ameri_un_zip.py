# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024
@author: ammar

Automated version: Processes all zip files in a folder
"""

import os
import zipfile
from pathlib import Path

def unzip_ameriflux_data(zip_file_path, extracted_folder):
    """Extracts specific AmeriFlux files containing 'HH' or 'BASE_HH' to the specified folder."""
    
    # Ensure the extraction folder exists
    os.makedirs(extracted_folder, exist_ok=True)
    
    try:
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            file_list = zip_ref.namelist()
            print(f"\nProcessing: {os.path.basename(zip_file_path)}")
            print(f"Files in ZIP: {file_list[:5]}...")  # Show first 5 files
            
            extracted = False
            for file_name in file_list:
                if 'HH' in file_name or 'BASE_HH' in file_name:
                    try:
                        extracted_file_path = os.path.join(extracted_folder, file_name)
                        
                        # Check if file already exists to avoid duplicates
                        if os.path.exists(extracted_file_path):
                            print(f"  File already exists: {file_name}")
                        else:
                            # Extract file
                            with zip_ref.open(file_name) as source, open(extracted_file_path, "wb") as target:
                                target.write(source.read())
                            print(f"  Extracted: {file_name}")
                            extracted = True
                            
                    except Exception as e:
                        print(f"  Error extracting {file_name}: {e}")
            
            if not extracted:
                print(f"  No matching files found in {os.path.basename(zip_file_path)}")
                
    except Exception as e:
        print(f"Error opening zip file {zip_file_path}: {e}")

def process_all_zip_files(input_folder, output_folder):
    """
    Process all zip files in the input folder and extract relevant files to output folder.
    
    Parameters:
    input_folder: Path to folder containing multiple zip files
    output_folder: Path where extracted files will be saved
    """
    
    # Convert to Path objects for easier handling
    input_path = Path(input_folder)
    output_path = Path(output_folder)
    
    # Check if input folder exists
    if not input_path.exists():
        print(f"Error: Input folder '{input_folder}' does not exist!")
        return
    
    # Create output folder if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all zip files in the input folder
    zip_files = list(input_path.glob("*.zip"))
    
    if not zip_files:
        print(f"No zip files found in {input_folder}")
        return
    
    print(f"Found {len(zip_files)} zip file(s) to process")
    print(f"Input folder: {input_folder}")
    print(f"Output folder: {output_folder}")
    print("-" * 50)
    
    # Process each zip file
    for i, zip_file in enumerate(zip_files, 1):
        print(f"\nProcessing file {i}/{len(zip_files)}")
        unzip_ameriflux_data(str(zip_file), str(output_path))
    
    print("\n" + "=" * 50)
    print(f"Processing complete! Processed {len(zip_files)} zip file(s)")
    print(f"Extracted files are saved in: {output_folder}")

# Example usage
if __name__ == "__main__":
    # Using the paths you provided
    input_folder = r"M:\Research\WUE_CUE\ameri_data"  # Folder containing multiple zip files
    output_folder = r"M:\Research\WUE_CUE\ameri_data\ameri_gaps"  # Output folder
    
    # Process all zip files
    process_all_zip_files(input_folder, output_folder)