# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:22:26 2024

@author: ammar
"""


import os
import datetime
import pandas as pd
import numpy as np
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt
# this code will unzip the amerifux file, name it and save it appropriately in a separate folder




def unzip_ameriflux_data(zip_file_path, extracted_folder):
    """Extracts specific AmeriFlux files containing 'HH' or 'BASE_HH' to the specified folder."""
    
    # Ensure the extraction folder exists
    os.makedirs(extracted_folder, exist_ok=True)
    
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        file_list = zip_ref.namelist()
        print("Files in ZIP:", file_list)  # Debugging step
        
        extracted = False
        for file_name in file_list:
            print(f"Checking file: {file_name}")  # Debugging step
            if 'HH' in file_name or 'BASE_HH' in file_name:
                try:
                    extracted_file_path = os.path.join(extracted_folder, file_name)
                    print(f"Attempting to extract to: {extracted_file_path}")  # Debugging step

                    # Manually extract file to ensure correct handling of network paths
                    with zip_ref.open(file_name) as source, open(extracted_file_path, "wb") as target:
                        target.write(source.read())

                    print(f"Extracted successfully: {extracted_file_path}")
                    extracted = True
                except Exception as e:
                    print(f"Error extracting {file_name}: {e}")
                    return
        
        if not extracted:
            print("No matching files found for extraction.")


 

# Example usage
zip_file_path = r'\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\ameri_data\\AMF_US-A10_BASE-BADM_4-5.zip'
extracted_folder = r'\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\ameri_data\\ameri_gaps'  
  
# call the function and then use the line below 

unzip_ameriflux_data(zip_file_path, extracted_folder)
  
