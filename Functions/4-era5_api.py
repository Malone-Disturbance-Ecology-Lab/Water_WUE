# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:30:27 2025

@author: ammara
"""

import cdsapi
import datetime
import os


def fetch_cds_data(area, year_range, month_range, day_range, time_range, output_file):
    # Set the API credentials as environment variables
    os.environ['CDSAPI_URL'] = 'https://cds.climate.copernicus.eu/api'
    os.environ['CDSAPI_KEY'] = '0a2b0eec-ed62-4a42-8fb8-83ca726ef382'  # Replace with your actual CDS API key
    
    # Initialize CDS API client
    client = cdsapi.Client()
    
    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [ # "2m_dewpoint_temperature", # "2m_temperature" # "surface_solar_radiation_downwards"
            
          # 'surface_net_solar_radiation',
           # 'surface_net_thermal_radiation'
             'surface_pressure'        
            
        ],
        #"year": [str(year) for year in range(year_range[0], year_range[1] + 1)], #for mutiple years
        "year": [str(year) for year in (range(year_range[0], year_range[1] + 1) if isinstance(year_range, tuple) else [year_range])],
        "month": [f"{month:02d}" for month in range(month_range[0], month_range[1] + 1)],
        "day": [f"{day:02d}" for day in range(day_range[0], day_range[1] + 1)],
        "time": [f"{hour:02d}:00" for hour in range(time_range[0], time_range[1] + 1)],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": area
    }
    
    client.retrieve(dataset, request).download(target=output_file)
    print(f"Data downloaded successfully to {output_file}")



# two years at a time

# Example usage


fetch_cds_data(
    # if you want excat ameriflux point, male middle points closer to ameriflux point 
    #area=[33.33, -79.25, 33.32, -79.24], # bigger absolute value of lat and long should come first
    area=[33.35,-79.20,	33.34,	-79.19],
    year_range=(2019,2020), # range of years  #2020-2022
    month_range=(1, 12), # range of months
    #month_range=(7, 12), # range of months
    day_range=(1, 31), # range of days
    time_range=(0, 23),
    output_file="output.nc"
)


##########################################################################################################


## ERA data for P and netrad

  #  area=[70.50, -149.90, 70.49,-149.88], # US-A03   2015, 2020
  #  area=[71.33,-156.62,71.32,-156.61], # US-A10   2012, 2020
  #   area=[70.47,-157.41,70.46,-157.40], # US-Atq   2003, 2006
  #   area=[71.28,-156.61,71.27,-156.60], #US-NGB    2013, 2023
  #  area=[38.01,-121.67,38.00,-121.66], #US-Dmg 
   # area=[37.62,-122.12,37.61,-122.11], # US-EDN
   # area=[36.81, -121.76, 36.80, -121.75]# US-EKH
 #   area=[36.86, -121.74, 36.85, -121.75], US-EKP
 #  area=[38.01,-121.77, 38.00, -121.76],US-Myb
 #   area=[38.20, -122.03, 38.19, -122.02], US-Srr
 #   area=[38.11, -121.65, 38.10, -121.64], US-Tw1
 #   area=[39.09, -75.44, 39.08, -75.43],US-StJ
 #   area=[28.71,-80.74,28.70,-80.73], US-KS3
 #area=[28.60, -80.73, 28.59, -80.72], US-KS4
 #area=[25.37,-81.08,	25.35,-81.07], US-Skr
#area=[25.35,-80.38, 25.34,-80.37], US-EvM
 # area=[25.20, -80.64, 25.19, -80.63], US-TaS
 #   area=[29.50, -90.44, 29.49, -90.43], US-LA1
#    area=[29.86, -90.29, 29.85, -90.28], US-LA2
#   area=[29.49,-89.92,29.48, -89.91], US-LA3
 # area=[42.74, -70.83, 42.73, -70.82], US-PhM
#    area=[40.82, -74.04, 40.81, -74.03], US-MRM
   # area=[33.32, -79.24, 33.31, -79.23], US-HB2
#    area=[40.77, -74.09, 40.76, -74.08],  US-HPY
#  area=[33.33, -79.25, 33.32,	-79.24], US-StS
# area=[33.35,-79.20,	33.34,	-79.19],  US-HB1
