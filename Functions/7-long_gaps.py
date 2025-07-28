# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 12:36:14 2025

@author: ammar
"""

import os
import pandas as pd
import numpy as np



import os
import pandas as pd

def fill_long_nans_batch(path, save_path):
    os.makedirs(save_path, exist_ok=True)

    csv_files = [f for f in os.listdir(path) if f.endswith('.csv')]

    for file in csv_files:
        try:
            input_path = os.path.join(path, file)
            df = pd.read_csv(input_path, parse_dates=['TIMESTAMP'])

            df_b = df.copy()

            if 'Month' not in df.columns:
                df['Month'] = df['TIMESTAMP'].dt.month
                df_b['Month'] = df['TIMESTAMP'].dt.month

            avg_by_hour_doy = df.groupby(['DoY', 'Hour'])[['LE', 'NEE']].mean().reset_index()

            for year in df['Year'].unique():
                df_year = df[df['Year'] == year]
                growing_season = df_year[df_year['Month'].isin([5, 6, 7, 8])]
                n_total = len(growing_season)

                for var in ['LE', 'NEE']:
                    n_available = growing_season[var].notna().sum()

                    if n_total > 0 and (n_available / n_total) >= 0.5:
                        for month in [4, 9]:
                            month_data = df_b[(df_b['Year'] == year) & (df_b['Month'] == month)]
                            idx = month_data.index
                            values = df_b.loc[idx, var]

                            is_nan = values.isna()
                            group_id = (is_nan != is_nan.shift()).cumsum()
                            group_sizes = is_nan.groupby(group_id).transform('sum')
                            long_nan_mask = (is_nan) & (group_sizes >= 48 * 7)

                            for i in idx[long_nan_mask].unique():
                                row = df_b.loc[i]
                                doy = row['DoY']
                                hour = row['Hour']
                                match = avg_by_hour_doy[
                                    (avg_by_hour_doy['DoY'] == doy) &
                                    (avg_by_hour_doy['Hour'] == hour)
                                ]
                                if not match.empty:
                                    df_b.at[i, var] = round(match[var].values[0], 3)

            for var in ['LE', 'NEE']:
                df_b[var] = df_b[var].apply(lambda x: round(x, 3) if pd.notna(x) else x)

            output_path = os.path.join(save_path, file)
            df_b.to_csv(output_path, index=False)
            print(f"✔ Processed and saved: {file}")

        except Exception as e:
            print(f"❌ Error processing {file}: {e}")



fill_long_nans_batch(path, save_path)



path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps'
save_path=r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps2'










######################################################################################################################

## dont need to loop through files below. It will go through all files 

################################################################################################################3
## Florida, 5 sites
file="gaps_blend_US-EvM.csv"

file="gaps_blend_US-KS3.csv"


file="gaps_blend_US-KS4.csv"

file="gaps_blend_US-Skr.csv"

file="gaps_blend_US-TaS.csv"

############################################################################################


## california, 7 sites
file="gaps_blend_US-Dmg.csv"


file="gaps_blend_US-EKH.csv"

file="gaps_blend_US-EDN.csv"


file="gaps_blend_US-EKP.csv"


file="gaps_blend_US-Myb.csv"


file="gaps_blend_US-Srr.csv"


file="gaps_blend_US-Tw1.csv"


##########################################################################################
# alaska, 4 sites

file="gaps_blend_US-A03.csv"


file="gaps_blend_US-A10.csv"


file="gaps_blend_US-NGB.csv"


file="gaps_blend_US-Atq.csv"


############################################################################################
# south carolina, 4 sites

file="gaps_blend_US-StS.csv"


file="gaps_blend_US-HB1.csv"

file="gaps_blend_US-HB2.csv"


file="gaps_blend_US-HB3.csv"


##########################################################################################
# delaware: 1 site

file="gaps_blend_US-StJ.csv"


############################################################################################
# louisana: 3 sites

file="gaps_blend_US-LA1.csv"


file="gaps_blend_US-LA2.csv"


file="gaps_blend_US-LA3.csv"


##############################################################################################
# New Jersey 2 sites

file="gaps_blend_US-HPY.csv"


file="gaps_blend_US-MRM.csv"

#############################################################################################
# Massachusetts 1 site 

file="gaps_blend_US-PHM.csv"


#############################################################################################




df,df_b=fill_long_nans(path, file, save_path)












