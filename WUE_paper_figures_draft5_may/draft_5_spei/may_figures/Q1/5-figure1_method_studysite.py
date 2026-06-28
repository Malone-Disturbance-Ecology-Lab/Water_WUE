# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 20:56:53 2026

@author: ammar
"""
import panda as pd

#Step 1 — Define your “analysis sites” (ground truth)
wue_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
df = pd.read_csv(wue_path)

analysis_sites = set(df['site_name'].dropna().unique())
print(len(analysis_sites))  # should be 70
#Step 2 — Filter metadata to ONLY those 70 sites
meta_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\site_metadata_with_salinity_SPEIinfo.csv"
meta = pd.read_csv(meta_path)

meta_filtered = meta[meta['site_name'].isin(analysis_sites)]

print(meta_filtered['site_name'].nunique())  # MUST be 70
#Step 3 — Identify excluded sites (important for sanity check)
excluded_sites = set(meta['site_name']) - analysis_sites
print("Excluded sites:", excluded_sites)


meta_filtered['Salinity_Category'].value_counts()

meta_filtered['climate'].value_counts()



#3. Ecosystem (IGBP → your categories)
meta_filtered['IGBP'].value_counts()


def classify_coast(row):
    if row['lat'] > 50:
        return 'Alaska'
    elif row['long'] > -100:
        return 'Atlantic'
    elif row['long'] < -120:
        return 'Pacific'
    else:
        return 'Gulf of Mexico'

meta_filtered['coast'] = meta_filtered.apply(classify_coast, axis=1)

meta_filtered['coast'].value_counts()