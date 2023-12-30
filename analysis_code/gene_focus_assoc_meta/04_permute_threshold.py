import glob
import pandas as pd
import os
import sys
import numpy as np


CNV_type = sys.argv[1]
WORKING_DIR = "/stanley/huang_lab/home/ychen/proj-CNV/08CNV_meta_analysis/02gene_assoc_meta/"

# Set the file pattern to match all relevant files
files = glob.glob(f"{WORKING_DIR}/output/permute/META_PGC_EAS.*_*_1.txt")

# Initialize an empty list to store P-values from all files
all_p_values = []

# Iterate over the list of files, reading data from each
for file in files:
    data = pd.read_csv(file, sep='\t')  # Assuming fields are tab-separated
    all_p_values.extend(data['P-value'].tolist())  # Add all P-values from the file

# Sort all the P-values
sorted_p_values = sorted(all_p_values)

# Select the most significant 5% of P-values
top_5_percent = sorted_p_values[:int(0.05 * len(sorted_p_values))]

# Find the 5th percentile P-value within this top 5%
percentile_5th_of_top_5_percent = np.percentile(top_5_percent, 5)

print("The p-value at the 5% ranking is:", percentile_5th_of_top_5_percent)
