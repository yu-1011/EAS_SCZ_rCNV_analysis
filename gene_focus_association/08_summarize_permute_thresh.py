import pandas as pd
import os
import sys

# Enter your directory path
directory = '/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/output/permute'
CNV_type = sys.argv[1]

# Store all extracted LOG10P values
all_log10p_values = []

# Iterate through all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(f"_{CNV_type}.step2_AFF.regenie"):
        # Construct the full path of the file
        file_path = os.path.join(directory, filename)
        
        # Read the file
        try:
            data = pd.read_csv(file_path, delim_whitespace=True)
        except Exception as e:
            print(f"Unable to read file {file_path}: {e}")
            continue
        
        # Sort the LOG10P column and get the top 5% of the data
        top_5_percent = data['LOG10P'].nlargest(int(len(data) * 0.05))
        
        # Add these values to the list
        all_log10p_values.extend(top_5_percent)

# Convert LOG10P values to p-values
p_values = [10 ** -value for value in all_log10p_values]

# Sort the p-values
p_values.sort()

# Calculate the position of the 5% ranking
rank_idx = int(len(p_values) * 0.05)

# Get the p-value at the 5% percentile
p_value_5_percentile = p_values[rank_idx]

print("The p-value at 5% ranking is:", p_value_5_percentile)
