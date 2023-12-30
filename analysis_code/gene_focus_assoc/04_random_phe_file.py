import pandas as pd
import numpy as np

# Set random seed for reproducibility
np.random.seed(42)

WORKING_DIR = "/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc"

# Read the file
file_path = f"{WORKING_DIR}/input/cnv_bed/all_cohorts.reformatted.phe"
df = pd.read_csv(file_path, sep=' ')

# Check column names
print(df.columns)

# Assuming the third column's name is "AFF"
for i in range(100):
    # Shuffle the data in the third column
    df["AFF"] = df["AFF"].sample(frac=1).reset_index(drop=True) 
    # Convert "AFF" column to integer type
    df["AFF"] = df["AFF"].fillna(0).astype(int) 
    # Save the new file
    output_path = f"{WORKING_DIR}/input/cnv_bed/random_phe_{i}.phe"
    df.to_csv(output_path, sep='\t', header=True, index=False)
