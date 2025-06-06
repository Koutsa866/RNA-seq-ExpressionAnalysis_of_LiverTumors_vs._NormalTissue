
# Importing necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import os
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### This code reads in a CSV file containing gene expression data, transposes it, and saves the transposed data to a new CSV file.
# Reading in CSV file through pandas. 
import pandas as pd
df = pd.read_csv('GSE62232_ExpData.csv', index_col=0) # Reading in CSV file
# Converting CSV file to transposed format and saving it as a new CSV file
df_transposed = df.transpose()
print(df_transposed.head())
# Ensure the output directory exists
output_dir = 'OutputFiles'
os.makedirs(output_dir, exist_ok=True)
# Save the transposed data
transposed_file_path = os.path.join(output_dir, 'transposed_expression.csv')
df_transposed.to_csv(transposed_file_path)
print(f'Transposed data saved to: {transposed_file_path}')
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### This code downloads and loads the GEO Series GSE62232, which contains gene expression data for breast cancer samples.
#Import necessary libraries
import pandas as pd
import GEOparse 
# Download and load the GEO Series
gse = GEOparse.get_GEO("GSE62232", destdir=".")
# Extract sample titles to determine condition
for gsm_name, gsm in gse.gsms.items():
    print(gsm_name, gsm.metadata['title'][0])
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Building metadata for GSE62232 samples
import GEOparse
import pandas as pd
import os
# Load the dataset
gse = GEOparse.get_GEO("GSE62232")

# Build metadata
sample_conditions = []
for gsm_name, gsm in gse.gsms.items():
    title = gsm.metadata['title'][0]
    if "N" in title and not "T" in title:
        condition = "Normal"
    else:
        condition = "Tumor"
    sample_conditions.append({"Sample": gsm_name, "Title": title, "Condition": condition})
# Convert to DataFrame and save
meta_df = pd.DataFrame(sample_conditions)
output_dir = 'OutputFiles'
os.makedirs(output_dir, exist_ok=True)
meta_df.to_csv("GSE62232_sample_metadata.csv", index=False)
print(meta_df.head())
# Align metadata and expression data!
import pandas as pd
import os
# Ensure the output directory exists
output_dir = 'OutputFiles'
os.makedirs(output_dir, exist_ok=True)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Merging expression data with metadata
import pandas as pd
import os
# Define file paths
output_dir = 'OutputFiles'
expr_file = os.path.join(output_dir, 'transposed_expression.csv')
meta_file = os.path.join(output_dir, 'GSE62232_sample_metadata.csv')

# Load expression data (samples as rows, genes as columns)
expr_data = pd.read_csv(expr_file, index_col=0)

# Load metadata (with columns: Sample, Title, Condition)
meta_data = pd.read_csv(meta_file)

# Rename 'Sample' to 'sample_name' to avoid confusion
meta_data.rename(columns={'Sample': 'sample_name'}, inplace=True)

# Set 'sample_name' as the index to match expression data index
meta_data.set_index('sample_name', inplace=True)

# Merge metadata (left) with expression data (right) on sample ID (index)
merged_data = meta_data.join(expr_data, how ='inner')
# Display merged data
print(merged_data.head())
# Save merged data
merged_file_path = os.path.join(output_dir, 'merged_expression_metadata.csv')
merged_data.to_csv(merged_file_path)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Playing around with the data. 
# 1-> (Bar plot of Normal vs Tumor samples)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the merged data directly from the CSV file
merged_data = pd.read_csv('OutputFiles/merged_expression_metadata.csv', index_col=0)

# Count the number of samples for each condition
condition_counts = merged_data['Condition'].value_counts().reset_index()
condition_counts.columns = ['Condition', 'Sample_Count']

# Display the counts
print(condition_counts)

# Plot the bar plot with solid bars
plt.figure(figsize=(8, 6))
sns.barplot(data=condition_counts, x='Condition', y='Sample_Count', palette='viridis')
plt.title("Number of Tumor vs Normal Samples")
plt.xlabel("Condition")
plt.ylabel("Number of Samples")
plt.show()
# Save the plot 
output_dir = 'OutputFiles'
os.makedirs(output_dir, exist_ok=True)
plt.savefig(os.path.join(output_dir, 'normal_vs_tumor_samples_bar_plot.png'))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 2-> PCA Plot (Scatterm Plot-Principal Component Analysis (PCA). Reduces the complexity of high-dimensional data (like gene expression profiles) into just 2 or 3 new axes that capture the most variation in your data.
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScalerimport 
import seaborn as sns
# Load the merged data directly from the CSV file
merged_data = pd.read_csv('OutputFiles/merged_expression_metadata.csv', index_col=0)
# Separate expression data from metadata
expr_data = df.drop(columns=["Sample", "Title", "Condition"])
meta_data = df[["Sample", "Condition"]]