## Statistical Analysis / Differential Expression Analysis
# This script has been moved into a new python script called "DEG_StatisticalAnalysis.py" for better organization and modularity.
# This script performs differential expression analysis on gene expression data, comparing tumor and normal samples.
# It calculates log2 fold changes and p-values for each gene, and saves the results in chunks to avoid memory issues.
# It also filters the results for significant differentially expressed genes (DEGs) based on p-value and fold change thresholds.
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import os

# STEP 1: Load the merged data
file_path = 'OutputFiles/merged_expression_metadata.csv'
if not os.path.exists(file_path):
    raise FileNotFoundError(f"CSV not found at: {file_path}")

df = pd.read_csv(file_path, index_col=0)

# STEP 2: Separate metadata and expression data
meta = df[['Title', 'Condition']]
expr = df.drop(columns=['Title', 'Condition'])

# STEP 3: Split expression data by condition: Tumor vs Normal
tumor_samples = meta[meta['Condition'] == 'Tumor'].index
normal_samples = meta[meta['Condition'] == 'Normal'].index

expr_tumor = expr.loc[tumor_samples]
expr_normal = expr.loc[normal_samples]

# === NEW: Chunked Processing === #
chunk_size = 1000
gene_chunks = np.array_split(expr.columns, int(np.ceil(len(expr.columns) / chunk_size)))

os.makedirs('OutputFiles/chunks', exist_ok=True)

for i, chunk in enumerate(gene_chunks):
    chunk_results = []

    for gene in chunk:
        tumor_vals = expr_tumor[gene].dropna()
        normal_vals = expr_normal[gene].dropna()

        # Avoid division by zero
        mean_tumor = np.mean(tumor_vals) if len(tumor_vals) > 0 else np.nan
        mean_normal = np.mean(normal_vals) if len(normal_vals) > 0 else np.nan

        log2_fc = np.log2(mean_tumor / mean_normal) if mean_normal and mean_normal > 0 else np.nan

        # Perform t-test only if both groups have >1 value
        if len(tumor_vals) > 1 and len(normal_vals) > 1:
            _, p_value = ttest_ind(tumor_vals, normal_vals, equal_var=False)
        else:
            p_value = np.nan

        chunk_results.append({
            'Gene': gene,
            'mean_tumor': mean_tumor,
            'mean_normal': mean_normal,
            'log2FC': log2_fc,
            'p_value': p_value
        })

    chunk_df = pd.DataFrame(chunk_results)
    chunk_df.to_csv(f'OutputFiles/chunks/DEG_chunk_{i+1}.csv', index=False)
    print(f"✅ Processed chunk {i+1}/{len(gene_chunks)}")

# === Combine All Chunks into One DataFrame === #
all_chunks = []
for i in range(len(gene_chunks)):
    chunk_file = f'OutputFiles/chunks/DEG_chunk_{i+1}.csv'
    if os.path.exists(chunk_file):
        all_chunks.append(pd.read_csv(chunk_file))
    else:
        print(f"⚠️ Missing chunk file: {chunk_file}")

if all_chunks:
    combined_df = pd.concat(all_chunks, ignore_index=True)
    combined_df = combined_df.sort_values('p_value')

    # STEP 6: Filter for significant DEGs
    deg_filtered = combined_df[(combined_df['p_value'] < 0.05) & (combined_df['log2FC'].abs() > 1)]

    # STEP 7: Save final results
    combined_df.to_csv('OutputFiles/DEG_full_results.csv', index=False)
    deg_filtered.to_csv('OutputFiles/DEG_significant_filtered.csv', index=False)

    # Final output
    print(" DEG analysis complete.")
    print(" Full results saved to: OutputFiles/DEG_full_results.csv")
    print(" Significant DEGs saved to: OutputFiles/DEG_significant_filtered.csv")
else:
    print(" No chunk results found. Please check earlier steps.")







#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 🧩 Script 1: prep_data.py
# Load and save tumor/normal matrices separately
df = pd.read_csv('OutputFiles/merged_expression_metadata.csv', index_col=0)
meta = df[['Title', 'Condition']]
expr = df.drop(columns=['Title', 'Condition'])
tumor = expr.loc[meta[meta['Condition'] == 'Tumor'].index]
normal = expr.loc[meta[meta['Condition'] == 'Normal'].index]
tumor.to_csv('OutputFiles/tumor_expr.csv')
normal.to_csv('OutputFiles/normal_expr.csv')

# 🧪 Script 2: deg_analysis_chunks.py
import pandas as pd, numpy as np
from scipy.stats import ttest_ind
import os

tumor = pd.read_csv('OutputFiles/tumor_expr.csv', index_col=0)
normal = pd.read_csv('OutputFiles/normal_expr.csv', index_col=0)

genes = tumor.columns
chunk_size = 500
results = []

for i in range(0, len(genes), chunk_size):
    chunk_genes = genes[i:i+chunk_size]
    chunk_result = []
    for gene in chunk_genes:
        t_vals, n_vals = tumor[gene], normal[gene]
        mean_t, mean_n = t_vals.mean(), n_vals.mean()
        log2fc = np.log2(mean_t / mean_n) if mean_n != 0 else np.nan
        _, p = ttest_ind(t_vals, n_vals, equal_var=False)
        chunk_result.append({
            'Gene': gene,
            'mean_tumor': mean_t,
            'mean_normal': mean_n,
            'log2FC': log2fc,
            'p_value': p
        })
    pd.DataFrame(chunk_result).to_csv(f'OutputFiles/chunks/DEG_chunk_{i//chunk_size + 1}.csv', index=False)
