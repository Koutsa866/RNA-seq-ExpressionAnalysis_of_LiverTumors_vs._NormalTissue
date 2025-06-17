# 🧬 RNA-seq Expression Analysis of Liver Tumors vs. Normal Tissue

This project investigates gene expression differences between hepatocellular carcinoma (HCC) and healthy liver tissue using public RNA-seq data. The goal is to identify differentially expressed genes (DEGs) and interpret the biological mechanisms driving liver tumorigenesis.

---

## 🎯 Research Purpose

Hepatocellular carcinoma is one of the leading causes of cancer-related death worldwide. Understanding transcriptional changes between healthy and tumor tissue can:

- Reveal potential biomarkers
- Highlight therapeutic targets
- Elucidate dysregulated pathways in cancer

---

## ✅ Dataset Used

- [GSE62232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62232)  
  *(Primary dataset used for differential expression and statistical analysis)*

---

## 🚀 Workflow Overview

### 🔹 1. Data Acquisition
- Imported the **GSE62232** dataset using the Python library [`GEOparse`](https://github.com/guma44/GEOparse)

### 🔹 2. Expression Matrix Preparation
- Transposed GSE62232 expression data to gene × sample format
- Extracted and saved **metadata** and **expression data** into separate CSV files
- Merged into a single `merged_expression_metadata.csv` for downstream analysis

### 🔹 3. Quality Control *(Skipped)*
- Raw-level QC was not needed since GSE62232 contains preprocessed data

### 🔹 4. Exploratory Data Analysis
- Bar plot of sample counts (Normal vs. Tumor)
- Principal Component Analysis (PCA) on expression profiles
- Summary statistics: gene-wise means and variance
- KMeans clustering of samples based on expression profiles

### 🔹 5. Differential Expression Analysis
- Conducted using **Welch’s t-test** (`scipy.stats.ttest_ind`) between tumor and normal groups
- Computed log2 fold change for each gene
- Implemented chunked processing (1,000 genes per batch) for performance
- Filtered DEGs by:
  - p-value < 0.05
  - |log2FC| > 1

---

## 📊 Methods Summary

**Statistical Approach:**
- Welch’s t-test assuming unequal variance
- Log2 Fold Change (log2FC)
- Significance cutoff: p < 0.05, |log2FC| > 1
- Chunked DEG testing for scalability on large datasets

---

## 📌 Outputs

- ✅ `merged_expression_metadata.csv`: Combined expression and sample metadata
- ✅ `chunks/`: DEG test results for each 1000-gene batch
- ✅ `DEG_full_results.csv`: All genes tested with statistics
- ✅ `DEG_significant_filtered.csv`: Filtered, statistically significant DEGs
- ✅ Figures:
  - PCA plot
  - Clustering heatmap
  - Sample distribution bar chart

---

## 🧭 Next Steps

- Apply False Discovery Rate (FDR) correction (Benjamini-Hochberg)
- Perform GO / pathway enrichment analysis on significant DEGs
- Visualize results via volcano plots and heatmaps
- Compare DEGs with public cancer databases (e.g., TCGA)
- Optionally integrate with methylation data in future research

---

## 📂 Project Structure

Needs Updating!!!


---

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
