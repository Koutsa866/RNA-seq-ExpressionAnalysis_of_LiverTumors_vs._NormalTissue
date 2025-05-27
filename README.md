# 🧬 RNA-seq Expression Analysis of Liver Tumors vs. Normal Tissue

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


## 🚀 Workflow Overview

1. Data Acquisition
   - Download FASTQ/SRA files using `sra-tools` or `GEOquery`

2. Quality Control
   - Assess with `FastQC`; aggregate reports using `MultiQC`

3. Transcript Quantification  
   - Use `Salmon` or `kallisto` for reference-based quantification

4. Import & Format Counts  
   - Summarize transcript-level counts to genes using `tximport` in R

5. Differential Expression Analysis
   - Use `DESeq2` to identify significantly regulated genes

6. Visualization & Interpretation
   - Generate PCA, heatmaps, and volcano plots
   - (Optional) Perform pathway or GO enrichment analysis

---

## ✅ Data sets utilized

**What genes are differentially expressed between hepatocellular carcinoma (HCC) and normal liver tissue?**
Using RNA-seq data from the following datasets, this project aims to  identify differentially expressed genes between tumor and non-tumor samples, and visualizes key patterns in expression.

1. https://www.ncbi.nlm.nih.gov/gds/ 

2. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520

3. https://www.ncbi.nlm.nih.gov/geo/query`/acc.cgi?acc=GSE136103

4. https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP174003&o=acc_s%3Aa

5. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520


---

## 📌 Planned Outputs

- Gene count matrices (TPM, raw counts, normalized)
- DE analysis tables (CSV)
- QC reports (FastQC/MultiQC)
- Visual summaries: PCA plots, volcano plots, heatmaps
