<p align="center">
  <img src="assets/pipeline_banner.png" alt="Integrated Multi-Omics Research Workflow" width="100%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/R-4.2%2B-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R" />
  <img src="https://img.shields.io/badge/Bioconductor-3.16-green?style=for-the-badge" alt="Bioconductor" />
  <img src="https://img.shields.io/badge/Focus-Multi--Omics%20%26%20Oncology-blueviolet?style=for-the-badge" alt="Focus" />
  <img src="https://img.shields.io/badge/Status-Maintained-success?style=for-the-badge" alt="Status" />
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge" alt="License" />
</p>

---

## 🔬 Overview & Scientific Scope

This repository presents an **end-to-end, reproducible computational pipeline** for multi-omics data integration, cancer genomics, co-expression network architecture, machine learning prognostic modeling, and translational biomarker discovery.

Developed and maintained by **Fatemeh Akhoundi** (Computational Biologist & Cancer Bioinformatician), this repository highlights production-grade analytical workflows implemented in **R**.

---

## 🧩 Analytical Modules

| Stage | Module Name | Key Methodologies | Core R Packages & Tools |
| :--- | :--- | :--- | :--- |
| **01** | **Data Acquisition & Harmonization** | TCGA/GEO/GTEx querying, GDC API, RNA-seq count processing, TMM/vst normalization | `TCGAbiolinks`, `GEOquery`, `SummarizedExperiment`, `DESeq2` |
| **02** | **Differential Expression & Functional Profiling** | Linear models for microarrays/RNA-seq, empirical Bayes moderation, GSEA, ORA | `limma`, `edgeR`, `clusterProfiler`, `GSVA`, `org.Hs.eg.db` |
| **03** | **Systems Biology & Co-expression Networks** | Scale-free topology fit, soft-thresholding, dynamic tree cut, module-trait eigengenes | `WGCNA`, `igraph`, `Cytoscape` (RCy3) |
| **04** | **Tumor Microenvironment & Immune Deconvolution** | Single-cell clustering, marker discovery, bulk immune infiltration estimation | `Seurat`, `SingleR`, `immunedeconv`, `ESTIMATE` |
| **05** | **Machine Learning & Prognostic Modeling** | LASSO feature selection, univariate/multivariate Cox regression, time-dependent ROC, Nomograms | `glmnet`, `survival`, `survminer`, `timeROC`, `rms` |

---

## 📁 Repository Structure
```text
├── assets/
│   └── pipeline_banner.png          # High-resolution pipeline workflow banner
├── scripts/
│   └── Cancer_MultiOmics_Pipeline_Demo.R  # Modular workflow demonstration script
├── .gitignore
├── LICENSE
└── README.md                        # Documentation & pipeline overview


---

## 🔬 Methodological Rationale & Key Computational Decisions

### 1. TCGA GDC RNA-Seq Data Harmonization (STAR - Counts)
* **Assay Selection (`unstranded`):** Following the current GDC data workflow (STAR algorithm), the resulting `SummarizedExperiment` objects provide six distinct assays (`unstranded`, `stranded_first`, `stranded_second`, `tpm_unstrand`, `fpkm_unstrand`, and `fpkm_uq_unstrand`). For standard differential count modeling (`limma-voom` / `edgeR`), the **`unstranded`** raw integer count matrix is explicitly extracted to prevent orientation bias and satisfy negative binomial distribution requirements.

### 2. Weighted Gene Co-expression Network Analysis (WGCNA)
* **Network Topology (`signed hybrid`):** Biological gene regulation strictly differentiates between positive co-activation and negative repression. The **`signed hybrid`** model is employed:
  - If correlation > 0: $a_{ij} = (\text{cor}(x_i, x_j))^\beta$
  - If correlation $\le$ 0: $a_{ij} = 0$
  This ensures negatively correlated genes receive zero connectivity weight, preventing biologically opposing regulators from clustering into identical functional modules.
* **Scale-Free Criterion:** The soft-thresholding power ($\beta$) is selected where the scale-free topology fit index ($R^2$) surpasses **0.85**, preserving power-law biological architecture ($P(k) \sim k^{-\gamma}$).

### 3. Prognostic Signature Construction & Dynamic Validation
* **High-Dimensional Regularization (LASSO-Cox):** Implements $L_1$-penalized log-partial likelihood optimization with multi-fold cross-validation (`cv.glmnet`) to effectively eliminate multi-collinearity among co-expressed biomarkers.
* **Dynamic Temporal Assessment (TimeROC):** Prognostic performance is validated through cumulative/dynamic Time-Dependent ROC curves across distinct temporal clinical endpoints (e.g., 1-, 3-, and 5-year OS), integrated with Kaplan-Meier survival stratification and multivariate Cox proportional hazards modeling.


