# TCGA LUAD — Multi-Modal Survival Analysis

End-to-end survival and biomarker analysis of lung adenocarcinoma using the TCGA 
Pan-Cancer Atlas 2018 dataset, integrating clinical and somatic mutation data.

## Overview

Lung adenocarcinoma (LUAD) is the most common subtype of non-small cell lung cancer. 
This project analyzes overall survival across clinical and molecular dimensions, 
combining classical biostatistical methods with genomic biomarker analysis — 
reflecting real-world evidence (RWE) methodology used in oncology research.

All data sourced from cBioPortal (publicly available). No PHI involved.

## Dataset

**Source:** [cBioPortal — TCGA LUAD Pan-Cancer Atlas 2018](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018)  
**Patients:** 505 with complete survival data (566 total, 61 excluded — missing OS)  
**Median OS:** 49.3 months  
**OS event rate:** 36.0% (182/505)  
**Data types:** Clinical, somatic mutations, TMB, MSI, copy number (arm-level)

## Analyses Completed

### 1. Overall Survival
Kaplan-Meier estimate for the full cohort. Median OS ~49 months with long-term 
survivors beyond 200 months.

### 2. Survival by Pathologic Stage
Strong stage-dependent separation (log-rank p=1.38e-12). Stage I patients show 
dramatically better long-term survival. Notable crossing of Stage III and IV curves 
early (~10 months), potentially reflecting more aggressive treatment in Stage III — 
to be explored in treatment analysis.

### 3. Tumor Mutational Burden (TMB)
TMB stratified by FDA clinical threshold (≥10 mut/Mb): 179 high, 326 low. 
No significant overall survival difference (p=0.249), consistent with TCGA data 
predating widespread immunotherapy use. Visual separation emerges after ~50 months, 
suggesting TMB may be more relevant as a predictive biomarker for immunotherapy 
response than as a general prognostic factor.

### 4. Microsatellite Instability (MSI)
Only 3/504 patients (0.6%) meet MSI-High threshold (≥3.5) — consistent with LUAD 
being predominantly microsatellite stable. MSI is not a clinically actionable 
stratification factor in this cohort.

### 5. Driver Mutation Prevalence
Somatic mutation data analyzed for 9 clinically relevant genes:

| Gene | Patients mutated | Prevalence |
|------|-----------------|------------|
| TP53 | 262 | 51.9% |
| KRAS | 155 | 30.7% |
| KEAP1 | 96 | 19.0% |
| STK11 | 75 | 14.9% |
| EGFR | 70 | 13.9% |
| BRAF | 41 | 8.1% |
| ALK | 39 | 7.7% |
| ROS1 | 31 | 6.1% |
| MET | 21 | 4.2% |

Prevalences are consistent with published LUAD literature, validating data quality.

### 6. Survival by Driver Mutation (EGFR and KRAS)
EGFR-mutated patients show a trend toward worse overall survival, with no patients 
surviving beyond 120 months vs. long-term survivors in the wild-type group — 
not reaching significance (p=0.197) likely due to limited power (n=70). 
This may reflect TCGA predating widespread third-generation TKI use (osimertinib).

KRAS survival curves cross multiple times (p=0.327), indicating heterogeneity 
within the KRAS-mutated group. KRAS comprises multiple subtypes (G12C, G12V, G12D) 
with distinct biology — to be explored in mutation co-occurrence analysis.

## Next Steps
- KRAS subtype analysis (G12C, G12V, G12D)
- Co-mutation analysis (KRAS + STK11, KRAS + KEAP1)
- Multivariable Cox model integrating stage + mutation status
- ML-based survival prediction integrating clinical + molecular features

## Stack

- Python, pandas — data ingestion and processing
- lifelines — Kaplan-Meier, Cox PH, log-rank tests
- scikit-learn — ML survival prediction
- Matplotlib, seaborn — visualizations
- Jupyter — documented analysis notebook

## Project Structure
```
tcga-luad-rwe/
├── data/
│   ├── raw/          # cBioPortal source files (not tracked in git)
│   └── processed/    # cleaned datasets (not tracked in git)
├── notebooks/
│   ├── 01_eda.ipynb  # main analysis notebook
│   └── figures/      # saved plots
├── src/
│   └── ingest.py     # data ingestion and cleaning
├── requirements.txt
├── .gitignore
└── README.md
```
## Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python src/ingest.py
jupyter notebook
```

---

**Author:** Raquel (Kely) Norel, PhD  
**Domain:** Oncology / Real-World Evidence  
**Status:** 🔄 In progress
