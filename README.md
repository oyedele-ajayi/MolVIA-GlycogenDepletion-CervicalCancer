# Transcriptomic Analysis of Glycogen Depletion in Cervical Cancer

This repository contains the complete R analysis pipeline for the manuscript:

**Title:** "Transcriptomic Evidence for Glycogen Depletion as the Molecular Basis of Acetowhitening in Cervical Dysplasia: The Glycogen Depletion Index"  
**Authors:** Oyedele Ajayi, Efosa Odigie  
**Journal:** [Journal Name - to be inserted upon acceptance]  
**DOI:** [Manuscript DOI - to be inserted upon acceptance]  
**Status:** Under Review

---

## Overview

This study characterizes glycogen metabolic reprogramming across the spectrum of cervical intraepithelial neoplasia (CIN) and invasive carcinoma. We introduce the Glycogen Depletion Index (GDI) as a novel composite molecular biomarker derived from transcriptomic data.

The analysis demonstrates that:
- `PPP1R3C` silencing (log₂FC = −3.43) is a primary molecular driver of glycogen synthesis impairment.
- The GDI achieves an AUC of 0.814 in the discovery dataset and is validated with an AUC of 0.895 in an independent, cross-platform dataset.
- Glycogen pathway reprogramming is grade-dependent, emerging at CIN3 and maximized in invasive carcinoma.

This repository provides all the necessary code to reproduce these findings from publicly available data.

## Datasets

All primary data are publicly available from the NCBI Gene Expression Omnibus (GEO). **No raw data files are stored in this repository.** The analysis scripts download all required data programmatically.

| Dataset   | Role       | Platform           | n   | GEO Accession                                                                |
|:----------|:-----------|:-------------------|:----|:-----------------------------------------------------------------------------|
| GSE63514  | Discovery  | GPL570 (HGU133Plus2) | 128 | [GSE63514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63514) |
| GSE9750   | Validation | GPL96 (HGU133A)    | 57  | [GSE9750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9750)   |

### Sample Composition

#### GSE63514 (Discovery)
| Group                      | n   |
|:---------------------------|:----|
| Normal cervical epithelium | 24  |
| CIN1                       | 14  |
| CIN2                       | 22  |
| CIN3                       | 40  |
| Invasive carcinoma         | 28  |
| **Total**                  | **128** |

#### GSE9750 (Validation)
| Group                      | n   |
|:---------------------------|:----|
| Normal cervical epithelium | 24  |
| Invasive carcinoma         | 33  |
| **Total**                  | **57**  |
*Note: 9 cervical cancer cell line samples were excluded from GSE9750 prior to analysis as detailed in the scripts.*

---

## How to Run the Analysis

### Requirements
- R version 4.2 or higher
- An active internet connection (to download data from GEO)
- See `session_info/R_session_info.txt` for a complete list of R packages and their versions.

### Instructions
1. **Clone the repository:**
   ```bash
   git clone https://github.com/oyedele-ajayi/MolVIA-GlycogenDepletion-CervicalCancer.git
