# Bio-STRATA Analysis

Custom Python scripts used for amplicon sequencing analysis in the Bio-STRATA study.

This repository contains separate preprocessing and SNP/MT-fraction analysis workflows for the R1 and R2 sequencing datasets from samples A, B, C, D, and F.

## Analysis Overview

The workflow consists of four scripts:

1. `01_adapter_trim_and_normalize.py`  
   Preprocesses R1 sequencing reads.

2. `02_snp_mt_ratio_analysis.py`  
   Performs SNP classification and WT/MT fraction calculation for R1.

3. `03_R2_adapter_trim_and_normalize.py`  
   Preprocesses R2 sequencing reads.

4. `04_R2_snp_mt_ratio_analysis.py`  
   Performs SNP classification and WT/MT fraction calculation for R2.

R1 and R2 read datasets are analyzed separately. After preprocessing, both workflows use the same SNP-calling definition on sequences normalized to the forward orientation.

## Requirements

- Python 3
- pandas
- openpyxl

Install the required packages with:

```bash
pip install pandas openpyxl
