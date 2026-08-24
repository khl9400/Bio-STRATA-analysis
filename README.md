# Bio-STRATA Analysis

Analysis scripts for amplicon sequencing data generated in the Bio-STRATA study.

## Overview

This repository contains custom Python scripts used to process amplicon sequencing reads and quantify wild-type (WT) and mutant (MT) synthetic mitochondrial DNA variants.

The analysis workflow consists of:

1. Adapter trimming
2. Insert-length filtering
3. Read orientation determination and normalization
4. SNP classification
5. Calculation of WT and MT read fractions

## Analysis workflow

### 1. Adapter trimming and orientation normalization

`01_adapter_trim_and_normalize.py`

This script:

* detects and removes sequencing adapter sequences using sequence similarity;
* retains the expected 120-bp amplicon insert;
* determines whether each insert is in the forward or reverse orientation;
* reverse-complements reverse-oriented reads so that all retained sequences are represented in the forward orientation.

The script automatically processes available sample files following the expected sample naming convention.

### 2. SNP and MT fraction analysis

`02_snp_mt_ratio_analysis.py`

Normalized sequences are classified according to the nucleotide present at SNP position 52 (1-based):

* **T:** WT
* **C:** MT
* **Other bases or insufficient sequence length:** Unknown

The MT fraction is calculated as:

`MT fraction (%) = MT reads / (WT reads + MT reads) × 100`

Unknown reads are excluded from the denominator.

## Requirements

* Python 3
* pandas
* openpyxl

Required packages can be installed using:

```bash
pip install pandas openpyxl
```

## Usage

Place the analysis scripts and input files in the same directory.

Run the preprocessing script first:

```bash
python 01_adapter_trim_and_normalize.py
```

Then run the SNP analysis:

```bash
python 02_snp_mt_ratio_analysis.py
```

The scripts automatically process the available Bio-STRATA samples and generate sample-specific analysis files.

## Output

The preprocessing step generates files containing:

* adapter-trimming results;
* read-orientation classification;
* normalized forward-oriented sequences;
* excluded or unclassified reads.

The SNP analysis generates:

* total read counts;
* WT read counts;
* MT read counts;
* Unknown read counts;
* WT fraction;
* MT fraction;
* read-level SNP classification results.

## Data availability

Raw sequencing data associated with this study will be deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession **[PRJNAxxxxxx]**.

## License

This repository is distributed under the MIT License.
