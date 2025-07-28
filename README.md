# FASTQ Quality Analysis

Python script for analyzing FASTQ file quality metrics.

The script processes FASTQ files and reports:
- Read length statistics
- Quality score distributions
- Quality classification (Q20/Q30 thresholds)
- Statistical plots

## Data

Sample FASTQ files `SRR562646_1.fastq` and `SRR562646_2.fastq` (paired-end reads) are from the European Nucleotide Archive (ENA): https://www.ebi.ac.uk/ena/browser/view/SRR562646

## Requirements

```
biopython
matplotlib
seaborn
numpy
```

Install: `pip install biopython matplotlib seaborn numpy`

## Usage

```bash
python analyze_fastq.py
```

The script processes the first 100,000 reads from `SRR562646_1.fastq` and generates:
- Terminal output with statistics
- `fastq_analysis.png` plot file

To analyze different files, modify the `fastq_file` variable in the script.
