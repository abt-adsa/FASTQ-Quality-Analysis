# FASTQ Quality Analysis Pipeline
A Python-based bioinformatics tool for analyzing paired-end FASTQ sequencing data with comprehensive quality metrics, statistical analysis, and publication-ready visualizations.

## Key Features
- **Paired-End Data Support**: Handles R1/R2 FASTQ files from paired-end sequencing
- **Quality Score Analysis**: Phred quality score assessment with Q20/Q30 thresholds
- **Statistical Reporting**: Comprehensive read length and quality metrics
- **Robust File Parsing**: Error handling for large genomic datasets
- **Progress Monitoring**: Real-time progress indicators for batch processing
- **Data Visualization**: Publication-ready histograms with statistical overlays
- **Performance Optimization**: Efficient processing of large FASTQ files
- **Command-Line Interface**: Flexible file input with argument support

## Technology Stack
- **Platform**: Python 3.10+
- **Bioinformatics**: BioPython for sequence parsing
- **Data Analysis**: NumPy for statistical computations
- **Visualization**: Matplotlib, Seaborn for scientific plotting
- **File Processing**: Robust FASTQ format handling

## Project Structure
```
├── analyze_fastq.py        # Main analysis script
├── requirements.txt        # Python dependencies
├── README.md              # Project documentation
└── fastq_analysis_results.png  # Generated visualization output
```

## Sample Data
The analysis was performed on SRR562646 paired-end sequencing data from the European Nucleotide Archive (ENA):
- **Data Source**: https://www.ebi.ac.uk/ena/browser/view/SRR562646
- **File Format**: Paired-end FASTQ files (R1/R2)
- **Note**: FASTQ files are not included in this repository due to large file size (~1.7GB each)

To reproduce the analysis:
1. Download the FASTQ files from the ENA link above
2. Place them in the project directory as `SRR562646_1.fastq` and `SRR562646_2.fastq`
3. Run the analysis script

## Technical Skills Demonstrated
- Bioinformatics data processing and FASTQ format handling
- Python programming with scientific computing libraries
- Statistical analysis and quality control in genomics
- Data visualization and scientific plotting
- Error handling and performance optimization for large datasets
- Command-line tool development and user interface design

---
*Developed for demonstrating bioinformatics, genomic data analysis, and scientific computing capabilities.*
