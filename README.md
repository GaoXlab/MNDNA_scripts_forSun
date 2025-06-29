# Paper Related Code and Data

This repository contains code from Sun et al., Cell Res (DOI: 10.1038/s41422-025-01122-7).

This project provides a bioinformatics analysis pipeline for preprocessing and analyzing the DNA remnants in red blood cells (rbcDNA, also referred to as MN-DNA in some original figures or code) whole-genome sequencing data. The preprocessing pipeline is implemented in shell scripts, and the subsequent analysis for human and mouse samples is organized in the `Human` and `Mouse` folders, respectively. The `Figure` folder contains the code for generating each figure in the article.

## 1. Directory structure
```text
├── Figure
├── Human
│   ├── 02.mndna_enriched_peak_calling
│   ├── 03.features_annotation
│   └── 04.mndna_model
│       ├── modelData
│       ├── results
│       └── script
└── Mouse
    ├── 01.H460added
    ├── 02.mndna_mouse_svcmodel
    ├── 03.MS
    ├── 04.scRNAseq_process
    └── 05.ChIPseq_process
```

## 2. Preprocessing the rbcDNA WGS data
The analysis pipeline for preprocessing the rbcDNA whole-genome sequencing data is `01.pipeline_preprocess.sh`. (The reference genomes for human (hg38) and mouse (mm10) used were restricted to primary chromosomes, excluding alternate contigs and unlocalized or unplaced sequences (e.g., _alt, _random, chrUn).) 
The pipeline reads the location of the FastQ files and the output directory from environment variables, as well as the type of genome. When using the pipeline, ensure that the specified locations contain the `sample_name_1.fq.gz` and `sample_name_2.fq.gz` files.
```bash

# ./00.build_primary_fa.sh GRCh38  # Generates GRCh38.fa
# ./00.build_primary_fa.sh mm10    # Generates mm10.fa

export SAMPLE_NAME=$SAMPLE_NAME;export SOURCE=$SOURCE_DIR;export OUTPUT_DIR=$OUTPUT_DIR;export GENOME_TYPE=$GENOME_TYPE;./01.pipeline_preprocess.sh

# For human genome-wide rbcDNA distribution, related to Fig1e,f
Rscript ./02.pipeline_QDNAseq.r ./ $SAMPLE_NAME.60m.nodup.q30.bam 1000 hg38
Rscript ./02.pipeline_QDNAseq.r ./ $SAMPLE_NAME.60m.nodup.q30.bam 100 hg38
# For mouse sample 
Rscript ./02.pipeline_QDNAseq.r ./ $SAMPLE_NAME.8m.mm10.nodup.q30.bam 10 mm10
```


## 3. Subsequent Analysis for Human Samples  

- **02.mndna_enriched_peak_calling**: Contains the scripts and processes for generating rbcDNA-enriched regions from 10 deep sequencing samples.   
- **03.features_annotation**: Contains the scripts and processes for the rbcDNA-enriched regions and tumor-associated rbcDNA regions annotation.  
- **04.mndna_model**: Contains the scripts and processes for building the model to detect colorectal cancer or pan-cancer patients using tumor-associated rbcDNA regions. The specific running process is described in the subfolders.  

## 4. Subsequent Analysis for Mouse Samples  

- **01.H460added**: Contains the script `mapgenome_bbsplit.sh` for mapping reads to multiple reference genomes simultaneously.
- **02.mndna_mouse_svcmodel**: Contains the scripts and processes for building the rbcDNA classification model for APC mice.
- **03.MS**: Contains the mass spectrometry data results.
- **04.scRNAseq_process**: Contains the script and processes for analyzing single-cell RNA sequencing data for APC mice.


## Software version and hardware requirements

- bwa-mem2@2.2.1
- samtools@1.10
- bedtools@2.27.1
- php@7.4
- deeptools@3.3.2
- Python@3.8
- R@4.2.2+