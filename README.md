# Paper Related Code and Data

This project provides a bioinformatics analysis pipeline for preprocessing and analyzing the MN-DNA whole-genome sequencing data. The preprocessing pipeline is implemented in shell scripts, and the subsequent analysis for human and mouse samples is organized in the `Human` and `Mouse` folders, respectively. The `Figure` folder contains the code for generating each figure in the article.

# 1. directory structure
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
    └── 04.scRNAseq_process
```

# 2. preprocessing the MN-DNA WGS data
The analysis pipeline for preprocessing the MN-DNA whole-genome sequencing data is 01.pipeline_preprocess.sh.  
The pipeline reads the location of the FastQ files and the output directory from environment variables, as well as the type of genome. When using the pipeline, ensure that the specified locations contain the `sample_name_1.fq.gz` and `sample_name_2.fq.gz` files.
```bash
export SAMPLE_NAME=$SAMPLE_NAME;export SOURCE=$SOURCE_DIR;export OUTPUT_DIR=$OUTPUT_DIR;export GENOME_TYPE=$GENOME_TYPE;./01.pipeline_preprocess.sh
```


# 3. Subsequent Analysis for Human Samples  

- **02.mndna_enriched_peak_calling**: Contains the scripts and processes for generating MN-DNA-enriched regions from 10 deep sequencing samples.   
- **03.features_annotation**: Contains the scripts and processes for the MN-DNA-enriched regions and tumor-associated MN-DNA regions annotation.  
- **04.mndna_model**: Contains the scripts and processes for building the model to detect colorectal cancer or pan-cancer patients using tumor-associated MN-DNA regions. The specific running process is described in the subfolders.  

# 4. Subsequent Analysis for Mouse Samples  

- **01.H460added**: Contains the script mapgenome_bbsplit.sh for mapping reads to multiple reference genomes simultaneously.
- **02.mndna_mouse_svcmodel**: Contains the scripts and processes for building the MN-DNA classification model for APC mice.
- **03.MS**: Contains the mass spectrometry data results.
- **04.scRNAseq_process**: Contains the script and processes for analyzing single-cell RNA sequencing data for APC mice.
