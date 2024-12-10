# normalized data
Rscript scripts/STEP1_normalization.R --path modelData --files trn.8m.nodup.q30.bam.10kb_readcounts.txt --output_name Input_raw_nor
Rscript scripts/STEP1_normalization.R --path modelData --files test.8m.nodup.q30.bam.10kb_readcounts.txt --output_name Input_raw_nor.test

# select features
Rscript scripts/STEP2_SelectedSig.R --path modelData --files Input_raw_nor.RData --samples_type1 A_Ctrl --samples_type2 A_Apc --output results

# run model
python3 scripts/STEP3_classification.py