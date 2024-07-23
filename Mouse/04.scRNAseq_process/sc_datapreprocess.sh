#!/bin/bash
#SBATCH -J BM4
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH -c 40
#SBATCH --mem=400G

module load cellranger/6.1.1
module load R/4.2.1
module load gcc/11.2.0
module load fftw/3.3.10

### 1.preprocess
### for only single-cell RNA-seq:
rawdata_path='Data/CleanData/'
sample="APC1 APC2 APC3 APC4 WT1 WT2 WT3 WT4"
for s in $sample
do
cellranger count --id=$s --fastqs=$rawdata_path/$s --sample=$s --transcriptome=refdata-gex-mm10-2020-A/ --localcores=40
done
### for single-cell RNA-seq and TCR-seq:
# cellranger multi --id=D01_multi_result --csv=D01_multi-config-template.csv --localcores=40

### 2.quality control (single sample)
R_script="Rscript/"
sample="APC1 APC2 APC3 APC4 WT1 WT2 WT3 WT4"
for s in $sample
do
subpath="${s}/outs/filtered_feature_bc_matrix/"
mkdir $s
Rscript ${R_script}/1.Seurat3_QualityControl.r -i 1.rawdata/${subpath} -w 2.cellanno/$s/ -o ${s} -s Mus > 2.cellanno/${s}/1.QualityControl.log
done


### 3.combined
combined="Final_20w_WTAPC"
mkdir $combined
cd $combined
ln -s ../../2.cellanno/APC1/*_raw.rds APC1_raw.rds
ln -s ../../2.cellanno/APC2/*_raw.rds APC2_raw.rds
ln -s ../../2.cellanno/APC3/*_raw.rds APC3_raw.rds
ln -s ../../2.cellanno/APC4/*_raw.rds APC4_raw.rds
ln -s ../../2.cellanno/WT1/*_raw.rds WT1_raw.rds
ln -s ../../2.cellanno/WT2/*_raw.rds WT2_raw.rds
ln -s ../../2.cellanno/WT3/*_raw.rds WT3_raw.rds
ln -s ../../2.cellanno/WT4/*_raw.rds WT4_raw.rds
cd ..


sample_names="APC1,APC2,APC3,APC4,WT1,WT2,WT3,WT4"
rds_files="APC1_raw.rds,APC2_raw.rds,APC3_raw.rds,APC4_raw.rds,WT1_raw.rds,WT2_raw.rds,WT3_raw.rds,WT4_raw.rds"
filter_fea="7500,7500,7500,7500,7500,7500,7500,7500"
filter_mt="10,10,10,10,10,10,10,10"
Rscript ${R_script}/2.IntegratedSamples.Seurat3_PCAselection.r -w 3.combination/$combined -l $sample_names -f $rds_files -u $filter_fea -m $filter_mt -o ${combined} > 3.combination/2.PCAselection.log
Rscript ${R_script}/CCA_3.Seurat3_tSNEorUMAP.r -w 3.combination/$combined -p 30 -f ${combined}_Origin_Integrated.rds -o ${combined} > 3.combination/3.CCA_p30.log
