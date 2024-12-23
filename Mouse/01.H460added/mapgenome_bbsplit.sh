#!/bin/bash

#SBATCH -J gz
#SBATCH -p amd-ep2-short,intel-sc3
#SBATCH -q normal
#SBATCH --mem=64G
#SBATCH -c 32

module load samtools/1.11

GENOME_DIR_hg38=/mnt/sfs-genome/hg38
GENOME_NAME_hg38=GRCh38
GENOME_DIR_mm10=/mnt/sfs-genome/mm10
GENOME_NAME_mm10=mm10
CORES=30

echo "seqID,Input,reads_mapped2hg38,reads_mapped2mm10" > "mapping.log"

sample="H460_1 H460_2 H460_MN1 H460_MN2 H460_g1 H460_g2"
for SAMPLE_NAME in $sample
do
./Software/bbmap/bbsplit.sh ref="${GENOME_DIR_hg38}/${GENOME_NAME_hg38}".fa,"${GENOME_DIR_mm10}/${GENOME_NAME_mm10}".fa in1=$SAMPLE_NAME.1.fastq.clipper in2=$SAMPLE_NAME.2.fastq.clipper basename=$SAMPLE_NAME%_#.fq

alignment_input=$(( $(wc -l "$SAMPLE_NAME".1.fastq.clipper | cut -d' ' -f1) / 4 ))

b0=$(cat $SAMPLE_NAME'GRCh38_1.fq'|wc -l);
b0=$[b0/4]

b1=$(cat $SAMPLE_NAME'mm10_1.fq'|wc -l);
b1=$[b1/4]

echo -en "$SAMPLE_NAME,$alignment_input,$b0,$b1\n" >> "mapping.log" 
done