if [ "$GENOME_TYPE" = "hg38" ]
then
  GENOME_DIR=/mnt/sfs-genome/hg38
  GENOME_NAME=GRCh38
  GENOME_BLACKLIST_NAME=hg38
  CORES=30
else
  GENOME_DIR=/mnt/sfs-genome/mm10
  GENOME_NAME=mm10
  GENOME_BLACKLIST_NAME=mm10
  CORES=30
fi

sample="WT_NR4A1 WT_IgG APC_NR4A1 APC_Igg"
for SAMPLE_NAME in $sample
do

if ! [ -f "$SAMPLE_NAME".2.fastq -a -f "$SAMPLE_NAME".1.fastq ]
then
  echo '[Running] decompressing'
  gzip -dkfc "./${SAMPLE_NAME}_1".fq.gz > "$SAMPLE_NAME".1.fq && mv "$SAMPLE_NAME".1.fq "$SAMPLE_NAME".1.fastq &
  gzip -dkfc "./${SAMPLE_NAME}_2".fq.gz > "$SAMPLE_NAME".2.fq && mv "$SAMPLE_NAME".2.fq "$SAMPLE_NAME".2.fastq &
  wait;
else
  echo "[Skip] decompress"
fi

  gzip -dkfc "./${s}_1".fq.gz > "$s".1.fq && mv "$s".1.fq "$s".1.fastq &
  gzip -dkfc "./${s}_2".fq.gz > "$s".2.fq && mv "$s".2.fq "$s".2.fastq &
wait;

if [ ! -f "$SAMPLE_NAME".2.fastq.clipper ]
then
  echo "[Running] trimmomatic adapter sequence"

  java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 30 ${SAMPLE_NAME}.1.fastq ${SAMPLE_NAME}.2.fastq "$SAMPLE_NAME.1.tmp.clipper.fastq" "$SAMPLE_NAME.1.tmp.unpaired.clipper.fastq" "$SAMPLE_NAME.2.tmp.clipper.fastq" "$SAMPLE_NAME.2.tmp.unpaired.clipper.fastq" ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 > logging/"$SAMPLE_NAME".trimmomatic.log 2>&1 \
  && mv "$SAMPLE_NAME.1.tmp.clipper.fastq" "$SAMPLE_NAME.1.fastq.clipper" \
  && mv "$SAMPLE_NAME.2.tmp.clipper.fastq" "$SAMPLE_NAME.2.fastq.clipper"

else
  echo "[Skip] trimmomatic"
fi

if [ ! -f "$SAMPLE_NAME".bowtie2.sam ]
then
  echo "[Running] bowtie2 mapping"
  bowtie2 -x $GENOME_DIR/$GENOME_NAME -q -k 2 -p "$CORES"  -1 "$SAMPLE_NAME".1.fastq.clipper -2 "$SAMPLE_NAME".2.fastq.clipper --no-discordant --no-mixed --no-unal -S "$SAMPLE_NAME".bowtie2.sam --un-conc "$SAMPLE_NAME".bowtie2.un

else
  echo "[Skip] mapping"
fi

if [ ! -f "$SAMPLE_NAME".bam ]
then
  echo "[Running] samtools filter "
  ##Filter reads: read unmapped (0x4);not primary alignment (0x100);read is PCR or optical duplicate (0x400)
  samtools view -Sb -F 1284 "$SAMPLE_NAME".bowtie2.sam |samtools sort - -o "$SAMPLE_NAME".bam \
  && samtools index "$SAMPLE_NAME".bam
else
  echo "[Skip] .bam"
fi

if [ ! -f "$SAMPLE_NAME".tdf ]
then
  echo "[Running] igvtools"
  igvtools count -z 10 -w 5 "$SAMPLE_NAME".bam "$SAMPLE_NAME".tdf $GENOME_DIR/$GENOME_NAME.chrom.sizes
else
  echo "[Skip] .tdf"
fi

if [ ! -f "$SAMPLE_NAME".tdf ]
then
  echo "[Running] picard markdup"
  java -jar picard.jar MarkDuplicates INPUT="$SAMPLE_NAME".bam OUTPUT="$SAMPLE_NAME".nodup.bam METRICS_FILE="$SAMPLE_NAME".dup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=true
else
  echo "[Skip] .tdf"
fi

done