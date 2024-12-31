name=$1
id=$2
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
BASE_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

if [ ! -f "${id}.10m.nodup.q30.bam" ]; then
    echo " BAM file not found."
    exit 1
fi

BAM_FILE="${id}.10m.nodup.q30.bam"
if [ ! -f $BASE_DIR/${name}/origin/"${id}".raw ]
then
  echo "Working ${id}"
  TMP_FILE="$BASE_DIR/${name}/origin/${id}.raw.tmp.$$"
  bedtools coverage -a "$BASE_DIR/${name}/sorted.tab.index" -b "${id}.10m.nodup.q30.bam" -counts -sorted -g $BASE_DIR/genome.txt | cut -f4 > $TMP_FILE
  echo "'${id}.uniq.nodup.bam'" | cat - $TMP_FILE > "$BASE_DIR/${name}/origin/${id}.raw"
  rm $TMP_FILE
else
  echo "Skip ${id}"
  exit 1
fi
# use samtools idxstats to get the total number of reads in the BAM file
# if the index file does not exist, build it first
[ ! -f "$BAM_FILE".bai ] && echo "Build Index" && samtools index $BAM_FILE -@ 6
# get total reads from BAM file excluding chrX, chrY, chrM
TOTAL_COUNT=$(samtools idxstats $BAM_FILE | grep -v '^[XYM]' | awk '{SUM+=$3} END {print SUM}')

# calculate CPM reads for each feature in the tab file, and write to a new file
echo $TOTAL_COUNT
cat $BASE_DIR/"${name}"/origin/"${id}".raw | awk -v total=$TOTAL_COUNT '{if (NR == 1) {print $0} else{cpm = ($1 / total) * 10000000 ;printf "%.1f\n", cpm} }' > $BASE_DIR/"${name}"/cleaned/"${id}".raw