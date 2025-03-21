#!/bin/bash
WORKING_DIR=`pwd`

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

#FILE_LOCATION_10M_BAMS=`cd $SCRIPT_DIR/../../10m;pwd`
MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

message() {
  local message="$1"
  echo $message
}

# build selected feature bed for rfe
#$SCRIPT_DIR/new_mode.sh "manu_${TYPE}_240103" $FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index
#cd $FILE_LOCATION_10M_BAMS
#ls *.10m.nodup.q30.bam 2>/dev/null | cut -f 1 -d . | xargs -n 1 -P 8 -I %1 $SCRIPT_DIR/tab.sh "manu_${TYPE}_240103" %1

# for TYPE in bc_ca crc_ca gc_ca lc_ca tc_ca
for TYPE in mc_bc_ca mc_crc_ca mc_gc_ca mc_lc_ca mc_tc_ca
do
  cd $WORKING_DIR

  echo Working $TYPE ...
  # make vector from 40 info files
  $SCRIPT_DIR/make_tab.sh $MODEL_DATA_DIR/"${TYPE}_p80_ids.txt" trim_q30_cpm train.tab

  echo "make 1000w tab files"
  # make 10k - 1000k tab
  seq 1 100 | xargs -n 1 -P 50 -I %1 php $SCRIPT_DIR/make_agg_tab_files.php %1
  # repeat 40 times for random ids
  # repeat 40 times for random ids
  if [ -f "all.${TYPE}.sample.info.1" ]; then
      message "Skip gen repeat 80 ids"
  else
      message "Gen repeat 80 ids"
      php $SCRIPT_DIR/repeat_p80.php $MODEL_DATA_DIR/${TYPE}.pos.ids.txt $MODEL_DATA_DIR/${TYPE}.neg.ids.txt 40 "all.${TYPE}.sample.info"
      message "Gen repeat 80 ids finished"
  fi
  # calc 10m feature scores and select top 1000 for 40 random repeat
  seq 1 40 | xargs -n 1 -I %1 $SCRIPT_DIR/fs.sh $TYPE %1

  # merge all 40 top 1000 feature scores
  php $SCRIPT_DIR/merge_p80.php $TYPE
  $SCRIPT_DIR/bed_select "all.${TYPE}.bed" "$FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index" 100
  cut -f 1-3 "$FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index" > temp.$TYPE.bed && mv temp.$TYPE.bed "$FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index"
  # clean up workspace and backup all random ids
#  mv "all.${TYPE}.sample.info".* $FEATURE_SELECTION_OUTPUT_DIR
  rm train.tab.*

  cd $FEATURE_SELECTION_OUTPUT_DIR
  tar zcf $OUTPUT_PREFIX.sample.info.tar.gz all.$TYPE.sample.info.*
  rm $FEATURE_SELECTION_OUTPUT_DIR/all.$TYPE.sample.info.*
done