#!/bin/bash
WORKING_DIR=`pwd`

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

#FILE_LOCATION_10M_BAMS=`cd $SCRIPT_DIR/../../10m;pwd`
MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

TYPE=$1
# build selected feature bed for rfe
$SCRIPT_DIR/new_mode.sh "manu_${TYPE}_240103" $FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index
cd $FILE_LOCATION_10M_BAMS
ls *.10m.nodup.q30.bam 2>/dev/null | cut -f 1 -d . | xargs -n 1 -P 8 -I %1 $SCRIPT_DIR/tab.sh "manu_${TYPE}_240103" %1

cd $WORKING_DIR
# output rfe tab for backup
$SCRIPT_DIR/make_all_tab.sh manu_crc_240103 all.manu_crc_240103.tab
# repeat 40 times for random ids
php $SCRIPT_DIR/repeat_p80.php $MODEL_DATA_DIR/${TYPE}.pos.ids.txt $MODEL_DATA_DIR/${TYPE}.neg.ids.txt 40 "all.${TYPE}.rfe.sample.info"
# make vector from 40 info files
seq 1 40 | xargs -n 1 -P 40 -I %1 php $SCRIPT_DIR/make_vector_from_info.php "all.${TYPE}.rfe.sample.info.%1" "manu_${TYPE}_240103" $FEATURE_SELECTION_OUTPUT_DIR/manu_${TYPE}_240103_p80.sorted.tab.index learn.vector.train.%1

cd $WORKING_DIR
# using rfe for feature selection
seq 1 40 | xargs -n 1 $SCRIPT_DIR/rfe.sh
# merge output.[1-40] to manu_crc_240103.bed
php $SCRIPT_DIR/merge_p80_rfe_result.php > $FEATURE_REDUCTION_OUTPUT_DIR/manu_"${TYPE}"_240103.bed
mv all.${TYPE}.rfe.sample.info.* $FEATURE_REDUCTION_OUTPUT_DIR
# clean up and backup all random ids
cd $FEATURE_REDUCTION_OUTPUT_DIR
tar zcf $TYPE.rfe.sample.info.tar.gz all.${TYPE}.rfe.sample.info.*