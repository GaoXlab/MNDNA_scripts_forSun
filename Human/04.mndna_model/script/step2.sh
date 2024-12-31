#!/bin/bash

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)
echo "Making train.tab"
OUTPUT_PREFIX=$1
TYPE=$OUTPUT_PREFIX
#$SCRIPT_DIR/make_all_tab.sh trim_q30_cpm train.tab
$SCRIPT_DIR/make_tab.sh $MODEL_DATA_DIR/"${TYPE}_p80_ids.txt" trim_q30_cpm train.tab


echo "make 1000w tab files"
# make 10k - 1000k tab
seq 1 100 | xargs -n 1 -P 50 -I %1 php $SCRIPT_DIR/make_agg_tab_files.php %1
# repeat 40 times for random ids
php $SCRIPT_DIR/repeat_p80.php $MODEL_DATA_DIR/${TYPE}.pos.ids.txt $MODEL_DATA_DIR/${TYPE}.neg.ids.txt 40 "all.${TYPE}.sample.info"
# calc 10m feature scores and select top 1000 for 40 random repeat
seq 1 40 | xargs -n 1 -I %1 $SCRIPT_DIR/fs.sh $TYPE %1

# merge all 40 top 1000 feature scores
php $SCRIPT_DIR/merge_p80_feature_result.php $TYPE> "$FEATURE_SELECTION_OUTPUT_DIR/manu_${OUTPUT_PREFIX}_240103_p80.sorted.tab.index"

# clean up workspace and backup all random ids
mv "all.${TYPE}.sample.info".* $FEATURE_SELECTION_OUTPUT_DIR
rm train.tab.*

cd $FEATURE_SELECTION_OUTPUT_DIR
tar zcf $OUTPUT_PREFIX.sample.info.tar.gz all.$TYPE.sample.info.*
rm $FEATURE_SELECTION_OUTPUT_DIR/all.$TYPE.sample.info.*



