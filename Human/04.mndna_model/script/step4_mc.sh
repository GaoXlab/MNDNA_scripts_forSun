#!/bin/bash
WORKING_DIR=`pwd`

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

# make vector from ids and final whitelist
php $SCRIPT_DIR/make_vector_mc.php "$MODEL_DATA_DIR/mc.train.ids.txt" "$MODEL_DATA_DIR/mc.test.ids.txt" manu_mc_240103 $FEATURE_REDUCTION_OUTPUT_DIR/manu_mc_240103.bed
# model prediction
python $SCRIPT_DIR/train_and_test_xgbm_mc.py

mv xgbm.* $FEATURE_CLASSIFICATION_DIR