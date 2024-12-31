#!/bin/bash
WORKING_DIR=`pwd`

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

TYPE=$1

# make vector from ids and final whitelist
php $SCRIPT_DIR/make_vector.php "$MODEL_DATA_DIR/${TYPE}.train.pos.ids.txt" "$MODEL_DATA_DIR/${TYPE}.train.neg.ids.txt" "$MODEL_DATA_DIR/${TYPE}.test.ids.txt" manu_${TYPE}_240103 $FEATURE_REDUCTION_OUTPUT_DIR/manu_${TYPE}_240103.bed
# model prediction
python $SCRIPT_DIR/train_and_test_xgbm.py

mv xgbm.* $FEATURE_CLASSIFICATION_DIR