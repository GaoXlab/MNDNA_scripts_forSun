#!/bin/bash

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/2_FeatureSelection; pwd)
FEATURE_REDUCTION_OUTPUT_DIR=$(cd $SCRIPT_DIR/../results/3_FeatureReduction; pwd)
FEATURE_CLASSIFICATION_DIR=$(cd $SCRIPT_DIR/../results/4_Classification; pwd)

MODEL_DATA_DIR=$(cd $SCRIPT_DIR/../modelData; pwd)

WORKING_DIR=`pwd`

TMP_DIR=$WORKING_DIR/tmp/$$

mkdir -p $TMP_DIR
cd $TMP_DIR
cp $WORKING_DIR/learn.vector.train.$1 $SCRIPT_DIR/trainIO.py $SCRIPT_DIR/train_and_test_xgbm_rand.py .
mv learn.vector.train.$1 learn.vector.train
cp $WORKING_DIR/learn.window.selection ./learn.window.selection
python train_and_test_xgbm_rand.py
cp output.RFE.txt $WORKING_DIR/output.bed.$1
rm -rf $TMP_DIR