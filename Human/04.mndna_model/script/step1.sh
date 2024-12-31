#!/bin/bash
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
FEATURE_SELECTION_OUTPUT_DIR=`realpath $SCRIPT_DIR/../results/2_FeatureSelection`
FEATURE_REDUCTION_OUTPUT_DIR=`realpath $SCRIPT_DIR/../results/3_FeatureReduction`
FEATURE_CLASSIFICATION_DIR=`realpath $SCRIPT_DIR/../results/4_Classification`
FILE_LOCATION_10M_BAMS=`cd $SCRIPT_DIR/../../10m;pwd`

cd $FILE_LOCATION_10M_BAMS
ls *.10m.nodup.q30.bam | cut -f 1 -d . | xargs -n 1 -P 8 -I %1 $SCRIPT_DIR/tab.sh  trim_q30_cpm %1
