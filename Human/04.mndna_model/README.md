# This file shows how to use the pipeline module to reproduce the results of the paper.

# 1. directory structure
```text
├── 10m
├── modelData
│   └── empty
│       ├── cleaned
│       └── origin
├── results
│   ├── 2_FeatureSelection
│   ├── 3_FeatureReduction
│   └── 4_Classification
│       ├── crc
│       └── panca
└── script 
```
You should put 10m data in the 10m directory, and the module data in the modelData directory. The results will be saved in the results directory.

# 2. run the pipeline
```bash
# Build 10k cpm data
./script/step1.sh

# crc pipeline
# 1. Feature selection from 10m features
./script/step2.sh crc
# 2. Feature reduction
./script/step3.sh crc
# 3. Classification
./script/step4.sh crc

# panca pipeline
# 1. Feature selection from 10m features
./script/step2.sh panca
# 2. Feature reduction
./script/step3.sh panca
# 3. Classification
./script/step4.sh panca

``` 
Feature selection results will be saved in the 2_FeatureSelection directory, feature reduction results will be saved in the 3_FeatureReduction directory, and classification results will be saved in the 4_Classification directory.
