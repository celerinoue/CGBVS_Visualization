#!/bin/sh

## set model
DATASET_NAME=kinase_chembl
MAX_LEN_SEQ=4128
MODEL=model_10


## build dataset
#sh ./script/A_processing/make_dataset/build_dataset_train.sh $DATASET_NAME $MAX_LEN_SEQ


## cross validation
sh ./script/A_processing/make_dataset/split_cv.sh $DATASET_NAME $MAX_LEN_SEQ $MODEL

## visualization
## データセットの中身 可視化のためのコードとか
