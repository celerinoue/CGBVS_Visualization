#!/bin/sh

# set model
DATA_DIR=$1
DATASET=$2
MAX_LEN_SEQ=$3
SEQ_NAME=$4

# mkdir
mkdir -p data/dataset/${DATASET}_test
mkdir -p log/${DATASET}/build_dataset
mkdir -p data/dataset/${DATASET}_test/${SEQ_NAME}


# run
kgcn-chem \
--assay_dir ${DATA_DIR}/${SEQ_NAME} \
-a 50 \
--output data/dataset/${DATASET}_test/${SEQ_NAME}/${DATASET}_test_${SEQ_NAME}.jbl \
--multimodal \
--no_pseudo_negative \
--max_len_seq $MAX_LEN_SEQ \
> log/${DATASET}/build_dataset/build_dataset_test_${SEQ_NAME}.log 2>&1


#
mv multimodal_data_index.csv \
data/dataset/${DATASET}_test/${SEQ_NAME}/multimodal_data_index.csv
