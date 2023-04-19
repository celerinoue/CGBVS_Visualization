## config
TRAIN_DATASET=kinase_chembl
MODEL=model_10-cv
GPU=0


## train
#<< COMMENTOUT
#sh script/B_run/train_cv_chembl.sh $TRAIN_DATASET $MODEL $GPU
#COMMENTOUT


## train_cv
sh ./script/B_run/train_cv_chembl.sh $TRAIN_DATASET $MODEL $GPU
