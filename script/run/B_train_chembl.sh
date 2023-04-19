## config
TRAIN_DATASET=kinase_chembl
MODEL="model_10"
GPU=1


## train
#<< COMMENTOUT
#sh script/run_py/train_cv_chembl.sh $TRAIN_DATASET $MODEL $GPU
#COMMENTOUT


## train_cv
#sh ./script/B_run/train_cv_chembl.sh $TRAIN_DATASET ${MODEL}-cv $GPU

##
sh script/B_run/train_chembl.sh $TRAIN_DATASET $MODEL $GPU

#kgcn train_cv --config setting/config/tmp/tmp_model_10.json --dataset data/dataset/kinase_chembl/kinase_chembl.jbl --gpu 1 > log/kinase_chembl/train/model_10/train_cv.log 2>&1
