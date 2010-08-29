#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results_data/merged

echo "Now merging EE datasets..."
hadd results_data/merged/dataset_eg_ee.root results_data/Data7TeV/dataset_eg/*datasetEE.root

echo "Now merging MM datasets..."
# to be implemented on Mu dataset

echo "Now merging EM datasets..."
hadd results_data/merged/dataset_eg_em.root results_data/Data7TeV/dataset_eg/*datasetEM.root
# missing mu dataset


