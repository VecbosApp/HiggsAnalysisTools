#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

# usage: ./mergeTrees.sh dataset ("eg", "mu")

dataset=$1

mkdir -p results_data/merged

if [ "$dataset" = "eg" ]; then

    echo "Now merging EE datasets..."
    hadd results_data/merged/dataset_eg_ee.root results_data/dataset_eg_Sep3rdReReco/*datasetEE.root results_data/PDElectron_11pbTo40pb/*datasetEE.root

    echo "Now merging EM datasets..."
    hadd results_data/merged/dataset_eg_em.root results_data/dataset_eg_Sep3rdReReco/*datasetEM.root results_data/PDElectron_11pbTo40pb/*datasetEM.root
elif [ "$dataset" = "mu" ]; then

    echo "Now merging MM datasets..."
    hadd results_data/merged/dataset_mu_mm.root results_data/Data7TeV/dataset_mu/*datasetMM.root

    echo "Now merging EM datasets..."
    hadd results_data/merged/dataset_mu_em.root results_data/Data7TeV/dataset_mu/*datasetEM.root
else
    echo "dataset can be only eg or mu. Nothing done."
fi
