#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results_data/merged

echo "Now merging EE datasets..."
hadd results_data/merged/dataset_ee.root results_data/Data7TeV/DoubleElectron/*datasetEE.root

echo "Now merging MM datasets..."
hadd results_data/merged/dataset_mm.root results_data/Data7TeV/DoubleMu/*datasetMM.root results_data/Data7TeV/SingleMu/*datasetMM.root

echo "Now merging EM datasets..."
hadd results_data/merged/dataset_em.root results_data/Data7TeV/MuEG/*datasetEM.root  results_data/Data7TeV/SingleMu/*datasetEM.root

echo "Now merging ME datasets..."
hadd results_data/merged/dataset_me.root results_data/Data7TeV/MuEG/*datasetME.root results_data/Data7TeV/SingleMu/*datasetME.root
