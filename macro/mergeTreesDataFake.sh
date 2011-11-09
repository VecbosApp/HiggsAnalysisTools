#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results_fake_ee/merged

echo "Now merging EE datasets..."
hadd results_fake_ee/merged/dataset_DoubleElectron_ee.root          results_fake_ee/DoubleElectron/*datasetEE.root
hadd results_fake_ee/merged/dataset_DoubleElectron_Run2011B_ee.root results_fake_ee/DoubleElectron_Run2011B/*datasetEE.root

hadd results_fake_ee/merged/dataset_SingleElectron_ee.root          results_fake_ee/SingleElectron/*datasetEE.root
hadd results_fake_ee/merged/dataset_SingleElectron_Run2011B_ee.root results_fake_ee/SingleElectron_Run2011B/*datasetEE.root

#echo "Now merging ME datasets..."
#hadd results_fake_ee/merged/dataset_MuEG_me.root results_fake_ee/MuEG/*datasetME.root 
#hadd results_fake_ee/merged/dataset_SingleMu_me.root results_fake_ee/SingleMu/*datasetME.root 
#hadd results_fake_ee/merged/dataset_SingleElectron_me.root results_fake_ee/SingleElectron/*datasetME.root

