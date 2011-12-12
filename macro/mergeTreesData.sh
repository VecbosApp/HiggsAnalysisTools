#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results_data/merged

echo "Now merging EE datasets..."
hadd results_data/merged/dataset_DoubleElectron_ee.root results_data/Data7TeVHWW/DoubleElectron/*datasetEE.root results_data/Data7TeVHWW/DoubleElectron_Run2011B/*datasetEE.root
hadd results_data/merged/dataset_SingleElectron_ee.root results_data/Data7TeVHWW/SingleElectron/*datasetEE.root results_data/Data7TeVHWW/SingleElectron_Run2011B/*datasetEE.root

echo "Now merging MM datasets..."
hadd results_data/merged/dataset_DoubleMu_mm.root results_data/Data7TeVHWW/DoubleMu/*datasetMM.root results_data/Data7TeVHWW/DoubleMu_Run2011B/*datasetMM.root 
hadd results_data/merged/dataset_SingleMu_mm.root results_data/Data7TeVHWW/SingleMu/*datasetMM.root results_data/Data7TeVHWW/SingleMu_Run2011B/*datasetMM.root

echo "Now merging EM datasets..."
hadd results_data/merged/dataset_MuEG_em.root results_data/Data7TeVHWW/MuEG/*datasetEM.root results_data/Data7TeVHWW/MuEG_Run2011B/*datasetEM.root 
hadd results_data/merged/dataset_SingleMu_em.root results_data/Data7TeVHWW/SingleMu/*datasetEM.root results_data/Data7TeVHWW/SingleMu_Run2011B/*datasetEM.root 
hadd results_data/merged/dataset_SingleElectron_em.root results_data/Data7TeVHWW/SingleElectron/*datasetEM.root results_data/Data7TeVHWW/SingleElectron_Run2011B/*datasetEM.root

echo "Now merging ME datasets..."
hadd results_data/merged/dataset_MuEG_me.root results_data/Data7TeVHWW/MuEG/*datasetME.root results_data/Data7TeVHWW/MuEG_Run2011B/*datasetME.root 
hadd results_data/merged/dataset_SingleMu_me.root results_data/Data7TeVHWW/SingleMu/*datasetME.root results_data/Data7TeVHWW/SingleMu_Run2011B/*datasetME.root 
hadd results_data/merged/dataset_SingleElectron_me.root results_data/Data7TeVHWW/SingleElectron/*datasetME.root results_data/Data7TeVHWW/SingleElectron_Run2011B/*datasetME.root

