#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

echo "Now merging EE datasets..."
hadd results/merged/H160_ee.root results/HiggsWW/GluGluToHToWWTo2L2Nu_M-160/*datasetEE.root
hadd results/merged/TTbar_ee.root results/TTbar/TTJets_TuneD6T/*datasetEE.root
hadd results/merged/SingleTop_sChannel_ee.root results/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tChannel_ee.root results/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tWChannel_ee.root results/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/Wenu_ee.root results/WPYTHIA/WToENu_TuneZ2/*datasetEE.root
hadd results/merged/Wmunu_ee.root results/WPYTHIA/WToMuNu_TuneZ2/*datasetEE.root
hadd results/merged/Wtaunu_ee.root results/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*datasetEE.root
hadd results/merged/Zee_Hi_ee.root results/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*datasetEE.root
hadd results/merged/Zmm_Hi_ee.root results/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*datasetEE.root
hadd results/merged/Ztautau_Hi_ee.root results/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetEE.root
hadd results/merged/Zee_Lo_ee.root results/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*datasetEE.root
hadd results/merged/Zmm_Lo_ee.root results/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*datasetEE.root
hadd results/merged/Ztautau_Lo_ee.root results/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*datasetEE.root
hadd results/merged/WW_ee.root results/DiBosons/WWTo2L2Nu_TuneZ2/*datasetEE.root
hadd results/merged/ggWW_ee.root results/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*datasetEE.root
hadd results/merged/WZ_ee.root results/DiBosons/WZTo3LNu_TuneZ2/*datasetEE.root
hadd results/merged/ZZ_ee.root results/DiBosons/ZZtoAnything_TuneZ2/*datasetEE.root
#hadd results/merged/Wgamma_ee.root results/DiBosons/Wgamma/*datasetEE.root

echo "Now merging MM datasets..."
hadd results/merged/H160_mm.root results/HiggsWW/GluGluToHToWWTo2L2Nu_M-160/*datasetMM.root
hadd results/merged/TTbar_mm.root results/TTbar/TTJets_TuneD6T/*datasetMM.root
hadd results/merged/SingleTop_sChannel_mm.root results/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tChannel_mm.root results/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tWChannel_mm.root results/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/Wenu_mm.root results/WPYTHIA/WToENu_TuneZ2/*datasetMM.root
hadd results/merged/Wmunu_mm.root results/WPYTHIA/WToMuNu_TuneZ2/*datasetMM.root
hadd results/merged/Wtaunu_mm.root results/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*datasetMM.root
hadd results/merged/Zee_Hi_mm.root results/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*datasetMM.root
hadd results/merged/Zmm_Hi_mm.root results/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*datasetMM.root
hadd results/merged/Ztautau_Hi_mm.root results/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetMM.root
hadd results/merged/Zee_Lo_mm.root results/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*datasetMM.root
hadd results/merged/Zmm_Lo_mm.root results/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*datasetMM.root
hadd results/merged/Ztautau_Lo_mm.root results/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*datasetMM.root
hadd results/merged/WW_mm.root results/DiBosons/WWTo2L2Nu_TuneZ2/*datasetMM.root
hadd results/merged/ggWW_mm.root results/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*datasetMM.root
hadd results/merged/WZ_mm.root results/DiBosons/WZTo3LNu_TuneZ2/*datasetMM.root
hadd results/merged/ZZ_mm.root results/DiBosons/ZZtoAnything_TuneZ2/*datasetMM.root
#hadd results/merged/Wgamma_mm.root results/DiBosons/Wgamma/*datasetMM.root

echo "Now merging EM datasets..."
hadd results/merged/H160_em.root results/HiggsWW/GluGluToHToWWTo2L2Nu_M-160/*datasetEM.root
hadd results/merged/TTbar_em.root results/TTbar/TTJets_TuneD6T/*datasetEM.root
hadd results/merged/SingleTop_sChannel_em.root results/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tChannel_em.root results/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tWChannel_em.root results/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/Wenu_em.root results/WPYTHIA/WToENu_TuneZ2/*datasetEM.root
hadd results/merged/Wmunu_em.root results/WPYTHIA/WToMuNu_TuneZ2/*datasetEM.root
hadd results/merged/Wtaunu_em.root results/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*datasetEM.root
hadd results/merged/Zee_Hi_em.root results/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*datasetEM.root
hadd results/merged/Zmm_Hi_em.root results/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*datasetEM.root
hadd results/merged/Ztautau_Hi_em.root results/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetEM.root
hadd results/merged/Zee_Lo_em.root results/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*datasetEM.root
hadd results/merged/Zmm_Lo_em.root results/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*datasetEM.root
hadd results/merged/Ztautau_Lo_em.root results/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*datasetEM.root
hadd results/merged/WW_em.root results/DiBosons/WWTo2L2Nu_TuneZ2/*datasetEM.root
hadd results/merged/ggWW_em.root results/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*datasetEM.root
hadd results/merged/WZ_em.root results/DiBosons/WZTo3LNu_TuneZ2/*datasetEM.root
hadd results/merged/ZZ_em.root results/DiBosons/ZZtoAnything_TuneZ2/*datasetEM.root
#hadd results/merged/Wgamma_em.root results/DiBosons/Wgamma/*datasetEM.root
