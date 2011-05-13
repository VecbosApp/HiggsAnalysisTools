#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

Hmass=$1

mkdir -p results/merged

echo "Now merging EE datasets for mass $Hmass ..."
hadd results/merged/H$Hmass\_ee.root results/Spring11_V2/GluGluToHToWWTo2L2Nu_M-$Hmass\_7TeV-powheg-pythia6/*datasetEE.root
hadd results/merged/TTbar_ee.root results/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetEE.root
hadd results/merged/SingleTop_sChannel_ee.root results/Spring11_V2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tChannel_ee.root results/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tWChannel_ee.root results/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/Wjets_ee.root results/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetEE.root
hadd results/merged/Zee_Hi_ee.root results/Spring11_V2/DYToEE_M-20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Zmm_Hi_ee.root results/Spring11_V2/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Ztautau_Hi_ee.root results/Spring11_V2/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Zee_Lo_ee.root results/Spring11_V2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Zmm_Lo_ee.root results/Spring11_V2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Ztautau_Lo_ee.root results/Spring11_V2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/WW_ee.root results/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/ggWW_ee.root results/Spring11_V2/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetEE.root
hadd results/merged/WZ_ee.root results/Spring11_V2/WZTo3LNu_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/ZZ_ee.root results/Spring11_V2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*datasetEE.root
hadd results/merged/Wgamma_ee.root results/Spring11_V2/PhotonVJets_7TeV-madgraph/*datasetEE.root

echo "Now merging MM datasets..."
hadd results/merged/H$Hmass\_mm.root results/Spring11_V2/GluGluToHToWWTo2L2Nu_M-$Hmass\_7TeV-powheg-pythia6/*datasetMM.root
hadd results/merged/TTbar_mm.root results/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetMM.root
hadd results/merged/SingleTop_sChannel_mm.root results/Spring11_V2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tChannel_mm.root results/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tWChannel_mm.root results/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/Wjets_mm.root results/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetMM.root
hadd results/merged/Zee_Hi_mm.root results/Spring11_V2/DYToEE_M-20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Zmm_Hi_mm.root results/Spring11_V2/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Ztautau_Hi_mm.root results/Spring11_V2/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Zee_Lo_mm.root results/Spring11_V2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Zmm_Lo_mm.root results/Spring11_V2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Ztautau_Lo_mm.root results/Spring11_V2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/WW_mm.root results/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/ggWW_mm.root results/Spring11_V2/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetMM.root
hadd results/merged/WZ_mm.root results/Spring11_V2/WZTo3LNu_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/ZZ_mm.root results/Spring11_V2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*datasetMM.root
hadd results/merged/Wgamma_mm.root results/Spring11_V2/PhotonVJets_7TeV-madgraph/*datasetMM.root

echo "Now merging EM datasets..."
hadd results/merged/H$Hmass\_em.root results/Spring11_V2/GluGluToHToWWTo2L2Nu_M-$Hmass\_7TeV-powheg-pythia6/*datasetEM.root
hadd results/merged/TTbar_em.root results/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetEM.root
hadd results/merged/SingleTop_sChannel_em.root results/Spring11_V2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tChannel_em.root results/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tWChannel_em.root results/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/Wjets_em.root results/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetEM.root
hadd results/merged/Zee_Hi_em.root results/Spring11_V2/DYToEE_M-20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Zmm_Hi_em.root results/Spring11_V2/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Ztautau_Hi_em.root results/Spring11_V2/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Zee_Lo_em.root results/Spring11_V2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Zmm_Lo_em.root results/Spring11_V2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Ztautau_Lo_em.root results/Spring11_V2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/WW_em.root results/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/ggWW_em.root results/Spring11_V2/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetEM.root
hadd results/merged/WZ_em.root results/Spring11_V2/WZTo3LNu_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/ZZ_em.root results/Spring11_V2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*datasetEM.root
hadd results/merged/Wgamma_em.root results/Spring11_V2/PhotonVJets_7TeV-madgraph/*datasetEM.root

echo "Now merging ME datasets..."
hadd results/merged/H$Hmass\_me.root results/Spring11_V2/GluGluToHToWWTo2L2Nu_M-$Hmass\_7TeV-powheg-pythia6/*datasetME.root
hadd results/merged/TTbar_me.root results/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetME.root
hadd results/merged/SingleTop_sChannel_me.root results/Spring11_V2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetME.root
hadd results/merged/SingleTop_tChannel_me.root results/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetME.root
hadd results/merged/SingleTop_tWChannel_me.root results/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetME.root
hadd results/merged/Wjets_me.root results/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetME.root
hadd results/merged/Zee_Hi_me.root results/Spring11_V2/DYToEE_M-20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Zmm_Hi_me.root results/Spring11_V2/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Ztautau_Hi_me.root results/Spring11_V2/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Zee_Lo_me.root results/Spring11_V2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Zmm_Lo_me.root results/Spring11_V2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Ztautau_Lo_me.root results/Spring11_V2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/WW_me.root results/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/ggWW_me.root results/Spring11_V2/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetME.root
hadd results/merged/WZ_me.root results/Spring11_V2/WZTo3LNu_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/ZZ_me.root results/Spring11_V2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*datasetME.root
hadd results/merged/Wgamma_me.root results/Spring11_V2/PhotonVJets_7TeV-madgraph/*datasetME.root

