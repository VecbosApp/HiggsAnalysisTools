#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

echo "Now merging EE datasets for mass $Hmass ..."
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/merged/ggH$i\2L2Nu_ee.root results/Summer11_V1/GluGluToHToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetEE.root
  hadd results/merged/ggH$i\LNuTauNu_ee.root results/Summer11_V1/GluGluToHToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetEE.root
  hadd results/merged/qqH$i\2L2Nu_ee.root results/Summer11_V1/VBF_HToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetEE.root
  hadd results/merged/qqH$i\LNuTauNu_ee.root results/Summer11_V1/VBF_HToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetEE.root
done
hadd results/merged/TTbar_ee.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetEE.root
#hadd results/merged/SingleTop_sChannel_ee.root results/Summer11_V1/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEE.root
#hadd results/merged/SingleTop_tChannel_ee.root results/Summer11_V1/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEE.root
#hadd results/merged/SingleTop_tWChannel_ee.root results/Summer11_V1/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEE.root
hadd results/merged/SingleTop_sChannel_ee.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetEE.root
hadd results/merged/SingleTop_tChannel_ee.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEE.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEE.root
hadd results/merged/SingleTop_tWChannel_ee.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEE.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEE.root
hadd results/merged/Wjets_ee.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetEE.root
hadd results/merged/Zee_Hi_ee.root results/Summer11_V1/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetEE.root
hadd results/merged/Zmm_Hi_ee.root results/Summer11_V1/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetEE.root
hadd results/merged/Ztautau_Hi_ee.root results/Summer11_V1/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetEE.root
hadd results/merged/Zee_Lo_ee.root results/Summer11_V1/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Zmm_Lo_ee.root results/Summer11_V1/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetEE.root
hadd results/merged/Ztautau_Lo_ee.root results/Summer11_V1/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/*datasetEE.root
hadd results/merged/WW_ee.root results/Summer11_V1/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/*datasetEE.root 
hadd results/merged/ggWW_ee.root results/Summer11_V1/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetEE.root
hadd results/merged/WZ_ee.root results/Summer11_V1/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/*datasetEE.root
hadd results/merged/ZZ_ee.root results/Summer11_V1/ZZ_TuneZ2_7TeV_pythia6_tauola/*datasetEE.root
hadd results/merged/Vgamma_ee.root results/Summer11_V1/GVJets_7TeV-madgraph/*datasetEE.root

echo "Now merging MM datasets..."
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/merged/ggH$i\2L2Nu_mm.root results/Summer11_V1/GluGluToHToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetMM.root
  hadd results/merged/ggH$i\LNuTauNu_mm.root results/Summer11_V1/GluGluToHToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetMM.root
  hadd results/merged/qqH$i\2L2Nu_mm.root results/Summer11_V1/VBF_HToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetMM.root
  hadd results/merged/qqH$i\LNuTauNu_mm.root results/Summer11_V1/VBF_HToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetMM.root
done
hadd results/merged/TTbar_mm.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetMM.root
#hadd results/merged/SingleTop_sChannel_mm.root results/Summer11_V1/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetMM.root
#hadd results/merged/SingleTop_tChannel_mm.root results/Summer11_V1/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetMM.root
#hadd results/merged/SingleTop_tWChannel_mm.root results/Summer11_V1/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetMM.root
hadd results/merged/SingleTop_sChannel_mm.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetMM.root
hadd results/merged/SingleTop_tChannel_mm.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetMM.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetMM.root
hadd results/merged/SingleTop_tWChannel_mm.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetMM.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetMM.root
hadd results/merged/Wjets_mm.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetMM.root
hadd results/merged/Zee_Hi_mm.root results/Summer11_V1/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetMM.root
hadd results/merged/Zmm_Hi_mm.root results/Summer11_V1/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetMM.root
hadd results/merged/Ztautau_Hi_mm.root results/Summer11_V1/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetMM.root
hadd results/merged/Zee_Lo_mm.root results/Summer11_V1/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Zmm_Lo_mm.root results/Summer11_V1/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetMM.root
hadd results/merged/Ztautau_Lo_mm.root results/Summer11_V1/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/*datasetMM.root
hadd results/merged/WW_mm.root results/Summer11_V1/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/*datasetMM.root
hadd results/merged/ggWW_mm.root results/Summer11_V1/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetMM.root
hadd results/merged/WZ_mm.root results/Summer11_V1/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/*datasetMM.root
hadd results/merged/ZZ_mm.root results/Summer11_V1/ZZ_TuneZ2_7TeV_pythia6_tauola/*datasetMM.root
hadd results/merged/Vgamma_mm.root results/Summer11_V1/GVJets_7TeV-madgraph/*datasetMM.root

echo "Now merging EM datasets..."
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/merged/ggH$i\2L2Nu_em.root results/Summer11_V1/GluGluToHToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetEM.root
  hadd results/merged/ggH$i\LNuTauNu_em.root results/Summer11_V1/GluGluToHToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetEM.root
  hadd results/merged/qqH$i\2L2Nu_em.root results/Summer11_V1/VBF_HToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetEM.root
  hadd results/merged/qqH$i\LNuTauNu_em.root results/Summer11_V1/VBF_HToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetEM.root
done
hadd results/merged/TTbar_em.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetEM.root
#hadd results/merged/SingleTop_sChannel_em.root results/Summer11_V1/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetEM.root
#hadd results/merged/SingleTop_tChannel_em.root results/Summer11_V1/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetEM.root
#hadd results/merged/SingleTop_tWChannel_em.root results/Summer11_V1/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetEM.root
hadd results/merged/SingleTop_sChannel_em.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetEM.root
hadd results/merged/SingleTop_tChannel_em.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEM.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetEM.root
hadd results/merged/SingleTop_tWChannel_em.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEM.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetEM.root
hadd results/merged/Wjets_em.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetEM.root
hadd results/merged/Zee_Hi_em.root results/Summer11_V1/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetEM.root
hadd results/merged/Zmm_Hi_em.root results/Summer11_V1/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetEM.root
hadd results/merged/Ztautau_Hi_em.root results/Summer11_V1/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetEM.root
hadd results/merged/Zee_Lo_em.root results/Summer11_V1/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Zmm_Lo_em.root results/Summer11_V1/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetEM.root
hadd results/merged/Ztautau_Lo_em.root results/Summer11_V1/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/*datasetEM.root
hadd results/merged/WW_em.root results/Summer11_V1/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/*datasetEM.root
hadd results/merged/ggWW_em.root results/Summer11_V1/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetEM.root
hadd results/merged/WZ_em.root results/Summer11_V1/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/*datasetEM.root
hadd results/merged/ZZ_em.root results/Summer11_V1/ZZ_TuneZ2_7TeV_pythia6_tauola/*datasetEM.root
hadd results/merged/Vgamma_em.root results/Summer11_V1/GVJets_7TeV-madgraph/*datasetEM.root

echo "Now merging ME datasets..."
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/merged/ggH$i\2L2Nu_me.root results/Summer11_V1/GluGluToHToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetME.root
  hadd results/merged/ggH$i\LNuTauNu_me.root results/Summer11_V1/GluGluToHToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetME.root
  hadd results/merged/qqH$i\2L2Nu_me.root results/Summer11_V1/VBF_HToWWTo2L2Nu_M-$i\_7TeV-powheg-pythia6/*datasetME.root
  hadd results/merged/qqH$i\LNuTauNu_me.root results/Summer11_V1/VBF_HToWWToLNuTauNu_M-$i\_7TeV-powheg-pythia6/*datasetME.root
done
hadd results/merged/TTbar_me.root results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*datasetME.root
#hadd results/merged/SingleTop_sChannel_me.root results/Summer11_V1/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*datasetME.root
#hadd results/merged/SingleTop_tChannel_me.root results/Summer11_V1/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*datasetME.root
#hadd results/merged/SingleTop_tWChannel_me.root results/Summer11_V1/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*datasetME.root
hadd results/merged/SingleTop_sChannel_me.root results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*datasetME.root
hadd results/merged/SingleTop_tChannel_me.root results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetME.root results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*datasetME.root
hadd results/merged/SingleTop_tWChannel_me.root results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetME.root results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*datasetME.root
hadd results/merged/Wjets_me.root results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*datasetME.root
hadd results/merged/Zee_Hi_me.root results/Summer11_V1/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetME.root
hadd results/merged/Zmm_Hi_me.root results/Summer11_V1/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/*datasetME.root
hadd results/merged/Ztautau_Hi_me.root results/Summer11_V1/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*datasetME.root
hadd results/merged/Zee_Lo_me.root results/Summer11_V1/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Zmm_Lo_me.root results/Summer11_V1/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*datasetME.root
hadd results/merged/Ztautau_Lo_me.root results/Summer11_V1/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola/*datasetME.root
hadd results/merged/WW_me.root results/Summer11_V1/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/*datasetME.root
hadd results/merged/ggWW_me.root results/Summer11_V1/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*datasetME.root
hadd results/merged/WZ_me.root results/Summer11_V1/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/*datasetME.root
hadd results/merged/ZZ_me.root results/Summer11_V1/ZZ_TuneZ2_7TeV_pythia6_tauola/*datasetME.root
hadd results/merged/Vgamma_me.root results/Summer11_V1/GVJets_7TeV-madgraph/*datasetME.root

