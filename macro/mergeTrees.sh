#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

echo "Now merging EE datasets..."
hadd results/merged/H130_ee.root results/HiggsWW/H130_2W_2lnu_gluonfusion_7TeV/*datasetEE.root
hadd results/merged/H160_ee.root results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*datasetEE.root
hadd results/merged/H190_ee.root results/HiggsWW/H190_2W_2lnu_gluonfusion_7TeV/*datasetEE.root
hadd results/merged/TTbar_ee.root results/TTbar/TTbarJets-madgraph/*datasetEE.root
hadd results/merged/SingleTop_sChannel_ee.root results/SingleTop/SingleTop_sChannel-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tChannel_ee.root results/SingleTop/SingleTop_tChannel-madgraph/*datasetEE.root
hadd results/merged/SingleTop_tWChannel_ee.root results/SingleTop/SingleTop_tWChannel-madgraph/*datasetEE.root
hadd results/merged/Wgamma_ee.root results/DiBosons/Wgamma/*datasetEE.root
hadd results/merged/WjetsMadgraph_ee.root results/WJetsMADGRAPH/WJets-madgraph/*datasetEE.root
hadd results/merged/WW_ee.root results/DiBosons/WW_2l_7TeV/*datasetEE.root
hadd results/merged/WZ_ee.root results/DiBosons/WZ_3l_7TeV/*datasetEE.root
hadd results/merged/ZjetsMadgraph_ee.root results/ZJetsMADGRAPH/ZJets-madgraph/*datasetEE.root
hadd results/merged/ZZ_ee.root results/DiBosons/ZZ_2l2nu/*datasetEE.root


echo "Now merging MM datasets..."
hadd results/merged/H130_mm.root results/HiggsWW/H130_2W_2lnu_gluonfusion_7TeV/*datasetMM.root
hadd results/merged/H160_mm.root results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*datasetMM.root
hadd results/merged/H190_mm.root results/HiggsWW/H190_2W_2lnu_gluonfusion_7TeV/*datasetMM.root
hadd results/merged/TTbar_mm.root results/TTbar/TTbarJets-madgraph/*datasetMM.root
hadd results/merged/SingleTop_sChannel_mm.root results/SingleTop/SingleTop_sChannel-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tChannel_mm.root results/SingleTop/SingleTop_tChannel-madgraph/*datasetMM.root
hadd results/merged/SingleTop_tWChannel_mm.root results/SingleTop/SingleTop_tWChannel-madgraph/*datasetMM.root
hadd results/merged/Wgamma_mm.root results/DiBosons/Wgamma/*datasetMM.root
hadd results/merged/WjetsMadgraph_mm.root results/WJetsMADGRAPH/WJets-madgraph/*datasetMM.root
hadd results/merged/WW_mm.root results/DiBosons/WW_2l_7TeV/*datasetMM.root
hadd results/merged/WZ_mm.root results/DiBosons/WZ_3l_7TeV/*datasetMM.root
hadd results/merged/ZjetsMadgraph_mm.root results/ZJetsMADGRAPH/ZJets-madgraph/*datasetMM.root
hadd results/merged/ZZ_mm.root results/DiBosons/ZZ_2l2nu/*datasetMM.root


echo "Now merging EM datasets..."
hadd results/merged/H130_em.root results/HiggsWW/H130_2W_2lnu_gluonfusion_7TeV/*datasetEM.root
hadd results/merged/H160_em.root results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*datasetEM.root
hadd results/merged/H190_em.root results/HiggsWW/H190_2W_2lnu_gluonfusion_7TeV/*datasetEM.root
hadd results/merged/TTbar_em.root results/TTbar/TTbarJets-madgraph/*datasetEM.root
hadd results/merged/SingleTop_sChannel_em.root results/SingleTop/SingleTop_sChannel-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tChannel_em.root results/SingleTop/SingleTop_tChannel-madgraph/*datasetEM.root
hadd results/merged/SingleTop_tWChannel_em.root results/SingleTop/SingleTop_tWChannel-madgraph/*datasetEM.root
hadd results/merged/Wgamma_em.root results/DiBosons/Wgamma/*datasetEM.root
hadd results/merged/WjetsMadgraph_em.root results/WJetsMADGRAPH/WJets-madgraph/*datasetEM.root
hadd results/merged/WW_em.root results/DiBosons/WW_2l_7TeV/*datasetEM.root
hadd results/merged/WZ_em.root results/DiBosons/WZ_3l_7TeV/*datasetEM.root
hadd results/merged/ZjetsMadgraph_em.root results/ZJetsMADGRAPH/ZJets-madgraph/*datasetEM.root
hadd results/merged/ZZ_em.root results/DiBosons/ZZ_2l2nu/*datasetEM.root


