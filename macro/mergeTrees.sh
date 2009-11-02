#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results/merged

echo "Now merging EE datasets..."
hadd results/merged/H120_ee.root results/H120/*datasetEE.root
hadd results/merged/H130_ee.root results/H130/*datasetEE.root
hadd results/merged/H140_ee.root results/H140/*datasetEE.root
hadd results/merged/H150_ee.root results/H150/*datasetEE.root
hadd results/merged/H155_ee.root results/H155/*datasetEE.root
hadd results/merged/H160_ee.root results/H160/*datasetEE.root
hadd results/merged/H165_ee.root results/H165/*datasetEE.root
hadd results/merged/H170_ee.root results/H170/*datasetEE.root
hadd results/merged/H175_ee.root results/H175/*datasetEE.root
hadd results/merged/H180_ee.root results/H180/*datasetEE.root
hadd results/merged/H190_ee.root results/H190/*datasetEE.root
hadd results/merged/H200_ee.root results/H200/*datasetEE.root
hadd results/merged/H210_ee.root results/H210/*datasetEE.root
hadd results/merged/H220_ee.root results/H220/*datasetEE.root
hadd results/merged/H230_ee.root results/H230/*datasetEE.root
hadd results/merged/H240_ee.root results/H240/*datasetEE.root
hadd results/merged/H250_ee.root results/H250/*datasetEE.root
hadd results/merged/H275_ee.root results/H275/*datasetEE.root
hadd results/merged/H300_ee.root results/H300/*datasetEE.root
hadd results/merged/H350_ee.root results/H350/*datasetEE.root
hadd results/merged/H400_ee.root results/H400/*datasetEE.root
hadd results/merged/H450_ee.root results/H450/*datasetEE.root
hadd results/merged/H500_ee.root results/H500/*datasetEE.root
hadd results/merged/H550_ee.root results/H550/*datasetEE.root
hadd results/merged/H600_ee.root results/H600/*datasetEE.root

hadd results/merged/PhotonJet_Pt0to15_ee.root results/PhotonJet_Pt0to15/*datasetEE.root
hadd results/merged/PhotonJet_Pt120to170_ee.root  results/PhotonJet_Pt120to170/*datasetEE.root
hadd results/merged/PhotonJet_Pt15to20_ee.root results/PhotonJet_Pt15to20/*datasetEE.root
hadd results/merged/PhotonJet_Pt170to300_ee.root results/PhotonJet_Pt170to300/*datasetEE.root
hadd results/merged/PhotonJet_Pt20to30_ee.root results/PhotonJet_Pt20to30/*datasetEE.root
hadd results/merged/PhotonJet_Pt300to500_ee.root results/PhotonJet_Pt300to500/*datasetEE.root
hadd results/merged/PhotonJet_Pt30to50_ee.root results/PhotonJet_Pt30to50/*datasetEE.root
hadd results/merged/PhotonJet_Pt500toInf_ee.root results/PhotonJet_Pt500toInf/*datasetEE.root
hadd results/merged/PhotonJet_Pt50to80_ee.root results/PhotonJet_Pt50to80/*datasetEE.root
hadd results/merged/PhotonJet_Pt80to120_ee.root results/PhotonJet_Pt80to120/*datasetEE.root
hadd results/merged/QCD_BCtoE_Pt20to30Ele10_ee.root results/QCD_BCtoE_Pt20to30Ele10/*datasetEE.root
hadd results/merged/QCD_BCtoE_Pt30to80Ele10_ee.root results/QCD_BCtoE_Pt30to80Ele10/*datasetEE.root
hadd results/merged/QCD_BCtoE_Pt80to170Ele10_ee.root results/QCD_BCtoE_Pt80to170Ele10/*datasetEE.root
hadd results/merged/QCD_EMEnriched_Pt20to30Ele10_ee.root results/QCD_EMEnriched_Pt20to30Ele10/*datasetEE.root
hadd results/merged/QCD_EMEnriched_Pt30to80Ele10_ee.root results/QCD_EMEnriched_Pt30to80Ele10/*datasetEE.root
hadd results/merged/QCD_EMEnriched_Pt80to170Ele10_ee.root results/QCD_EMEnriched_Pt80to170Ele10/*datasetEE.root
hadd results/merged/TTbar_ee.root results/TTbar/*datasetEE.root
hadd results/merged/TTbarJetsMadgraph_ee.root results/TTbarJetsMadgraph/*datasetEE.root
hadd results/merged/Wgamma_ee.root results/Wgamma/*datasetEE.root
hadd results/merged/WjetsMadgraph_ee.root results/WjetsMadgraph/*datasetEE.root
hadd results/merged/WW_ee.root results/WW/*datasetEE.root
hadd results/merged/WZ_ee.root results/WZ/*datasetEE.root
hadd results/merged/ZjetsMadgraph_ee.root results/ZjetsMadgraph/*datasetEE.root
hadd results/merged/ZZ_ee.root results/ZZ/*datasetEE.root


echo "Now merging MM datasets..."
hadd results/merged/H120_mm.root results/H120/*datasetMM.root
hadd results/merged/H130_mm.root results/H130/*datasetMM.root
hadd results/merged/H140_mm.root results/H140/*datasetMM.root
hadd results/merged/H150_mm.root results/H150/*datasetMM.root
hadd results/merged/H155_mm.root results/H155/*datasetMM.root
hadd results/merged/H160_mm.root results/H160/*datasetMM.root
hadd results/merged/H165_mm.root results/H165/*datasetMM.root
hadd results/merged/H170_mm.root results/H170/*datasetMM.root
hadd results/merged/H175_mm.root results/H175/*datasetMM.root
hadd results/merged/H180_mm.root results/H180/*datasetMM.root
hadd results/merged/H190_mm.root results/H190/*datasetMM.root
hadd results/merged/H200_mm.root results/H200/*datasetMM.root
hadd results/merged/H210_mm.root results/H210/*datasetMM.root
hadd results/merged/H220_mm.root results/H220/*datasetMM.root
hadd results/merged/H230_mm.root results/H230/*datasetMM.root
hadd results/merged/H240_mm.root results/H240/*datasetMM.root
hadd results/merged/H250_mm.root results/H250/*datasetMM.root
hadd results/merged/H275_mm.root results/H275/*datasetMM.root
hadd results/merged/H300_mm.root results/H300/*datasetMM.root
hadd results/merged/H350_mm.root results/H350/*datasetMM.root
hadd results/merged/H400_mm.root results/H400/*datasetMM.root
hadd results/merged/H450_mm.root results/H450/*datasetMM.root
hadd results/merged/H500_mm.root results/H500/*datasetMM.root
hadd results/merged/H550_mm.root results/H550/*datasetMM.root
hadd results/merged/H600_mm.root results/H600/*datasetMM.root

hadd results/merged/InclusiveMu15_mm.root results/InclusiveMu15/*datasetMM.root
hadd results/merged/TTbar_mm.root results/TTbar/*datasetMM.root
hadd results/merged/TTbarJetsMadgraph_mm.root results/TTbarJetsMadgraph/*datasetMM.root
hadd results/merged/Wgamma_mm.root results/Wgamma/*datasetMM.root
hadd results/merged/WjetsMadgraph_mm.root results/WjetsMadgraph/*datasetMM.root
hadd results/merged/WW_mm.root results/WW/*datasetMM.root
hadd results/merged/WZ_mm.root results/WZ/*datasetMM.root
hadd results/merged/ZjetsMadgraph_mm.root results/ZjetsMadgraph/*datasetMM.root
hadd results/merged/ZZ_mm.root results/ZZ/*datasetMM.root

echo "Now merging EM datasets..."
hadd results/merged/H120_em.root results/H120/*datasetEM.root
hadd results/merged/H130_em.root results/H130/*datasetEM.root
hadd results/merged/H140_em.root results/H140/*datasetEM.root
hadd results/merged/H150_em.root results/H150/*datasetEM.root
hadd results/merged/H155_em.root results/H155/*datasetEM.root
hadd results/merged/H160_em.root results/H160/*datasetEM.root
hadd results/merged/H165_em.root results/H165/*datasetEM.root
hadd results/merged/H170_em.root results/H170/*datasetEM.root
hadd results/merged/H175_em.root results/H175/*datasetEM.root
hadd results/merged/H180_em.root results/H180/*datasetEM.root
hadd results/merged/H190_em.root results/H190/*datasetEM.root
hadd results/merged/H200_em.root results/H200/*datasetEM.root
hadd results/merged/H210_em.root results/H210/*datasetEM.root
hadd results/merged/H220_em.root results/H220/*datasetEM.root
hadd results/merged/H230_em.root results/H230/*datasetEM.root
hadd results/merged/H240_em.root results/H240/*datasetEM.root
hadd results/merged/H250_em.root results/H250/*datasetEM.root
hadd results/merged/H275_em.root results/H275/*datasetEM.root
hadd results/merged/H300_em.root results/H300/*datasetEM.root
hadd results/merged/H350_em.root results/H350/*datasetEM.root
hadd results/merged/H400_em.root results/H400/*datasetEM.root
hadd results/merged/H450_em.root results/H450/*datasetEM.root
hadd results/merged/H500_em.root results/H500/*datasetEM.root
hadd results/merged/H550_em.root results/H550/*datasetEM.root
hadd results/merged/H600_em.root results/H600/*datasetEM.root

hadd results/merged/InclusiveMu15_em.root results/InclusiveMu15/*datasetEM.root
hadd results/merged/PhotonJet_Pt0to15_em.root results/PhotonJet_Pt0to15/*datasetEM.root
hadd results/merged/PhotonJet_Pt120to170_em.root  results/PhotonJet_Pt120to170/*datasetEM.root
hadd results/merged/PhotonJet_Pt15to20_em.root results/PhotonJet_Pt15to20/*datasetEM.root
hadd results/merged/PhotonJet_Pt170to300_em.root results/PhotonJet_Pt170to300/*datasetEM.root
hadd results/merged/PhotonJet_Pt20to30_em.root results/PhotonJet_Pt20to30/*datasetEM.root
hadd results/merged/PhotonJet_Pt300to500_em.root results/PhotonJet_Pt300to500/*datasetEM.root
hadd results/merged/PhotonJet_Pt30to50_em.root results/PhotonJet_Pt30to50/*datasetEM.root
hadd results/merged/PhotonJet_Pt500toInf_em.root results/PhotonJet_Pt500toInf/*datasetEM.root
hadd results/merged/PhotonJet_Pt50to80_em.root results/PhotonJet_Pt50to80/*datasetEM.root
hadd results/merged/PhotonJet_Pt80to120_em.root results/PhotonJet_Pt80to120/*datasetEM.root
hadd results/merged/QCD_BCtoE_Pt20to30Ele10_em.root results/QCD_BCtoE_Pt20to30Ele10/*datasetEM.root
hadd results/merged/QCD_BCtoE_Pt30to80Ele10_em.root results/QCD_BCtoE_Pt30to80Ele10/*datasetEM.root
hadd results/merged/QCD_BCtoE_Pt80to170Ele10_em.root results/QCD_BCtoE_Pt80to170Ele10/*datasetEM.root
hadd results/merged/QCD_EMEnriched_Pt20to30Ele10_em.root results/QCD_EMEnriched_Pt20to30Ele10/*datasetEM.root
hadd results/merged/QCD_EMEnriched_Pt30to80Ele10_em.root results/QCD_EMEnriched_Pt30to80Ele10/*datasetEM.root
hadd results/merged/QCD_EMEnriched_Pt80to170Ele10_em.root results/QCD_EMEnriched_Pt80to170Ele10/*datasetEM.root
hadd results/merged/TTbar_em.root results/TTbar/*datasetEM.root
hadd results/merged/TTbarJetsMadgraph_em.root results/TTbarJetsMadgraph/*datasetEM.root
hadd results/merged/Wgamma_em.root results/Wgamma/*datasetEM.root
hadd results/merged/WjetsMadgraph_em.root results/WjetsMadgraph/*datasetEM.root
hadd results/merged/WW_em.root results/WW/*datasetEM.root
hadd results/merged/WZ_em.root results/WZ/*datasetEM.root
hadd results/merged/ZjetsMadgraph_em.root results/ZjetsMadgraph/*datasetEM.root
hadd results/merged/ZZ_em.root results/ZZ/*datasetEM.root
