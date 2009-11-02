#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

echo "Adding weights for ee datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H120_ee.root", 0.000391868);
addWeights("results/merged/H130_ee.root", 0.000735827);
addWeights("results/merged/H140_ee.root", 0.00105772);
addWeights("results/merged/H150_ee.root", 0.00130271);
addWeights("results/merged/H155_ee.root", 0.00144067);
addWeights("results/merged/H160_ee.root", 0.00151712);
addWeights("results/merged/H165_ee.root", 0.00162467);
addWeights("results/merged/H170_ee.root", 0.00139067);
addWeights("results/merged/H175_ee.root", 0.00142666);
addWeights("results/merged/H180_ee.root", 0.00119266);
addWeights("results/merged/H190_ee.root", 0.000885987);
addWeights("results/merged/H200_ee.root", 0.000604428);
addWeights("results/merged/H210_ee.root", 0.000667591);
addWeights("results/merged/H220_ee.root", 0.000604428);
addWeights("results/merged/H230_ee.root", 0.000566099);
addWeights("results/merged/H240_ee.root", 0.000521067);
addWeights("results/merged/H250_ee.root", 0.000475435);
addWeights("results/merged/H275_ee.root", 0.000384236);
addWeights("results/merged/H300_ee.root", 0.000383913);
addWeights("results/merged/H350_ee.root", 0.000404639);
addWeights("results/merged/H400_ee.root", 0.000290823);
addWeights("results/merged/H450_ee.root", 0.000188096);
addWeights("results/merged/H500_ee.root", 0.000119614);
addWeights("results/merged/H550_ee.root", 8.35818e-05);
addWeights("results/merged/H600_ee.root", 5.42998e-05);

addWeights("results/merged/PhotonJet_Pt0to15_ee.root", 96084.9);
addWeights("results/merged/PhotonJet_Pt120to170_ee.root", 0.156946);
addWeights("results/merged/PhotonJet_Pt15to20_ee.root", 172.343);
addWeights("results/merged/PhotonJet_Pt170to300_ee.root", 0.0433001);
addWeights("results/merged/PhotonJet_Pt20to30_ee.root", 79.9604);
addWeights("results/merged/PhotonJet_Pt300to500_ee.root", 0.00408739);
addWeights("results/merged/PhotonJet_Pt30to50_ee.root", 21.9547);
addWeights("results/merged/PhotonJet_Pt500toInf_ee.root", 0.000316517);
addWeights("results/merged/PhotonJet_Pt50to80_ee.root", 4.58895);
addWeights("results/merged/PhotonJet_Pt80to120_ee.root", 0.786361);
addWeights("results/merged/QCD_BCtoE_Pt20to30Ele10_ee.root", 7.77958);
addWeights("results/merged/QCD_BCtoE_Pt30to80Ele10_ee.root", 11.7589);
addWeights("results/merged/QCD_BCtoE_Pt80to170Ele10_ee.root", 2.1881);
addWeights("results/merged/QCD_EMEnriched_Pt20to30Ele10_ee.root", 9.45263);
addWeights("results/merged/QCD_EMEnriched_Pt30to80Ele10_ee.root", 11.1153);
addWeights("results/merged/QCD_EMEnriched_Pt80to170Ele10_ee.root", 5.30233);
addWeights("results/merged/TTbar_ee.root", 0.291393);
addWeights("results/merged/TTbarJetsMadgraph_ee.root", 0.259565);
addWeights("results/merged/Wgamma_ee.root", 11.5203);
addWeights("results/merged/WjetsMadgraph_ee.root", 5.15638);
addWeights("results/merged/WW_ee.root", 0.000578358);
addWeights("results/merged/WZ_ee.root", 0.00034796);
addWeights("results/merged/ZjetsMadgraph_ee.root", 2.47585); 
addWeights("results/merged/ZZ_ee.root", 0.00023806);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H120_mm.root", 0.000391868);
addWeights("results/merged/H130_mm.root", 0.000735827);
addWeights("results/merged/H140_mm.root", 0.00105772);
addWeights("results/merged/H150_mm.root", 0.00130271);
addWeights("results/merged/H155_mm.root", 0.00144067);
addWeights("results/merged/H160_mm.root", 0.00151712);
addWeights("results/merged/H165_mm.root", 0.00162467);
addWeights("results/merged/H170_mm.root", 0.00139067);
addWeights("results/merged/H175_mm.root", 0.00142666);
addWeights("results/merged/H180_mm.root", 0.00119266);
addWeights("results/merged/H190_mm.root", 0.000885987);
addWeights("results/merged/H200_mm.root", 0.000604428);
addWeights("results/merged/H210_mm.root", 0.000667591);
addWeights("results/merged/H220_mm.root", 0.000604428);
addWeights("results/merged/H230_mm.root", 0.000566099);
addWeights("results/merged/H240_mm.root", 0.000521067);
addWeights("results/merged/H250_mm.root", 0.000475435);
addWeights("results/merged/H275_mm.root", 0.000384236);
addWeights("results/merged/H300_mm.root", 0.000383913);
addWeights("results/merged/H350_mm.root", 0.000404639);
addWeights("results/merged/H400_mm.root", 0.000290823);
addWeights("results/merged/H450_mm.root", 0.000188096);
addWeights("results/merged/H500_mm.root", 0.000119614);
addWeights("results/merged/H550_mm.root", 8.35818e-05);
addWeights("results/merged/H600_mm.root", 5.42998e-05);

addWeights("results/merged/InclusiveMu15_mm.root", 8.92846);
addWeights("results/merged/TTbar_mm.root", 0.291393);
addWeights("results/merged/TTbarJetsMadgraph_mm.root", 0.259565);
addWeights("results/merged/Wgamma_mm.root", 11.5203);
addWeights("results/merged/WjetsMadgraph_mm.root", 5.15638);
addWeights("results/merged/WW_mm.root", 0.000578358);
addWeights("results/merged/WZ_mm.root", 0.00034796);
addWeights("results/merged/ZjetsMadgraph_mm.root", 2.47585);
addWeights("results/merged/ZZ_mm.root", 0.00023806);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H120_em.root", 0.000391868);
addWeights("results/merged/H130_em.root", 0.000735827);
addWeights("results/merged/H140_em.root", 0.00105772);
addWeights("results/merged/H150_em.root", 0.00130271);
addWeights("results/merged/H155_em.root", 0.00144067);
addWeights("results/merged/H160_em.root", 0.00151712);
addWeights("results/merged/H165_em.root", 0.00162467);
addWeights("results/merged/H170_em.root", 0.00139067);
addWeights("results/merged/H175_em.root", 0.00142666);
addWeights("results/merged/H180_em.root", 0.00119266);
addWeights("results/merged/H190_em.root", 0.000885987);
addWeights("results/merged/H200_em.root", 0.000604428);
addWeights("results/merged/H210_em.root", 0.000667591);
addWeights("results/merged/H220_em.root", 0.000604428);
addWeights("results/merged/H230_em.root", 0.000566099);
addWeights("results/merged/H240_em.root", 0.000521067);
addWeights("results/merged/H250_em.root", 0.000475435);
addWeights("results/merged/H275_em.root", 0.000384236);
addWeights("results/merged/H300_em.root", 0.000383913);
addWeights("results/merged/H350_em.root", 0.000404639);
addWeights("results/merged/H400_em.root", 0.000290823);
addWeights("results/merged/H450_em.root", 0.000188096);
addWeights("results/merged/H500_em.root", 0.000119614);
addWeights("results/merged/H550_em.root", 8.35818e-05);
addWeights("results/merged/H600_em.root", 5.42998e-05);

addWeights("results/merged/InclusiveMu15_em.root", 8.92846);
addWeights("results/merged/PhotonJet_Pt0to15_em.root", 96084.9);
addWeights("results/merged/PhotonJet_Pt120to170_em.root", 0.156946);
addWeights("results/merged/PhotonJet_Pt15to20_em.root", 172.343);
addWeights("results/merged/PhotonJet_Pt170to300_em.root", 0.0433001);
addWeights("results/merged/PhotonJet_Pt20to30_em.root", 79.9604);
addWeights("results/merged/PhotonJet_Pt300to500_em.root", 0.00408739);
addWeights("results/merged/PhotonJet_Pt30to50_em.root", 21.9547);
addWeights("results/merged/PhotonJet_Pt500toInf_em.root", 0.000316517);
addWeights("results/merged/PhotonJet_Pt50to80_em.root", 4.58895);
addWeights("results/merged/PhotonJet_Pt80to120_em.root", 0.786361);
addWeights("results/merged/QCD_BCtoE_Pt20to30Ele10_em.root", 7.77958);
addWeights("results/merged/QCD_BCtoE_Pt30to80Ele10_em.root", 11.7589);
addWeights("results/merged/QCD_BCtoE_Pt80to170Ele10_em.root", 2.1881);
addWeights("results/merged/QCD_EMEnriched_Pt20to30Ele10_em.root", 9.45263);
addWeights("results/merged/QCD_EMEnriched_Pt30to80Ele10_em.root", 11.1153);
addWeights("results/merged/QCD_EMEnriched_Pt80to170Ele10_em.root", 5.30233);
addWeights("results/merged/TTbar_em.root", 0.291393);
addWeights("results/merged/TTbarJetsMadgraph_em.root", 0.259565);
addWeights("results/merged/Wgamma_em.root", 11.5203);
addWeights("results/merged/WjetsMadgraph_em.root", 5.15638);
addWeights("results/merged/WW_em.root", 0.000578358);
addWeights("results/merged/WZ_em.root", 0.00034796);
addWeights("results/merged/ZjetsMadgraph_em.root", 2.47585);
addWeights("results/merged/ZZ_em.root", 0.00023806);

.q

EOF

echo "done weighting."
