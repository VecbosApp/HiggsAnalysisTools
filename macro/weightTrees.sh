#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

# PS. weights are compute for 1000 pb-1

echo "Adding weights for ee datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H130_ee.root", 0.001 * 0.00452859);
addWeights("results/merged/H160_ee.root", 0.001 * 0.00818843);
addWeights("results/merged/H190_ee.root", 0.001 * 0.00486346);
addWeights("results/merged/TTbar_ee.root", 0.001 * 0.122721);
addWeights("results/merged/SingleTop_sChannel_ee.root", 0.001 * 0.00437115);
addWeights("results/merged/SingleTop_tChannel_ee.root", 0.001 * 0.0437332);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 0.001 * 0.644887);
addWeights("results/merged/Wgamma_ee.root", 0.001 * 0.390098);
addWeights("results/merged/WjetsMadgraph_ee.root", 0.001 * 3.20548);
addWeights("results/merged/WW_ee.root", 0.001 * 0.0203685);
addWeights("results/merged/WZ_ee.root", 0.001 * 0.00549947);
addWeights("results/merged/ZjetsMadgraph_ee.root", 0.001 * 2.94516); 
addWeights("results/merged/ZZ_ee.root", 0.001 * 0.00231513);
addWeights("results/merged/W1Jets_Pt0to100-alpgen_ee.root", 0.001 * 5.08688);
addWeights("results/merged/W1Jets_Pt100to300-alpgen_ee.root", 0.001 * 0.490897);
addWeights("results/merged/W1Jets_Pt300to800-alpgen_ee.root", 0.001 * 0.00686286);
addWeights("results/merged/W1Jets_Pt800to1600-alpgen_ee.root", 0.001 * 3.40704e-05);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H130_mm.root", 0.001 * 0.00452859);
addWeights("results/merged/H160_mm.root", 0.001 * 0.00818843);
addWeights("results/merged/H190_mm.root", 0.001 * 0.00486346);
addWeights("results/merged/TTbar_mm.root", 0.001 * 0.122721);
addWeights("results/merged/SingleTop_sChannel_mm.root", 0.001 * 0.00437115);
addWeights("results/merged/SingleTop_tChannel_mm.root", 0.001 * 0.0437332);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 0.001 * 0.644887);
addWeights("results/merged/Wgamma_mm.root", 0.001 * 0.390098);
addWeights("results/merged/WjetsMadgraph_mm.root", 0.001 * 3.20548);
addWeights("results/merged/WW_mm.root", 0.001 * 0.0203685);
addWeights("results/merged/WZ_mm.root", 0.001 * 0.00549947);
addWeights("results/merged/ZjetsMadgraph_mm.root", 0.001 * 2.94516); 
addWeights("results/merged/ZZ_mm.root", 0.001 * 0.00231513);
addWeights("results/merged/W1Jets_Pt0to100-alpgen_mm.root", 0.001 * 5.08688);
addWeights("results/merged/W1Jets_Pt100to300-alpgen_mm.root", 0.001 * 0.490897);
addWeights("results/merged/W1Jets_Pt300to800-alpgen_mm.root", 0.001 * 0.00686286);
addWeights("results/merged/W1Jets_Pt800to1600-alpgen_mm.root", 0.001 * 3.40704e-05);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H130_em.root", 0.001 * 0.00452859);
addWeights("results/merged/H160_em.root", 0.001 * 0.00818843);
addWeights("results/merged/H190_em.root", 0.001 * 0.00486346);
addWeights("results/merged/TTbar_em.root", 0.001 * 0.122721);
addWeights("results/merged/SingleTop_sChannel_em.root", 0.001 * 0.00437115);
addWeights("results/merged/SingleTop_tChannel_em.root", 0.001 * 0.0437332);
addWeights("results/merged/SingleTop_tWChannel_em.root", 0.001 * 0.644887);
addWeights("results/merged/Wgamma_em.root", 0.001 * 0.390098);
addWeights("results/merged/WjetsMadgraph_em.root", 0.001 * 3.20548);
addWeights("results/merged/WW_em.root", 0.001 * 0.0203685);
addWeights("results/merged/WZ_em.root", 0.001 * 0.00549947);
addWeights("results/merged/ZjetsMadgraph_em.root", 0.001 * 2.94516); 
addWeights("results/merged/ZZ_em.root", 0.001 * 0.00231513);
addWeights("results/merged/W1Jets_Pt0to100-alpgen_em.root", 0.001 * 5.08688);
addWeights("results/merged/W1Jets_Pt100to300-alpgen_em.root", 0.001 * 0.490897);
addWeights("results/merged/W1Jets_Pt300to800-alpgen_em.root", 0.001 * 0.00686286);
addWeights("results/merged/W1Jets_Pt800to1600-alpgen_em.root", 0.001 * 3.40704e-05);

.q

EOF

echo "done weighting."
