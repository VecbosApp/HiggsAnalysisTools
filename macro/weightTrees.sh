#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

# PS. weights are compute for 100 pb-1

echo "Adding weights for ee datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H130_ee.root", 0.00452859);
addWeights("results/merged/H160_ee.root", 0.00818843);
addWeights("results/merged/H190_ee.root", 0.00486346);
addWeights("results/merged/TTbar_ee.root", 0.128565);
addWeights("results/merged/SingleTop_sChannel_ee.root", 0.00477608);
addWeights("results/merged/SingleTop_tChannel_ee.root", 0.04265);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 0.644887);
addWeights("results/merged/Wgamma_ee.root", 0.390098);
addWeights("results/merged/WjetsMadgraph_ee.root", 2.97267);
addWeights("results/merged/WW_ee.root", 0.0203685);
addWeights("results/merged/WZ_ee.root", 0.00549947);
addWeights("results/merged/ZjetsMadgraph_ee.root", 4.65717); 
addWeights("results/merged/ZZ_ee.root", 0.00231513);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H130_mm.root", 0.00452859);
addWeights("results/merged/H160_mm.root", 0.00818843);
addWeights("results/merged/H190_mm.root", 0.00486346);
addWeights("results/merged/TTbar_mm.root", 0.128565);
addWeights("results/merged/SingleTop_sChannel_mm.root", 0.00477608);
addWeights("results/merged/SingleTop_tChannel_mm.root", 0.04265);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 0.644887);
addWeights("results/merged/Wgamma_mm.root", 0.390098);
addWeights("results/merged/WjetsMadgraph_mm.root", 2.97267);
addWeights("results/merged/WW_mm.root", 0.0203685);
addWeights("results/merged/WZ_mm.root", 0.00549947);
addWeights("results/merged/ZjetsMadgraph_mm.root", 4.65717); 
addWeights("results/merged/ZZ_mm.root", 0.00231513);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

gSystem->Load("addWeightsToTree_cc.so");

addWeights("results/merged/H130_em.root", 0.00452859);
addWeights("results/merged/H160_em.root", 0.00818843);
addWeights("results/merged/H190_em.root", 0.00486346);
addWeights("results/merged/TTbar_em.root", 0.128565);
addWeights("results/merged/SingleTop_sChannel_em.root", 0.00477608);
addWeights("results/merged/SingleTop_tChannel_em.root", 0.04265);
addWeights("results/merged/SingleTop_tWChannel_em.root", 0.644887);
addWeights("results/merged/Wgamma_em.root", 0.390098);
addWeights("results/merged/WjetsMadgraph_em.root", 2.97267);
addWeights("results/merged/WW_em.root", 0.0203685);
addWeights("results/merged/WZ_em.root", 0.00549947);
addWeights("results/merged/ZjetsMadgraph_em.root", 4.65717); 
addWeights("results/merged/ZZ_em.root", 0.00231513);

.q

EOF

echo "done weighting."
