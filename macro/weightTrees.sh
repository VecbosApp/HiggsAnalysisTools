#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

# PS. weights are compute for 1000 pb-1
mass=$1

echo "Adding weights for ee datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_ee.root", 2.69559e-06);
addWeights("results/merged/TTbar_ee.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_ee.root", 3.45355e-06);
addWeights("results/merged/SingleTop_tChannel_ee.root", 4.32393e-05);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 6.93873e-06);
addWeights("results/merged/Wjets_ee.root", 0.000899277);
addWeights("results/merged/Zjets_ee.root", 0.0103643); 
addWeights("results/merged/WW_ee.root", 4.09406e-05);
addWeights("results/merged/WZ_ee.root", 5.44947e-06);
addWeights("results/merged/ZZ_ee.root", 2.79175e-06);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_mm.root", 2.69559e-06);
addWeights("results/merged/TTbar_mm.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_mm.root", 3.45355e-06);
addWeights("results/merged/SingleTop_tChannel_mm.root", 4.32393e-05);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 6.93873e-06);
addWeights("results/merged/Wjets_mm.root", 0.000899277);
addWeights("results/merged/Zjets_mm.root", 0.0103643); 
addWeights("results/merged/WW_mm.root", 4.09406e-05);
addWeights("results/merged/WZ_mm.root", 5.44947e-06);
addWeights("results/merged/ZZ_mm.root", 2.79175e-06);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_em.root", 2.69559e-06);
addWeights("results/merged/TTbar_em.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_em.root", 3.45355e-06);
addWeights("results/merged/SingleTop_tChannel_em.root", 4.32393e-05);
addWeights("results/merged/SingleTop_tWChannel_em.root", 6.93873e-06);
addWeights("results/merged/Wjets_em.root", 0.000899277);
addWeights("results/merged/Zjets_em.root", 0.0103643); 
addWeights("results/merged/WW_em.root", 4.09406e-05);
addWeights("results/merged/WZ_em.root", 5.44947e-06);
addWeights("results/merged/ZZ_em.root", 2.79175e-06);

.q

EOF

echo "done weighting."
