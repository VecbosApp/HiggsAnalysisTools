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

addWeights("results/merged/H130_ee.root", 4.11709e-06);
addWeights("results/merged/TTbar_ee.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_ee.root", 2.54689e-06);
addWeights("results/merged/SingleTop_tChannel_ee.root", 2.25913e-05);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 6.93873e-06);
addWeights("results/merged/Zee_ee.root", 0.000813084); 
addWeights("results/merged/Zmm_ee.root", 0.000877336); 
addWeights("results/merged/Ztautau_ee.root", 0.000901408); 
addWeights("results/merged/Wenu_ee.root", 0.00276752);
addWeights("results/merged/Wmunu_ee.root", 0.00145284);
addWeights("results/merged/Wtaunu_ee.root", 0.00197614);
addWeights("results/merged/WW_ee.root", 4.09406e-05);
addWeights("results/merged/WZ_ee.root", 5.44947e-06);
addWeights("results/merged/ZZ_ee.root", 1.89505e-06);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H130_mm.root", 4.11709e-06);
addWeights("results/merged/TTbar_mm.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_mm.root", 2.54689e-06);
addWeights("results/merged/SingleTop_tChannel_mm.root", 2.25913e-05);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 6.93873e-06);
addWeights("results/merged/Zee_mm.root", 0.000813084); 
addWeights("results/merged/Zmm_mm.root", 0.000877336); 
addWeights("results/merged/Ztautau_mm.root", 0.000901408); 
addWeights("results/merged/Wenu_mm.root", 0.00276752);
addWeights("results/merged/Wmunu_mm.root", 0.00145284);
addWeights("results/merged/Wtaunu_mm.root", 0.00197614);
addWeights("results/merged/WW_mm.root", 4.09406e-05);
addWeights("results/merged/WZ_mm.root", 5.44947e-06);
addWeights("results/merged/ZZ_mm.root", 1.89505e-06);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H130_em.root", 4.11709e-06);
addWeights("results/merged/TTbar_em.root", 8.20589e-05);
addWeights("results/merged/SingleTop_sChannel_em.root", 2.54689e-06);
addWeights("results/merged/SingleTop_tChannel_em.root", 2.25913e-05);
addWeights("results/merged/SingleTop_tWChannel_em.root", 6.93873e-06);
addWeights("results/merged/Zee_em.root", 0.000813084); 
addWeights("results/merged/Zmm_em.root", 0.000877336); 
addWeights("results/merged/Ztautau_em.root", 0.000901408); 
addWeights("results/merged/Wenu_em.root", 0.00276752);
addWeights("results/merged/Wmunu_em.root", 0.00145284);
addWeights("results/merged/Wtaunu_em.root", 0.00197614);
addWeights("results/merged/WW_em.root", 4.09406e-05);
addWeights("results/merged/WZ_em.root", 5.44947e-06);
addWeights("results/merged/ZZ_em.root", 1.89505e-06);

.q

EOF

echo "done weighting."
