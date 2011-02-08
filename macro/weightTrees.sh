#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTrees.sh

# PS. weights are compute for 1000 pb-1
lumi=$1

echo "Adding weights for ee datasets for " $lumi " pb-1..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_ee.root", 2.695590e-06*$lumi, 0);
addWeights("results/merged/WW_ee.root", 4.094063e-05*$lumi, 0);
addWeights("results/merged/ggWW_ee.root", 1.538923e-05*$lumi, 0);
addWeights("results/merged/WZ_ee.root", 5.449473e-06*$lumi, 0);
addWeights("results/merged/ZZ_ee.root", 2.791750e-06*$lumi, 0);
addWeights("results/merged/Wenu_ee.root", 4.737725e-03*$lumi, 0);
addWeights("results/merged/Wmunu_ee.root", 1.452839e-03*$lumi, 0);
addWeights("results/merged/Wtaunu_ee.root", 1.976142e-03*$lumi, 0);
addWeights("results/merged/Zee_Hi_ee.root", 8.334209e-04*$lumi, 0);
addWeights("results/merged/Zmm_Hi_ee.root", 8.773362e-04*$lumi, 0);
addWeights("results/merged/Ztautau_Hi_ee.root", 9.014079e-04*$lumi, 0);
addWeights("results/merged/Zee_Lo_ee.root", 4.540642e-01*$lumi, 0);
addWeights("results/merged/Zmm_Lo_ee.root", 1.795360e-02*$lumi, 0);
addWeights("results/merged/Ztautau_Lo_ee.root", 4.565363e-04*$lumi, 0);
addWeights("results/merged/TTbar_ee.root", 1.313869e-04*$lumi, 0);
addWeights("results/merged/SingleTop_sChannel_ee.root", 3.065486e-06*$lumi, 0);
addWeights("results/merged/SingleTop_tChannel_ee.root", 4.323927e-05*$lumi, 0);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 6.938728e-06*$lumi, 0);
.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_mm.root", 2.695590e-06*$lumi, 1);
addWeights("results/merged/WW_mm.root", 4.094063e-05*$lumi, 1);
addWeights("results/merged/ggWW_mm.root", 1.538923e-05*$lumi, 1);
addWeights("results/merged/WZ_mm.root", 5.449473e-06*$lumi, 1);
addWeights("results/merged/ZZ_mm.root", 2.791750e-06*$lumi, 1);
addWeights("results/merged/Wenu_mm.root", 4.737725e-03*$lumi, 1);
addWeights("results/merged/Wmunu_mm.root", 1.452839e-03*$lumi, 1);
addWeights("results/merged/Wtaunu_mm.root", 1.976142e-03*$lumi, 1);
addWeights("results/merged/Zee_Hi_mm.root", 8.334209e-04*$lumi, 1);
addWeights("results/merged/Zmm_Hi_mm.root", 8.773362e-04*$lumi, 1);
addWeights("results/merged/Ztautau_Hi_mm.root", 9.014079e-04*$lumi, 1);
addWeights("results/merged/Zee_Lo_mm.root", 4.540642e-01*$lumi, 1);
addWeights("results/merged/Zmm_Lo_mm.root", 1.795360e-02*$lumi, 1);
addWeights("results/merged/Ztautau_Lo_mm.root", 4.565363e-04*$lumi, 1);
addWeights("results/merged/TTbar_mm.root", 1.313869e-04*$lumi, 1);
addWeights("results/merged/SingleTop_sChannel_mm.root", 3.065486e-06*$lumi, 1);
addWeights("results/merged/SingleTop_tChannel_mm.root", 4.323927e-05*$lumi, 1);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 6.938728e-06*$lumi, 1);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results/merged/H160_em.root", 2.695590e-06*$lumi, 2);
addWeights("results/merged/WW_em.root", 4.094063e-05*$lumi, 2);
addWeights("results/merged/ggWW_em.root", 1.538923e-05*$lumi, 2);
addWeights("results/merged/WZ_em.root", 5.449473e-06*$lumi, 2);
addWeights("results/merged/ZZ_em.root", 2.791750e-06*$lumi, 2);
addWeights("results/merged/Wenu_em.root", 4.737725e-03*$lumi, 2);
addWeights("results/merged/Wmunu_em.root", 1.452839e-03*$lumi, 2);
addWeights("results/merged/Wtaunu_em.root", 1.976142e-03*$lumi, 2);
addWeights("results/merged/Zee_Hi_em.root", 8.334209e-04*$lumi, 2);
addWeights("results/merged/Zmm_Hi_em.root", 8.773362e-04*$lumi, 2);
addWeights("results/merged/Ztautau_Hi_em.root", 9.014079e-04*$lumi, 2);
addWeights("results/merged/Zee_Lo_em.root", 4.540642e-01*$lumi, 2);
addWeights("results/merged/Zmm_Lo_em.root", 1.795360e-02*$lumi, 2);
addWeights("results/merged/Ztautau_Lo_em.root", 4.565363e-04*$lumi, 2);
addWeights("results/merged/TTbar_em.root", 1.313869e-04*$lumi, 2);
addWeights("results/merged/SingleTop_sChannel_em.root", 3.065486e-06*$lumi, 2);
addWeights("results/merged/SingleTop_tChannel_em.root", 4.323927e-05*$lumi, 2);
addWeights("results/merged/SingleTop_tWChannel_em.root", 6.938728e-06*$lumi, 2);

.Q

EOF

echo "done weighting."
