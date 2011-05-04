#! /bin/sh


lumiEE=$1
lumiMM=$2
lumiEM=$3
echo "===> THIS WEIGHTING IS FOR MH = 140 <====" 
echo "Adding weights for ee datasets for " $lumiEE " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wjets_ee.root", 0.00534013*$lumiEE, 0);
addWeights("results/merged/Zjets_Hi_ee.root", 0.00297612*$lumiEE, 0);
addWeights("results/merged/Zjets_Lo_ee.root", 0.00175734*$lumiEE, 0);
addWeights("results/merged/SingleTop_sChannel_ee.root", 0*$lumiEE, 0);
addWeights("esults/merged/SingleTop_tChannel_ee.root", 4.20658e-05*$lumiEE, 0);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 2.10968e-05*$lumiEE, 0);
addWeights("results/merged/TTbar_ee.root", 9.11884e-05*$lumiEE, 0);
addWeights("results/merged/H140_ee.root", 2.6057e-06*$lumiEE, 0);
addWeights("results/merged/Wgamma_ee.root", 0.000152149*$lumiEE, 0);
addWeights("results/merged/WW_ee.root", 3.356e-05*$lumiEE, 0);
addWeights("results/merged/ggWW_ee.root", 1.36384e-06*$lumiEE, 0);
addWeights("results/merged/WZ_ee.root", 5.27635e-06*$lumiEE, 0);
addWeights("results/merged/ZZ_ee.root", 3.12428e-06*$lumiEE, 0);
.q

EOF

echo "Adding weights for mm datasets for " $lumiMM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wjets_mm.root", 0.00534013*$lumiMM, 1);
addWeights("results/merged/Zjets_Hi_mm.root", 0.00297612*$lumiMM, 1);
addWeights("results/merged/Zjets_Lo_mm.root", 0.00175734*$lumiMM, 1);
addWeights("results/merged/SingleTop_sChannel_mm.root", 0*$lumiMM, 1);
addWeights("esults/merged/SingleTop_tChannel_mm.root", 4.20658e-05*$lumiMM, 1);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 2.10968e-05*$lumiMM, 1);
addWeights("results/merged/TTbar_mm.root", 9.11884e-05*$lumiMM, 1);
addWeights("results/merged/H140_mm.root", 2.6057e-06*$lumiMM, 1);
addWeights("results/merged/Wgamma_mm.root", 0.000152149*$lumiMM, 1);
addWeights("results/merged/WW_mm.root", 3.356e-05*$lumiMM, 1);
addWeights("results/merged/ggWW_mm.root", 1.36384e-06*$lumiMM, 1);
addWeights("results/merged/WZ_mm.root", 5.27635e-06*$lumiMM, 1);
addWeights("results/merged/ZZ_mm.root", 3.12428e-06*$lumiMM, 1);
.q

EOF

echo "Adding weights for em datasets for " $lumiEM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wjets_em.root", 0.00534013*$lumiEM, 2);
addWeights("results/merged/Zjets_Hi_em.root", 0.00297612*$lumiEM, 2);
addWeights("results/merged/Zjets_Lo_em.root", 0.00175734*$lumiEM, 2);
addWeights("results/merged/SingleTop_sChannel_em.root", 0*$lumiEM, 2);
addWeights("esults/merged/SingleTop_tChannel_em.root", 4.20658e-05*$lumiEM, 2);
addWeights("results/merged/SingleTop_tWChannel_em.root", 2.10968e-05*$lumiEM, 2);
addWeights("results/merged/TTbar_em.root", 9.11884e-05*$lumiEM, 2);
addWeights("results/merged/H140_em.root", 2.6057e-06*$lumiEM, 2);
addWeights("results/merged/Wgamma_em.root", 0.000152149*$lumiEM, 2);
addWeights("results/merged/WW_em.root", 3.356e-05*$lumiEM, 2);
addWeights("results/merged/ggWW_em.root", 1.36384e-06*$lumiEM, 2);
addWeights("results/merged/WZ_em.root", 5.27635e-06*$lumiEM, 2);
addWeights("results/merged/ZZ_em.root", 3.12428e-06*$lumiEM, 2);
.q

EOF

echo "Adding weights for me datasets for " $lumiEM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wjets_me.root", 0.00534013*$lumiEM, 2);
addWeights("results/merged/Zjets_Hi_me.root", 0.00297612*$lumiEM, 2);
addWeights("results/merged/Zjets_Lo_me.root", 0.00175734*$lumiEM, 2);
addWeights("results/merged/SingleTop_sChannel_me.root", 0*$lumiEM, 2);
addWeights("esults/merged/SingleTop_tChannel_me.root", 4.20658e-05*$lumiEM, 2);
addWeights("results/merged/SingleTop_tWChannel_me.root", 2.10968e-05*$lumiEM, 2);
addWeights("results/merged/TTbar_me.root", 9.11884e-05*$lumiEM, 2);
addWeights("results/merged/H140_me.root", 2.6057e-06*$lumiEM, 2);
addWeights("results/merged/Wgamma_me.root", 0.000152149*$lumiEM, 2);
addWeights("results/merged/WW_me.root", 3.356e-05*$lumiEM, 2);
addWeights("results/merged/ggWW_me.root", 1.36384e-06*$lumiEM, 2);
addWeights("results/merged/WZ_me.root", 5.27635e-06*$lumiEM, 2);
addWeights("results/merged/ZZ_me.root", 3.12428e-06*$lumiEM, 2);
.q

EOF

echo "done weighting."
