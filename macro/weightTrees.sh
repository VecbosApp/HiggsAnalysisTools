#! /bin/sh


lumiEE=$1
lumiMM=$2
lumiEM=$3
echo "===> THIS WEIGHTING IS FOR MH = 160 <====" 
echo "Adding weights for em datasets for " $lumiEM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wjets_em.root", 0.0053428*$lumiEM, 2);
addWeights("results/merged/Zee_Lo_em.root", 0.00070545*$lumiEM, 2);
addWeights("results/merged/Zmm_Lo_em.root", 0.00169952*$lumiEM, 2);
addWeights("results/merged/Ztautau_Lo_em.root", 0*$lumiEM, 2);
addWeights("results/merged/Zee_Hi_em.root", 0.000824758*$lumiEM, 2);
addWeights("results/merged/Zmm_Hi_em.root", 0.000742182*$lumiEM, 2);
addWeights("results/merged/Ztautau_Hi_em.root", 0.000809742*$lumiEM, 2);
addWeights("results/merged/SingleTop_sChannel_em.root", 0*$lumiEM, 2);
addWeights("esults/merged/SingleTop_tChannel_em.root", 5.44977e-05*$lumiEM, 2);
addWeights("results/merged/SingleTop_tWChannel_em.root", 2.41229e-05*$lumiEM, 2);
addWeights("results/merged/TTbar_em.root", 0.000155294*$lumiEM, 2);
addWeights("results/merged/H160_em.root", 3.74881e-06*$lumiEM, 2);
addWeights("results/merged/Wgamma_em.root", 0.00331462*$lumiEM, 2);
addWeights("results/merged/WW_em.root", 4.09406e-05*$lumiEM, 2);
addWeights("results/merged/ggWW_em.root", 1.53892e-05*$lumiEM, 2);
addWeights("results/merged/WZ_em.root", 9.9907e-06*$lumiEM, 2);
addWeights("results/merged/ZZ_em.root", 3.55322e-06*$lumiEM, 2);
.q

EOF


echo "done weighting."
