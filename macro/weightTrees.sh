#! /bin/sh


lumi=$1
echo "===> THIS WEIGHTING IS FOR MH = 120 <====" 
echo "Adding weights for ee datasets for " $lumi " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wenu_ee.root", 0.00218859*$lumi, 0);
addWeights("results/merged/Wmunu_ee.root", 0.00153712*$lumi, 0);
addWeights("results/merged/Wtaunu_ee.root", 0.00162241*$lumi, 0);
addWeights("results/merged/Zee_Lo_ee.root", 0.000542764*$lumi, 0);
addWeights("results/merged/Zmm_Lo_ee.root", 0.00130758*$lumi, 0);
addWeights("results/merged/Ztautau_Lo_ee.root", 0*$lumi, 0);
addWeights("results/merged/Zee_Hi_ee.root", 0.000824758*$lumi, 0);
addWeights("results/merged/Zmm_Hi_ee.root", 0.000742182*$lumi, 0);
addWeights("results/merged/Ztautau_Hi_ee.root", 0.000809742*$lumi, 0);
addWeights("results/merged/SingleTop_sChannel_ee.root", 0*$lumi, 0);
addWeights("esults/merged/SingleTop_tChannel_ee.root", 5.44977e-05*$lumi, 0);
addWeights("results/merged/SingleTop_tWChannel_ee.root", 2.41229e-05*$lumi, 0);
addWeights("results/merged/TTbar_ee.root", 0.000155294*$lumi, 0);
addWeights("results/merged/H120_ee.root", 2.34177e-06*$lumi, 0);
addWeights("results/merged/WW_ee.root", 4.09406e-05*$lumi, 0);
addWeights("results/merged/ggWW_ee.root", 1.53892e-05*$lumi, 0);
addWeights("results/merged/WZ_ee.root", 0*$lumi, 0);
addWeights("results/merged/ZZ_ee.root", 0*$lumi, 0);
.q

EOF

echo "Adding weights for mm datasets for " $lumi " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wenu_mm.root", 0.00218859*$lumi, 1);
addWeights("results/merged/Wmunu_mm.root", 0.00153712*$lumi, 1);
addWeights("results/merged/Wtaunu_mm.root", 0.00162241*$lumi, 1);
addWeights("results/merged/Zee_Lo_mm.root", 0.000542764*$lumi, 1);
addWeights("results/merged/Zmm_Lo_mm.root", 0.00130758*$lumi, 1);
addWeights("results/merged/Ztautau_Lo_mm.root", 0*$lumi, 1);
addWeights("results/merged/Zee_Hi_mm.root", 0.000824758*$lumi, 1);
addWeights("results/merged/Zmm_Hi_mm.root", 0.000742182*$lumi, 1);
addWeights("results/merged/Ztautau_Hi_mm.root", 0.000809742*$lumi, 1);
addWeights("results/merged/SingleTop_sChannel_mm.root", 0*$lumi, 1);
addWeights("esults/merged/SingleTop_tChannel_mm.root", 5.44977e-05*$lumi, 1);
addWeights("results/merged/SingleTop_tWChannel_mm.root", 2.41229e-05*$lumi, 1);
addWeights("results/merged/TTbar_mm.root", 0.000155294*$lumi, 1);
addWeights("results/merged/H120_mm.root", 2.34177e-06*$lumi, 1);
addWeights("results/merged/WW_mm.root", 4.09406e-05*$lumi, 1);
addWeights("results/merged/ggWW_mm.root", 1.53892e-05*$lumi, 1);
addWeights("results/merged/WZ_mm.root", 0*$lumi, 1);
addWeights("results/merged/ZZ_mm.root", 0*$lumi, 1);
.q

EOF

echo "Adding weights for em datasets for " $lumi " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wenu_em.root", 0.00218859*$lumi, 2);
addWeights("results/merged/Wmunu_em.root", 0.00153712*$lumi, 2);
addWeights("results/merged/Wtaunu_em.root", 0.00162241*$lumi, 2);
addWeights("results/merged/Zee_Lo_em.root", 0.000542764*$lumi, 2);
addWeights("results/merged/Zmm_Lo_em.root", 0.00130758*$lumi, 2);
addWeights("results/merged/Ztautau_Lo_em.root", 0*$lumi, 2);
addWeights("results/merged/Zee_Hi_em.root", 0.000824758*$lumi, 2);
addWeights("results/merged/Zmm_Hi_em.root", 0.000742182*$lumi, 2);
addWeights("results/merged/Ztautau_Hi_em.root", 0.000809742*$lumi, 2);
addWeights("results/merged/SingleTop_sChannel_em.root", 0*$lumi, 2);
addWeights("esults/merged/SingleTop_tChannel_em.root", 5.44977e-05*$lumi, 2);
addWeights("results/merged/SingleTop_tWChannel_em.root", 2.41229e-05*$lumi, 2);
addWeights("results/merged/TTbar_em.root", 0.000155294*$lumi, 2);
addWeights("results/merged/H120_em.root", 2.34177e-06*$lumi, 2);
addWeights("results/merged/WW_em.root", 4.09406e-05*$lumi, 2);
addWeights("results/merged/ggWW_em.root", 1.53892e-05*$lumi, 2);
addWeights("results/merged/WZ_em.root", 0*$lumi, 2);
addWeights("results/merged/ZZ_em.root", 0*$lumi, 2);
.q

EOF

echo "Adding weights for me datasets for " $lumi " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/Wenu_me.root", 0.00218859*$lumi, 2);
addWeights("results/merged/Wmunu_me.root", 0.00153712*$lumi, 2);
addWeights("results/merged/Wtaunu_me.root", 0.00162241*$lumi, 2);
addWeights("results/merged/Zee_Lo_me.root", 0.000542764*$lumi, 2);
addWeights("results/merged/Zmm_Lo_me.root", 0.00130758*$lumi, 2);
addWeights("results/merged/Ztautau_Lo_me.root", 0*$lumi, 2);
addWeights("results/merged/Zee_Hi_me.root", 0.000824758*$lumi, 2);
addWeights("results/merged/Zmm_Hi_me.root", 0.000742182*$lumi, 2);
addWeights("results/merged/Ztautau_Hi_me.root", 0.000809742*$lumi, 2);
addWeights("results/merged/SingleTop_sChannel_me.root", 0*$lumi, 2);
addWeights("esults/merged/SingleTop_tChannel_me.root", 5.44977e-05*$lumi, 2);
addWeights("results/merged/SingleTop_tWChannel_me.root", 2.41229e-05*$lumi, 2);
addWeights("results/merged/TTbar_me.root", 0.000155294*$lumi, 2);
addWeights("results/merged/H120_me.root", 2.34177e-06*$lumi, 2);
addWeights("results/merged/WW_me.root", 4.09406e-05*$lumi, 2);
addWeights("results/merged/ggWW_me.root", 1.53892e-05*$lumi, 2);
addWeights("results/merged/WZ_me.root", 0*$lumi, 2);
addWeights("results/merged/ZZ_me.root", 0*$lumi, 2);
.q

EOF

echo "done weighting."
