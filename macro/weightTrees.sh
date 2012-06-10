#! /bin/sh


mkdir -p results/merged_skim
lumiEE=$1
lumiMM=$2
lumiEM=$3
echo "Adding weights for ee datasets for " $lumiEE " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/ggH1102LAndTau2Nu_ee.root", 4.26205e-07*$lumiEE, 9110 ,1, 0);
addWeights("results/merged/qqH1102LAndTau2Nu_ee.root", 3.02383e-08*$lumiEE, 8110 ,1, 0);
addWeights("results/merged/ggH1152LAndTau2Nu_ee.root", 8.35876e-07*$lumiEE, 9115 ,1, 0);
addWeights("results/merged/qqH1152LAndTau2Nu_ee.root", 0*$lumiEE, 8115 ,1, 0);
addWeights("results/merged/ggH1202LAndTau2Nu_ee.root", 1.06704e-06*$lumiEE, 9120 ,1, 0);
addWeights("results/merged/qqH1202LAndTau2Nu_ee.root", 0*$lumiEE, 8120 ,1, 0);
addWeights("results/merged/ggH1252LAndTau2Nu_ee.root", 1.4755e-06*$lumiEE, 9125 ,1, 0);
addWeights("results/merged/qqH1252LAndTau2Nu_ee.root", 1.41577e-07*$lumiEE, 8125 ,1, 0);
addWeights("results/merged/ggH1302LAndTau2Nu_ee.root", 1.92872e-06*$lumiEE, 9130 ,1, 0);
addWeights("results/merged/qqH1302LAndTau2Nu_ee.root", 1.59205e-07*$lumiEE, 8130 ,1, 0);
addWeights("results/merged/ggH1352LAndTau2Nu_ee.root", 2.36788e-06*$lumiEE, 9135 ,1, 0);
addWeights("results/merged/qqH1352LAndTau2Nu_ee.root", 0*$lumiEE, 8135 ,1, 0);
addWeights("results/merged/ggH1402LAndTau2Nu_ee.root", 2.75125e-06*$lumiEE, 9140 ,1, 0);
addWeights("results/merged/qqH1402LAndTau2Nu_ee.root", 2.40548e-07*$lumiEE, 8140 ,1, 0);
addWeights("results/merged/ggH1452LAndTau2Nu_ee.root", 3.07371e-06*$lumiEE, 9145 ,1, 0);
addWeights("results/merged/qqH1452LAndTau2Nu_ee.root", 2.75451e-07*$lumiEE, 8145 ,1, 0);
addWeights("results/merged/ggH1502LAndTau2Nu_ee.root", 3.33417e-06*$lumiEE, 9150 ,1, 0);
addWeights("results/merged/qqH1502LAndTau2Nu_ee.root", 3.05939e-07*$lumiEE, 8150 ,1, 0);
addWeights("results/merged/ggH1552LAndTau2Nu_ee.root", 3.55826e-06*$lumiEE, 9155 ,1, 0);
addWeights("results/merged/qqH1552LAndTau2Nu_ee.root", 3.33338e-07*$lumiEE, 8155 ,1, 0);
addWeights("results/merged/ggH1602LAndTau2Nu_ee.root", 3.79728e-06*$lumiEE, 9160 ,1, 0);
addWeights("results/merged/qqH1602LAndTau2Nu_ee.root", 4.37525e-07*$lumiEE, 8160 ,1, 0);
addWeights("results/merged/ggH1702LAndTau2Nu_ee.root", 3.41403e-06*$lumiEE, 9170 ,1, 0);
addWeights("results/merged/qqH1702LAndTau2Nu_ee.root", 3.61367e-07*$lumiEE, 8170 ,1, 0);
addWeights("results/merged/ggH1802LAndTau2Nu_ee.root", 2.89425e-06*$lumiEE, 9180 ,1, 0);
addWeights("results/merged/qqH1802LAndTau2Nu_ee.root", 0*$lumiEE, 8180 ,1, 0);
addWeights("results/merged/ggH1902LAndTau2Nu_ee.root", 2.16421e-06*$lumiEE, 9190 ,1, 0);
addWeights("results/merged/qqH1902LAndTau2Nu_ee.root", 2.51662e-07*$lumiEE, 8190 ,1, 0);
addWeights("results/merged/ggH2002LAndTau2Nu_ee.root", 3.69593e-06*$lumiEE, 9200 ,1, 0);
addWeights("results/merged/qqH2002LAndTau2Nu_ee.root", 2.19108e-07*$lumiEE, 8200 ,1, 0);
addWeights("results/merged/ggH2502LAndTau2Nu_ee.root", 1.178e-06*$lumiEE, 9250 ,1, 0);
addWeights("results/merged/qqH2502LAndTau2Nu_ee.root", 0*$lumiEE, 8250 ,1, 0);
addWeights("results/merged/ggH3002LAndTau2Nu_ee.root", 8.75005e-07*$lumiEE, 9300 ,1, 0);
addWeights("results/merged/qqH3002LAndTau2Nu_ee.root", 0*$lumiEE, 8300 ,1, 0);
addWeights("results/merged/ggH3502LAndTau2Nu_ee.root", 2.41703e-06*$lumiEE, 9350 ,1, 0);
addWeights("results/merged/qqH3502LAndTau2Nu_ee.root", 1.51389e-07*$lumiEE, 8350 ,1, 0);
addWeights("results/merged/ggH4002LAndTau2Nu_ee.root", 5.99725e-07*$lumiEE, 9400 ,1, 0);
addWeights("results/merged/qqH4002LAndTau2Nu_ee.root", 5.19202e-08*$lumiEE, 8400 ,1, 0);
addWeights("results/merged/ggH4502LAndTau2Nu_ee.root", 3.86613e-07*$lumiEE, 9450 ,1, 0);
addWeights("results/merged/qqH4502LAndTau2Nu_ee.root", 4.65678e-08*$lumiEE, 8450 ,1, 0);
addWeights("results/merged/ggH5002LAndTau2Nu_ee.root", 2.47147e-07*$lumiEE, 9500 ,1, 0);
addWeights("results/merged/qqH5002LAndTau2Nu_ee.root", 0*$lumiEE, 8500 ,1, 0);
addWeights("results/merged/ggH5502LAndTau2Nu_ee.root", 1.58151e-07*$lumiEE, 9550 ,1, 0);
addWeights("results/merged/qqH5502LAndTau2Nu_ee.root", 0*$lumiEE, 8550 ,1, 0);
addWeights("results/merged/ggH6002LAndTau2Nu_ee.root", 1.03169e-07*$lumiEE, 9600 ,1, 0);
addWeights("results/merged/qqH6002LAndTau2Nu_ee.root", 1.90246e-08*$lumiEE, 8600 ,1, 0);
addWeights("results/merged/ggH7002LAndTau2Nu_ee.root", 4.62052e-08*$lumiEE, 9700 ,1, 0);
addWeights("results/merged/qqH7002LAndTau2Nu_ee.root", 1.54801e-08*$lumiEE, 8700 ,1, 0);
addWeights("results/merged/ggH8002LAndTau2Nu_ee.root", 2.27653e-08*$lumiEE, 9800 ,1, 0);
addWeights("results/merged/qqH8002LAndTau2Nu_ee.root", 0*$lumiEE, 8800 ,1, 0);
addWeights("results/merged/ggH9002LAndTau2Nu_ee.root", 1.21167e-08*$lumiEE, 9900 ,1, 0);
addWeights("results/merged/qqH9002LAndTau2Nu_ee.root", 6.82463e-09*$lumiEE, 8900 ,1, 0);
addWeights("results/merged/ggH10002LAndTau2Nu_ee.root", 6.88254e-09*$lumiEE, 10000 ,1, 0);
addWeights("results/merged/qqH10002LAndTau2Nu_ee.root", 0*$lumiEE, 9000 ,1, 0);
addWeights("results/merged/Wjets_ee.root", 0.00233401*$lumiEE, 80 ,1, 0);
addWeights("results/merged/Zll_Lo_ee.root", 0.000123242*$lumiEE, 36 ,1, 0);
addWeights("results/merged/Zll_Hi_ee.root", 0.000117914*$lumiEE, 37 ,1, 0);
addWeights("results/merged/SingleT_tWChannel_ee.root", 2.81079e-05*$lumiEE, 19 ,1, 0);
addWeights("results/merged/SingleTbar_tWChannel_ee.root", 2.26509e-05*$lumiEE, 20 ,1, 0);
addWeights("results/merged/TTbar_ee.root", 3.36812e-05*$lumiEE, 10 ,1, 0);
addWeights("results/merged/WW_ee.root", 3.00652e-06*$lumiEE, 0 ,1, 0);
addWeights("results/merged/ggWW_ee.root", 1.66249e-06*$lumiEE, 1 ,1, 0);
addWeights("results/merged/WZ_ee.root", 2.48505e-07*$lumiEE, 74 ,1, 0);
addWeights("results/merged/ZZ_ee.root", 7.32849e-07*$lumiEE, 71 ,1, 0);
addWeights("results/merged/WW_pythia_ee.root", 1.1645e-05*$lumiEE, 6 ,1, 0);
addWeights("results/merged/WGToLNuG_ee.root", 0.000113575*$lumiEE, 85 ,1, 0);
addWeights("results/merged/WGstarToLNu2E_ee.root", 1.86651e-05*$lumiEE, 82 ,1, 0);
addWeights("results/merged/WGstarToLNu2Mu_ee.root", 6.38057e-06*$lumiEE, 83 ,1, 0);
addWeights("results/merged/WGstarToLNu2Tau_ee.root", 6.7204e-06*$lumiEE, 84 ,1, 0);
.q

EOF

echo "Adding weights for mm datasets for " $lumiMM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/ggH1102LAndTau2Nu_mm.root", 4.26205e-07*$lumiMM, 9110 ,0, 0);
addWeights("results/merged/qqH1102LAndTau2Nu_mm.root", 3.02383e-08*$lumiMM, 8110 ,0, 0);
addWeights("results/merged/ggH1152LAndTau2Nu_mm.root", 8.35876e-07*$lumiMM, 9115 ,0, 0);
addWeights("results/merged/qqH1152LAndTau2Nu_mm.root", 0*$lumiMM, 8115 ,0, 0);
addWeights("results/merged/ggH1202LAndTau2Nu_mm.root", 1.06704e-06*$lumiMM, 9120 ,0, 0);
addWeights("results/merged/qqH1202LAndTau2Nu_mm.root", 0*$lumiMM, 8120 ,0, 0);
addWeights("results/merged/ggH1252LAndTau2Nu_mm.root", 1.4755e-06*$lumiMM, 9125 ,0, 0);
addWeights("results/merged/qqH1252LAndTau2Nu_mm.root", 1.41577e-07*$lumiMM, 8125 ,0, 0);
addWeights("results/merged/ggH1302LAndTau2Nu_mm.root", 1.92872e-06*$lumiMM, 9130 ,0, 0);
addWeights("results/merged/qqH1302LAndTau2Nu_mm.root", 1.59205e-07*$lumiMM, 8130 ,0, 0);
addWeights("results/merged/ggH1352LAndTau2Nu_mm.root", 2.36788e-06*$lumiMM, 9135 ,0, 0);
addWeights("results/merged/qqH1352LAndTau2Nu_mm.root", 0*$lumiMM, 8135 ,0, 0);
addWeights("results/merged/ggH1402LAndTau2Nu_mm.root", 2.75125e-06*$lumiMM, 9140 ,0, 0);
addWeights("results/merged/qqH1402LAndTau2Nu_mm.root", 2.40548e-07*$lumiMM, 8140 ,0, 0);
addWeights("results/merged/ggH1452LAndTau2Nu_mm.root", 3.07371e-06*$lumiMM, 9145 ,0, 0);
addWeights("results/merged/qqH1452LAndTau2Nu_mm.root", 2.75451e-07*$lumiMM, 8145 ,0, 0);
addWeights("results/merged/ggH1502LAndTau2Nu_mm.root", 3.33417e-06*$lumiMM, 9150 ,0, 0);
addWeights("results/merged/qqH1502LAndTau2Nu_mm.root", 3.05939e-07*$lumiMM, 8150 ,0, 0);
addWeights("results/merged/ggH1552LAndTau2Nu_mm.root", 3.55826e-06*$lumiMM, 9155 ,0, 0);
addWeights("results/merged/qqH1552LAndTau2Nu_mm.root", 3.33338e-07*$lumiMM, 8155 ,0, 0);
addWeights("results/merged/ggH1602LAndTau2Nu_mm.root", 3.79728e-06*$lumiMM, 9160 ,0, 0);
addWeights("results/merged/qqH1602LAndTau2Nu_mm.root", 4.37525e-07*$lumiMM, 8160 ,0, 0);
addWeights("results/merged/ggH1702LAndTau2Nu_mm.root", 3.41403e-06*$lumiMM, 9170 ,0, 0);
addWeights("results/merged/qqH1702LAndTau2Nu_mm.root", 3.61367e-07*$lumiMM, 8170 ,0, 0);
addWeights("results/merged/ggH1802LAndTau2Nu_mm.root", 2.89425e-06*$lumiMM, 9180 ,0, 0);
addWeights("results/merged/qqH1802LAndTau2Nu_mm.root", 0*$lumiMM, 8180 ,0, 0);
addWeights("results/merged/ggH1902LAndTau2Nu_mm.root", 2.16421e-06*$lumiMM, 9190 ,0, 0);
addWeights("results/merged/qqH1902LAndTau2Nu_mm.root", 2.51662e-07*$lumiMM, 8190 ,0, 0);
addWeights("results/merged/ggH2002LAndTau2Nu_mm.root", 3.69593e-06*$lumiMM, 9200 ,0, 0);
addWeights("results/merged/qqH2002LAndTau2Nu_mm.root", 2.19108e-07*$lumiMM, 8200 ,0, 0);
addWeights("results/merged/ggH2502LAndTau2Nu_mm.root", 1.178e-06*$lumiMM, 9250 ,0, 0);
addWeights("results/merged/qqH2502LAndTau2Nu_mm.root", 0*$lumiMM, 8250 ,0, 0);
addWeights("results/merged/ggH3002LAndTau2Nu_mm.root", 8.75005e-07*$lumiMM, 9300 ,0, 0);
addWeights("results/merged/qqH3002LAndTau2Nu_mm.root", 0*$lumiMM, 8300 ,0, 0);
addWeights("results/merged/ggH3502LAndTau2Nu_mm.root", 2.41703e-06*$lumiMM, 9350 ,0, 0);
addWeights("results/merged/qqH3502LAndTau2Nu_mm.root", 1.51389e-07*$lumiMM, 8350 ,0, 0);
addWeights("results/merged/ggH4002LAndTau2Nu_mm.root", 5.99725e-07*$lumiMM, 9400 ,0, 0);
addWeights("results/merged/qqH4002LAndTau2Nu_mm.root", 5.19202e-08*$lumiMM, 8400 ,0, 0);
addWeights("results/merged/ggH4502LAndTau2Nu_mm.root", 3.86613e-07*$lumiMM, 9450 ,0, 0);
addWeights("results/merged/qqH4502LAndTau2Nu_mm.root", 4.65678e-08*$lumiMM, 8450 ,0, 0);
addWeights("results/merged/ggH5002LAndTau2Nu_mm.root", 2.47147e-07*$lumiMM, 9500 ,0, 0);
addWeights("results/merged/qqH5002LAndTau2Nu_mm.root", 0*$lumiMM, 8500 ,0, 0);
addWeights("results/merged/ggH5502LAndTau2Nu_mm.root", 1.58151e-07*$lumiMM, 9550 ,0, 0);
addWeights("results/merged/qqH5502LAndTau2Nu_mm.root", 0*$lumiMM, 8550 ,0, 0);
addWeights("results/merged/ggH6002LAndTau2Nu_mm.root", 1.03169e-07*$lumiMM, 9600 ,0, 0);
addWeights("results/merged/qqH6002LAndTau2Nu_mm.root", 1.90246e-08*$lumiMM, 8600 ,0, 0);
addWeights("results/merged/ggH7002LAndTau2Nu_mm.root", 4.62052e-08*$lumiMM, 9700 ,0, 0);
addWeights("results/merged/qqH7002LAndTau2Nu_mm.root", 1.54801e-08*$lumiMM, 8700 ,0, 0);
addWeights("results/merged/ggH8002LAndTau2Nu_mm.root", 2.27653e-08*$lumiMM, 9800 ,0, 0);
addWeights("results/merged/qqH8002LAndTau2Nu_mm.root", 0*$lumiMM, 8800 ,0, 0);
addWeights("results/merged/ggH9002LAndTau2Nu_mm.root", 1.21167e-08*$lumiMM, 9900 ,0, 0);
addWeights("results/merged/qqH9002LAndTau2Nu_mm.root", 6.82463e-09*$lumiMM, 8900 ,0, 0);
addWeights("results/merged/ggH10002LAndTau2Nu_mm.root", 6.88254e-09*$lumiMM, 10000 ,0, 0);
addWeights("results/merged/qqH10002LAndTau2Nu_mm.root", 0*$lumiMM, 9000 ,0, 0);
addWeights("results/merged/Wjets_mm.root", 0.00233401*$lumiMM, 80 ,0, 0);
addWeights("results/merged/Zll_Lo_mm.root", 0.000123242*$lumiMM, 36 ,0, 0);
addWeights("results/merged/Zll_Hi_mm.root", 0.000117914*$lumiMM, 37 ,0, 0);
addWeights("results/merged/SingleT_tWChannel_mm.root", 2.81079e-05*$lumiMM, 19 ,0, 0);
addWeights("results/merged/SingleTbar_tWChannel_mm.root", 2.26509e-05*$lumiMM, 20 ,0, 0);
addWeights("results/merged/TTbar_mm.root", 3.36812e-05*$lumiMM, 10 ,0, 0);
addWeights("results/merged/WW_mm.root", 3.00652e-06*$lumiMM, 0 ,0, 0);
addWeights("results/merged/ggWW_mm.root", 1.66249e-06*$lumiMM, 1 ,0, 0);
addWeights("results/merged/WZ_mm.root", 2.48505e-07*$lumiMM, 74 ,0, 0);
addWeights("results/merged/ZZ_mm.root", 7.32849e-07*$lumiMM, 71 ,0, 0);
addWeights("results/merged/WW_pythia_mm.root", 1.1645e-05*$lumiMM, 6 ,0, 0);
addWeights("results/merged/WGToLNuG_mm.root", 0.000113575*$lumiMM, 85 ,0, 0);
addWeights("results/merged/WGstarToLNu2E_mm.root", 1.86651e-05*$lumiMM, 82 ,0, 0);
addWeights("results/merged/WGstarToLNu2Mu_mm.root", 6.38057e-06*$lumiMM, 83 ,0, 0);
addWeights("results/merged/WGstarToLNu2Tau_mm.root", 6.7204e-06*$lumiMM, 84 ,0, 0);
.q

EOF

echo "Adding weights for em datasets for " $lumiEM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/ggH1102LAndTau2Nu_em.root", 4.26205e-07*$lumiEM, 9110 ,2, 0);
addWeights("results/merged/qqH1102LAndTau2Nu_em.root", 3.02383e-08*$lumiEM, 8110 ,2, 0);
addWeights("results/merged/ggH1152LAndTau2Nu_em.root", 8.35876e-07*$lumiEM, 9115 ,2, 0);
addWeights("results/merged/qqH1152LAndTau2Nu_em.root", 0*$lumiEM, 8115 ,2, 0);
addWeights("results/merged/ggH1202LAndTau2Nu_em.root", 1.06704e-06*$lumiEM, 9120 ,2, 0);
addWeights("results/merged/qqH1202LAndTau2Nu_em.root", 0*$lumiEM, 8120 ,2, 0);
addWeights("results/merged/ggH1252LAndTau2Nu_em.root", 1.4755e-06*$lumiEM, 9125 ,2, 0);
addWeights("results/merged/qqH1252LAndTau2Nu_em.root", 1.41577e-07*$lumiEM, 8125 ,2, 0);
addWeights("results/merged/ggH1302LAndTau2Nu_em.root", 1.92872e-06*$lumiEM, 9130 ,2, 0);
addWeights("results/merged/qqH1302LAndTau2Nu_em.root", 1.59205e-07*$lumiEM, 8130 ,2, 0);
addWeights("results/merged/ggH1352LAndTau2Nu_em.root", 2.36788e-06*$lumiEM, 9135 ,2, 0);
addWeights("results/merged/qqH1352LAndTau2Nu_em.root", 0*$lumiEM, 8135 ,2, 0);
addWeights("results/merged/ggH1402LAndTau2Nu_em.root", 2.75125e-06*$lumiEM, 9140 ,2, 0);
addWeights("results/merged/qqH1402LAndTau2Nu_em.root", 2.40548e-07*$lumiEM, 8140 ,2, 0);
addWeights("results/merged/ggH1452LAndTau2Nu_em.root", 3.07371e-06*$lumiEM, 9145 ,2, 0);
addWeights("results/merged/qqH1452LAndTau2Nu_em.root", 2.75451e-07*$lumiEM, 8145 ,2, 0);
addWeights("results/merged/ggH1502LAndTau2Nu_em.root", 3.33417e-06*$lumiEM, 9150 ,2, 0);
addWeights("results/merged/qqH1502LAndTau2Nu_em.root", 3.05939e-07*$lumiEM, 8150 ,2, 0);
addWeights("results/merged/ggH1552LAndTau2Nu_em.root", 3.55826e-06*$lumiEM, 9155 ,2, 0);
addWeights("results/merged/qqH1552LAndTau2Nu_em.root", 3.33338e-07*$lumiEM, 8155 ,2, 0);
addWeights("results/merged/ggH1602LAndTau2Nu_em.root", 3.79728e-06*$lumiEM, 9160 ,2, 0);
addWeights("results/merged/qqH1602LAndTau2Nu_em.root", 4.37525e-07*$lumiEM, 8160 ,2, 0);
addWeights("results/merged/ggH1702LAndTau2Nu_em.root", 3.41403e-06*$lumiEM, 9170 ,2, 0);
addWeights("results/merged/qqH1702LAndTau2Nu_em.root", 3.61367e-07*$lumiEM, 8170 ,2, 0);
addWeights("results/merged/ggH1802LAndTau2Nu_em.root", 2.89425e-06*$lumiEM, 9180 ,2, 0);
addWeights("results/merged/qqH1802LAndTau2Nu_em.root", 0*$lumiEM, 8180 ,2, 0);
addWeights("results/merged/ggH1902LAndTau2Nu_em.root", 2.16421e-06*$lumiEM, 9190 ,2, 0);
addWeights("results/merged/qqH1902LAndTau2Nu_em.root", 2.51662e-07*$lumiEM, 8190 ,2, 0);
addWeights("results/merged/ggH2002LAndTau2Nu_em.root", 3.69593e-06*$lumiEM, 9200 ,2, 0);
addWeights("results/merged/qqH2002LAndTau2Nu_em.root", 2.19108e-07*$lumiEM, 8200 ,2, 0);
addWeights("results/merged/ggH2502LAndTau2Nu_em.root", 1.178e-06*$lumiEM, 9250 ,2, 0);
addWeights("results/merged/qqH2502LAndTau2Nu_em.root", 0*$lumiEM, 8250 ,2, 0);
addWeights("results/merged/ggH3002LAndTau2Nu_em.root", 8.75005e-07*$lumiEM, 9300 ,2, 0);
addWeights("results/merged/qqH3002LAndTau2Nu_em.root", 0*$lumiEM, 8300 ,2, 0);
addWeights("results/merged/ggH3502LAndTau2Nu_em.root", 2.41703e-06*$lumiEM, 9350 ,2, 0);
addWeights("results/merged/qqH3502LAndTau2Nu_em.root", 1.51389e-07*$lumiEM, 8350 ,2, 0);
addWeights("results/merged/ggH4002LAndTau2Nu_em.root", 5.99725e-07*$lumiEM, 9400 ,2, 0);
addWeights("results/merged/qqH4002LAndTau2Nu_em.root", 5.19202e-08*$lumiEM, 8400 ,2, 0);
addWeights("results/merged/ggH4502LAndTau2Nu_em.root", 3.86613e-07*$lumiEM, 9450 ,2, 0);
addWeights("results/merged/qqH4502LAndTau2Nu_em.root", 4.65678e-08*$lumiEM, 8450 ,2, 0);
addWeights("results/merged/ggH5002LAndTau2Nu_em.root", 2.47147e-07*$lumiEM, 9500 ,2, 0);
addWeights("results/merged/qqH5002LAndTau2Nu_em.root", 0*$lumiEM, 8500 ,2, 0);
addWeights("results/merged/ggH5502LAndTau2Nu_em.root", 1.58151e-07*$lumiEM, 9550 ,2, 0);
addWeights("results/merged/qqH5502LAndTau2Nu_em.root", 0*$lumiEM, 8550 ,2, 0);
addWeights("results/merged/ggH6002LAndTau2Nu_em.root", 1.03169e-07*$lumiEM, 9600 ,2, 0);
addWeights("results/merged/qqH6002LAndTau2Nu_em.root", 1.90246e-08*$lumiEM, 8600 ,2, 0);
addWeights("results/merged/ggH7002LAndTau2Nu_em.root", 4.62052e-08*$lumiEM, 9700 ,2, 0);
addWeights("results/merged/qqH7002LAndTau2Nu_em.root", 1.54801e-08*$lumiEM, 8700 ,2, 0);
addWeights("results/merged/ggH8002LAndTau2Nu_em.root", 2.27653e-08*$lumiEM, 9800 ,2, 0);
addWeights("results/merged/qqH8002LAndTau2Nu_em.root", 0*$lumiEM, 8800 ,2, 0);
addWeights("results/merged/ggH9002LAndTau2Nu_em.root", 1.21167e-08*$lumiEM, 9900 ,2, 0);
addWeights("results/merged/qqH9002LAndTau2Nu_em.root", 6.82463e-09*$lumiEM, 8900 ,2, 0);
addWeights("results/merged/ggH10002LAndTau2Nu_em.root", 6.88254e-09*$lumiEM, 10000 ,2, 0);
addWeights("results/merged/qqH10002LAndTau2Nu_em.root", 0*$lumiEM, 9000 ,2, 0);
addWeights("results/merged/Wjets_em.root", 0.00233401*$lumiEM, 80 ,2, 0);
addWeights("results/merged/Zll_Lo_em.root", 0.000123242*$lumiEM, 36 ,2, 0);
addWeights("results/merged/Zll_Hi_em.root", 0.000117914*$lumiEM, 37 ,2, 0);
addWeights("results/merged/SingleT_tWChannel_em.root", 2.81079e-05*$lumiEM, 19 ,2, 0);
addWeights("results/merged/SingleTbar_tWChannel_em.root", 2.26509e-05*$lumiEM, 20 ,2, 0);
addWeights("results/merged/TTbar_em.root", 3.36812e-05*$lumiEM, 10 ,2, 0);
addWeights("results/merged/WW_em.root", 3.00652e-06*$lumiEM, 0 ,2, 0);
addWeights("results/merged/ggWW_em.root", 1.66249e-06*$lumiEM, 1 ,2, 0);
addWeights("results/merged/WZ_em.root", 2.48505e-07*$lumiEM, 74 ,2, 0);
addWeights("results/merged/ZZ_em.root", 7.32849e-07*$lumiEM, 71 ,2, 0);
addWeights("results/merged/WW_pythia_em.root", 1.1645e-05*$lumiEM, 6 ,2, 0);
addWeights("results/merged/WGToLNuG_em.root", 0.000113575*$lumiEM, 85 ,2, 0);
addWeights("results/merged/WGstarToLNu2E_em.root", 1.86651e-05*$lumiEM, 82 ,2, 0);
addWeights("results/merged/WGstarToLNu2Mu_em.root", 6.38057e-06*$lumiEM, 83 ,2, 0);
addWeights("results/merged/WGstarToLNu2Tau_em.root", 6.7204e-06*$lumiEM, 84 ,2, 0);
.q

EOF

echo "Adding weights for me datasets for " $lumiEM " pb-1..."
root -l -b <<EOF
.L addWeightsToTree.cc+
addWeights("results/merged/ggH1102LAndTau2Nu_me.root", 4.26205e-07*$lumiEM, 9110 ,3, 0);
addWeights("results/merged/qqH1102LAndTau2Nu_me.root", 3.02383e-08*$lumiEM, 8110 ,3, 0);
addWeights("results/merged/ggH1152LAndTau2Nu_me.root", 8.35876e-07*$lumiEM, 9115 ,3, 0);
addWeights("results/merged/qqH1152LAndTau2Nu_me.root", 0*$lumiEM, 8115 ,3, 0);
addWeights("results/merged/ggH1202LAndTau2Nu_me.root", 1.06704e-06*$lumiEM, 9120 ,3, 0);
addWeights("results/merged/qqH1202LAndTau2Nu_me.root", 0*$lumiEM, 8120 ,3, 0);
addWeights("results/merged/ggH1252LAndTau2Nu_me.root", 1.4755e-06*$lumiEM, 9125 ,3, 0);
addWeights("results/merged/qqH1252LAndTau2Nu_me.root", 1.41577e-07*$lumiEM, 8125 ,3, 0);
addWeights("results/merged/ggH1302LAndTau2Nu_me.root", 1.92872e-06*$lumiEM, 9130 ,3, 0);
addWeights("results/merged/qqH1302LAndTau2Nu_me.root", 1.59205e-07*$lumiEM, 8130 ,3, 0);
addWeights("results/merged/ggH1352LAndTau2Nu_me.root", 2.36788e-06*$lumiEM, 9135 ,3, 0);
addWeights("results/merged/qqH1352LAndTau2Nu_me.root", 0*$lumiEM, 8135 ,3, 0);
addWeights("results/merged/ggH1402LAndTau2Nu_me.root", 2.75125e-06*$lumiEM, 9140 ,3, 0);
addWeights("results/merged/qqH1402LAndTau2Nu_me.root", 2.40548e-07*$lumiEM, 8140 ,3, 0);
addWeights("results/merged/ggH1452LAndTau2Nu_me.root", 3.07371e-06*$lumiEM, 9145 ,3, 0);
addWeights("results/merged/qqH1452LAndTau2Nu_me.root", 2.75451e-07*$lumiEM, 8145 ,3, 0);
addWeights("results/merged/ggH1502LAndTau2Nu_me.root", 3.33417e-06*$lumiEM, 9150 ,3, 0);
addWeights("results/merged/qqH1502LAndTau2Nu_me.root", 3.05939e-07*$lumiEM, 8150 ,3, 0);
addWeights("results/merged/ggH1552LAndTau2Nu_me.root", 3.55826e-06*$lumiEM, 9155 ,3, 0);
addWeights("results/merged/qqH1552LAndTau2Nu_me.root", 3.33338e-07*$lumiEM, 8155 ,3, 0);
addWeights("results/merged/ggH1602LAndTau2Nu_me.root", 3.79728e-06*$lumiEM, 9160 ,3, 0);
addWeights("results/merged/qqH1602LAndTau2Nu_me.root", 4.37525e-07*$lumiEM, 8160 ,3, 0);
addWeights("results/merged/ggH1702LAndTau2Nu_me.root", 3.41403e-06*$lumiEM, 9170 ,3, 0);
addWeights("results/merged/qqH1702LAndTau2Nu_me.root", 3.61367e-07*$lumiEM, 8170 ,3, 0);
addWeights("results/merged/ggH1802LAndTau2Nu_me.root", 2.89425e-06*$lumiEM, 9180 ,3, 0);
addWeights("results/merged/qqH1802LAndTau2Nu_me.root", 0*$lumiEM, 8180 ,3, 0);
addWeights("results/merged/ggH1902LAndTau2Nu_me.root", 2.16421e-06*$lumiEM, 9190 ,3, 0);
addWeights("results/merged/qqH1902LAndTau2Nu_me.root", 2.51662e-07*$lumiEM, 8190 ,3, 0);
addWeights("results/merged/ggH2002LAndTau2Nu_me.root", 3.69593e-06*$lumiEM, 9200 ,3, 0);
addWeights("results/merged/qqH2002LAndTau2Nu_me.root", 2.19108e-07*$lumiEM, 8200 ,3, 0);
addWeights("results/merged/ggH2502LAndTau2Nu_me.root", 1.178e-06*$lumiEM, 9250 ,3, 0);
addWeights("results/merged/qqH2502LAndTau2Nu_me.root", 0*$lumiEM, 8250 ,3, 0);
addWeights("results/merged/ggH3002LAndTau2Nu_me.root", 8.75005e-07*$lumiEM, 9300 ,3, 0);
addWeights("results/merged/qqH3002LAndTau2Nu_me.root", 0*$lumiEM, 8300 ,3, 0);
addWeights("results/merged/ggH3502LAndTau2Nu_me.root", 2.41703e-06*$lumiEM, 9350 ,3, 0);
addWeights("results/merged/qqH3502LAndTau2Nu_me.root", 1.51389e-07*$lumiEM, 8350 ,3, 0);
addWeights("results/merged/ggH4002LAndTau2Nu_me.root", 5.99725e-07*$lumiEM, 9400 ,3, 0);
addWeights("results/merged/qqH4002LAndTau2Nu_me.root", 5.19202e-08*$lumiEM, 8400 ,3, 0);
addWeights("results/merged/ggH4502LAndTau2Nu_me.root", 3.86613e-07*$lumiEM, 9450 ,3, 0);
addWeights("results/merged/qqH4502LAndTau2Nu_me.root", 4.65678e-08*$lumiEM, 8450 ,3, 0);
addWeights("results/merged/ggH5002LAndTau2Nu_me.root", 2.47147e-07*$lumiEM, 9500 ,3, 0);
addWeights("results/merged/qqH5002LAndTau2Nu_me.root", 0*$lumiEM, 8500 ,3, 0);
addWeights("results/merged/ggH5502LAndTau2Nu_me.root", 1.58151e-07*$lumiEM, 9550 ,3, 0);
addWeights("results/merged/qqH5502LAndTau2Nu_me.root", 0*$lumiEM, 8550 ,3, 0);
addWeights("results/merged/ggH6002LAndTau2Nu_me.root", 1.03169e-07*$lumiEM, 9600 ,3, 0);
addWeights("results/merged/qqH6002LAndTau2Nu_me.root", 1.90246e-08*$lumiEM, 8600 ,3, 0);
addWeights("results/merged/ggH7002LAndTau2Nu_me.root", 4.62052e-08*$lumiEM, 9700 ,3, 0);
addWeights("results/merged/qqH7002LAndTau2Nu_me.root", 1.54801e-08*$lumiEM, 8700 ,3, 0);
addWeights("results/merged/ggH8002LAndTau2Nu_me.root", 2.27653e-08*$lumiEM, 9800 ,3, 0);
addWeights("results/merged/qqH8002LAndTau2Nu_me.root", 0*$lumiEM, 8800 ,3, 0);
addWeights("results/merged/ggH9002LAndTau2Nu_me.root", 1.21167e-08*$lumiEM, 9900 ,3, 0);
addWeights("results/merged/qqH9002LAndTau2Nu_me.root", 6.82463e-09*$lumiEM, 8900 ,3, 0);
addWeights("results/merged/ggH10002LAndTau2Nu_me.root", 6.88254e-09*$lumiEM, 10000 ,3, 0);
addWeights("results/merged/qqH10002LAndTau2Nu_me.root", 0*$lumiEM, 9000 ,3, 0);
addWeights("results/merged/Wjets_me.root", 0.00233401*$lumiEM, 80 ,3, 0);
addWeights("results/merged/Zll_Lo_me.root", 0.000123242*$lumiEM, 36 ,3, 0);
addWeights("results/merged/Zll_Hi_me.root", 0.000117914*$lumiEM, 37 ,3, 0);
addWeights("results/merged/SingleT_tWChannel_me.root", 2.81079e-05*$lumiEM, 19 ,3, 0);
addWeights("results/merged/SingleTbar_tWChannel_me.root", 2.26509e-05*$lumiEM, 20 ,3, 0);
addWeights("results/merged/TTbar_me.root", 3.36812e-05*$lumiEM, 10 ,3, 0);
addWeights("results/merged/WW_me.root", 3.00652e-06*$lumiEM, 0 ,3, 0);
addWeights("results/merged/ggWW_me.root", 1.66249e-06*$lumiEM, 1 ,3, 0);
addWeights("results/merged/WZ_me.root", 2.48505e-07*$lumiEM, 74 ,3, 0);
addWeights("results/merged/ZZ_me.root", 7.32849e-07*$lumiEM, 71 ,3, 0);
addWeights("results/merged/WW_pythia_me.root", 1.1645e-05*$lumiEM, 6 ,3, 0);
addWeights("results/merged/WGToLNuG_me.root", 0.000113575*$lumiEM, 85 ,3, 0);
addWeights("results/merged/WGstarToLNu2E_me.root", 1.86651e-05*$lumiEM, 82 ,3, 0);
addWeights("results/merged/WGstarToLNu2Mu_me.root", 6.38057e-06*$lumiEM, 83 ,3, 0);
addWeights("results/merged/WGstarToLNu2Tau_me.root", 6.7204e-06*$lumiEM, 84 ,3, 0);
.q

EOF

echo "done weighting."
