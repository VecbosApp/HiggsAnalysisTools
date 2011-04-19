#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTreesData.sh

echo "Adding weights for ee datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results_data/merged/dataset_eg_ee.root", 1.0, 0);

.q

EOF

echo "Adding weights for mm datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results_data/merged/dataset_mu_mm.root", 1.0, 1);

.q

EOF


echo "Adding weights for em datasets..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addWeights("results_data/merged/dataset_mu_em.root", 1.0, 2);
addWeights("results_data/merged/dataset_eg_em.root", 1.0, 2);

.q

EOF

echo "done weighting."