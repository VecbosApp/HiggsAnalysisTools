#! /bin/sh 
# this script replaces every merged tree with the same tree with one more branch, containing the event weight for that sample
# the event weight is evaluated with the total number of generated events, cross section and eventual prescale and the wanted luminosity
# the values can be evaluated with the program weights.cc

# usage: ./weightTreesData.sh

mkdir results_data/merged_skim

echo "Adding weights..."
root -l -b <<EOF

.L addWeightsToTree.cc+

addFRWeights();

addWeights("results_data/merged/dataset_fake_ee.root", 1.0, -1, 1);
//addWeights("results_data/merged/dataset_fake_mm.root", 1.0, -1, 0);
//addWeights("results_data/merged/dataset_fake_em.root", 1.0, -1, 2);
addWeights("results_data/merged/dataset_fake_me.root", 1.0, -1, 3);

.q

EOF

echo "done weighting."
