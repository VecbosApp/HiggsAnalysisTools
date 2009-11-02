#! /bin/sh 
# ./mergeTrees.sh expects results of single sample files in results/merged/
# it creates a merged root file for each single specie of the fit in results/datasets

mkdir -p results/datasets_trees

echo "Now merging species for ee..."
# signal is always a species per se
cp results/merged/H*_ee.root results/datasets_trees

# WW is a species per se
cp results/merged/WW_ee.root results/datasets_trees

# ttbar is a species per se
cp results/merged/TTbarJetsMadgraph_ee.root results/datasets_trees
cp results/merged/TTbar_ee.root results/datasets_trees

# merging all the samples with a Z boson
hadd results/datasets_trees/Zspecies_ee.root  results/merged/ZjetsMadgraph_ee.root results/merged/WZ_ee.root results/merged/ZZ_ee.root

# merging all other backgrounds
hadd results/datasets_trees/others_ee.root results/merged/WjetsMadgraph_ee.root results/merged/Wgamma_ee.root results/merged/QCD_*_ee.root results/merged/PhotonJet_*_ee.root


echo "Now merging species for mm..."
# signal is always a species per se
cp results/merged/H*_mm.root results/datasets_trees

# WW is a species per se
cp results/merged/WW_mm.root results/datasets_trees

# ttbar is a species per se
cp results/merged/TTbarJetsMadgraph_mm.root results/datasets_trees
cp results/merged/TTbar_mm.root results/datasets_trees

# merging all the samples with a Z boson
hadd results/datasets_trees/Zspecies_mm.root  results/merged/ZjetsMadgraph_mm.root results/merged/WZ_mm.root results/merged/ZZ_mm.root

# merging all other backgrounds
hadd results/datasets_trees/others_mm.root results/merged/WjetsMadgraph_mm.root results/merged/Wgamma_mm.root results/merged/InclusiveMu15_mm.root



echo "Now merging species for em..."
# signal is always a species per se
cp results/merged/H*_em.root results/datasets_trees

# WW is a species per se
cp results/merged/WW_em.root results/datasets_trees

# ttbar is a species per se
cp results/merged/TTbarJetsMadgraph_em.root results/datasets_trees
cp results/merged/TTbar_em.root results/datasets_trees

# merging all other backgrounds
hadd results/datasets_trees/others_em.root results/merged/ZjetsMadgraph_em.root results/merged/WZ_em.root results/merged/ZZ_em.root results/merged/WjetsMadgraph_em.root results/merged/Wgamma_em.root results/merged/QCD_*_em.root results/merged/PhotonJet_*_em.root results/merged/InclusiveMu15_em.root
