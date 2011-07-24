#! /bin/sh 
# ./mergeTrees.sh expects results of single sample files in results/merged/
# it creates a merged root file for each single specie of the fit in results/datasets

mkdir -p results/datasets_trees

echo "Now merging species for ee..."
# signal is always a species per se
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/datasets_trees/H$i\_ee.root results/merged/ggH$i\2L2Nu_ee.root results/merged/ggH$i\LNuTauNu_ee.root results/merged/qqH$i\2L2Nu_ee.root results/merged/qqH$i\LNuTauNu_ee.root
done

# WW is a species per se
hadd results/datasets_trees/WW_ee.root results/merged/WW_ee.root results/merged/ggWW_ee.root

# Wjets is a species per se
cp results/merged/Wjets_ee.root results/datasets_trees/Wjets_ee.root

# Zjets is a species per se
hadd results/datasets_trees/Zjets_ee.root results/merged/Zee_*_ee.root results/merged/Zmm_*_ee.root results/merged/Ztautau_*_ee.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_ee.root results/merged/TTbar_ee.root results/merged/SingleTop_sChannel_ee.root results/merged/SingleTop_tChannel_ee.root results/merged/SingleTop_tWChannel_ee.root 

# merging all other backgrounds
hadd results/datasets_trees/others_ee.root results/merged/WZ_ee.root results/merged/ZZ_ee.root


echo "Now merging species for mm..."
# signal is always a species per se
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/datasets_trees/H$i\_mm.root results/merged/ggH$i\2L2Nu_mm.root results/merged/ggH$i\LNuTauNu_mm.root results/merged/qqH$i\2L2Nu_mm.root results/merged/qqH$i\LNuTauNu_mm.root
done

# WW is a species per se
hadd results/datasets_trees/WW_mm.root results/merged/WW_mm.root results/merged/ggWW_mm.root

# Wjets is a species per se
cp results/merged/Wjets_mm.root results/datasets_trees/Wjets_mm.root

# Zjets is a species per se
hadd results/datasets_trees/Zjets_mm.root results/merged/Zee_*_mm.root results/merged/Zmm_*_mm.root results/merged/Ztautau_*_mm.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_mm.root results/merged/TTbar_mm.root results/merged/SingleTop_sChannel_mm.root results/merged/SingleTop_tChannel_mm.root results/merged/SingleTop_tWChannel_mm.root 

# merging all other backgrounds
hadd results/datasets_trees/others_mm.root results/merged/WZ_mm.root results/merged/ZZ_mm.root


echo "Now merging species for em..."
# signal is always a species per se
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/datasets_trees/H$i\_em.root results/merged/ggH$i\2L2Nu_em.root results/merged/ggH$i\LNuTauNu_em.root results/merged/qqH$i\2L2Nu_em.root results/merged/qqH$i\LNuTauNu_em.root
done

# WW is a species per se
hadd results/datasets_trees/WW_em.root results/merged/WW_em.root results/merged/ggWW_em.root

# Wjets is a species per se
cp results/merged/Wjets_em.root results/datasets_trees/Wjets_em.root

# Zjets is a species per se
hadd results/datasets_trees/Zjets_em.root results/merged/Zee_*_em.root results/merged/Zmm_*_em.root results/merged/Ztautau_*_em.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_em.root results/merged/TTbar_em.root results/merged/SingleTop_sChannel_em.root results/merged/SingleTop_tChannel_em.root results/merged/SingleTop_tWChannel_em.root 

# merging all other backgrounds
hadd results/datasets_trees/others_em.root results/merged/WZ_em.root results/merged/ZZ_em.root



echo "Now merging species for me..."
# signal is always a species per se
for i in 120 130 140 150 160 170 180 190 200 250 300 350 400 450 500 550 600
do
  hadd results/datasets_trees/H$i\_me.root results/merged/ggH$i\2L2Nu_me.root results/merged/ggH$i\LNuTauNu_me.root results/merged/qqH$i\2L2Nu_me.root results/merged/qqH$i\LNuTauNu_me.root
done

# WW is a species per se
hadd results/datasets_trees/WW_me.root results/merged/WW_me.root results/merged/ggWW_me.root

# Wjets is a species per se
cp results/merged/Wjets_me.root results/datasets_trees/Wjets_me.root

# Zjets is a species per se
hadd results/datasets_trees/Zjets_me.root results/merged/Zee_*_me.root results/merged/Zmm_*_me.root results/merged/Ztautau_*_me.root

# merging ttbar and single t in a species
hadd results/datasets_trees/top_me.root results/merged/TTbar_me.root results/merged/SingleTop_sChannel_me.root results/merged/SingleTop_tChannel_me.root results/merged/SingleTop_tWChannel_me.root 

# merging all other backgrounds
hadd results/datasets_trees/others_me.root results/merged/WZ_me.root results/merged/ZZ_me.root


