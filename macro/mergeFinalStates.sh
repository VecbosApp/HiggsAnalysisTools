#! /bin/sh
# ./mergeFinalStates.sh merges the trees for ee/emu/mm final states in one for a combined fit
# the variable "finalstate in the tree will distinguish the three"

Hmass=$1

# usage: ./mergeFinalStates.sh Hmass

hadd -f results/datasets_trees/H$Hmass\_ll.root results/datasets_trees/H$Hmass\_ee.root results/datasets_trees/H$Hmass\_mm.root results/datasets_trees/H$Hmass\_em.root results/datasets_trees/H$Hmass\_me.root
hadd -f results/datasets_trees/WW_ll.root results/datasets_trees/WW_ee.root results/datasets_trees/WW_mm.root results/datasets_trees/WW_em.root results/datasets_trees/WW_me.root
hadd -f results/datasets_trees/top_ll.root results/datasets_trees/top_ee.root results/datasets_trees/top_mm.root results/datasets_trees/top_em.root results/datasets_trees/top_me.root
hadd -f results/datasets_trees/Wjets_ll.root results/datasets_trees/Wjets_ee.root results/datasets_trees/Wjets_mm.root results/datasets_trees/Wjets_em.root results/datasets_trees/Wjets_me.root
hadd -f results/datasets_trees/Zjets_ll.root results/datasets_trees/Zjets_ee.root results/datasets_trees/Zjets_mm.root results/datasets_trees/Zjets_em.root results/datasets_trees/Zjets_me.root
hadd -f results/datasets_trees/others_ll.root results/datasets_trees/others_ee.root results/datasets_trees/others_mm.root results/datasets_trees/others_em.root results/datasets_trees/others_me.root
hadd -f results_data/merged/dataset_ll.root results_data/merged/dataset_ee.root results_data/merged/dataset_mm.root results_data/merged/dataset_em.root results_data/merged/dataset_me.root
