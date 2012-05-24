#define addWeightsToTreeLooseLoose_cxx
#include "addWeightsToTreeLooseLoose.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void addWeightsToTreeLooseLoose::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L addWeightsToTreeLooseLoose.C
  //      Root > addWeightsToTreeLooseLoose t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  float WWSel0j, WWSel1j;

  TFile *fileNew = TFile::Open("/cmsrm/pc24_2/emanuele/data/Higgs5.2.X/Data_May18JSON_V1/datasets_trees/dataset_looseloose_wwbits.root","recreate");
  TTree *treeNew = new TTree("latino","tree with only selected events");

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("WWSel0j", &WWSel0j, "WWSel0j/F");
    theTreeNew->Branch("WWSel1j", &WWSel1j, "WWSel1j/F");
  }

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000 == 0) std::cout << ">>> Weighting event # " << jentry << " / " << nentries << " entries" << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    WWSel0j = trigger==1. && pfmet>20. && mll>12 && zveto==1 && mpmet>20. && bveto_mu==1 && nextra==0 && (bveto_ip==1 &&  (njet != 1  || nbjet==0) && ((njet<2 || njet>3) || (jetbjpb1<=1.05 && jetbjpb2<=1.05)))  && ptll>45. &&   ( !sameflav || ( (njet==0 || dymva0>0.60) && (njet==1 || dymva1>0.30) && ( njet==0 || njet==1 || (pfmet > (40.0+nvtx/2.0))) ) ) && njet==0;
    WWSel1j = trigger==1. && pfmet>20. && mll>12 && zveto==1 && mpmet>20. && bveto_mu==1 && nextra==0 && (bveto_ip==1 &&  (njet != 1  || nbjet==0) && ((njet<2 || njet>3) || (jetbjpb1<=1.05 && jetbjpb2<=1.05)))  && ptll>45. &&   ( !sameflav || ( (njet==0 || dymva0>0.60) && (njet==1 || dymva1>0.30) && ( njet==0 || njet==1 || (pfmet > (40.0+nvtx/2.0))) ) ) && njet==1;

    treeNew->Fill();
  }

  fileNew->cd();
  treeNew->Write();
  fileNew->Close();

}
