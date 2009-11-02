#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

using namespace std;

void addWeights(const char* filename, float weight) {

  cout << "Adding weight branch to file " << filename << " with weight " << weight << endl;

  TFile *file = 0;
  TTree *tree = 0;

  file = TFile::Open(filename,"Update");
  if( file ) {
    file->cd();
    tree = (TTree*)file->Get("T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if ( tree ) {
  int nentries = tree->GetEntries();

  TBranch *b_weight = tree->Branch("weight", &weight,  "weight/F");


  // add also a branch with jet category (1 for njets=0, -1 for njets=1: useful for the fit)

  int njets;
  tree->SetBranchAddress("njets", &njets);

  float jetcat = 0;
  TBranch *b_jetCat = tree->Branch("jetcat", &jetcat,  "jetcat/F");

  for(int i=0; i<nentries; i++) {
    if (i%1000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentries << " entries" << std::endl;
    tree->GetEntry(i);
    if (njets==0) jetcat = 1;
    else if(njets==1) jetcat = -1;
    else jetcat = -2;

    b_weight->Fill();
    b_jetCat->Fill();
  }
  
  tree->Write("", TObject::kOverwrite);
  file->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }

}
