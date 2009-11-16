#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

using namespace std;

void addWeights(const char* filename, float weight) {

  cout << "Adding weight branch to file " << filename << " with weight " << weight << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if ( treeOrig ) {
  int nentriesOrig = treeOrig->GetEntries();

  TFile *fileNew = TFile::Open(filename,"recreate");
  TTree *treeNew = new TTree("T1","tree with only selected events");

  // add also a branch with jet category (1 for njets=0, -1 for njets=1: useful for the fit)
  // and a branch with float final selection bool (for roofit)
  Int_t           njets;
  Float_t         met;
  Float_t         deltaPhi;
  Float_t         eleInvMass;
  Float_t         maxPtEle;
  Float_t         bTagImpPar;
  Char_t          finalSelection;
  treeOrig->SetBranchAddress("njets", &njets);
  treeOrig->SetBranchAddress("met", &met);
  treeOrig->SetBranchAddress("deltaPhi", &deltaPhi);
  treeOrig->SetBranchAddress("eleInvMass", &eleInvMass);
  treeOrig->SetBranchAddress("maxPtEle", &maxPtEle);
  treeOrig->SetBranchAddress("bTagImpPar", &bTagImpPar);
  treeOrig->SetBranchAddress("finalSelection", &finalSelection);

  // copy branches
  treeNew->Branch("met",                 &met,                 "met/F");  
  treeNew->Branch("deltaPhi",            &deltaPhi,            "deltaPhi/F");  
  treeNew->Branch("eleInvMass",          &eleInvMass,          "eleInvMass/F");  
  treeNew->Branch("maxPtEle",            &maxPtEle,            "maxPtEle/F");  
  treeNew->Branch("bTagImpPar",          &bTagImpPar,          "bTagImpPar/F");

  float jetcat = 0;
  treeNew->Branch("jetcat", &jetcat,  "jetcat/F");
  float event = -1;
  treeNew->Branch("event", &event, "event/F");
  treeNew->Branch("weight", &weight,  "weight/F");

  int j =0;

  for(int i=0; i<nentriesOrig; i++) {
    if (i%1000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl;
    treeOrig->GetEntry(i);

    if (njets==0) jetcat = 1;
    else if(njets==1) jetcat = -1;
    else jetcat = -2;

    // consider only events with 0 or 1 jet
    // and fit variables within fit range
    if (njets<=1 && 
        met>=0 && met<=200 &&
        deltaPhi>=0 && deltaPhi<=180 &&
        maxPtEle>=20 && maxPtEle<=200 &&
        eleInvMass>=12 && eleInvMass<=150 &&
        bTagImpPar>=-1001 && bTagImpPar<=2 &&
        weight>=0 && weight<=10000 &&
        finalSelection) {
      event = (float)j;
      treeNew->Fill();
      j++;
    }

  }
  
  fileNew->cd();
  treeNew->Write();
  fileNew->Close();

  fileOrig->cd();
  fileOrig->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }

}
