#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void countEvents() {
  
  char nametree[400];
  sprintf(nametree,"PRESELECTION_EVENT_COUNTER");
  
  cout << "nametree = " << nametree << endl;
  
  TChain *chains[6];
  for(int isample=0; isample<6; isample++) {
    chains[isample] = new TChain(nametree);
  }
  
  chains[0]->Add("results/WJetsMADGRAPH/WJets-madgraph/*Counters.root");
  
  chains[1]->Add("results/WJetsALPGEN/W1Jets_Pt0to100-alpgen/*Counters.root");
  chains[2]->Add("results/WJetsALPGEN/W1Jets_Pt100to300-alpgen/*Counters.root");
  chains[3]->Add("results/WJetsALPGEN/W1Jets_Pt300to800-alpgen/*Counters.root");
  chains[4]->Add("results/WJetsALPGEN/W1Jets_Pt800to1600-alpgen/*Counters.root");
  
  chains[5]->Add("results/DiBosons/WW_2l_7TeV/*Counters.root");
  
  cout << "chains added. " << endl;

  std::vector<std::string> sampleName;

  sampleName.push_back("WJetsMADGRAPH");
  sampleName.push_back("WJetsALPGEN_0-100");
  sampleName.push_back("WJetsALPGEN_100-300");
  sampleName.push_back("WJetsALPGEN_300-800");
  sampleName.push_back("WJetsALPGEN_800-1600");
  sampleName.push_back("WW");

  float nEv[6];
  
  for(int isample=0; isample<6; isample++) {
    nEv[isample] = 0.0;
  }
  
  for(int isample=0; isample<6; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;
    
    Int_t           nCuts;
    Float_t         nSel[10];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;
      
      nEv[isample] += nSel[0];
    }
  }
  
  for(int isample=0; isample<6; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
  }
  
}

