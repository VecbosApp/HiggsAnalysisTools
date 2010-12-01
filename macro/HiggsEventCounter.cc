#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void countEvents() {

  char nametree[200];
  sprintf(nametree,"PRESELECTION_EVENT_COUNTER");

  cout << "nametree = " << nametree << endl;

  TChain *chains[17];
  for(int isample=0; isample<17; isample++) {
    chains[isample] = new TChain(nametree);
  }

  chains[0]->Add("results/WPYTHIA/WToENu_TuneZ2/*Counters.root");
  chains[1]->Add("results/WPYTHIA/WToMuNu_TuneZ2/*Counters.root");
  chains[2]->Add("results/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*Counters.root");

  chains[3]->Add("results/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*Counters.root");
  chains[4]->Add("results/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*Counters.root");
  chains[5]->Add("results/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*Counters.root");

  chains[6]->Add("results/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*Counters.root");
  chains[7]->Add("results/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*Counters.root");
  chains[8]->Add("results/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*Counters.root");

  chains[9]->Add("results/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root");
  chains[10]->Add("results/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root");
  chains[11]->Add("results/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root");

  chains[12]->Add("results/TTbar/TTJets_TuneD6T/*Counters.root");

  chains[13]->Add("results/HiggsWW/GluGluToHToWWTo2L2Nu_M-130/*Counters.root");

  //  chains[34]->Add("results/DiBosons/Wgamma/*Counters.root");
  chains[14]->Add("results/DiBosons/WWTo2L2Nu_TuneZ2/*Counters.root");
  chains[15]->Add("results/DiBosons/WZTo3LNu_TuneZ2/*Counters.root");
  chains[16]->Add("results/DiBosons/ZZtoAnything_TuneZ2/*Counters.root");

  cout << "chains added. " << endl;

  std::vector<std::string> sampleName;

  sampleName.push_back("W->enu");
  sampleName.push_back("W->munu");
  sampleName.push_back("W->taunu");

  sampleName.push_back("DY(ee) 10<m<20 GeV");
  sampleName.push_back("DY(mm) 10<m<20 GeV");
  sampleName.push_back("DY(tautau) 10<m<20 GeV");

  sampleName.push_back("DY(ee) m>20 GeV");
  sampleName.push_back("DY(mm) m>20 GeV");
  sampleName.push_back("DY(tautau) m>20 GeV");

  sampleName.push_back("single top s");
  sampleName.push_back("single top t");
  sampleName.push_back("single top tW");

  sampleName.push_back("ttbar");

  sampleName.push_back("Higgs");

  //  sampleName.push_back("Wgamma");
  sampleName.push_back("WW");
  sampleName.push_back("WZ");
  sampleName.push_back("ZZ");

  float nEv[17];
  for(int isample=0; isample<17; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<17; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;

    Int_t           nCuts;
    Float_t         nSel[20];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    // loop over files (>1 if VecBos in batch, splitted in many jobs)
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;

      nEv[isample] += nSel[0];
    }
  }

  for(int isample=0; isample<17; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
  }
  
}

