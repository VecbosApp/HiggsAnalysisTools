#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

// arrays filled from counters
float nEv_init[9];          // mass
float nEv_endWW[9];
float nEv_end0j[9];
float nEv_end1j[9];

// arrays filled with final results
float numAtHiggs_0j[9];     // mass
float errAtHiggs_0j[9];     // mass
float numAtHiggs_1j[9];     // mass
float errAtHiggs_1j[9];     // mass

void countEvents(int mass);

void estimateWjets() {

  // W+jets number of events at WW level measured from data  
  float numAtWW = 7.45;
  float errAtWW = 1.27;     
  

  // -------------------------------------------------------------------
  // for 1 mass only: taking MC files with kine variables to estimate the W+jets amount in the WW control region
  char file_mcTree[1000];
  sprintf(file_mcTree,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V12/OptimMH120/datasets_trees/Wjets_ee.root");   
  // sprintf(file_mcTree,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V12/OptimMH120/datasets_trees/Wjets_me.root");
  cout << "reading file " << file_mcTree << endl;

  TFile *fileWjets = TFile::Open(file_mcTree);
  TTree *treeWjets = (TTree*)fileWjets->Get("T1");
    
  TH1F *WjetsH = new TH1F("WjetsH","",50,0,180);
  
  treeWjets->Project("WjetsH","deltaPhi","(WWSel)*weight");
  float nWjetsTotal = WjetsH->Integral();

  treeWjets->Project("WjetsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  float nWjetsControl = WjetsH->Integral();

  std::cout <<"in total = " << nWjetsTotal << " events, in the control region = " << nWjetsControl << std::endl;

  float ratioInOut   = nWjetsControl/nWjetsTotal;
  float numWWcontrol = numAtWW*ratioInOut;
  float errWWcontrol = errAtWW*ratioInOut;   


  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int myMass=0; myMass<9; myMass++) {
    
    int mass = 120 + myMass*10; 
    std::cout << "analyzing mass " << mass << std::endl;

    countEvents(myMass);

    float ratioEndWW_0j = nEv_end0j[myMass]/nEv_endWW[myMass];
    numAtHiggs_0j[myMass] = numAtWW*ratioEndWW_0j;
    errAtHiggs_0j[myMass] = errAtWW*ratioEndWW_0j;   
    // cout << nEv_end0j[myMass] << " " << nEv_endWW[myMass] << " " << ratioEndWW_0j  << endl;

    float ratioEndWW_1j = nEv_end1j[myMass]/nEv_endWW[myMass];
    numAtHiggs_1j[myMass] = numAtWW*ratioEndWW_1j;
    errAtHiggs_1j[myMass] = errAtWW*ratioEndWW_1j;   
    cout << nEv_end1j[myMass] << " " << nEv_endWW[myMass] << " " << ratioEndWW_1j  << endl;
  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "summary: -------------------------" << std::endl;
  std::cout << "in the WW control region = "   << numWWcontrol << " +- " << errWWcontrol << std::endl; 
  std::cout << std::endl;
  std::cout << "at the end of the Higgs, 0 jets, selection: " << std::endl;
  // for (int myMass=0; myMass<9; myMass++) 
  //  cout << "mass = " << 120 + myMass*10 << ": " << numAtHiggs_0j[myMass] << " +- " << errAtHiggs_0j[myMass] << std::endl;
  // std::cout << std::endl;
  for (int myMass=0; myMass<9; myMass++) 
    cout << "mass = " << 120 + myMass*10 << ": " << numAtHiggs_1j[myMass] << " +- " << errAtHiggs_1j[myMass] << std::endl;
  std::cout << std::endl;

}

void countEvents(int myMass) {

  // taking the EE or ME trees for the wanted mass
  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EE");  
  // sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_ME");
  TChain *theChain = new TChain(nametree);
  
  int mass = 120 + myMass*10; 
  std::cout << "in countEvents: analyzing mass " << mass << std::endl;
    
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V12/OptimMH%d/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root",mass);  
  theChain->Add(file_mc);
  cout << "reading tree " << nametree << " from file " << file_mc << endl;    
  
  // number of events at the wanted step of the selection
  nEv_init[myMass]  = 0.0;
  nEv_endWW[myMass] = 0.0;
  nEv_end0j[myMass] = 0.0;
  nEv_end1j[myMass] = 0.0;

  // reading the tree
  Int_t    nCuts;
  Float_t  nSel[25];   //[nCuts]                                      
  TBranch  *b_nCuts;   
  TBranch  *b_nSel;    
  theChain->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
  theChain->SetBranchAddress("nSel", nSel, &b_nSel);
  
  Long64_t nentries = theChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = theChain->GetEntry(jentry);   
    nbytes    += nb;
    nEv_init[myMass]  += nSel[0];
    nEv_endWW[myMass] += nSel[16];
    nEv_end0j[myMass] += nSel[22];
    nEv_end1j[myMass] += nSel[23];
  }
}

