#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>

enum { ee=0, mm=1, em=2, me=3 };

// arrays filled from counters
float nEv_init[4];          // channel
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

// arrays filled with final results
float numAtHiggs_0j[4];     // channel
float errAtHiggs_0j[4];     // channel
float numAtHiggs_1j[4];     // channel
float errAtHiggs_1j[4];     // channel

float numAtHiggsMC_0j[4];     // channel
float errAtHiggsMC_0j[4];     // channel
float numAtHiggsMC_1j[4];     // channel
float errAtHiggsMC_1j[4];     // channel

void countEvents(int mass, const char* channel);
float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);

void estimateWjets() {

  // W+jets number of events at WW level measured from data  
  float numAtWW[4];
  float errAtWW[4];

  numAtWW[ee] = 3.8;
  errAtWW[ee] = quadrSum(0.5, 0.6); // stat (+) syst

  numAtWW[me] = 8.8;
  errAtWW[me] = quadrSum(0.8, 1.1); // stat (+) syst
  
  numAtWW[mm] = 6.33;
  errAtWW[mm] = quadrSum(1.28, 3.15); // stat (+) syst

  numAtWW[em] = 6.75;
  errAtWW[em] = quadrSum(1.28, 3.38); // stat (+) syst

  // -------------------------------------------------------------------
  // for 1 mass only: taking MC files with kine variables to estimate the W+jets amount in the WW control region
  char file_mcTree[1000];
  sprintf(file_mcTree,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH120/datasets_trees/Wjets_ll.root");
  cout << "reading file " << file_mcTree << endl;

  TFile *fileWjets = TFile::Open(file_mcTree);
  TTree *treeWjets = (TTree*)fileWjets->Get("T1");
    
  TH1F *WjetsHEE = new TH1F("WjetsHEE","",50,0,180);
  TH1F *WjetsHMM = new TH1F("WjetsHMM","",50,0,180);
  TH1F *WjetsHEM = new TH1F("WjetsHEM","",50,0,180);
  TH1F *WjetsHME = new TH1F("WjetsHME","",50,0,180);

  TH1F *WjetsHiggsSelHEE = new TH1F("WjetsHiggsSelHEE","",50,0,180);
  TH1F *WjetsHiggsSelHMM = new TH1F("WjetsHiggsSelHMM","",50,0,180);
  TH1F *WjetsHiggsSelHEM = new TH1F("WjetsHiggsSelHEM","",50,0,180);
  TH1F *WjetsHiggsSelHME = new TH1F("WjetsHiggsSelHME","",50,0,180);

  TH1F *WjetsHiggsSel1jHEE = new TH1F("WjetsHiggsSel1jHEE","",50,0,180);
  TH1F *WjetsHiggsSel1jHMM = new TH1F("WjetsHiggsSel1jHMM","",50,0,180);
  TH1F *WjetsHiggsSel1jHEM = new TH1F("WjetsHiggsSel1jHEM","",50,0,180);
  TH1F *WjetsHiggsSel1jHME = new TH1F("WjetsHiggsSel1jHME","",50,0,180);

  TH1F *WjetsHEEOut = new TH1F("WjetsHEEOut","",50,0,180);
  TH1F *WjetsHMMOut = new TH1F("WjetsHMMOut","",50,0,180);
  TH1F *WjetsHEMOut = new TH1F("WjetsHEMOut","",50,0,180);
  TH1F *WjetsHMEOut = new TH1F("WjetsHMEOut","",50,0,180);
  
  treeWjets->Project("WjetsHEE","deltaPhi","(WWSel && finalstate==0)*weight*puweight");
  treeWjets->Project("WjetsHMM","deltaPhi","(WWSel && finalstate==1)*weight*puweight");
  treeWjets->Project("WjetsHEM","deltaPhi","(WWSel && finalstate==2)*weight*puweight");
  treeWjets->Project("WjetsHME","deltaPhi","(WWSel && finalstate==3)*weight*puweight");
  
  float nWjetsTotalEE = WjetsHEE->Integral();
  float nWjetsTotalMM = WjetsHMM->Integral();
  float nWjetsTotalEM = WjetsHEM->Integral();
  float nWjetsTotalME = WjetsHME->Integral();

  treeWjets->Project("WjetsHEEOut","deltaPhi","(eleInvMass>100 && WWSel && finalstate==0)*weight*puweight");
  treeWjets->Project("WjetsHMMOut","deltaPhi","(eleInvMass>100 && WWSel && finalstate==1)*weight*puweight");
  treeWjets->Project("WjetsHEMOut","deltaPhi","(eleInvMass>100 && WWSel && finalstate==2)*weight*puweight");
  treeWjets->Project("WjetsHMEOut","deltaPhi","(eleInvMass>100 && WWSel && finalstate==3)*weight*puweight");
  float nWjetsControlEE = WjetsHEEOut->Integral();
  float nWjetsControlEE_err = yieldErrPoisson(nWjetsControlEE,WjetsHEEOut->GetEntries());
  float nWjetsControlMM = WjetsHMMOut->Integral();
  float nWjetsControlMM_err = yieldErrPoisson(nWjetsControlMM,WjetsHMMOut->GetEntries());
  float nWjetsControlEM = WjetsHEMOut->Integral();
  float nWjetsControlEM_err = yieldErrPoisson(nWjetsControlEM,WjetsHEMOut->GetEntries());
  float nWjetsControlME = WjetsHMEOut->Integral();
  float nWjetsControlME_err = yieldErrPoisson(nWjetsControlME,WjetsHMEOut->GetEntries());

  float nWjetsControl_Tot = nWjetsControlEE+nWjetsControlMM+nWjetsControlEM+nWjetsControlME;
  float nWjetsControl_Tot_err = quadrSum(nWjetsControlEE_err,nWjetsControlMM_err,nWjetsControlEM_err,nWjetsControlME_err);

  treeWjets->Project("WjetsHiggsSelHEE","deltaPhi","(finalSelection && njets==0 && finalstate==0)*weight*puweight");
  treeWjets->Project("WjetsHiggsSelHMM","deltaPhi","(finalSelection && njets==0 && finalstate==1)*weight*puweight");
  treeWjets->Project("WjetsHiggsSelHEM","deltaPhi","(finalSelection && njets==0 && finalstate==2)*weight*puweight");
  treeWjets->Project("WjetsHiggsSelHME","deltaPhi","(finalSelection && njets==0 && finalstate==3)*weight*puweight");

  treeWjets->Project("WjetsHiggsSel1jHEE","deltaPhi","(finalSelection && njets==1 && finalstate==0)*weight*puweight");
  treeWjets->Project("WjetsHiggsSel1jHMM","deltaPhi","(finalSelection && njets==1 && finalstate==1)*weight*puweight");
  treeWjets->Project("WjetsHiggsSel1jHEM","deltaPhi","(finalSelection && njets==1 && finalstate==2)*weight*puweight");
  treeWjets->Project("WjetsHiggsSel1jHME","deltaPhi","(finalSelection && njets==1 && finalstate==3)*weight*puweight");
  
  float nWjetsHiggsSelTotal[4], nWjetsHiggsSelTotal_err[4];

  nWjetsHiggsSelTotal[ee] = WjetsHiggsSelHEE->Integral();
  nWjetsHiggsSelTotal_err[ee] = yieldErrPoisson(nWjetsHiggsSelTotal[ee],WjetsHiggsSelHEE->GetEntries());

  nWjetsHiggsSelTotal[mm] = WjetsHiggsSelHMM->Integral();
  nWjetsHiggsSelTotal_err[mm] = yieldErrPoisson(nWjetsHiggsSelTotal[mm],WjetsHiggsSelHMM->GetEntries());

  nWjetsHiggsSelTotal[em] = WjetsHiggsSelHEM->Integral();
  nWjetsHiggsSelTotal_err[em] = yieldErrPoisson(nWjetsHiggsSelTotal[em],WjetsHiggsSelHEM->GetEntries());

  nWjetsHiggsSelTotal[me] = WjetsHiggsSelHME->Integral();
  nWjetsHiggsSelTotal_err[me] = yieldErrPoisson(nWjetsHiggsSelTotal[me],WjetsHiggsSelHME->GetEntries());

  float nWjetsHiggsSel1jTotal[4], nWjetsHiggsSel1jTotal_err[4];

  nWjetsHiggsSel1jTotal[ee] = WjetsHiggsSel1jHEE->Integral();
  nWjetsHiggsSel1jTotal_err[ee] = yieldErrPoisson(nWjetsHiggsSel1jTotal[ee],WjetsHiggsSel1jHEE->GetEntries());

  nWjetsHiggsSel1jTotal[mm] = WjetsHiggsSel1jHMM->Integral();
  nWjetsHiggsSel1jTotal_err[mm] = yieldErrPoisson(nWjetsHiggsSel1jTotal[mm],WjetsHiggsSel1jHMM->GetEntries());

  nWjetsHiggsSel1jTotal[em] = WjetsHiggsSel1jHEM->Integral();
  nWjetsHiggsSel1jTotal_err[em] = yieldErrPoisson(nWjetsHiggsSel1jTotal[em],WjetsHiggsSel1jHEM->GetEntries());

  nWjetsHiggsSel1jTotal[me] = WjetsHiggsSel1jHME->Integral();
  nWjetsHiggsSel1jTotal_err[me] = yieldErrPoisson(nWjetsHiggsSel1jTotal[me],WjetsHiggsSel1jHME->GetEntries());
  
  std::cout <<"MC: in total = " << nWjetsTotalEE+nWjetsTotalMM+nWjetsTotalEM+nWjetsTotalME
            << " events, in the control region = " << nWjetsControl_Tot << " +/- " << nWjetsControl_Tot_err << std::endl;

  float ratioInOut[4];
  float numWWcontrol[4], errWWcontrol[4];

  ratioInOut[ee]   = nWjetsControlEE/nWjetsTotalEE;
  numWWcontrol[ee] = numAtWW[ee]*ratioInOut[ee];
  errWWcontrol[ee] = errAtWW[ee]*ratioInOut[ee];   

  ratioInOut[mm]   = nWjetsControlMM/nWjetsTotalMM;
  numWWcontrol[mm] = numAtWW[mm]*ratioInOut[mm];
  errWWcontrol[mm] = errAtWW[mm]*ratioInOut[mm];   

  ratioInOut[em]   = nWjetsControlEM/nWjetsTotalEM;
  numWWcontrol[em] = numAtWW[em]*ratioInOut[em];
  errWWcontrol[em] = errAtWW[em]*ratioInOut[em];   

  ratioInOut[me]   = nWjetsControlME/nWjetsTotalME;
  numWWcontrol[me] = numAtWW[me]*ratioInOut[me];
  errWWcontrol[me] = errAtWW[me]*ratioInOut[me];   

  std::cout << "summary WW control region:: -------------------------" << std::endl;
  for(int icha=0;icha<4;icha++) {
    char channelName[2];
    if(icha==ee) sprintf(channelName,"EE");
    if(icha==mm) sprintf(channelName,"MM");
    if(icha==em) sprintf(channelName,"EM");
    if(icha==me) sprintf(channelName,"ME");
    std::cout.precision(3);
    std::cout << channelName << ": data = " << numWWcontrol[icha] << " +/- " << errWWcontrol[icha] << std::endl;
   
  }

  float nWjetsControlData_Tot = numWWcontrol[ee] + numWWcontrol[mm] + numWWcontrol[em] + numWWcontrol[me];
  float nWjetsControlData_Tot_err = quadrSum(errWWcontrol[ee],errWWcontrol[mm],errWWcontrol[em],errWWcontrol[me]);

  std::cout << "W+jets TOTAL in WW control region in MC = " << nWjetsControl_Tot << " +/- " << nWjetsControl_Tot_err
            << "\t data = " << nWjetsControlData_Tot << " +/- " << nWjetsControlData_Tot_err << std::endl;
 
  std::cout << "----------------------------------------------------" << std::endl;

  ofstream textfile;
  textfile.open("WjetsYieldsData.txt", ios_base::app);
  textfile.precision(2);

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];
    //    std::cout << "analyzing mass " << mass << std::endl;

    countEvents(mass,"EE");
    countEvents(mass,"MM");
    countEvents(mass,"EM");
    countEvents(mass,"ME");

    for(int icha=0;icha<4;icha++) {
      float ratioEndWW_0j = (nEv_endWW[icha]==0) ? 0. : nEv_end0j[icha]/nEv_endWW[icha];
      numAtHiggs_0j[icha] = numAtWW[icha]*ratioEndWW_0j;
      errAtHiggs_0j[icha] = errAtWW[icha]*ratioEndWW_0j;   

      numAtHiggsMC_0j[icha] = nWjetsHiggsSelTotal[icha];
      errAtHiggsMC_0j[icha] = nWjetsHiggsSelTotal_err[icha];
      
      float ratioEndWW_1j = (nEv_endWW[icha]==0) ? 0. : nEv_end1j[icha]/nEv_endWW[icha];
      numAtHiggs_1j[icha] = numAtWW[icha]*ratioEndWW_1j;
      errAtHiggs_1j[icha] = errAtWW[icha]*ratioEndWW_1j;   

      numAtHiggsMC_1j[icha] = nWjetsHiggsSel1jTotal[icha];
      errAtHiggsMC_1j[icha] = nWjetsHiggsSel1jTotal_err[icha];

      char channelName[2];
      if(icha==ee) sprintf(channelName,"EE");
      if(icha==mm) sprintf(channelName,"MM");
      if(icha==em) sprintf(channelName,"EM");
      if(icha==me) sprintf(channelName,"ME");
      
      textfile << channelName << ": Higgs Mass = " << mass 
               << "\tdata 0 jet = " << numAtHiggs_0j[icha] << " +/- " << errAtHiggs_0j[icha] 
               << "\tdata 1 jet = " << numAtHiggs_1j[icha] << " +/- " << errAtHiggs_1j[icha] 
               << "\tMC 0 jet = " << numAtHiggsMC_0j[icha] << " +/- " << errAtHiggsMC_0j[icha] 
               << "\tMC 1 jet = " << numAtHiggsMC_1j[icha] << " +/- " << errAtHiggsMC_1j[icha]
               << std::endl;
    }

    float numAtHiggs_0j_Tot = numAtHiggs_0j[ee] + numAtHiggs_0j[mm] + numAtHiggs_0j[em] + numAtHiggs_0j[me];
    float errAtHiggs_0j_Tot = quadrSum(errAtHiggs_0j[ee],errAtHiggs_0j[mm],errAtHiggs_0j[em],errAtHiggs_0j[me]); 

    float numAtHiggsMC_0j_Tot = numAtHiggsMC_0j[ee] + numAtHiggsMC_0j[mm] + numAtHiggsMC_0j[em] + numAtHiggsMC_0j[me];
    float errAtHiggsMC_0j_Tot = quadrSum(errAtHiggsMC_0j[ee],errAtHiggsMC_0j[mm],errAtHiggsMC_0j[em],errAtHiggsMC_0j[me]); 
    
    float numAtHiggs_1j_Tot = numAtHiggs_1j[ee] + numAtHiggs_1j[mm] + numAtHiggs_1j[em] + numAtHiggs_1j[me];
    float errAtHiggs_1j_Tot = quadrSum(errAtHiggs_1j[ee],errAtHiggs_1j[mm],errAtHiggs_1j[em],errAtHiggs_1j[me]); 

    float numAtHiggsMC_1j_Tot = numAtHiggsMC_1j[ee] + numAtHiggsMC_1j[mm] + numAtHiggsMC_1j[em] + numAtHiggsMC_1j[me];
    float errAtHiggsMC_1j_Tot = quadrSum(errAtHiggsMC_1j[ee],errAtHiggsMC_1j[mm],errAtHiggsMC_1j[em],errAtHiggsMC_1j[me]); 

    textfile << "\t===>> TOTAL: Higgs Mass = " << mass 
             << "\tdata 0 jet = " << numAtHiggs_0j_Tot << " +/- " << errAtHiggs_0j_Tot 
             << "\tdata 1 jet = " << numAtHiggs_1j_Tot << " +/- " << errAtHiggs_1j_Tot 
             << "\tMC 0 jet = " << numAtHiggsMC_0j_Tot << " +/- " << errAtHiggsMC_0j_Tot 
             << "\tMC 1 jet = " << numAtHiggsMC_1j_Tot << " +/- " << errAtHiggsMC_1j_Tot 
             << std::endl;

  }

  std::cout << "Full W+jets yields in data in:  WjetsYieldsData.txt " << std::endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

float yieldErrPoisson(float nEst1, float n1, float nEst2, float n2, float nEst3, float n3, float nEst4, float n4, float nEst5, float n5, float nEst6, float n6) {

  float sum=0;
  if(n1>0) sum += pow(nEst1,2)/n1;
  if(n2>0) sum += pow(nEst2,2)/n2;
  if(n3>0) sum += pow(nEst3,2)/n3;
  if(n4>0) sum += pow(nEst4,2)/n4;
  if(n5>0) sum += pow(nEst5,2)/n5;
  if(n6>0) sum += pow(nEst6,2)/n6;
  
  return sqrt(sum);
}

void countEvents(int mass, const char *channel) {

  // taking the EE or ME trees for the wanted mass
  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_%s",channel);  
  TChain *theChain = new TChain(nametree);
  
  //  std::cout << "in countEvents: analyzing mass " << mass << std::endl;
    
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc23_2/crovelli/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH%d/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root",mass);  
  theChain->Add(file_mc);
  //  cout << "reading tree " << nametree << " from file " << file_mc << endl;    

  int theCha=-1;
  if(TString(channel).Contains("EE")) theCha=ee;
  if(TString(channel).Contains("MM")) theCha=mm;
  if(TString(channel).Contains("EM")) theCha=em;
  if(TString(channel).Contains("ME")) theCha=me;
  
  // number of events at the wanted step of the selection
  nEv_init[theCha]  = 0.0;
  nEv_endWW[theCha] = 0.0;
  nEv_end0j[theCha] = 0.0;
  nEv_end1j[theCha] = 0.0;

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
    nEv_init[theCha]  += nSel[0];
    nEv_endWW[theCha] += nSel[16];
    nEv_end0j[theCha] += nSel[22];
    nEv_end1j[theCha] += nSel[23];
  }
}

