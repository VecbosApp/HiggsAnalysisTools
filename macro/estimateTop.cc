#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass, const char* channel);

void estimateTop() {

  // constants
  float eff_2b = 0.4983;
  float eff_2b_err = 0.03; 
  float eff_2b_softmu = 0.1677; // MC 

  TFile *fileData = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeTop= (TTree*)fileTop->Get("T1");

  TH1F *topHEE = new TH1F("topHEE","",50,0,180);
  TH1F *topHMM = new TH1F("topHMM","",50,0,180);
  TH1F *topHEM = new TH1F("topHEM","",50,0,180);
  TH1F *topHME = new TH1F("topHME","",50,0,180);
  TH1F *btagHDataEE = new TH1F("btagHDataEE","",50,0,180);
  TH1F *btagHDataMM = new TH1F("btagHDataMM","",50,0,180);
  TH1F *btagHDataEM = new TH1F("btagHDataEM","",50,0,180);
  TH1F *btagHDataME = new TH1F("btagHDataME","",50,0,180);

  treeTop->Project("topHEE","deltaPhi","(WWSel && finalstate==0)*weight*puweight");
  treeTop->Project("topHMM","deltaPhi","(WWSel && finalstate==1)*weight*puweight");
  treeTop->Project("topHEM","deltaPhi","(WWSel && finalstate==2)*weight*puweight");
  treeTop->Project("topHME","deltaPhi","(WWSel && finalstate==3)*weight*puweight");
  treeData->Project("btagHDataEE","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==0");
  treeData->Project("btagHDataMM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==1");
  treeData->Project("btagHDataEM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==2");
  treeData->Project("btagHDataME","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==3");

  ///// TOP ESTIMATION ///////
  float nTopMC[4];
  float nTopMC_err[4];
  
  nTopMC[ee] = topHEE->Integral();
  nTopMC_err[ee] = yieldErrPoisson(nTopMC[ee],topHEE->GetEntries());
  nTopMC[mm] = topHMM->Integral();
  nTopMC_err[mm] = yieldErrPoisson(nTopMC[mm],topHMM->GetEntries());
  nTopMC[em] = topHEM->Integral();
  nTopMC_err[em] = yieldErrPoisson(nTopMC[em],topHEM->GetEntries());
  nTopMC[me] = topHME->Integral();
  nTopMC_err[me] = yieldErrPoisson(nTopMC[me],topHME->GetEntries());
  
  // top estimation from data (0-jet bin method)
  float nTopData[4];
  float nTopData_err[4];
  // EE
  float nBTagTag_data = btagHDataEE->Integral();
  float nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  float nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  float nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  float nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EE = " << nBTagTag_data << std::endl;

  nTopData[ee] = nTopSoftMuVeto_data;
  nTopData_err[ee] = nTopSoftMuVeto_data_err;

  // MM
  nBTagTag_data = btagHDataMM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "MM = " << nBTagTag_data << std::endl;

  nTopData[mm] = nTopSoftMuVeto_data;
  nTopData_err[mm] = nTopSoftMuVeto_data_err;

  // EM
  nBTagTag_data = btagHDataEM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EM = " << nBTagTag_data << std::endl;

  nTopData[em] = nTopSoftMuVeto_data;
  nTopData_err[em] = nTopSoftMuVeto_data_err;

  // ME
  nBTagTag_data = btagHDataME->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "ME = " << nBTagTag_data << std::endl;

  nTopData[me] = nTopSoftMuVeto_data;
  nTopData_err[me] = nTopSoftMuVeto_data_err;

  std::cout << "Number of Top events from data at W+W- level: " << std::endl
            << "*\tEE = " << (150./126.) * nTopData[ee] << " +/- " << nTopData_err[ee] << std::endl
            << "*\tMM = " << (150./126.) * nTopData[mm] << " +/- " << nTopData_err[mm] << std::endl
            << "*\tEM = " << (150./126.) * nTopData[em] << " +/- " << nTopData_err[em] << std::endl
            << "*\tME = " << (150./126.) * nTopData[me] << " +/- " << nTopData_err[me] << std::endl;

  std::cout << "Number of Top events from MC at W+W- level: " << std::endl
            << "*\tEE = " << (150./126.) * nTopMC[ee] << " +/- " << nTopMC_err[ee] << std::endl
            << "*\tMM = " << (150./126.) * nTopMC[mm] << " +/- " << nTopMC_err[mm] << std::endl
            << "*\tEM = " << (150./126.) * nTopMC[em] << " +/- " << nTopMC_err[em] << std::endl
            << "*\tME = " << (150./126.) * nTopMC[me] << " +/- " << nTopMC_err[me] << std::endl;

  ofstream textfile;
  textfile.open("TopYieldsData.txt", ios_base::app);
  textfile.precision(2);

  ofstream tablefile1;
  tablefile1.open("TopYieldsData_ForTable_0j.txt", ios_base::app);
  tablefile1.precision(2);

  ofstream tablefile2;
  tablefile2.open("TopYieldsData_ForTable_1j.txt", ios_base::app);
  tablefile2.precision(2);

  ofstream tablefile3;
  tablefile3.open("TopYieldsMC_ForTable_0j.txt", ios_base::app);
  tablefile3.precision(2);

  ofstream tablefile4;
  tablefile4.open("TopYieldsMC_ForTable_1j.txt", ios_base::app);
  tablefile4.precision(2);

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    countEvents(mass,"EE");
    countEvents(mass,"MM");
    countEvents(mass,"EM");
    countEvents(mass,"ME");

    float nTopData_HiggsSel_0j[4], nTopData_HiggsSel_1j[4];
    float nTopData_HiggsSel_0j_err[4], nTopData_HiggsSel_1j_err[4];

    float nTopMC_HiggsSel_0j[4],  nTopMC_HiggsSel_1j[4];
    float nTopMC_HiggsSel_0j_err[4], nTopMC_HiggsSel_1j_err[4];

    for(int icha=0;icha<4;icha++) {
      float eff_0j = nEv_end0j[icha]/nEv_endWW[icha];
      float eff_0j_err = sqrt(eff_0j*(1-eff_0j)/nEv_endWW[icha]);
      float effErrRel = (eff_0j==0) ? 0. : eff_0j_err/eff_0j;
    
      nTopData_HiggsSel_0j[icha] = nTopData[icha] * eff_0j;
      float topDataErrRel = (nTopData[icha]==0) ? 0. : nTopData_err[icha]/nTopData[icha];
      nTopData_HiggsSel_0j_err[icha] = nTopData_HiggsSel_0j[icha] * quadrSum(topDataErrRel,effErrRel);

      nTopMC_HiggsSel_0j[icha] = nTopMC[icha] * eff_0j;
      float topMCErrRel = (nTopMC[icha]==0) ? 0. : nTopMC_err[icha]/nTopMC[icha];
      nTopMC_HiggsSel_0j_err[icha] = nTopMC_HiggsSel_0j[icha] * quadrSum(topMCErrRel,effErrRel);

      float eff_1j = nEv_end1j[icha]/nEv_endWW[icha];
      float eff_1j_err = sqrt(eff_1j*(1-eff_1j)/nEv_endWW[icha]);
      effErrRel = (eff_1j==0) ? 0. : eff_1j_err/eff_1j;
    
      nTopData_HiggsSel_1j[icha] = nTopData[icha] * eff_1j;
      topDataErrRel = (nTopData[icha]==0) ? 0. : nTopData_err[icha]/nTopData[icha];
      nTopData_HiggsSel_1j_err[icha] = nTopData_HiggsSel_1j[icha] * quadrSum(topDataErrRel,effErrRel);

      nTopMC_HiggsSel_1j[icha] = nTopMC[icha] * eff_1j;
      topMCErrRel = (nTopMC[icha]==0) ? 0. : nTopMC_err[icha]/nTopMC[icha];
      nTopMC_HiggsSel_1j_err[icha] = nTopMC_HiggsSel_1j[icha] * quadrSum(topMCErrRel,effErrRel);

      char channelName[2];
      if(icha==ee) sprintf(channelName,"EE");
      if(icha==mm) sprintf(channelName,"MM");
      if(icha==em) sprintf(channelName,"EM");
      if(icha==me) sprintf(channelName,"ME");
      
      textfile << channelName << ": Higgs Mass = " << mass 
               << "\tdata 0 jet = " << nTopData_HiggsSel_0j[icha] << " +/- " << nTopData_HiggsSel_0j_err[icha] 
               << "\tdata 1 jet = " << nTopData_HiggsSel_1j[icha] << " +/- " << nTopData_HiggsSel_1j_err[icha] 
               << "\tMC 0 jet = " << nTopMC_HiggsSel_0j[icha] << " +/- " << nTopMC_HiggsSel_0j_err[icha] 
               << "\tMC 1 jet = " << nTopMC_HiggsSel_1j[icha] << " +/- " << nTopMC_HiggsSel_1j_err[icha]
               << std::endl;
    }

    // summary table for limits
    if (i==0) { 
      tablefile1 << "# zero jets bin data" << endl;
      tablefile1 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile1 << mass 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_0j[1] << " +/- " <<  nTopData_HiggsSel_0j_err[1] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_0j[3] << " +/- " <<  nTopData_HiggsSel_0j_err[3] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_0j[2] << " +/- " <<  nTopData_HiggsSel_0j_err[2] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_0j[0] << " +/- " <<  nTopData_HiggsSel_0j_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile2 << "#one jet bin data" << endl;
      tablefile2 << "#\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile2 << mass 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_1j[1] << " +/- " <<  nTopData_HiggsSel_1j_err[1] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_1j[3] << " +/- " <<  nTopData_HiggsSel_1j_err[3] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_1j[2] << " +/- " <<  nTopData_HiggsSel_1j_err[2] 
	       << " " << "\t" << (150./126.) * nTopData_HiggsSel_1j[0] << " +/- " <<  nTopData_HiggsSel_1j_err[0] 
	       << std::endl;

    if (i==0) { 
      tablefile3 << "# zero jets bin MC" << endl;
      tablefile3 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile3 << mass 
	       << " " << "\t" << nTopMC_HiggsSel_0j[1] << " +/- " <<  nTopMC_HiggsSel_0j_err[1] 
	       << " " << "\t" << nTopMC_HiggsSel_0j[3] << " +/- " <<  nTopMC_HiggsSel_0j_err[3] 
	       << " " << "\t" << nTopMC_HiggsSel_0j[2] << " +/- " <<  nTopMC_HiggsSel_0j_err[2] 
	       << " " << "\t" << nTopMC_HiggsSel_0j[0] << " +/- " <<  nTopMC_HiggsSel_0j_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile4 << "#one jet bin MC" << endl;
      tablefile4 << "#\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile4 << mass 
	       << " " << "\t" << nTopMC_HiggsSel_1j[1] << " +/- " <<  nTopMC_HiggsSel_1j_err[1] 
	       << " " << "\t" << nTopMC_HiggsSel_1j[3] << " +/- " <<  nTopMC_HiggsSel_1j_err[3] 
	       << " " << "\t" << nTopMC_HiggsSel_1j[2] << " +/- " <<  nTopMC_HiggsSel_1j_err[2] 
	       << " " << "\t" << nTopMC_HiggsSel_1j[0] << " +/- " <<  nTopMC_HiggsSel_1j_err[0] 
	       << std::endl;

    float nTopData_HiggsSel_0j_Tot = nTopData_HiggsSel_0j[ee] + nTopData_HiggsSel_0j[mm] + nTopData_HiggsSel_0j[em] + nTopData_HiggsSel_0j[me];
    float nTopData_HiggsSel_0j_Tot_err = quadrSum(nTopData_HiggsSel_0j_err[ee],nTopData_HiggsSel_0j_err[mm],nTopData_HiggsSel_0j_err[em],nTopData_HiggsSel_0j_err[me]);

    float nTopMC_HiggsSel_0j_Tot = nTopMC_HiggsSel_0j[ee] + nTopMC_HiggsSel_0j[mm] + nTopMC_HiggsSel_0j[em] + nTopMC_HiggsSel_0j[me];
    float nTopMC_HiggsSel_0j_Tot_err = quadrSum(nTopMC_HiggsSel_0j_err[ee],nTopMC_HiggsSel_0j_err[mm],nTopMC_HiggsSel_0j_err[em],nTopMC_HiggsSel_0j_err[me]);

    float nTopData_HiggsSel_1j_Tot = nTopData_HiggsSel_1j[ee] + nTopData_HiggsSel_1j[mm] + nTopData_HiggsSel_1j[em] + nTopData_HiggsSel_1j[me];
    float nTopData_HiggsSel_1j_Tot_err = quadrSum(nTopData_HiggsSel_1j_err[ee],nTopData_HiggsSel_1j_err[mm],nTopData_HiggsSel_1j_err[em],nTopData_HiggsSel_1j_err[me]);

    float nTopMC_HiggsSel_1j_Tot = nTopMC_HiggsSel_1j[ee] + nTopMC_HiggsSel_1j[mm] + nTopMC_HiggsSel_1j[em] + nTopMC_HiggsSel_1j[me];
    float nTopMC_HiggsSel_1j_Tot_err = quadrSum(nTopMC_HiggsSel_1j_err[ee],nTopMC_HiggsSel_1j_err[mm],nTopMC_HiggsSel_1j_err[em],nTopMC_HiggsSel_1j_err[me]);
    
    textfile.precision(2);
    textfile << "\t===>> TOTAL: Higgs Mass = " << mass 
             << "\tdata 0 jet = " << nTopData_HiggsSel_0j_Tot << " +/- " << nTopData_HiggsSel_0j_Tot_err 
             << "\tdata 1 jet = " << nTopData_HiggsSel_1j_Tot << " +/- " << nTopData_HiggsSel_1j_Tot_err
             << "\tMC 0 jet = " << nTopMC_HiggsSel_0j_Tot << " +/- " << nTopMC_HiggsSel_0j_Tot_err 
             << "\tMC 1 jet = " << nTopMC_HiggsSel_1j_Tot << " +/- " << nTopMC_HiggsSel_1j_Tot_err
             << std::endl;

  }

  std::cout << "Full top yields in data in:  TopYieldsData.txt " << std::endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr) {
  float val = ntag * (1-eff2b) / eff2b;
  float err = ntag * eff2berr / pow(eff2b,2);
  return std::make_pair(val,err);
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
  
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH%d/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root",mass);  
  theChain->Add(file_mc);
  //  cout << "reading tree " << nametree << " from file " << file_mc << endl;    
  
  int theCha=-1;
  if(TString(channel).Contains("EE")) theCha=ee;
  if(TString(channel).Contains("MM")) theCha=mm;
  if(TString(channel).Contains("EM")) theCha=em;
  if(TString(channel).Contains("ME")) theCha=me;

  // number of events at the wanted step of the selection
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
    nEv_endWW[theCha] += nSel[16];
    nEv_end0j[theCha] += nSel[22];
    nEv_end1j[theCha] += nSel[23];
  }
}

