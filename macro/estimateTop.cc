#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include "massDependentCuts.cc"

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float usedLumi = 1.091;
float wantedLumi = 1.091;
float scaleFactorLumi = wantedLumi/usedLumi;

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass, const char* channel);

void estimateTop(int njets) {

  // constants
  float eff_2b[2] = { 0.520, 0.66 };
  float eff_2b_err[2] = { 0.080, 0.020 }; 

  char njcut[30];
  sprintf(njcut, "njets==%d", njets);
  char wwselcut[30];
  if(njets==0) sprintf(wwselcut,"WWSel");
  else if(njets==1) sprintf(wwselcut,"WWSel1j");
  else {
    std::cout << "Jet bin must be 0/1" << std::endl;
    return;
  }

  // scale factors for the backgrounds
  float WWDataOverMC[2] = {1.2, 1.3} ; // estimation 1.1 fb-1 
  float DYDataOverMC[2] = {3.0, 8.0}; // estimation 1.1 fb-1, 0j 
  float WjDataTot[4][2]; // [icha][jetbin]
  WjDataTot[ee][0] = 10.6;
  WjDataTot[mm][0] = 12.1;
  WjDataTot[em][0] = 72.9;
  WjDataTot[me][0] = 39.7;
  WjDataTot[ee][1] = 5.0;
  WjDataTot[mm][1] = 6.1;
  WjDataTot[em][1] = 26.4;
  WjDataTot[me][1] = 12.8;
  

  TFile *fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileDY = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileWj = TFile::Open("results/datasets_trees/Wjets_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeTop= (TTree*)fileTop->Get("T1");
  TTree *treeWW= (TTree*)fileWW->Get("T1");
  TTree *treeDY= (TTree*)fileDY->Get("T1");
  TTree *treeWj= (TTree*)fileWj->Get("T1");

  // these used  for the channel-split estimates (not now)
  TH1F *topHEE = new TH1F("topHEE","",50,0,180);
  TH1F *topHMM = new TH1F("topHMM","",50,0,180);
  TH1F *topHEM = new TH1F("topHEM","",50,0,180);
  TH1F *topHME = new TH1F("topHME","",50,0,180);
  TH1F *btagHDataEE = new TH1F("btagHDataEE","",50,0,180);
  TH1F *btagHDataMM = new TH1F("btagHDataMM","",50,0,180);
  TH1F *btagHDataEM = new TH1F("btagHDataEM","",50,0,180);
  TH1F *btagHDataME = new TH1F("btagHDataME","",50,0,180);

  // backgrounds in the tagged region
  TH1F *btagWWHLL = new TH1F("btagWWHLL","",50,0,180);
  TH1F *btagDYHLL = new TH1F("btagDYHLL","",50,0,180);
  TH1F *CJVWjHEE = new TH1F("CJVWjHEE","",50,0,180);
  TH1F *CJVWjHMM = new TH1F("CJVWjHMM","",50,0,180);
  TH1F *CJVWjHEM = new TH1F("CJVWjHEM","",50,0,180);
  TH1F *CJVWjHME = new TH1F("CJVWjHME","",50,0,180);
  TH1F *btagWjHEE = new TH1F("btagWjHEE","",50,0,180);
  TH1F *btagWjHMM = new TH1F("btagWjHMM","",50,0,180);
  TH1F *btagWjHEM = new TH1F("btagWjHEM","",50,0,180);
  TH1F *btagWjHME = new TH1F("btagWjHME","",50,0,180);

  // WW background after CJV (used with the mistag rate)
  TH1F *btagWWAllHLL = new TH1F("btagWWAllHLL","",50,0,180);

  treeTop->Project("topHEE","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==1)*baseW*puW")).Data());
  treeTop->Project("topHMM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==0)*baseW*puW")).Data());
  treeTop->Project("topHEM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==2)*baseW*puW")).Data());
  treeTop->Project("topHME","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==3)*baseW*puW")).Data());

  if(njets==0) {
    treeData->Project("btagHDataEE","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && !bveto && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && !bveto && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && !bveto && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && !bveto && channel==3")).Data());
  } else {
    treeData->Project("btagHDataEE","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && channel==3")).Data());
  }

  treeWW->Project("btagWWHLL","dphill",(TString("(finalLeptons && met>20 && ((channel<2 && mpmet>40 && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15)) || (channel>=2 && mpmet>20)) && ")+TString(njcut)+TString( " && mll>12 && zveto && !bveto)*baseW*puW")).Data());
  treeDY->Project("btagDYHLL","dphill",(TString("(finalLeptons && met>20 && ((channel<2 && mpmet>40 && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15)) || (channel>=2 && mpmet>20)) && ")+TString(njcut)+TString( " && mll>12 && zveto && !bveto)*baseW*puW")).Data());

  treeWW->Project("btagWWAllHLL","dphill",(TString("(finalLeptons && met>20 && ((channel<2 && mpmet>40 && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15)) || (channel>=2 && mpmet>20)) && ")+TString(njcut)+TString( " && mll>12 && zveto)*baseW*puW")).Data());

  treeWj->Project("btagWjHEE","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && !bveto && channel==1")).Data());
  treeWj->Project("btagWjHMM","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && !bveto && channel==0")).Data());
  treeWj->Project("btagWjHEM","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && !bveto && channel==2")).Data());
  treeWj->Project("btagWjHME","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && !bveto && channel==3")).Data());

  treeWj->Project("CJVWjHEE","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && channel==1")).Data());
  treeWj->Project("CJVWjHMM","dphill",(TString("finalLeptons && met>20 && mpmet>40 && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && mll>12 && zveto && channel==0")).Data());
  treeWj->Project("CJVWjHEM","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && channel==2")).Data());
  treeWj->Project("CJVWjHME","dphill",(TString("finalLeptons && met>20 && mpmet>20 && ")+TString(njcut)+TString( " && mll>12 && channel==3")).Data());

  // backgrounds in the tagged region (inclusive)
  float effBtagWj[4] = { (CJVWjHEE->Integral()>0) ? btagWjHEE->Integral()/CJVWjHEE->Integral() : 0.0, 
                         (CJVWjHMM->Integral()>0) ? btagWjHMM->Integral()/CJVWjHMM->Integral() : 0.0,
                         (CJVWjHEM->Integral()>0) ? btagWjHEM->Integral()/CJVWjHEM->Integral() : 0.0,
                         (CJVWjHME->Integral()>0) ? btagWjHME->Integral()/CJVWjHME->Integral() : 0.0 };

  float WW_tot = wantedLumi * btagWWHLL->Integral() * WWDataOverMC[njets];
  float WW_tot_err = wantedLumi * yieldErrPoisson(WW_tot,btagWWHLL->Integral());
  float DY_tot = wantedLumi * btagDYHLL->Integral() * DYDataOverMC[njets];
  float DY_tot_err = wantedLumi * yieldErrPoisson(DY_tot,btagDYHLL->Integral());
  float Wjets_tot = WjDataTot[ee][njets] * effBtagWj[ee] + WjDataTot[mm][njets] * effBtagWj[mm] +
    WjDataTot[em][njets] * effBtagWj[em] + WjDataTot[me][njets] * effBtagWj[me];
  float Wjets_tot_err = 0.40 * Wjets_tot; // approximation

  // as x-check use the WW(all) x mistag rate from data
  float WW_befTag = wantedLumi * btagWWAllHLL->Integral() * WWDataOverMC[njets];
  float WW_befTag_err =  wantedLumi * yieldErrPoisson(WW_befTag,btagWWAllHLL->Integral());
  float mistagR_MC = WW_tot / WW_befTag;
  float mistagR_SF = 1.1; // from BTV-11-001
  float mistagR_SF_err = 0.11;
  float WW_tot_2 = WW_befTag * mistagR_MC * mistagR_SF;
  float WW_tot_2_err = WW_tot_2 * quadrSum(WW_befTag_err/WW_befTag,mistagR_SF_err/mistagR_SF);

  float tagBkg_tot = WW_tot_2 + DY_tot + Wjets_tot;
  float tagBkg_tot_err = quadrSum(WW_tot_2_err,DY_tot_err,Wjets_tot_err);

  std::cout << "--- Background estimations: ---" << std::endl;
  std::cout << "WW = " << WW_tot << " +/-" << WW_tot_err << std::endl;
  std::cout << "WW mistag rate from sim * SF = " << WW_tot / WW_befTag * mistagR_SF << std::endl;
  std::cout << "WW (with data mistag) = " << WW_tot_2 << " +/- " << WW_tot_2_err << std::endl;
  std::cout << "DY = " << DY_tot << " +/-" << DY_tot_err << std::endl;
  std::cout << "Wjets = " << Wjets_tot << " +/-" << Wjets_tot_err << std::endl;
  std::cout << "Tot background to tagged events = " << tagBkg_tot << " +/- " << tagBkg_tot_err << std::endl;
  std::cout << "---- end backgrounds ---\n\n" << std::endl;

  ///// TOP ESTIMATION ///////
  float nTopMC[4];
  float nTopMC_err[4];
  
  nTopMC[ee] = wantedLumi * topHEE->Integral();
  nTopMC_err[ee] = wantedLumi * yieldErrPoisson(nTopMC[ee],topHEE->GetEntries());
  nTopMC[mm] = wantedLumi * topHMM->Integral();
  nTopMC_err[mm] = wantedLumi * yieldErrPoisson(nTopMC[mm],topHMM->GetEntries());
  nTopMC[em] = wantedLumi * topHEM->Integral();
  nTopMC_err[em] = wantedLumi * yieldErrPoisson(nTopMC[em],topHEM->GetEntries());
  nTopMC[me] = wantedLumi * topHME->Integral();
  nTopMC_err[me] = wantedLumi * yieldErrPoisson(nTopMC[me],topHME->GetEntries());
  
  float nTopMC_tot = nTopMC[ee] + nTopMC[mm] + nTopMC[em] + nTopMC[me];
  float nTopMC_tot_err = quadrSum(nTopMC_err[ee],nTopMC_err[mm],nTopMC_err[em],nTopMC_err[me]);

  /// TOP CHANNEL FRACTIONS FROM MC ///
  // they are used in case of low stat to us.e the summed WW estimation in the 4 sub-channels
  std::cout << "Using channel fractions from MC!" << std::endl;
  float frac[4];
  for(int icha=0; icha<4; icha++) {
    frac[icha] = nTopMC[icha]/nTopMC_tot;
    std::cout << "Fraction of Top in channel " << icha << " = " << frac[icha] << std::endl;
  }

  // top estimation from data (0-jet bin method)
  float nTopData[4];
  float nTopData_err[4];
  //  float eff_2b_softmu = 1 - pow(1-eff_1b_softmu,2);

  // EE
  float nBTagTag_data = btagHDataEE->Integral();
  float nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).first;
  float nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).second; 
  float nTopSoftMuVeto_data = nTopBTagVeto_data; //* (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's). Now included in eff_2b
  float nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EE = " << nBTagTag_data << std::endl;

  nTopData[ee] = nTopSoftMuVeto_data;
  nTopData_err[ee] = nTopSoftMuVeto_data_err;

  // MM
  nBTagTag_data = btagHDataMM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; //  * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "MM = " << nBTagTag_data << std::endl;

  nTopData[mm] = nTopSoftMuVeto_data;
  nTopData_err[mm] = nTopSoftMuVeto_data_err;

  // EM
  nBTagTag_data = btagHDataEM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; // * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EM = " << nBTagTag_data << std::endl;

  nTopData[em] = nTopSoftMuVeto_data;
  nTopData_err[em] = nTopSoftMuVeto_data_err;

  // ME
  nBTagTag_data = btagHDataME->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; // * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "ME = " << nBTagTag_data << std::endl;

  float nBTagTag_data_tot =  btagHDataEE->Integral() +  btagHDataMM->Integral() +  btagHDataEM->Integral() +  btagHDataME->Integral();
  std::cout << "Tagged events in data: " << "TOT = " << nBTagTag_data_tot << std::endl;

  nTopData[me] = nTopSoftMuVeto_data;
  nTopData_err[me] = nTopSoftMuVeto_data_err;


  // here summing up the separate ones
  //   float nTopData_tot = nTopData[ee] + nTopData[mm] + nTopData[em] + nTopData[me];
  //   float nTopData_tot_err = quadrSum(nTopData_err[ee],nTopData_err[mm],nTopData_err[em],nTopData_err[me]);

  // LL
  nBTagTag_data = nBTagTag_data_tot - tagBkg_tot;
  std::cout << "number of tagged events after bkg subtraction = " << nBTagTag_data << std::endl;
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b[njets], eff_2b_err[njets])).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; // * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  float nTopData_tot = nTopSoftMuVeto_data;
  float nTopData_tot_err = nTopSoftMuVeto_data_err;

  std::cout << "Number of Top events from data at W+W- level: " << std::endl
    // commented because the bkg estimation is done only for the cumulative sample
//             << "*\tEE = " << nTopData[ee] << " +/- " << nTopData_err[ee] << std::endl
//             << "*\tMM = " << nTopData[mm] << " +/- " << nTopData_err[mm] << std::endl
//             << "*\tEM = " << nTopData[em] << " +/- " << nTopData_err[em] << std::endl
//             << "*\tME = " << nTopData[me] << " +/- " << nTopData_err[me] << std::endl
            << "*\tTOT = " << wantedLumi/usedLumi * nTopData_tot << " +/- " << nTopData_tot_err << std::endl;

  std::cout << "Number of Top events from MC at W+W- level: " << std::endl
//             << "*\tEE = " << nTopMC[ee] << " +/- " << nTopMC_err[ee] << std::endl
//             << "*\tMM = " << nTopMC[mm] << " +/- " << nTopMC_err[mm] << std::endl
//             << "*\tEM = " << nTopMC[em] << " +/- " << nTopMC_err[em] << std::endl
//             << "*\tME = " << nTopMC[me] << " +/- " << nTopMC_err[me] << std::endl
            << "*\tTOT = " << wantedLumi/usedLumi * nTopMC_tot << " +/- " << nTopMC_tot_err << std::endl;

  ofstream textfile;
  textfile.open("TopYieldsData.txt", ios_base::trunc);
  textfile.precision(2);

  ofstream tablefile1;
  tablefile1.open("TopYieldsData_ForTable_0j.txt", ios_base::trunc);
  tablefile1.precision(2);

  ofstream tablefile2;
  tablefile2.open("TopYieldsData_ForTable_1j.txt", ios_base::trunc);
  tablefile2.precision(2);

  ofstream tablefile3;
  tablefile3.open("TopYieldsMC_ForTable_0j.txt", ios_base::trunc);
  tablefile3.precision(2);

  ofstream tablefile4;
  tablefile4.open("TopYieldsMC_ForTable_1j.txt", ios_base::trunc);
  tablefile4.precision(2);

  // these are for limits
  ofstream cardfile[4]; //[cha]
  for(int icha=0; icha<4; icha++) {
    char fileName[2];
    if(icha==ee) sprintf(fileName,"TopCard_ee_%dj.txt",njets);
    if(icha==mm) sprintf(fileName,"TopCard_mm_%dj.txt",njets);
    if(icha==em) sprintf(fileName,"TopCard_em_%dj.txt",njets);
    if(icha==me) sprintf(fileName,"TopCard_me_%dj.txt",njets);
    cardfile[icha].open(fileName, ios_base::trunc);
    cardfile[icha].precision(2);
  }
    

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    TString higgsMassDependentCut = higgsCuts(mass,true);

    float eff[4], eff_err[4];
    TString HCut[4];

    // calculate 0 jet mass dependent effciencies
    HCut[ee] = TString("(")+TString("WWSel")+TString(" && ")+higgsMassDependentCut+TString(" && channel==1)*baseW*puW*kfW");
    HCut[mm] = TString("(")+TString("WWSel")+TString(" && ")+higgsMassDependentCut+TString(" && channel==0)*baseW*puW*kfW");
    HCut[em] = TString("(")+TString("WWSel")+TString(" && ")+higgsMassDependentCut+TString(" && channel==2)*baseW*puW*kfW");
    HCut[me] = TString("(")+TString("WWSel")+TString(" && ")+higgsMassDependentCut+TString(" && channel==3)*baseW*puW*kfW");
    
    treeTop->Project("topHEE","dphill",HCut[ee]);
    treeTop->Project("topHMM","dphill",HCut[mm]);
    treeTop->Project("topHEM","dphill",HCut[em]);
    treeTop->Project("topHME","dphill",HCut[me]);
    
    std::vector<TH1F*> TopFin;
    TopFin.push_back(topHEE);
    TopFin.push_back(topHMM);
    TopFin.push_back(topHEM);
    TopFin.push_back(topHME);

    for(int icha=0; icha<4; icha++) {
      eff[icha] = TopFin[icha]->Integral() / nTopMC[icha];
      eff_err[icha] = eff[icha] * quadrSum(yieldErrPoisson(TopFin[icha]->Integral(),TopFin[icha]->GetEntries())/TopFin[icha]->Integral(),
                                           nTopMC_err[icha]/nTopMC[icha] );
    }

    float nTopData_HiggsSel[4], nTopData_HiggsSel_err[4];
    float nTopMC_HiggsSel[4], nTopMC_HiggsSel_err[4];

    for(int icha=0;icha<4;icha++) {
      float effErrRel = (eff[icha]==0) ? 0. : eff_err[icha]/eff[icha];
    
      // this is the correct esztimation for when we have sufficient stat
      // nTopData_HiggsSel_0j[icha] = nTopData[icha] * eff_0j;
      // float topDataErrRel = (nTopData[icha]==0) ? 0. : nTopData_err[icha]/nTopData[icha];
      // nTopData_HiggsSel_0j_err[icha] = nTopData_HiggsSel_0j[icha] * quadrSum(topDataErrRel,effErrRel);
      nTopData_HiggsSel[icha] = nTopData_tot * frac[icha] * eff[icha];
      nTopData_HiggsSel_err[icha] = nTopData_HiggsSel[icha] * quadrSum(nTopData_tot_err/nTopData_tot,eff_err[icha]/eff[icha]);

      nTopMC_HiggsSel[icha] = nTopMC[icha] * eff[icha];
      float topMCErrRel = (nTopMC[icha]==0) ? 0. : nTopMC_err[icha]/nTopMC[icha];
      nTopMC_HiggsSel_err[icha] = nTopMC_HiggsSel[icha] * quadrSum(topMCErrRel,effErrRel);

      char channelName[2];
      if(icha==ee) sprintf(channelName,"EE");
      if(icha==mm) sprintf(channelName,"MM");
      if(icha==em) sprintf(channelName,"EM");
      if(icha==me) sprintf(channelName,"ME");
      
      // for Giovanni
      float alpha = nTopData_HiggsSel[icha] / nBTagTag_data_tot;
      float alpha_err = alpha * nTopData_HiggsSel_err[icha] / nTopData_HiggsSel[icha];

      cardfile[icha] << mass 
                     << "\t" << nBTagTag_data_tot << "\t" << alpha
                     << "\t" <<  alpha_err 
                     << std::endl;

      ///////////////

      textfile << channelName << ": Higgs Mass = " << mass 
               << "\tdata jet = " << scaleFactorLumi * nTopData_HiggsSel[icha] << " +/- " << nTopData_HiggsSel_err[icha] 
               << "\tMC jet = " << scaleFactorLumi * nTopMC_HiggsSel[icha] << " +/- " << nTopMC_HiggsSel_err[icha] 
               << std::endl;
    }

    // summary table for limits
    if (i==0) { 
      tablefile1 << "# " << njets << " jets bin data" << endl;
      tablefile1 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile1 << mass 
	       << " " << "\t\t" << scaleFactorLumi * nTopData_HiggsSel[1] << " +/- " <<  nTopData_HiggsSel_err[1] 
	       << " " << "\t\t" << scaleFactorLumi * nTopData_HiggsSel[3] << " +/- " <<  nTopData_HiggsSel_err[3] 
	       << " " << "\t\t" << scaleFactorLumi * nTopData_HiggsSel[2] << " +/- " <<  nTopData_HiggsSel_err[2] 
	       << " " << "\t\t" << scaleFactorLumi * nTopData_HiggsSel[0] << " +/- " <<  nTopData_HiggsSel_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile3 << "# " << njets << " jets bin MC" << endl;
      tablefile3 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile3 << mass 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[1] << " +/- " <<  nTopMC_HiggsSel_err[1] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[3] << " +/- " <<  nTopMC_HiggsSel_err[3] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[2] << " +/- " <<  nTopMC_HiggsSel_err[2] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[0] << " +/- " <<  nTopMC_HiggsSel_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile4 << "# " << njets << " jet bin MC" << endl;
      tablefile4 << "#\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile4 << mass 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[1] << " +/- " <<  nTopMC_HiggsSel_err[1] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[3] << " +/- " <<  nTopMC_HiggsSel_err[3] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[2] << " +/- " <<  nTopMC_HiggsSel_err[2] 
	       << " " << "\t\t" << scaleFactorLumi * nTopMC_HiggsSel[0] << " +/- " <<  nTopMC_HiggsSel_err[0] 
	       << std::endl;

    float nTopData_HiggsSel_Tot = nTopData_HiggsSel[ee] + nTopData_HiggsSel[mm] + nTopData_HiggsSel[em] + nTopData_HiggsSel[me];
    float nTopData_HiggsSel_Tot_err = quadrSum(nTopData_HiggsSel_err[ee],nTopData_HiggsSel_err[mm],nTopData_HiggsSel_err[em],nTopData_HiggsSel_err[me]);

    float nTopMC_HiggsSel_Tot = nTopMC_HiggsSel[ee] + nTopMC_HiggsSel[mm] + nTopMC_HiggsSel[em] + nTopMC_HiggsSel[me];
    float nTopMC_HiggsSel_Tot_err = quadrSum(nTopMC_HiggsSel_err[ee],nTopMC_HiggsSel_err[mm],nTopMC_HiggsSel_err[em],nTopMC_HiggsSel_err[me]);

    textfile.precision(2);
    textfile << "\t===>> TOTAL: Higgs Mass = " << mass 
             << "\tdata = " << scaleFactorLumi * nTopData_HiggsSel_Tot << " +/- " << nTopData_HiggsSel_Tot_err 
             << "\tMC = " << scaleFactorLumi * nTopMC_HiggsSel_Tot << " +/- " << nTopMC_HiggsSel_Tot_err 
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

  // assume that the final selection efficiency is the same for all the top samples and use average of it
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.2.X/MC2011_Merged_V5/OptimMH%d/Spring11_V5HWW/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.2.X/MC2011_Merged_V5/OptimMH%d/Spring11_V5HWW/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.2.X/MC2011_Merged_V5/OptimMH%d/Spring11_V5HWW/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.2.X/MC2011_Merged_V5/OptimMH%d/Spring11_V5HWW/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root",mass);  
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
    nEv_end0j[theCha] += nSel[23];
    nEv_end1j[theCha] += nSel[24];
  }
}

