#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include "massDependentCuts.cc"

enum { ee=0, mm=1, em=2, me=3 };

float wantedLumi = 4.63;

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
std::pair<float,float> nWjetsTag(float nnotag, float nnotagerr, float r, float rerr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
std::pair<float,float> estimateTopVetoEff(int njet, bool fromData=true);
std::pair<float,float> estimateTopVetoEff2(int njet, float x, bool effFromData);
std::pair<float,float> estimateTopVetoEffBkgSub(int njets, bool effFromData=true);
void estimateTop(int njets);
void closureTest(int njets);

void printLatex() {
  std::cout << "=============================================================" << std::endl;
  std::cout << "@@@@ ESTIMATION OF THE TOP BACKGROUND IN THE ZERO-JET BIN @@@" << std::endl;
  std::cout << "=============================================================" << std::endl;
  estimateTop(0);
  std::cout << "DONE." << std::endl;
  std::cout << "=============================================================" << std::endl << std::endl << std::endl;

  std::cout << "=============================================================" << std::endl;
  std::cout << "@@@@ ESTIMATION OF THE TOP BACKGROUND IN THE ONE-JET BIN @@@" << std::endl;
  std::cout << "=============================================================" << std::endl;
  estimateTop(1);
  std::cout << "DONE." << std::endl;
  std::cout << "=============================================================" << std::endl;
}

void estimateTop(int njets) {

  // constants
  std::pair<float,float> efftag = estimateTopVetoEffBkgSub(njets);
  float eff_2b = efftag.first;
  float eff_2b_err = efftag.second;
  //  std::cout << "===> DISCLAIMER: TAKING THE EFF 2B PRECOMPUTED!" << std::endl; for 0-jet: Matt eff estimation  
  // float eff_2b = 0.49;
  // float eff_2b_err = 0.01;

  char njcut[30];
  sprintf(njcut, "njet==%d", njets);
  char wwselcut[30];
  if(njets==0) sprintf(wwselcut,"WWSel");
  else if(njets==1) sprintf(wwselcut,"WWSel1j");
  else {
    std::cout << "Jet bin must be 0/1" << std::endl;
    return;
  }

  // scale factors for the backgrounds
  float WWDataOverMC[2] = {1.00, 1.00} ; // use one for first iteration
  float DYDataOverMC[2] = {6.0, 7.2};  // estimation 4.63 fb-1
  float WjDataTot[4][2];     // [icha][jetbin]
  // cut based eleID, all LP (scenario1)
//   WjDataTot[ee][0] = 19.5; // updated 2.12 fb-1
//   WjDataTot[mm][0] = 19.1 ; // updated 2.12 fb-1
//   WjDataTot[em][0] = 87.1; // updated 2.12 fb-1
//   WjDataTot[me][0] = 64.0; // updated 2.12 fb-1
//   WjDataTot[ee][1] = 8.1; // updated 2.12 fb-1
//   WjDataTot[mm][1] = 10.4; // updated 2.12 fb-1
//   WjDataTot[em][1] = 32.7; // updated 2.12 fb-1
//   WjDataTot[me][1] = 21.4; // updated 2.12 fb-1

  WjDataTot[ee][0] = 11.4; // updated 4.63 fb-1
  WjDataTot[mm][0] = 8.2 ; // updated 4.63 fb-1
  WjDataTot[em][0] = 52.9; // updated 4.63 fb-1
  WjDataTot[me][0] = 22.9; // updated 4.63 fb-1
  WjDataTot[ee][1] = 5.1; // updated 4.63 fb-1
  WjDataTot[mm][1] = 8.3; // updated 4.63 fb-1
  WjDataTot[em][1] = 35.6; // updated 4.63 fb-1
  WjDataTot[me][1] = 18.3; // updated 4.63 fb-1
  
  TFile *fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  TFile *fileTop  = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWW   = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileDY   = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileOthers   = TFile::Open("results/datasets_trees/others_ll.root");

  TTree *treeData = (TTree*)fileData->Get("latino");
  TTree *treeTop  = (TTree*)fileTop->Get("latino");
  TTree *treeWW   = (TTree*)fileWW->Get("latino");
  TTree *treeDY   = (TTree*)fileDY->Get("latino");
  TTree *treeOthers = (TTree*)fileOthers->Get("latino");

  // these used  for the channel-split estimates (not now)
  TH1F *topHEE = new TH1F("topHEE","",50,0,2*TMath::Pi());
  TH1F *topHMM = new TH1F("topHMM","",50,0,2*TMath::Pi());
  TH1F *topHEM = new TH1F("topHEM","",50,0,2*TMath::Pi());
  TH1F *topHME = new TH1F("topHME","",50,0,2*TMath::Pi());
  TH1F *btagHDataEE = new TH1F("btagHDataEE","",50,0,2*TMath::Pi());
  TH1F *btagHDataMM = new TH1F("btagHDataMM","",50,0,2*TMath::Pi());
  TH1F *btagHDataEM = new TH1F("btagHDataEM","",50,0,2*TMath::Pi());
  TH1F *btagHDataME = new TH1F("btagHDataME","",50,0,2*TMath::Pi());

  // WW,DY background after CJV (used with the mistag rate)
  TH1F *WWAllHLL = new TH1F("WWAllHLL","",50,0,2*TMath::Pi());
  TH1F *DYAllHLL = new TH1F("DYAllHLL","",50,0,2*TMath::Pi());
  TH1F *OthersAllHLL = new TH1F("OthersAllHLL","",50,0,2*TMath::Pi());

  treeTop->Project("topHEE","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==1)*baseW*puW*kfW*effW")).Data());
  treeTop->Project("topHMM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==0)*baseW*puW*kfW*effW")).Data());
  treeTop->Project("topHEM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==2)*baseW*puW*kfW*effW")).Data());
  treeTop->Project("topHME","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==3)*baseW*puW*kfW*effW")).Data());

  TString btagLevelCut("(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav))");

  if(njets==0) {
    treeData->Project("btagHDataEE","dphill",(btagLevelCut + TString(" && abs(dphilljet)<165 && nextra==0 && ")+TString(njcut)+TString( " && !bveto && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(btagLevelCut + TString(" && abs(dphilljet)<165 && nextra==0 && ")+TString(njcut)+TString( " && !bveto && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(btagLevelCut + TString(" && nextra==0 && ")+TString(njcut)+TString( " && !bveto && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(btagLevelCut + TString(" && nextra==0 && ")+TString(njcut)+TString( " && !bveto && channel==3")).Data());
  } else {
    treeData->Project("btagHDataEE","dphill",(btagLevelCut + TString(" && abs(dphilljet)<165 && nextra==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(btagLevelCut + TString(" && abs(dphilljet)<165 && nextra==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(btagLevelCut + TString(" && nextra==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(btagLevelCut + TString(" && nextra==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==3")).Data());
  }

  treeWW->Project("WWAllHLL","dphill",(TString("(") + btagLevelCut + TString("  && (abs(dphilljet)<165 || channel>1) && nextra==0 && ")+TString(njcut)+TString( ") *baseW*puW*kfW*effW")).Data());
  treeDY->Project("DYAllHLL","dphill",(TString("(") + btagLevelCut + TString("  && (abs(dphilljet)<165 || channel>1) && nextra==0 && ")+TString(njcut)+TString( ") *baseW*puW*kfW*effW")).Data());
  treeOthers->Project("OthersAllHLL","dphill",(TString("(") + btagLevelCut + TString("  && (abs(dphilljet)<165 || channel>1) && nextra==0 && ")+TString(njcut)+TString( ") *baseW*puW*kfW*effW")).Data());

  // for Wjets WW and DY use the Z events in data
  float mistagSig[2];
  TH1F *ZDataHisto = new TH1F("ZDataHisto","",50,0,2*TMath::Pi());
  treeData->Project("ZDataHisto","dphill","finalLeptons && abs(mll-91.1876)<7.5 && pfmet<30 && njet==0");
  float mistagDenom = ZDataHisto->Integral();
  treeData->Project("ZDataHisto","dphill","finalLeptons && abs(mll-91.1876)<7.5 && pfmet<30 && njet==0 && !bveto");
  float mistagNum = ZDataHisto->Integral();
  mistagSig[0] = mistagNum/mistagDenom;

  treeData->Project("ZDataHisto","dphill","finalLeptons && abs(mll-91.1876)<7.5 && pfmet<30 && njet==1");
  mistagDenom = ZDataHisto->Integral();
  treeData->Project("ZDataHisto","dphill","finalLeptons && abs(mll-91.1876)<7.5 && pfmet<30 && njet==1 && !bveto");
  mistagNum = ZDataHisto->Integral();
  mistagSig[1] = mistagNum/mistagDenom;

  std::cout << "Measured mistag rate from DY events in " << njets << " jet bin = " << mistagSig[njets] << std::endl;

  float WW_tot = wantedLumi * WWAllHLL->Integral() * WWDataOverMC[njets] * mistagSig[njets];
  float WW_tot_err = wantedLumi * yieldErrPoisson(WW_tot,WWAllHLL->GetEntries()) * mistagSig[njets];
  float DY_tot = wantedLumi * DYAllHLL->Integral() * DYDataOverMC[njets] * mistagSig[njets];
  float DY_tot_err = wantedLumi * yieldErrPoisson(DY_tot,DYAllHLL->GetEntries()) * mistagSig[njets];
  float Others_tot = wantedLumi * OthersAllHLL->Integral() * mistagSig[njets];
  float Others_tot_err = wantedLumi * yieldErrPoisson(Others_tot,OthersAllHLL->GetEntries()) * mistagSig[njets];

  float Wjets_tot_notag = WjDataTot[ee][njets] + WjDataTot[mm][njets] + WjDataTot[em][njets] + WjDataTot[me][njets];
  // approximation: assuming 40% uncertainty on W+jets and 10% on mistag rate
  float Wjets_tot = (nWjetsTag(Wjets_tot_notag, 0.40 * Wjets_tot_notag, mistagSig[njets], 0.10 * mistagSig[njets])).first;
  float Wjets_tot_err = (nWjetsTag(Wjets_tot_notag, 0.40 * Wjets_tot_notag, mistagSig[njets], 0.10 * mistagSig[njets])).second;

  float tagBkg_tot = WW_tot + DY_tot + Wjets_tot + Others_tot;
  float tagBkg_tot_err = quadrSum(WW_tot_err,DY_tot_err,Wjets_tot_err,Others_tot_err);

  std::cout << "--- Background estimations: ---" << std::endl;
  std::cout << "WW = " << WW_tot << " +/-" << WW_tot_err << std::endl;
  std::cout << "DY = " << DY_tot << " +/-" << DY_tot_err << std::endl;
  std::cout << "Wjets = " << Wjets_tot << " +/-" << Wjets_tot_err << std::endl;
  std::cout << "Others = " << Others_tot << " +/-" << Others_tot_err << std::endl;
  std::cout << "Tot background to tagged events = " << tagBkg_tot << " +/- " << tagBkg_tot_err << std::endl;
  std::cout << "---- end backgrounds ---\n\n" << std::endl;

  ///// TOP ESTIMATION ///////
  float nTopMC[4];
  float nTopMC_err[4];
  float nTopEntries[4];

  nTopMC[ee] = wantedLumi * topHEE->Integral();
  nTopMC_err[ee] = wantedLumi * yieldErrPoisson(nTopMC[ee],topHEE->GetEntries());
  nTopEntries[ee] = topHEE->GetEntries();
  nTopMC[mm] = wantedLumi * topHMM->Integral();
  nTopMC_err[mm] = wantedLumi * yieldErrPoisson(nTopMC[mm],topHMM->GetEntries());
  nTopEntries[mm] = topHMM->GetEntries();
  nTopMC[em] = wantedLumi * topHEM->Integral();
  nTopMC_err[em] = wantedLumi * yieldErrPoisson(nTopMC[em],topHEM->GetEntries());
  nTopEntries[em] = topHEM->GetEntries();
  nTopMC[me] = wantedLumi * topHME->Integral();
  nTopMC_err[me] = wantedLumi * yieldErrPoisson(nTopMC[me],topHME->GetEntries());
  nTopEntries[me] = topHME->GetEntries();

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
  float nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  float nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  float nTopSoftMuVeto_data = nTopBTagVeto_data; //* (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's). Now included in eff_2b
  float nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EE = " << nBTagTag_data << std::endl;

  nTopData[ee] = nTopSoftMuVeto_data;
  nTopData_err[ee] = nTopSoftMuVeto_data_err;

  // MM
  nBTagTag_data = btagHDataMM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; //  * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "MM = " << nBTagTag_data << std::endl;

  nTopData[mm] = nTopSoftMuVeto_data;
  nTopData_err[mm] = nTopSoftMuVeto_data_err;

  // EM
  nBTagTag_data = btagHDataEM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data; // * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err; // * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EM = " << nBTagTag_data << std::endl;

  nTopData[em] = nTopSoftMuVeto_data;
  nTopData_err[em] = nTopSoftMuVeto_data_err;

  // ME
  nBTagTag_data = btagHDataME->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
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
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
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
            << "*\tTOT = " << nTopData_tot << " +/- " << nTopData_tot_err << std::endl;

  std::cout << "Number of Top events from MC at W+W- level: " << std::endl
//             << "*\tEE = " << nTopMC[ee] << " +/- " << nTopMC_err[ee] << std::endl
//             << "*\tMM = " << nTopMC[mm] << " +/- " << nTopMC_err[mm] << std::endl
//             << "*\tEM = " << nTopMC[em] << " +/- " << nTopMC_err[em] << std::endl
//             << "*\tME = " << nTopMC[me] << " +/- " << nTopMC_err[me] << std::endl
            << "*\tTOT = " << nTopMC_tot << " +/- " << nTopMC_tot_err << std::endl;

  ofstream tablefileData;
  char nameFile[100];
  sprintf(nameFile,"TopYieldsData_ForTable_%dj.txt",njets);
  tablefileData.open(nameFile, ios_base::trunc);
  tablefileData.setf(ios::fixed,ios::floatfield);
  tablefileData.precision(2);

  tablefileData << "\\begin{table}[p]" << endl;
  tablefileData << "\\begin{small}" << endl;
  tablefileData << "\\begin{center}" << endl;
  tablefileData << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  tablefileData << "\\hline" << endl;
  tablefileData << "$m_{H}$ [GeV] \t & $\\mu\\mu$ \t & \t $\\mu$ e \t & \t e $\\mu$ \t & \t ee \t & \t $\\ell\\ell$ \\\\" << endl;
  tablefileData << "\\hline" << endl;

  ofstream tablefileMC;
  sprintf(nameFile,"TopYieldsMC_ForTable_%dj.txt",njets);
  tablefileMC.open(nameFile, ios_base::trunc);
  tablefileMC.setf(ios::fixed,ios::floatfield);
  tablefileMC.precision(2);

  tablefileMC << "\\begin{table}[p]" << endl;
  tablefileMC << "\\begin{small}" << endl;
  tablefileMC << "\\begin{center}" << endl;
  tablefileMC << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  tablefileMC << "\\hline" << endl;
  tablefileMC << "$m_{H}$ [GeV] \t & $\\mu\\mu$ \t & \t $\\mu$ e \t & \t e $\\mu$ \t & \t ee \t & \t $\\ell\\ell$ \\\\" << endl;
  tablefileMC << "\\hline" << endl;

  // these are for limits
  ofstream cardfile[4]; //[cha]
  for(int icha=0; icha<4; icha++) {
    char fileName[2];
    if(icha==ee) sprintf(fileName,"TopCard_ee_%dj.txt",njets);
    if(icha==mm) sprintf(fileName,"TopCard_mm_%dj.txt",njets);
    if(icha==em) sprintf(fileName,"TopCard_em_%dj.txt",njets);
    if(icha==me) sprintf(fileName,"TopCard_me_%dj.txt",njets);
    cardfile[icha].open(fileName, ios_base::trunc);
    cardfile[icha].setf(ios::fixed,ios::floatfield);
    cardfile[icha].precision(5);
  }
    

  int masses[19] = {110,115,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<19; i++) {
    
    int mass = masses[i];

    TString higgsMassDependentCut = higgsCuts(mass,true);

    float eff[4], eff_err[4];
    TString HCut[4];

    // calculate 0/1 jet mass dependent effciencies
    HCut[ee] = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && channel==1)*baseW*puW*effW*kfW");
    HCut[mm] = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && channel==0)*baseW*puW*effW*kfW");
    HCut[em] = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && channel==2)*baseW*puW*effW*kfW");
    HCut[me] = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && channel==3)*baseW*puW*effW*kfW");
    
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
      eff[icha] = wantedLumi * TopFin[icha]->Integral() / nTopMC[icha];
      eff_err[icha] = sqrt(eff[icha] * (1.0 - eff[icha]) / nTopEntries[icha]);
//       eff_err[icha] = eff[icha] * quadrSum(yieldErrPoisson(wantedLumi * TopFin[icha]->Integral(),
//                                                            TopFin[icha]->GetEntries())/(wantedLumi * TopFin[icha]->Integral()),
//                                            nTopMC_err[icha]/nTopMC[icha] );
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

    }

    // summary table for limits
    float nTopData_HiggsSel_tot = 0;
    float nTopData_HiggsSel_tot_err = 0;
    float nTopMC_HiggsSel_tot = 0;
    float nTopMC_HiggsSel_tot_err = 0;
    for(int icha=0; icha<4; icha++) {
      nTopData_HiggsSel_tot += nTopData_HiggsSel[icha];
      nTopData_HiggsSel_tot_err += pow(nTopData_HiggsSel_err[icha],2);

      nTopMC_HiggsSel_tot += nTopMC_HiggsSel[icha];
      nTopMC_HiggsSel_tot_err += pow(nTopMC_HiggsSel_err[icha],2);
    }
    nTopData_HiggsSel_tot_err = sqrt(nTopData_HiggsSel_tot_err);
    nTopMC_HiggsSel_tot_err = sqrt(nTopMC_HiggsSel_tot_err);

    tablefileData << mass 
                  << "& \t\t" << nTopData_HiggsSel[1] << " $\\pm$ " <<  nTopData_HiggsSel_err[1] << "\t & \t"
                  << "\t\t" << nTopData_HiggsSel[3] << " $\\pm$ " <<  nTopData_HiggsSel_err[3] << "\t & \t"
                  << "\t\t" << nTopData_HiggsSel[2] << " $\\pm$ " <<  nTopData_HiggsSel_err[2] << "\t & \t"
                  << "\t\t" << nTopData_HiggsSel[0] << " $\\pm$ " <<  nTopData_HiggsSel_err[0] << "\t & \t"
                  << "\t\t" << nTopData_HiggsSel_tot << " $\\pm$ " << nTopData_HiggsSel_tot_err << "\t \\\\"
                  << std::endl;
    
    tablefileMC << mass 
                << "& \t\t" << nTopMC_HiggsSel[1] << " $\\pm$ " <<  nTopMC_HiggsSel_err[1] << "\t & \t"
                << "\t\t" << nTopMC_HiggsSel[3] << " $\\pm$ " <<  nTopMC_HiggsSel_err[3] << "\t & \t"
                << "\t\t" << nTopMC_HiggsSel[2] << " $\\pm$ " <<  nTopMC_HiggsSel_err[2] << "\t & \t"
                << "\t\t" << nTopMC_HiggsSel[0] << " $\\pm$ " <<  nTopMC_HiggsSel_err[0] << "\t & \t"
                << "\t\t" << nTopMC_HiggsSel_tot << " $\\pm$ " << nTopMC_HiggsSel_tot_err << "\t \\\\"
                << std::endl;
  }

  tablefileData << "\\hline" << endl;
  tablefileData << "\\end{tabular}" << endl;
  tablefileData << "\\caption{Top yields at Higgs selection level from data control sample in the " << njets << " jet bin.}" << std::endl;
  tablefileData << "\\end{center}" << endl;
  tablefileData << "\\end{small}" << endl;
  tablefileData << "\\end{table}" << endl;

  tablefileMC << "\\hline" << endl;
  tablefileMC << "\\end{tabular}" << endl;
  tablefileMC << "\\caption{Top yields at Higgs selection level as estimated in simulation in the " << njets << " jet bin.}" << std::endl;
  tablefileMC << "\\end{center}" << endl;
  tablefileMC << "\\end{small}" << endl;
  tablefileMC << "\\end{table}" << endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr) {
  float val = ntag * (1-eff2b) / eff2b;
  float err = ntag * eff2berr / pow(eff2b,2);
  return std::make_pair(val,err);
}

std::pair<float,float> nWjetsTag(float nnotag, float nnotagerr, float r, float rerr) {

  // W jets number is anti-tagged. Extrapolate to the tagged region
  float val = nnotag * r/(1-r);
  float err = quadrSum(r/(1-r)*nnotagerr, nnotag/pow(r-1,2)*rerr);
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

std::pair<float,float> estimateTopVetoEff(int njets, bool effFromData) {

  TFile *fileData = 0;
  if(effFromData) fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  else fileData = TFile::Open("results/datasets_trees/top_ll.root");

  TTree *treeData = (TTree*)fileData->Get("latino");


  // 0-jet bin method
  if(njets==0) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    if(effFromData) treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1");
    else treeData->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1)*baseW*puW*effW");
    float Ncontrol = histo1->Integral();

    if(effFromData) treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0)");
    else treeData->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0))*baseW*puW*effW");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    // export to the 0-jet bin using the fraction of ttbar/single-t
    TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
    TTree *treeTop = (TTree*)fileTop->Get("latino");

    TH1F *histo2 = new TH1F("histo2","",50,0,2*TMath::Pi());

    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0)*baseW*puW*effW");
    float top_pretopveto = histo2->Integral();
    
    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0 && dataset==10)*baseW*puW*effW");
    float ttbar_pretopveto = histo2->Integral();

    float fttbar = ttbar_pretopveto / top_pretopveto;
    float fsinglet = 1.0 - fttbar;

    float fttbar_over_singlet_err = 0.17; // generator uncertainty
    float fttbar_err = 1.0/fttbar * fttbar_over_singlet_err;
    float fsinglet_err = 1.0/fsinglet * fttbar_over_singlet_err;

    // do the weighted average of the efficiency in ttbar and single-t
    float eff2b_tt = 1 - pow(1 - eff_softtoptag, 2);
    float eff2b_tt_err = 2 * (1-eff2b_tt) * eff_softtoptag_err;
    
    std::cout << "N^{control}_{top-tag} = " << Ncontrol_toptag << std::endl;
    std::cout << "N^{control} = " << Ncontrol << std::endl;
    std::cout << "eff^{soft}_{top-tag} = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;
    std::cout << "eff_2b^{soft}_{top-tag} = " << eff2b_tt << " +/- " << eff2b_tt_err << std::endl;
    std::cout << "f_{ttbar} = " << fttbar << " +/- " << fttbar_err << std::endl;
    std::cout << "f_{singlet} = " << fsinglet << " +/- " << fsinglet_err << std::endl;
    
    float eff_0j = fttbar * eff2b_tt + fsinglet * eff_softtoptag;
    float eff_0j_err = quadrSum( pow(fttbar_err,2)*pow(eff2b_tt,2) + pow(eff2b_tt_err,2)*pow(fttbar,2),
                                 pow(fsinglet_err,2)*pow(eff_softtoptag,2) + pow(eff_softtoptag_err,2)*pow(fsinglet,2) );
    
    std::cout << "===> BIN-0: eff_0j = " << eff_0j << " +/- " << eff_0j_err << std::endl;

    return std::make_pair(eff_0j,eff_0j_err);

  } else if(njets==1) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    std::cout << "N^{control 2 jets} = " << Ncontrol << std::endl;
    std::cout << "N^{leading jet tagged} = " << Ncontrol_toptag << std::endl;
    std::cout << "===> BIN-1: eff_1j = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;

    return std::make_pair(eff_softtoptag,eff_softtoptag_err); 
  }

  std::cout << "ERROR: njet must be 0 or 1" << std::endl;
  return std::make_pair(0,0);

}

std::pair<float,float> estimateTopVetoEff2(int njets, float x, bool effFromData) {

  TFile *fileData = 0;
  if(effFromData) fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  else fileData = TFile::Open("results/datasets_trees/top_ll.root");

  TTree *treeData = (TTree*)fileData->Get("latino");


  // 0-jet bin method
  if(njets==0) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0)");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    // export to the 0-jet bin using the fraction of ttbar/single-t
    TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
    TTree *treeTop = (TTree*)fileTop->Get("latino");

    TH1F *histo2 = new TH1F("histo2","",50,0,2*TMath::Pi());

    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0)*baseW*puW*effW");
    float top_pretopveto = histo2->Integral();
    
    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0 && dataset==10)*baseW*puW*effW");
    float ttbar_pretopveto = histo2->Integral();

    float fttbar = ttbar_pretopveto / top_pretopveto;
    float fsinglet = 1.0 - fttbar;

    float fttbar_over_singlet_err = 0.17; // generator uncertainty
    float fttbar_err = 1.0/fttbar * fttbar_over_singlet_err;
    float fsinglet_err = 1.0/fsinglet * fttbar_over_singlet_err;

    // do the weighted average of the efficiency in ttbar and single-t
    float eff2b_tt = 1 - pow(1 - eff_softtoptag, 2);
    float eff2b_tt_err = 2 * (1-eff2b_tt) * eff_softtoptag_err;
    
    std::cout << "N^{control}_{top-tag} = " << Ncontrol_toptag << std::endl;
    std::cout << "N^{control} = " << Ncontrol << std::endl;
    std::cout << "eff^{soft}_{top-tag} = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;
    std::cout << "eff_2b^{soft}_{top-tag} = " << eff2b_tt << " +/- " << eff2b_tt_err << std::endl;
    std::cout << "f_{ttbar} = " << fttbar << " +/- " << fttbar_err << std::endl;
    std::cout << "f_{singlet} = " << fsinglet << " +/- " << fsinglet_err << std::endl;
    
    float eff_0j = (fttbar+(fsinglet)*x) * eff2b_tt + fsinglet * (1-x) * eff_softtoptag;
    float eff_0j_err = quadrSum( pow(fttbar_err,2)*pow(eff2b_tt,2) + pow(eff2b_tt_err,2)*pow(fttbar,2),
                                 pow(fsinglet_err,2)*pow(eff_softtoptag,2) + pow(eff_softtoptag_err,2)*pow(fsinglet,2) );
    
    std::cout << "===> BIN-0: eff_0j = " << eff_0j << " +/- " << eff_0j_err << std::endl;

    return std::make_pair(eff_0j,eff_0j_err);

  } else if(njets==1) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    std::cout << "N^{control 2 jets} = " << Ncontrol << std::endl;
    std::cout << "N^{leading jet tagged} = " << Ncontrol_toptag << std::endl;
    std::cout << "===> BIN-1: eff_1j = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;

    return std::make_pair(eff_softtoptag,eff_softtoptag_err); 
  }

  std::cout << "ERROR: njet must be 0 or 1" << std::endl;
  return std::make_pair(0,0);

}

std::pair<float,float> estimateTopVetoEffBkgSub(int njets, bool effFromData) {

  TFile *fileData = 0;
  if(effFromData) fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  else fileData = TFile::Open("results/datasets_trees/top_ll.root");

  TTree *treeData = (TTree*)fileData->Get("latino");


  // 0-jet bin method
  if(njets==0) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    if(effFromData) treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1");
    else treeData->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1)*baseW*puW*effW");
    float Ncontrol_all = histo1->Integral();

    if(effFromData) treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0)");
    else treeData->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0))*baseW*puW*effW");
    float Ncontrol_toptag_all = histo1->Integral();

    // export to the 0-jet bin using the fraction of ttbar/single-t
    TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
    TTree *treeTop = (TTree*)fileTop->Get("latino");

    TH1F *histo2 = new TH1F("histo2","",50,0,2*TMath::Pi());

    // subtract the background to the numerator and denominator
    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && dataset!=10)*baseW*puW*effW");
    float Ncontrol_bkg = histo2->Integral();
    treeTop->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0) && dataset!=10)*baseW*puW*effW");
    float Ncontrol_toptag_bkg = histo2->Integral();

    // if data, subtract the other backgrounds also
    if(effFromData) {
      TFile *fileWW = TFile::Open("results/datasets_trees/WW_ll.root");
      TTree *treeWW = (TTree*)fileWW->Get("latino");
      TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");
      TTree *treeZjets = (TTree*)fileZjets->Get("latino");

      treeWW->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && dataset!=10)*baseW*puW*effW");
      float Ncontrol_WW = histo2->Integral();
      treeWW->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0) && dataset!=10)*baseW*puW*effW");
      float Ncontrol_toptag_WW = histo2->Integral();

      treeZjets->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && dataset!=10)*baseW*puW*effW");
      float Ncontrol_DY = histo2->Integral();
      treeZjets->Project("histo2","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0) && dataset!=10)*baseW*puW*effW");
      float Ncontrol_toptag_DY = histo2->Integral();

      Ncontrol_bkg += (Ncontrol_WW + Ncontrol_DY);
      Ncontrol_toptag_bkg += (Ncontrol_toptag_WW + Ncontrol_toptag_DY);
    }

    // subtract the background
    float Ncontrol = Ncontrol_all - Ncontrol_bkg;
    float Ncontrol_toptag = Ncontrol_toptag_all - Ncontrol_toptag_bkg;

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    TH1F *histo3 = new TH1F("histo3","",50,0,2*TMath::Pi());
    treeTop->Project("histo3","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0)*baseW*puW*effW");
    float top_pretopveto = histo3->Integral();
    
    treeTop->Project("histo3","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0 && dataset==10)*baseW*puW*effW");
    float ttbar_pretopveto = histo3->Integral();

    float fttbar = ttbar_pretopveto / top_pretopveto;
    float fsinglet = 1.0 - fttbar;

    float fttbar_over_singlet_err = 0.17; // generator uncertainty
    float fttbar_err = 1.0/fttbar * fttbar_over_singlet_err;
    float fsinglet_err = 1.0/fsinglet * fttbar_over_singlet_err;

    // do the weighted average of the efficiency in ttbar and single-t
    float eff2b_tt = 1 - pow(1 - eff_softtoptag, 2);
    float eff2b_tt_err = 2 * (1-eff2b_tt) * eff_softtoptag_err;
    
    std::cout << "N^{control}_{top-tag} (raw) = " << Ncontrol_toptag_all << std::endl;
    std::cout << "N^{control} (raw) = " << Ncontrol_all << std::endl;
    std::cout << "N^{control}_{top-tag} (subtr.) = " << Ncontrol_toptag << std::endl;
    std::cout << "N^{control} = " << Ncontrol << std::endl;
    std::cout << "N^{control}_{top-tag}^{bkg} = " << Ncontrol_toptag_bkg << std::endl;
    std::cout << "N^{control}^{bkg} = " << Ncontrol_bkg << std::endl;
    std::cout << "eff^{soft}_{top-tag} = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;
    std::cout << "eff_2b^{soft}_{top-tag} = " << eff2b_tt << " +/- " << eff2b_tt_err << std::endl;
    std::cout << "f_{ttbar} = " << fttbar << " +/- " << fttbar_err << std::endl;
    std::cout << "f_{singlet} = " << fsinglet << " +/- " << fsinglet_err << std::endl;
    
    float eff_0j = fttbar * eff2b_tt + fsinglet * eff_softtoptag;
    float eff_0j_err = quadrSum( pow(fttbar_err,2)*pow(eff2b_tt,2) + pow(eff2b_tt_err,2)*pow(fttbar,2),
                                 pow(fsinglet_err,2)*pow(eff_softtoptag,2) + pow(eff_softtoptag_err,2)*pow(fsinglet,2) );
    
    std::cout << "===> BIN-0: eff_0j = " << eff_0j << " +/- " << eff_0j_err << std::endl;

    return std::make_pair(eff_0j,eff_0j_err);

  } else if(njets==1) {

    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==2 && subleadingJetBTagTrackCount>2.1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    std::cout << "N^{control 2 jets} = " << Ncontrol << std::endl;
    std::cout << "N^{leading jet tagged} = " << Ncontrol_toptag << std::endl;
    std::cout << "===> BIN-1: eff_1j = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;

    return std::make_pair(eff_softtoptag,eff_softtoptag_err); 
  }

  std::cout << "ERROR: njet must be 0 or 1" << std::endl;
  return std::make_pair(0,0);

}


void closureTest(int njets) {

  std::cout << "====> PERFORMING TOP CLOSURE TEST ON " << njets << " JET BIN" << std::endl;

  // export to the 0-jet bin using the fraction of ttbar/single-t
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
  TTree *treeTop = (TTree*)fileTop->Get("latino");
  
  if(njets==0) {
    // MC truth efficiency
    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());
    treeTop->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0 && bveto)*baseW*puW*effW");
    float num_eff2b_MC = histo1->Integral();
    treeTop->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==0)*baseW*puW*effW");
    float denom_eff2b_MC = histo1->Integral();
    float eff2b_MC = num_eff2b_MC / denom_eff2b_MC;
    float eff2b_MC_err = sqrt(eff2b_MC*(1-eff2b_MC)/histo1->GetEntries());

    // 0-jet bin method
    std::pair<float,float> eff2b = estimateTopVetoEffBkgSub(0,false);
    float eff2b_est = eff2b.first;
    float eff2b_est_err = eff2b.second;
    
    std::cout << "=== CLOSURE TEST OUTPUT:===" << std::endl;
    std::cout << "MC truth eff2b = " << eff2b_MC << " +/- " << eff2b_MC_err << std::endl;
    std::cout << "estimated eff2b = " << eff2b_est << " +/- " << eff2b_est_err << std::endl;
    std::cout << "closure (estimate-MC)/estimate = " << (eff2b_est-eff2b_MC)/eff2b_MC << " +/- " << eff2b_est/eff2b_MC * quadrSum(eff2b_MC_err/eff2b_MC, eff2b_est_err/eff2b_est) << std::endl;

  } else if(njets==1) {
    // MC truth efficiency
    TH1F *histo1 = new TH1F("histo1","",50,0,2*TMath::Pi());

    treeTop->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1)*baseW*puW*effW");
    float Ncontrol = histo1->Integral();

    treeTop->Project("histo1","dphill","(step[9] && ptll>45 && pt1>20 && ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) && (mll>20 || !sameflav) && (abs(dphilljet)<165 || channel>1) && nextra==0 && njet==1 && leadingJetBTagTrackCount>2.1)*baseW*puW*effW");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    // 1-jet bin method
    std::pair<float,float> effHighestPt = estimateTopVetoEff(1,false);
    float effHighestPt_est = effHighestPt.first;
    float effHighestPt_est_err = effHighestPt.second;
    
    std::cout << "=== CLOSURE TEST OUTPUT:===" << std::endl;
    std::cout << "MC truth highest pT jet tagging (1jet) = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;
    std::cout << "estimated highest pT jet tagging (1jet) = " << effHighestPt_est << " +/- " << effHighestPt_est_err << std::endl;
    std::cout << "closure (estimate-MC)/estimate = " << (effHighestPt_est - eff_softtoptag)/effHighestPt_est << " +/- " 
              << effHighestPt_est/eff_softtoptag * quadrSum(effHighestPt_est_err/effHighestPt_est, eff_softtoptag_err/eff_softtoptag) << std::endl;
    
  } else {
    std::cout << "ONLY 0/1 jet bin implemented. Doing nothing. Bye." << std::endl;
    return;
  }

}
