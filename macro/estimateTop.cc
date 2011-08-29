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

float wantedLumi = 1.54;

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
std::pair<float,float> estimateTopVetoEff(int njets);
void estimateTop(int njets);

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
  std::pair<float,float> efftag = estimateTopVetoEff(njets);
  float eff_2b = efftag.first;
  float eff_2b_err = efftag.second;

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
  float WWDataOverMC[2] = {1.02, 1.09} ; // estimation 1.1 fb-1 
  float DYDataOverMC[2] = {2.2, 2.4};  // estimation 1.54 fb-1
  float WjDataTot[4][2];     // [icha][jetbin]
  WjDataTot[ee][0] = 22.30; // updated LP
  WjDataTot[mm][0] = 13.5 ; // updated LP
  WjDataTot[em][0] = 63.4; // updated LP
  WjDataTot[me][0] = 50.4; // updated LP
  WjDataTot[ee][1] = 8.1; // updated LP
  WjDataTot[mm][1] = 8.3; // updated LP
  WjDataTot[em][1] = 22.9; // updated LP
  WjDataTot[me][1] = 13.8; // updated LP
  
  TFile *fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  TFile *fileTop  = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWW   = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileDY   = TFile::Open("results/datasets_trees/Zjets_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeTop  = (TTree*)fileTop->Get("T1");
  TTree *treeWW   = (TTree*)fileWW->Get("T1");
  TTree *treeDY   = (TTree*)fileDY->Get("T1");

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

  // WW background after CJV (used with the mistag rate)
  TH1F *btagWWAllHLL = new TH1F("btagWWAllHLL","",50,0,180);

  treeTop->Project("topHEE","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==1)*baseW*puW*effW")).Data());
  treeTop->Project("topHMM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==0)*baseW*puW*effW")).Data());
  treeTop->Project("topHEM","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==2)*baseW*puW*effW")).Data());
  treeTop->Project("topHME","dphill",(TString("(")+TString(wwselcut)+TString(" && channel==3)*baseW*puW*effW")).Data());

  if(njets==0) {
    treeData->Project("btagHDataEE","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto && channel==3")).Data());

    treeWW->Project("btagWWHLL","dphill",(TString("(step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto)*baseW*puW*effW")).Data());
    treeDY->Project("btagDYHLL","dphill",(TString("(step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && !bveto)*baseW*puW*effW")).Data());
    
  } else {
    treeData->Project("btagHDataEE","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==1")).Data());
    treeData->Project("btagHDataMM","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==0")).Data());
    treeData->Project("btagHDataEM","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==2")).Data());
    treeData->Project("btagHDataME","dphill",(TString("step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0 && channel==3")).Data());

    treeWW->Project("btagWWHLL","dphill",(TString("(step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0)*baseW*puW*effW")).Data());
    treeDY->Project("btagDYHLL","dphill",(TString("(step[9] && nExtraLep==0 && ")+TString(njcut)+TString( " && leadingJetBTagTrackCount>2.1 && subleadingJetsMaxBTagTrackCount<=2.1 && nSoftMu==0)*baseW*puW*effW")).Data());

  }

  treeWW->Project("btagWWAllHLL","dphill",(TString("(step[9] && nExtraLep==0 && ")+TString(njcut)+TString( ") *baseW*puW*effW")).Data());

  // backgrounds in the tagged region
  // for W+jets estimated on data
  float effBtagWj[2] = { 0.046, 0.016 };
  float WW_tot = wantedLumi * btagWWHLL->Integral() * WWDataOverMC[njets];
  float WW_tot_err = wantedLumi * yieldErrPoisson(WW_tot,btagWWHLL->Integral());
  float DY_tot = wantedLumi * btagDYHLL->Integral() * DYDataOverMC[njets];
  float DY_tot_err = wantedLumi * yieldErrPoisson(DY_tot,btagDYHLL->Integral());
  float Wjets_tot = effBtagWj[njets] * (WjDataTot[ee][njets] + WjDataTot[mm][njets] +
                                        WjDataTot[em][njets] + WjDataTot[me][njets]);
  float Wjets_tot_err = 0.40 * Wjets_tot; // approximation

  // as x-check use the WW(all) x mistag rate from data
  float WW_befTag = wantedLumi * btagWWAllHLL->Integral() * WWDataOverMC[njets];
  float WW_befTag_err =  wantedLumi * yieldErrPoisson(WW_befTag,btagWWAllHLL->Integral());
  float mistagR_MC = WW_tot / WW_befTag;
  float mistagR_SF = 1.1; // from BTV-11-001
  float mistagR_SF_err = 0.11;
  float WW_tot_2 = WW_befTag * mistagR_MC * mistagR_SF;
  float WW_tot_2_err = WW_tot_2 * quadrSum(WW_befTag_err/WW_befTag,mistagR_SF_err/mistagR_SF);

  float tagBkg_tot = WW_tot + DY_tot + Wjets_tot;
  float tagBkg_tot_err = quadrSum(WW_tot_err,DY_tot_err,Wjets_tot_err);

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
    

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
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

std::pair<float,float> estimateTopVetoEff(int njets) {

  TFile *fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  TTree *treeData = (TTree*)fileData->Get("T1");


  // 0-jet bin method
  if(njets==0) {

    TH1F *histo1 = new TH1F("histo1","",50,0,180);

    treeData->Project("histo1","dphill","step[9] && nExtraLep==0 && njets==1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && nExtraLep==0 && njets==1 && leadingJetBTagTrackCount>2.1 && (subleadingJetsMaxBTagTrackCount>2.1 || nSoftMuNoJets>0)");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    // export to the 0-jet bin using the fraction of ttbar/single-t
    TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
    TTree *treeTop = (TTree*)fileTop->Get("T1");

    TH1F *histo2 = new TH1F("histo2","",50,0,180);

    treeTop->Project("histo2","dphill","(step[9] && nExtraLep==0 && njets==0)*baseW*puW*effW");
    float top_pretopveto = histo2->Integral();
    
    treeTop->Project("histo2","dphill","(step[9] && nExtraLep==0 && njets==0 && process==11)*baseW*puW*effW");
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

    TH1F *histo1 = new TH1F("histo1","",50,0,180);

    treeData->Project("histo1","dphill","step[9] && nExtraLep==0 && njets==2 && subleadingJetBTagTrackCount>2.1");
    float Ncontrol = histo1->Integral();

    treeData->Project("histo1","dphill","step[9] && nExtraLep==0 && njets==2 && subleadingJetBTagTrackCount>2.1 && leadingJetBTagTrackCount>2.1");
    float Ncontrol_toptag = histo1->Integral();

    float eff_softtoptag = Ncontrol_toptag / Ncontrol;
    float eff_softtoptag_err = sqrt(eff_softtoptag * (1-eff_softtoptag)/Ncontrol);

    std::cout << "N^{control 2 jets} = " << Ncontrol << std::endl;
    std::cout << "N^{leading jet tagged} = " << Ncontrol_toptag << std::endl;
    std::cout << "===> BIN-1: eff_1j = " << eff_softtoptag << " +/- " << eff_softtoptag_err << std::endl;

    return std::make_pair(eff_softtoptag,eff_softtoptag_err); 
  }

  std::cout << "ERROR: njets must be 0 or 1" << std::endl;
  return std::make_pair(0,0);

}
