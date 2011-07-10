#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float Rmm, Rmm_err, Ree, Ree_err;

float wantedLumi = 1091;
float usedLumi = 1091;
float scaleFactorLumi = wantedLumi/usedLumi;

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass, const char* channel);
void estimateRWWRegion();

void estimateWW() {

  // constants: backgrounds estimated in all region.
  float Wjets0j[4] = { 10.6, 15.23 * wantedLumi/710., 25.4, 25.1 };
  float Wjets0j_err[4] = { 1.4, quadrSum(2.81,4.99), 1.7, 1.5 };

  float R1j0j_wj = 38./132.; // from MC
  float Wjets1j[4], Wjets1j_err[4];

  for(int icha=0; icha<4; icha++) {
    Wjets1j[icha] = Wjets0j[icha] * R1j0j_wj;
    Wjets1j_err[icha] = Wjets0j_err[icha] * R1j0j_wj;
  }

  float DY0j[4] = { 0.74, 4.6, 0.0, 0.0 };
  float DY0j_err[4] = { 1.4, 2.7, 0., 0. };

  float DY1j[4] = { 3.7, 8.4, 0., 0.};
  float DY1j_err[4] = { 1.3, 2.3, 0., 0.};

  float Wjets0j_tot(0), Wjets0j_tot_err(0), Wjets1j_tot(0), Wjets1j_tot_err(0);
  float DY0j_tot(0), DY0j_tot_err(0), DY1j_tot(0), DY1j_tot_err(0);
  for(int icha=0; icha<4; icha++) {
    Wjets0j_tot += Wjets0j[icha];
    Wjets0j_tot_err += pow(Wjets0j_err[icha],2);
    Wjets1j_tot += Wjets1j[icha];
    Wjets1j_tot_err += pow(Wjets1j_err[icha],2);

    DY0j_tot += DY0j[icha];
    DY0j_tot_err += pow(DY0j_err[icha],2);
    DY1j_tot += DY1j[icha];
    DY1j_tot_err += pow(DY1j_err[icha],2);
  }

  Wjets0j_tot_err = sqrt(Wjets0j_tot_err);
  Wjets1j_tot_err = sqrt(Wjets1j_tot_err);
  DY0j_tot_err = sqrt(DY0j_tot_err);
  DY1j_tot_err = sqrt(DY1j_tot_err);
    
  float top0j_tot = 67.6; 
  float top0j_tot_err = 16.3;

  float top1j_tot = 122;
  float top1j_tot_err = 13.0;

  float dibosons0j_tot = 13.8/1000. * usedLumi ;
  float dibosons0j_tot_err = 0.4/1000. * usedLumi;

  float dibosons1j_tot = 9.2/1000. * usedLumi ;
  float dibosons1j_tot_err = 0.35/1000. * usedLumi;
  
  TFile *fileData = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWjets = TFile::Open("results/datasets_trees/Wjets_ll.root");
  TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileDiBosons = TFile::Open("results/datasets_trees/others_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeWW = (TTree*)fileWW->Get("T1");
  TTree *treeTop= (TTree*)fileTop->Get("T1");
  TTree *treeWjets = (TTree*)fileWjets->Get("T1");
  TTree *treeZjets = (TTree*)fileZjets->Get("T1");
  TTree *treeDiBosons = (TTree*)fileDiBosons->Get("T1");

  TH1F *dataH = new TH1F("dataH","",50,0,180);
  TH1F *WWH = new TH1F("WWH","",50,0,180); // all
  TH1F *topH = new TH1F("topH","",50,0,180);
  TH1F *WjetsH = new TH1F("WjetsH","",50,0,180);
  TH1F *DYH = new TH1F("DYH","",50,0,180);
  TH1F *DiBosonsH = new TH1F("DiBosonsH","",50,0,180);

  std::vector<TH1F*> histosMC;
  float integralData0jTot;
  std::vector<float> integrals0jTot, error0jTot, ratio0j, ratio0j_err;
  histosMC.push_back(WWH);
  histosMC.push_back(topH);
  histosMC.push_back(WjetsH);
  histosMC.push_back(DYH);
  histosMC.push_back(DiBosonsH);

  std::vector<string> namesMC;
  namesMC.push_back("WW");
  namesMC.push_back("top");
  namesMC.push_back("wjets");
  namesMC.push_back("dy");
  namesMC.push_back("dibosons");

  enum { ww=0, top=1, wjets=2, dy=3, dib=4 };

  // --- 0 jet ---
  // estimate the values in the full region
  treeData->Project("dataH","deltaPhi","WWSel");
  treeWW->Project("WWH","deltaPhi","(WWSel)*weight*puweight");
  treeTop->Project("topH","deltaPhi","(WWSel)*weight*puweight");
  treeWjets->Project("WjetsH","deltaPhi","(WWSel)*weight*puweight");
  treeZjets->Project("DYH","deltaPhi","(WWSel)*weight*puweight");
  treeDiBosons->Project("DiBosonsH","deltaPhi","(WWSel)*weight*puweight");

  integralData0jTot = dataH->Integral();
  for(int imc=0; imc<(int)histosMC.size(); imc++) {
    integrals0jTot.push_back(histosMC[imc]->Integral());
    error0jTot.push_back(yieldErrPoisson(histosMC[imc]->Integral(), histosMC[imc]->GetEntries()));
  }

  // fix DY 0 jet to Matt's
  std::cout << "Taking DY 0j from external input" << std::endl;
  integrals0jTot[dy] = 2.299;
  error0jTot[dy] = 0.2; // put 10%

  // estimate the yields in the control region
  float integralData0jCtrl;
  std::vector<float> integrals0jCtrl, error0jCtrl;
  treeData->Project("dataH","deltaPhi","eleInvMass>100 && WWSel");
  treeWW->Project("WWH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeTop->Project("topH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeWjets->Project("WjetsH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeZjets->Project("DYH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeDiBosons->Project("DiBosonsH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");

  integralData0jCtrl = dataH->Integral();
  std::cout << "Ratios All/Control, 0jet" << std::endl; 
  for(int imc=0; imc<(int)histosMC.size(); imc++) {
    integrals0jCtrl.push_back(histosMC[imc]->Integral());
    error0jCtrl.push_back(yieldErrPoisson(histosMC[imc]->Integral(), histosMC[imc]->GetEntries()));
    ratio0j.push_back(integrals0jCtrl[imc]/integrals0jTot[imc]);
    ratio0j_err.push_back(ratio0j[imc] * quadrSum(error0jCtrl[imc]/integrals0jCtrl[imc], error0jTot[imc]/integrals0jTot[imc]));
    std::cout << namesMC[imc] << " (all,mll>100 GeV,fraction) = \t (" << integrals0jTot[imc] << ", " << integrals0jCtrl[imc] << ", " 
              << ratio0j[imc] << " +/- " << ratio0j_err[imc] << ")" << std::endl; 
  }
  
  // fix DY 0 jet to Matt's
  integrals0jCtrl[dy] = 0.445;
  error0jCtrl[dy] = 0.04; // put 10%
  ratio0j[dy] = integrals0jCtrl[dy] / integrals0jTot[dy];
  ratio0j_err[dy] = ratio0j[dy] * quadrSum(error0jCtrl[dy]/integrals0jCtrl[dy],error0jTot[dy]/integrals0jTot[dy]);

  // try to fix Wjets
  ratio0j[wjets] = 0.05;
  ratio0j_err[wjets] = 0.02;

  // --- 1 jet ---
  // estimate the values in the full region
  std::vector<float> integrals1jTot, error1jTot, ratio1j, ratio1j_err;
  treeData->Project("dataH","deltaPhi","WWSel1j");
  treeWW->Project("WWH","deltaPhi","(WWSel1j)*weight*puweight");
  treeTop->Project("topH","deltaPhi","(WWSel1j)*weight*puweight");
  treeWjets->Project("WjetsH","deltaPhi","(WWSel1j)*weight*puweight");
  treeZjets->Project("DYH","deltaPhi","(WWSel1j)*weight*puweight");
  treeDiBosons->Project("DiBosonsH","deltaPhi","(WWSel1j)*weight*puweight");

  float integralData1jTot = dataH->Integral();
  for(int imc=0; imc<(int)histosMC.size(); imc++) {
    integrals1jTot.push_back(histosMC[imc]->Integral());
    error1jTot.push_back(yieldErrPoisson(histosMC[imc]->Integral(), histosMC[imc]->GetEntries()));
  }

  // estimate the yields in the control region
  std::vector<float> integrals1jCtrl, error1jCtrl;
  treeData->Project("dataH","deltaPhi","eleInvMass>100 && WWSel1j");
  treeWW->Project("WWH","deltaPhi","(eleInvMass>100 && WWSel1j)*weight*puweight");
  treeTop->Project("topH","deltaPhi","(eleInvMass>100 && WWSel1j)*weight*puweight");
  treeWjets->Project("WjetsH","deltaPhi","(eleInvMass>100 && WWSel1j)*weight*puweight");
  treeZjets->Project("DYH","deltaPhi","(eleInvMass>100 && WWSel1j)*weight*puweight");
  treeDiBosons->Project("DiBosonsH","deltaPhi","(eleInvMass>100 && WWSel1j)*weight*puweight");

  float integralData1jCtrl = dataH->Integral();
  std::cout << "Ratios All/Control, 1 jet" << std::endl; 
  for(int imc=0; imc<(int)histosMC.size(); imc++) {
    integrals1jCtrl.push_back(histosMC[imc]->Integral());
    error1jCtrl.push_back(yieldErrPoisson(histosMC[imc]->Integral(), histosMC[imc]->GetEntries()));
    ratio1j.push_back(integrals1jCtrl[imc]/integrals1jTot[imc]);
    ratio1j_err.push_back(ratio1j[imc] * quadrSum(error1jCtrl[imc]/integrals1jCtrl[imc], error1jTot[imc]/integrals1jTot[imc]));
    std::cout << namesMC[imc] << " (all,mll>100 GeV,fraction) = \t (" << integrals1jTot[imc] << ", " << integrals1jCtrl[imc] << ", " 
              << ratio1j[imc] << " +/- " << ratio1j_err[imc] << ")" << std::endl; 
  }

  // try to fix Wjets
  ratio1j[wjets] = 0.05;
  ratio1j_err[wjets] = 0.02;

  /// EXTRAPOLATE THE BKGS TO THE WW CONTROL REGION
  // --- 0 jet ---
  float nWjets_0j_Out = Wjets0j_tot * ratio0j[wjets];
  float nWjets_0j_Out_err = nWjets_0j_Out * quadrSum(ratio0j_err[wjets]/ratio0j[wjets], Wjets0j_tot_err/Wjets0j_tot);
  
  std::cout << DY0j_tot << " " << ratio0j[dy] << std::endl;
  float nDY_0j_Out = DY0j_tot * ratio0j[dy];
  float nDY_0j_Out_err = nDY_0j_Out * quadrSum(ratio0j_err[dy]/ratio0j[dy], (DY0j_tot==0) ? 0 : DY0j_tot_err/DY0j_tot);

  float nTop_0j_Out = top0j_tot * ratio0j[top];
  float nTop_0j_Out_err = nTop_0j_Out * quadrSum(ratio0j_err[top]/ratio0j[top], top0j_tot_err/top0j_tot);

  float nDiBosons_0j_Out = dibosons0j_tot * ratio0j[dib];
  float nDiBosons_0j_Out_err = nDiBosons_0j_Out * quadrSum(ratio0j_err[dib]/ratio0j[dib], dibosons0j_tot_err/dibosons0j_tot);

  float nBkg_0j_Out = nWjets_0j_Out + nDY_0j_Out + nTop_0j_Out + nDiBosons_0j_Out;
  float nBkg_0j_Out_err = quadrSum(nWjets_0j_Out_err,nDY_0j_Out_err,nTop_0j_Out_err,nDiBosons_0j_Out_err);

  // --- 1 jet ---
  float nWjets_1j_Out = Wjets1j_tot * ratio1j[wjets];
  float nWjets_1j_Out_err = nWjets_1j_Out * quadrSum(ratio1j_err[wjets]/ratio1j[wjets], Wjets1j_tot_err/Wjets1j_tot);

  float nDY_1j_Out = DY1j_tot * ratio1j[dy];
  float nDY_1j_Out_err = nDY_1j_Out * quadrSum(ratio1j_err[dy]/ratio1j[dy], (DY0j_tot==0) ? 0 : DY0j_tot_err/DY0j_tot);

  float nTop_1j_Out = top1j_tot * ratio1j[top];
  float nTop_1j_Out_err = nTop_1j_Out * quadrSum(ratio1j_err[top]/ratio1j[top], top1j_tot_err/top1j_tot);

  float nDiBosons_1j_Out = dibosons1j_tot * ratio1j[dib];
  float nDiBosons_1j_Out_err = nDiBosons_1j_Out * quadrSum(ratio1j_err[dib]/ratio1j[dib], dibosons1j_tot_err/dibosons1j_tot);

  float nBkg_1j_Out = nWjets_1j_Out + nDY_1j_Out + nTop_1j_Out + nDiBosons_1j_Out;
  float nBkg_1j_Out_err = quadrSum(nWjets_1j_Out_err,nDY_1j_Out_err,nTop_1j_Out_err,nDiBosons_1j_Out_err);

  ///// SUMMARY BKG ////
  // --- 0 jet ---
  std::cout << "---->  BACKGROUND SUMMARY 0 jet  <-------" << std::endl;
  std::cout << "bkg in data:" << std::endl;   
  std::cout.precision(3);
  std::cout << "W+jets =\t" << nWjets_0j_Out << " +/- " << nWjets_0j_Out_err << std::endl;
  std::cout << "top =\t\t" << nTop_0j_Out << " +/- " << nTop_0j_Out_err << std::endl;
  std::cout << "DY =\t\t" << nDY_0j_Out << " +/- " << nDY_0j_Out_err << std::endl;
  std::cout << "WZ,ZZ =\t\t" << nDiBosons_0j_Out << " +/- " << nDiBosons_0j_Out_err << std::endl;
  std::cout << "---------------------------------------------------------------" << std::endl;
  std::cout << "TOTAL:\n\t\t" << nBkg_0j_Out << " +/- " << nBkg_0j_Out_err << std::endl;
  std::cout << "---------------------------------------------------------------\n\n\n" << std::endl;

  // --- 1 jet ---
  std::cout << "---->  BACKGROUND SUMMARY 1 jet  <-------" << std::endl;
  std::cout << "bkg in data:" << std::endl;   
  std::cout.precision(3);
  std::cout << "W+jets =\t" << nWjets_1j_Out << " +/- " << nWjets_1j_Out_err << std::endl;
  std::cout << "top =\t\t" << nTop_1j_Out << " +/- " << nTop_1j_Out_err << std::endl;
  std::cout << "DY =\t\t" << nDY_1j_Out << " +/- " << nDY_1j_Out_err << std::endl;
  std::cout << "WZ,ZZ =\t\t" << nDiBosons_1j_Out << " +/- " << nDiBosons_1j_Out_err << std::endl;
  std::cout << "---------------------------------------------------------------" << std::endl;
  std::cout << "TOTAL:\n\t\t" << nBkg_1j_Out << " +/- " << nBkg_1j_Out_err << std::endl;
  std::cout << "---------------------------------------------------------------\n\n\n" << std::endl;
  /////////////////////


  // WW ESTIMATION AT THE MASS INDEPENDENT LEVEL
  // rescaling to the observed number of WW by Boris (naive)
  integralData0jCtrl *= 1.0;
  integralData1jCtrl *= 1.0;
  float nWW_0j_OutData_tot = integralData0jCtrl - nBkg_0j_Out;
  float nWW_0j_OutData_tot_err = quadrSum(sqrt(integralData0jCtrl),nBkg_0j_Out_err);

  float nWW_0j_Data_tot = nWW_0j_OutData_tot / ratio0j[ww];
  float nWW_0j_Data_tot_err = nWW_0j_Data_tot * quadrSum(nWW_0j_OutData_tot_err/nWW_0j_OutData_tot, ratio0j_err[ww]/ratio0j[ww]);
  
  float sf_0j = nWW_0j_Data_tot/integrals0jTot[ww];
  float sf_0j_err = sf_0j * quadrSum(nWW_0j_OutData_tot_err/nWW_0j_OutData_tot, error0jTot[ww]/integrals0jTot[ww]);

  float nWW_1j_OutData_tot = integralData1jCtrl - nBkg_1j_Out;
  float nWW_1j_OutData_tot_err = quadrSum(sqrt(integralData1jCtrl),nBkg_1j_Out_err);

  float nWW_1j_Data_tot = nWW_1j_OutData_tot / ratio1j[ww];
  float nWW_1j_Data_tot_err = nWW_1j_Data_tot * quadrSum(nWW_1j_OutData_tot_err/nWW_1j_OutData_tot, ratio1j_err[ww]/ratio1j[ww]);

  float sf_1j = nWW_1j_Data_tot/integrals1jTot[ww];
  float sf_1j_err = sf_1j * quadrSum(nWW_1j_OutData_tot_err/nWW_1j_OutData_tot, error1jTot[ww]/integrals1jTot[ww]);

  std::cout << "               \t\t Data selected \t\t All bkg \t\t WW \t\t MC \t\t SF(data/MC)" << std::endl;
  std::cout << " ---------- 0 jet ----------- " << std::endl;
  std::cout << "control region \t\t " 
            << integralData0jCtrl << " \t\t " 
            << nBkg_0j_Out << " +/- " << nBkg_0j_Out_err << " \t\t " 
            << nWW_0j_OutData_tot << " +/- " << nWW_0j_OutData_tot_err << " \t\t "
            << integrals0jCtrl[ww] << " +/- " << error0jCtrl[ww]
            << std::endl;
  std::cout << "all region \t\t " 
            << integralData0jTot << " \t\t " 
            << " not used \t\t " 
            << nWW_0j_Data_tot << " +/- " << nWW_0j_Data_tot_err << " \t\t "
            << integrals0jTot[ww] << " +/- " << error0jTot[ww] << " \t\t "
            << sf_0j << " +/- " << sf_0j_err
            << std::endl;
  std::cout << " ---------------------------- " << std::endl;
  std::cout << " ---------- 1 jet ----------- " << std::endl;
  std::cout << "control region \t\t " 
            << integralData1jCtrl << " \t\t " 
            << nBkg_1j_Out << " +/- " << nBkg_1j_Out_err << " \t\t " 
            << nWW_1j_OutData_tot << " +/- " << nWW_1j_OutData_tot_err << " \t\t "
            << integrals1jCtrl[ww] << " +/- " << error1jCtrl[ww]
            << std::endl;
  std::cout << "all region \t\t " 
            << integralData1jTot << " \t\t " 
            << " not used \t\t " 
            << nWW_1j_Data_tot << " +/- " << nWW_1j_Data_tot_err << " \t\t "
            << integrals1jTot[ww] << " +/- " << error1jTot[ww] << " \t\t "
            << sf_1j << " +/- " << sf_1j_err
            << std::endl;
  std::cout << " ---------------------------- " << std::endl;


  /// WW CHANNEL FRACTIONS FROM MC 0j ///
  // they are used in case of low stat to use the summed WW estimation in the 4 sub-channels
  TH1F *WWEEH = new TH1F("WWEEH","",50,0,180); // all
  TH1F *WWMMH = new TH1F("WWMMH","",50,0,180); // all
  TH1F *WWEMH = new TH1F("WWEMH","",50,0,180); // all
  TH1F *WWMEH = new TH1F("WWMEH","",50,0,180); // all
  treeWW->Project("WWEEH","deltaPhi","(WWSel && finalstate==0)*weight*puweight");
  treeWW->Project("WWMMH","deltaPhi","(WWSel && finalstate==1)*weight*puweight");
  treeWW->Project("WWEMH","deltaPhi","(WWSel && finalstate==2)*weight*puweight");
  treeWW->Project("WWMEH","deltaPhi","(WWSel && finalstate==3)*weight*puweight");

  float nWW_0j_MC[4] = { WWEEH->Integral(), WWMMH->Integral(), WWEMH->Integral(), WWMEH->Integral() };
  float nWW_0j_MC_err[4] = { yieldErrPoisson(nWW_0j_MC[ee], WWEEH->GetEntries()), yieldErrPoisson(nWW_0j_MC[mm], WWMMH->GetEntries()),
                             yieldErrPoisson(nWW_0j_MC[em], WWEMH->GetEntries()), yieldErrPoisson(nWW_0j_MC[me], WWMEH->GetEntries()) };
 
  std::cout << "Using channel fractions from MC!" << std::endl;
  float frac[4];
  for(int icha=0; icha<4; icha++) {
    frac[icha] = nWW_0j_MC[icha]/integrals0jTot[ww];
    std::cout << "Fraction of WW in channel " << icha << " = " << frac[icha] << std::endl;
  }

  /// 1 jet entries at the level of WW selection from the tree (counter not existent)
  treeWW->Project("WWEEH","deltaPhi","(WWSel1j && finalstate==0)*weight*puweight");
  treeWW->Project("WWMMH","deltaPhi","(WWSel1j && finalstate==1)*weight*puweight");
  treeWW->Project("WWEMH","deltaPhi","(WWSel1j && finalstate==2)*weight*puweight");
  treeWW->Project("WWMEH","deltaPhi","(WWSel1j && finalstate==3)*weight*puweight");
  float nWW_1j_MC[4] = { WWEEH->Integral(), WWMMH->Integral(), WWEMH->Integral(), WWMEH->Integral() };
  float nWW_1j_MC_err[4] = { yieldErrPoisson(nWW_1j_MC[ee], WWEEH->GetEntries()), yieldErrPoisson(nWW_1j_MC[mm], WWMMH->GetEntries()),
                             yieldErrPoisson(nWW_1j_MC[em], WWEMH->GetEntries()), yieldErrPoisson(nWW_1j_MC[me], WWMEH->GetEntries()) };
  float nEv_endWW1j[4] = { WWEEH->GetEntries(), WWMMH->GetEntries(), WWEMH->GetEntries(), WWMEH->GetEntries() };


  /////////////////////

  ofstream textfile;
  textfile.open("WWYieldsData.txt", ios_base::trunc);
  textfile.precision(2);

  ofstream tablefile1;
  tablefile1.open("WWYieldsData_ForTable_0j.txt", ios_base::trunc);
  tablefile1.precision(2);

  ofstream tablefile2;
  tablefile2.open("WWYieldsData_ForTable_1j.txt", ios_base::trunc);
  tablefile2.precision(2);

  ofstream tablefile3;
  tablefile3.open("WWYieldsMC_ForTable_0j.txt", ios_base::trunc);
  tablefile3.precision(2);

  ofstream tablefile4;
  tablefile4.open("WWYieldsMC_ForTable_1j.txt", ios_base::trunc);
  tablefile4.precision(2);

  // these are for limits
  ofstream cardfile[4][2]; //[cha][jetbin]
  for(int icha=0; icha<4; icha++) {
    for(int j=0; j<2; j++) {
      char fileName[2];
      if(icha==ee) sprintf(fileName,"WWCard_ee_%dj.txt",j);
      if(icha==mm) sprintf(fileName,"WWCard_mm_%dj.txt",j);
      if(icha==em) sprintf(fileName,"WWCard_em_%dj.txt",j);
      if(icha==me) sprintf(fileName,"WWCard_me_%dj.txt",j);
      cardfile[icha][j].open(fileName, ios_base::trunc);
      cardfile[icha][j].precision(3);
    }
  }
    

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    countEvents(mass,"EE");
    countEvents(mass,"MM");
    countEvents(mass,"EM");
    countEvents(mass,"ME");

    float nWWData_HiggsSel_0j[4], nWWData_HiggsSel_0j_err[4];
    float nWWMC_HiggsSel_0j[4], nWWMC_HiggsSel_0j_err[4];

    float nWWData_HiggsSel_1j[4], nWWData_HiggsSel_1j_err[4];
    float nWWMC_HiggsSel_1j[4], nWWMC_HiggsSel_1j_err[4];

    for(int icha=0;icha<4;icha++) {
      float eff_0j = (nEv_endWW[icha]==0) ? 0. : nEv_end0j[icha]/nEv_endWW[icha];
      float eff_0j_err = (nEv_endWW[icha]==0) ? 0. : sqrt(eff_0j*(1-eff_0j)/nEv_endWW[icha]);
      
      nWWData_HiggsSel_0j[icha] = nWW_0j_Data_tot * frac[icha] * eff_0j;
      nWWData_HiggsSel_0j_err[icha] = nWWData_HiggsSel_0j[icha] * quadrSum(nWW_0j_Data_tot_err/nWW_0j_Data_tot,eff_0j_err/eff_0j);

      nWWMC_HiggsSel_0j[icha] = nWW_0j_MC[icha] * eff_0j;
      nWWMC_HiggsSel_0j_err[icha] = nWWMC_HiggsSel_0j[icha] * quadrSum(nWW_0j_MC_err[icha]/nWW_0j_MC[icha],eff_0j_err/eff_0j);

      float eff_1j = nEv_end1j[icha]/nEv_endWW1j[icha];
      float eff_1j_err = sqrt(eff_1j*(1-eff_1j)/nEv_endWW1j[icha]);

      nWWData_HiggsSel_1j[icha] = nWW_1j_Data_tot * frac[icha] * eff_1j;
      nWWData_HiggsSel_1j_err[icha] =  nWWData_HiggsSel_1j[icha] * quadrSum(nWW_1j_Data_tot_err/nWW_1j_Data_tot,eff_1j_err/eff_1j);      

      nWWMC_HiggsSel_1j[icha] = nWW_1j_MC[icha] * eff_1j;
      nWWMC_HiggsSel_1j_err[icha] = nWWMC_HiggsSel_1j[icha] * quadrSum(nWW_1j_MC_err[icha]/nWW_1j_MC[icha],eff_1j_err/eff_1j);

      char channelName[2];
      if(icha==ee) sprintf(channelName,"EE");
      if(icha==mm) sprintf(channelName,"MM");
      if(icha==em) sprintf(channelName,"EM");
      if(icha==me) sprintf(channelName,"ME");

      // for Giovanni
      float alpha_0j = scaleFactorLumi * frac[icha] / ratio0j[ww] * eff_0j;
      float alpha_0j_err = alpha_0j * quadrSum(ratio0j_err[ww]/ratio0j[ww],eff_0j_err/eff_0j);

      float alpha_1j = scaleFactorLumi * frac[icha] / ratio1j[ww] * eff_1j;
      float alpha_1j_err = alpha_1j * quadrSum(ratio0j_err[ww]/ratio0j[ww],eff_1j_err/eff_1j);

      cardfile[icha][0] << mass 
                        << "\t" << nWW_0j_OutData_tot << "\t" << alpha_0j
                        << "\t" <<  alpha_0j_err 
                        << std::endl;

      cardfile[icha][1] << mass 
                        << "\t" << nWW_1j_OutData_tot << "\t" << alpha_1j
                        << "\t" <<  alpha_1j_err 
                        << std::endl;

      /////////////// 

      textfile << channelName << ": Higgs Mass = " << mass
               << "\tdata 0 jet = " << nWWData_HiggsSel_0j[icha] << " +/- " << nWWData_HiggsSel_0j_err[icha]
               << "\tMC 0 jet = " << nWWMC_HiggsSel_0j[icha] << " +/- " << nWWMC_HiggsSel_0j_err[icha]
               << "\tdata 1 jet = " << nWWData_HiggsSel_1j[icha] << " +/- " << nWWData_HiggsSel_1j_err[icha]
               << "\tMC 1 jet = " << nWWMC_HiggsSel_1j[icha] << " +/- " << nWWMC_HiggsSel_1j_err[icha]
               << std::endl;
    }


    // summary table for limits                                                                                          
    if (i==0) {
      tablefile1 << "zero jets bin data" << endl;
      tablefile1 << "\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile1 << mass
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_0j[1] << " +/- " << nWWData_HiggsSel_0j_err[1]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_0j[3] << " +/- " << nWWData_HiggsSel_0j_err[3]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_0j[2] << " +/- " << nWWData_HiggsSel_0j_err[2]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_0j[0] << " +/- " << nWWData_HiggsSel_0j_err[0]
               << std::endl;

    if (i==0) {
      tablefile2 << "one jet bin data" << endl;
      tablefile2 << "\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile2 << mass
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_1j[1] << " +/- " << nWWData_HiggsSel_1j_err[1]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_1j[3] << " +/- " << nWWData_HiggsSel_1j_err[3]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_1j[2] << " +/- " << nWWData_HiggsSel_1j_err[2]
               << " " << "\t" << scaleFactorLumi * nWWData_HiggsSel_1j[0] << " +/- " << nWWData_HiggsSel_1j_err[0]
               << std::endl;

    if (i==0) {
      tablefile3 << "zero jets bin MC" << endl;
      tablefile3 << "\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile3 << mass
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_0j[1] << " +/- " << nWWMC_HiggsSel_0j_err[1]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_0j[3] << " +/- " << nWWMC_HiggsSel_0j_err[3]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_0j[2] << " +/- " << nWWMC_HiggsSel_0j_err[2]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_0j[0] << " +/- " << nWWMC_HiggsSel_0j_err[0]
               << std::endl;

    if (i==0) {
      tablefile4 << "one jet bin MC" << endl;
      tablefile4 << "\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile4 << mass
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_1j[1] << " +/- " << nWWMC_HiggsSel_1j_err[1]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_1j[3] << " +/- " << nWWMC_HiggsSel_1j_err[3]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_1j[2] << " +/- " << nWWMC_HiggsSel_1j_err[2]
               << " " << "\t" << scaleFactorLumi * nWWMC_HiggsSel_1j[0] << " +/- " << nWWMC_HiggsSel_1j_err[0]
               << std::endl;


    float nWWData_HiggsSel_0j_Tot = nWWData_HiggsSel_0j[ee] + nWWData_HiggsSel_0j[mm] + nWWData_HiggsSel_0j[em] + nWWData_HiggsSel_0j[me];
    float nWWData_HiggsSel_0j_Tot_err = quadrSum(nWWData_HiggsSel_0j_err[ee],nWWData_HiggsSel_0j_err[mm],nWWData_HiggsSel_0j_err[em],nWWData_HiggsSel_0j_err[me]);

    float nWWMC_HiggsSel_0j_Tot = nWWMC_HiggsSel_0j[ee] + nWWMC_HiggsSel_0j[mm] + nWWMC_HiggsSel_0j[em] + nWWMC_HiggsSel_0j[me];
    float nWWMC_HiggsSel_0j_Tot_err = quadrSum(nWWMC_HiggsSel_0j_err[ee],nWWMC_HiggsSel_0j_err[mm],nWWMC_HiggsSel_0j_err[em],nWWMC_HiggsSel_0j_err[me]);

    float nWWData_HiggsSel_1j_Tot = nWWData_HiggsSel_1j[ee] + nWWData_HiggsSel_1j[mm] + nWWData_HiggsSel_1j[em] + nWWData_HiggsSel_1j[me];
    float nWWData_HiggsSel_1j_Tot_err = quadrSum(nWWData_HiggsSel_1j_err[ee],nWWData_HiggsSel_1j_err[mm],nWWData_HiggsSel_1j_err[em],nWWData_HiggsSel_1j_err[me]);

    float nWWMC_HiggsSel_1j_Tot = nWWMC_HiggsSel_1j[ee] + nWWMC_HiggsSel_1j[mm] + nWWMC_HiggsSel_1j[em] + nWWMC_HiggsSel_1j[me];
    float nWWMC_HiggsSel_1j_Tot_err = quadrSum(nWWMC_HiggsSel_1j_err[ee],nWWMC_HiggsSel_1j_err[mm],nWWMC_HiggsSel_1j_err[em],nWWMC_HiggsSel_1j_err[me]);
    
    
    textfile << "\t===>> TOTAL: Higgs Mass = " << mass 
             << "\tdata 0 jet = " << nWWData_HiggsSel_0j_Tot << " +/- " << nWWData_HiggsSel_0j_Tot_err 
             << "\tdata 1 jet = " << nWWData_HiggsSel_1j_Tot << " +/- " << nWWData_HiggsSel_1j_Tot_err
             << "\tMC 0 jet = " << nWWMC_HiggsSel_0j_Tot << " +/- " << nWWMC_HiggsSel_0j_Tot_err 
             << "\tMC 1 jet = " << nWWMC_HiggsSel_1j_Tot << " +/- " << nWWMC_HiggsSel_1j_Tot_err
             << std::endl;

  }

  std::cout << "Full WW yields in data in:  WWYieldsData.txt " << std::endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr) {
  float val = ntag * (1-eff2b) / eff2b;
  float err = ntag * eff2berr / pow(eff2b,2);
  return std::make_pair(val,err);
}

std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK) {
  float val = R * (nDYin - 0.5 * nemu * K);
  float err = sqrt(pow(sigmaR,2) * pow(nDYin - 0.5 * nemu * K , 2) + 
                   1./4. * pow(nemu,2) * R*R * sigmaK*sigmaK + 
                   nDYin * R*R +
                   nemu * 1./4. * K*K * R*R); 
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
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.2.X/MC2011_Merged_V5/OptimMH%d/Spring11_V5HWW/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola-WWFilter/*Counters.root",mass);  
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

void estimateRWWRegion() {
  TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TTree *treeZjets = (TTree*)fileZjets->Get("T1");
  TH1F *ZeejetsH_IN = new TH1F("ZeejetsH_IN","",50,0,180);
  TH1F *ZmmjetsH_IN = new TH1F("ZmmjetsH_IN","",50,0,180);
  TH1F *ZeejetsH_OUT = new TH1F("ZeejetsH_OUT","",50,0,180);
  TH1F *ZmmjetsH_OUT = new TH1F("ZmmjetsH_OUT","",50,0,180);

  treeZjets->Project("ZeejetsH_IN","deltaPhi","(finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0)*weight*puweight");
  treeZjets->Project("ZmmjetsH_IN","deltaPhi","(finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1)*weight*puweight");
  treeZjets->Project("ZeejetsH_OUT","deltaPhi","(finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && eleInvMass>100 && finalstate==0)*weight*puweight");
  treeZjets->Project("ZmmjetsH_OUT","deltaPhi","(finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && eleInvMass>100 && finalstate==1)*weight*puweight");

  float nZmmjetsH_OUT = ZmmjetsH_OUT->Integral();
  float nZmmjetsH_OUT_err = yieldErrPoisson(nZmmjetsH_OUT,ZmmjetsH_OUT->GetEntries());
  float nZeejetsH_OUT = ZeejetsH_OUT->Integral();
  float nZeejetsH_OUT_err = yieldErrPoisson(nZeejetsH_OUT,ZeejetsH_OUT->GetEntries());

  float nZmmjetsH_IN = ZmmjetsH_IN->Integral();
  float nZmmjetsH_IN_err = yieldErrPoisson(nZmmjetsH_IN,ZmmjetsH_IN->GetEntries());
  float nZeejetsH_IN = ZeejetsH_IN->Integral();
  float nZeejetsH_IN_err = yieldErrPoisson(nZeejetsH_IN,ZeejetsH_IN->GetEntries());

  Rmm = nZmmjetsH_OUT/nZmmjetsH_IN;
  Rmm_err = Rmm * quadrSum(nZmmjetsH_OUT_err/nZmmjetsH_OUT, nZmmjetsH_IN_err/nZmmjetsH_IN);

  Ree = nZeejetsH_OUT/nZeejetsH_IN;
  Ree_err = Ree * quadrSum(nZeejetsH_OUT_err/nZeejetsH_OUT, nZeejetsH_IN_err/nZeejetsH_IN);

}
