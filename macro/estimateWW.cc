#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
float nVeto(float ntag, float eff1b);

void estimateWW() {

  // constants
  float eff_2b = 0.57;
  float eff_2b_err = 0.02; 
  float eff_softmu = 0.36; // i.e. e1b soft mu = 0.20 w/o error

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

  std::vector<TFile*> bakgrounds;
  bakgrounds.push_back(fileTop);
  bakgrounds.push_back(fileWjets);
  bakgrounds.push_back(fileZjets);
  bakgrounds.push_back(fileDiBosons);

  std::vector<float> backgroundsUnc; // dummy values so far
  backgroundsUnc.push_back(0.5);
  backgroundsUnc.push_back(0.5);
  backgroundsUnc.push_back(0.1);
  backgroundsUnc.push_back(0.1);
  
  TH1F *dataH = new TH1F("dataH","",50,0,180);
  TH1F *WWH = new TH1F("WWH","",50,0,180);
  TH1F *topH = new TH1F("topH","",50,0,180);
  TH1F *btagHData = new TH1F("btagHData","",50,0,180);
  TH1F *WjetsH = new TH1F("WjetsH","",50,0,180);
  TH1F *ZjetsH = new TH1F("ZjetsH","",50,0,180);
  TH1F *DiBosonsH = new TH1F("DiBosonsH","",50,0,180);

  treeData->Project("dataH","deltaPhi","eleInvMass>100 && WWSel");
  treeWW->Project("WWH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeTop->Project("topH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeData->Project("btagHData","deltaPhi","eleInvMass>100 && jetVeto && bTagTrackCount>2.1");
  treeWjets->Project("WjetsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeZjets->Project("ZjetsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeDiBosons->Project("DiBosonsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");

  float nDataOut = dataH->Integral();
  float nWWOut = WWH->Integral();
  float nTopOut = topH->Integral();
  
  // top estimation from data (0-jet bin method)
  float nBTagTagOut_data = btagHData->Integral();
  float nTopOutBTagVeto_data = nVeto(nBTagTagOut_data, eff_2b);
  float nTopOutBTagVeto_data_err = nVeto(nBTagTagOut_data, eff_2b_err); // wrong for the moment
  float nTopOutSoftMuVeto_data = nTopOutBTagVeto_data * pow((1-eff_softmu),2); // efficiency of passing the soft muon veto (both the b's).
  float nTopOutSoftMuVeto_data_err = nTopOutBTagVeto_data_err; // for the moment...

  std::cout << "Top out from MC = " << nTopOut << " tagged events = " << nBTagTagOut_data
            << "   Top out from data after btag veto = " << nTopOutBTagVeto_data << " +/- " << nTopOutBTagVeto_data_err 
            << "   Top out from data after soft mu veto =  " << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;

  float nWjetsOut = WjetsH->Integral();
  float nZjetsOut = ZjetsH->Integral();
  float nDiBosonsOut = DiBosonsH->Integral();

  float nWWOut_Ne = WWH->GetEntries();
  float nTopOut_Ne = topH->GetEntries();
  float nWjetsOut_Ne = WjetsH->GetEntries();
  float nZjetsOut_Ne = ZjetsH->GetEntries();
  float nDiBosonsOut_Ne = DiBosonsH->GetEntries();

  float bkgOut_Ne = nTopOut_Ne + nWjetsOut_Ne + nZjetsOut_Ne + nDiBosonsOut_Ne;

  float bkgOut = nTopOut + nWjetsOut + nZjetsOut + nDiBosonsOut;
  float nWWDataOut = nDataOut - bkgOut;

  float nDataOut_err = sqrt(nDataOut);
  float nBkgOut_err_stat = sqrt(bkgOut_Ne)/bkgOut_Ne * bkgOut;
  float nWWDataOut_err = quadrSum(nDataOut_err,nBkgOut_err_stat);

  float dataOmc = nWWDataOut/nWWOut;
  float nWWOut_err = sqrt(nWWOut_Ne)/nWWOut_Ne * nWWOut;
  float dataOmc_err = dataOmc * quadrSum(nWWDataOut_err/nWWDataOut,nWWOut_err/nWWOut);

  std::cout << "Scale factor data / MC = " << dataOmc << " +/- " << dataOmc_err << std::endl;
  
  std::map<int,float> nWWIn;
  nWWIn.insert(std::make_pair(120,2.18));
  nWWIn.insert(std::make_pair(130,2.27));
  nWWIn.insert(std::make_pair(140,1.88));
  nWWIn.insert(std::make_pair(150,0.89));
  nWWIn.insert(std::make_pair(160,0.81));
  nWWIn.insert(std::make_pair(170,0.81));
  nWWIn.insert(std::make_pair(180,0.81));
  nWWIn.insert(std::make_pair(190,0.81));
  nWWIn.insert(std::make_pair(200,0.81));

  std::cout << "Data Tot = " << nDataOut << " " 
            << nTopOut << " " << nWjetsOut << " " << nZjetsOut << " " << nDiBosonsOut << " "
            << "Bkg tot = " << bkgOut
            << " nWW data out = " << nWWDataOut 
            << " nWW MC out = " << nWWOut << std::endl; 
  
//   for(int i=0; i<(int)nWWIn.size(); i++) {
//     int mass = (nWWIn.at(i)).first;
//     char header[20];
//     sprintf(header,"NWW[%d] (data) = ", mass);
//     std::cout << header << nWWIn[mass] * dataOmc << " +/- " << nWWIn[mass] * dataOmc_err << " (stat)" << std::endl;
//   }

  std::cout <<"NWW[120] corr. factor = " << dataOmc << " +/- " << dataOmc_err << " (stat)" << std::endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

float nVeto(float ntag, float eff2b) {
  return ntag * (1-eff2b) / eff2b;
}
