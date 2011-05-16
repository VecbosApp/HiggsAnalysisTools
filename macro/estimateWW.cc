#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK);

void estimateWW() {

  // constants
  float eff_2b = 0.4983;
  float eff_2b_err = 0.03; 
  float eff_2b_softmu = 0.1677; // MC 

  float eff_softmu_Z = 0.867;

  float Rmm = 0.187398;
  float Rmm_err = 0.00355062;
  float kmm = 0.592288;
  float kmm_err = 0.0139356;
  float Ree = 0.152089;
  float Ree_err = 0.00338336;
  float kee = 0.422092;
  float kee_err = 0.00874687;


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
  TH1F *ZeejetsH = new TH1F("ZeejetsH","",50,0,180);
  TH1F *ZmmjetsH = new TH1F("ZmmjetsH","",50,0,180);
  TH1F *ZemjetsH = new TH1F("ZemjetsH","",50,0,180);
  TH1F *neeInH = new TH1F("neeInH","",50,0,180);
  TH1F *nmmInH = new TH1F("nmmInH","",50,0,180);
  TH1F *nemInH = new TH1F("nemInH","",50,0,180);
  TH1F *DiBosonsH = new TH1F("DiBosonsH","",50,0,180);

  treeData->Project("dataH","deltaPhi","eleInvMass>100 && WWSel");
  treeWW->Project("WWH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeTop->Project("topH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeData->Project("btagHData","deltaPhi","eleInvMass>100 && jetVeto && bTagTrackCount>2.1");
  treeWjets->Project("WjetsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");
  treeZjets->Project("ZeejetsH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==0)*weight");
  treeZjets->Project("ZmmjetsH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==1)*weight");
  treeZjets->Project("ZemjetsH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==2)*weight");
  treeData->Project("neeInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nmmInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nemInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==2"); // missing softmu... not available in red trees... hopefully small contrib here
  treeDiBosons->Project("DiBosonsH","deltaPhi","(eleInvMass>100 && WWSel)*weight");

  float nDataOut = dataH->Integral();
  float nWWOut = WWH->Integral();
  float nTopOut = topH->Integral();
  
  // top estimation from data (0-jet bin method)
  float nBTagTagOut_data = btagHData->Integral();
  float nTopOutBTagVeto_data = (nVeto(nBTagTagOut_data, eff_2b, eff_2b_err)).first;
  float nTopOutBTagVeto_data_err = (nVeto(nBTagTagOut_data, eff_2b, eff_2b_err)).second; 
  float nTopOutSoftMuVeto_data = nTopOutBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  float nTopOutSoftMuVeto_data_err = nTopOutBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "TOP ESTIMATION..." << std::endl;
  std::cout << "Using eff_2b = " << eff_2b << " +/- " << eff_2b_err << std::endl;
  std::cout << "Using eff_2b_softmu = " << eff_2b_softmu << std::endl;
  std::cout << "Tagged events = " << nBTagTagOut_data
            << "   Top out from data after btag veto = " << nTopOutBTagVeto_data << " +/- " << nTopOutBTagVeto_data_err 
            << "   Top out from data after soft mu veto =  " << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;

  std::cout << "Top from MC = " << nTopOut << "\tTop from data = " << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;
  std::cout << "END TOP ESTIMATION." << std::endl;

  float nWjetsOut = WjetsH->Integral();

  std::cout << "DY ESTIMATION..." << std::endl;
  float nZeejetsOut = ZeejetsH->Integral();
  float nZmmjetsOut = ZmmjetsH->Integral();
  float nZemjetsOut = ZemjetsH->Integral();

  float neeIn = neeInH->Integral() * eff_softmu_Z;
  float nmmIn = nmmInH->Integral() * eff_softmu_Z;
  float nemIn = nemInH->Integral() * eff_softmu_Z;


  float neeExp = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).first;
  float neeExp_err = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).second;
  float nmmExp = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).first;
  float nmmExp_err = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).second;

  std::cout << "neeIn = " << neeIn << "\tnmmIn = " << nmmIn << "\tnemIn = " << nemIn << std::endl;
  std::cout << "nEE MC = " << nZeejetsOut << "\tData = " << neeExp << " +/- " << neeExp_err << std::endl;
  std::cout << "nMM MC = " << nZmmjetsOut << "\tData = " << nmmExp << " +/- " << nmmExp_err << std::endl;
  std::cout << "END DY ESTIMATION." << std::endl;

  float nDiBosonsOut = DiBosonsH->Integral();

  float nWWOut_Ne = WWH->GetEntries();
  float nTopOut_Ne = topH->GetEntries();
  float nWjetsOut_Ne = WjetsH->GetEntries();
  float nZjetsOut_Ne = ZeejetsH->GetEntries() + ZmmjetsH->GetEntries() + ZemjetsH->GetEntries();
  float nDiBosonsOut_Ne = DiBosonsH->GetEntries();

  float bkgOut_Ne = nTopOut_Ne + nWjetsOut_Ne + nZjetsOut_Ne + nDiBosonsOut_Ne;

  float bkgOut = nTopOut + nWjetsOut + nZeejetsOut + nZmmjetsOut + nZemjetsOut + nDiBosonsOut;
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
            << nTopOut << " " << nWjetsOut << " " << nZeejetsOut + nZmmjetsOut + nZemjetsOut << " " << nDiBosonsOut << " "
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
