#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
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

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK, float nZV);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);

void estimateDY(int njets, bool useDataRk) {

  // constants
  //  float eff_softmu_Z = 0.867;
  //  float eff_softmu_Z = 1;

  char njcut[30];
  sprintf(njcut, "njets==%d", njets);
  char wwselcut[30];
  if(njets==0) sprintf(wwselcut,"WWSel");
  else if(njets==1) sprintf(wwselcut,"WWSel1j");
  else {
    std::cout << "Jet bin must be 0/1" << std::endl;
    return;
  }

  TFile *fileData  = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileOthers = TFile::Open("results/datasets_trees/others_ll.root");

  TTree *treeData  = (TTree*)fileData->Get("T1");
  TTree *treeZjets = (TTree*)fileZjets->Get("T1");
  TTree *treeOthers = (TTree*)fileOthers->Get("T1");

  TH1F *ZeejetsH = new TH1F("ZeejetsH","",50,0,180);
  TH1F *ZmmjetsH = new TH1F("ZmmjetsH","",50,0,180);
  TH1F *ZemjetsH = new TH1F("ZemjetsH","",50,0,180);
  TH1F *ZmejetsH = new TH1F("ZmejetsH","",50,0,180);

  // to subtract WW,ZZ
  TH1F *OthersEEH = new TH1F("OthersEEH","",50,0,180);
  TH1F *OthersMMH = new TH1F("OthersMMH","",50,0,180);

  TH1F *neeInH = new TH1F("neeInH","",50,0,180);
  TH1F *nmmInH = new TH1F("nmmInH","",50,0,180);
  TH1F *nemInH = new TH1F("nemInH","",50,0,180);

  // to estimate R and k
  TH1F *neeLooseInH = new TH1F("neeLooseInH","",50,0,180);
  TH1F *neeLooseOutH = new TH1F("neeLooseOutH","",50,0,180);
  TH1F *nmmLooseInH = new TH1F("nmmLooseInH","",50,0,180);
  TH1F *nmmLooseOutH = new TH1F("nmmLooseOutH","",50,0,180);

  // to estimate R and k from data
  TH1F *neeLooseInHData = new TH1F("neeLooseInHData","",50,0,180);
  TH1F *neeLooseOutHData = new TH1F("neeLooseOutHData","",50,0,180);
  TH1F *nmmLooseInHData = new TH1F("nmmLooseInHData","",50,0,180);
  TH1F *nmmLooseOutHData = new TH1F("nmmLooseOutHData","",50,0,180);  
  TH1F *nemLooseInHData = new TH1F("nemLooseInHData","",50,0,180);
  TH1F *nemLooseOutHData = new TH1F("nemLooseOutHData","",50,0,180);  

  treeZjets->Project("ZeejetsH","deltaPhi",(TString("(")+TString(wwselcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && finalstate==0)*weight*puweight")).Data());
  treeZjets->Project("ZmmjetsH","deltaPhi",(TString("(")+TString(wwselcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && finalstate==1)*weight*puweight")).Data());
  treeZjets->Project("ZemjetsH","deltaPhi",(TString("(")+TString(wwselcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && finalstate==2)*weight*puweight")).Data());
  treeZjets->Project("ZmejetsH","deltaPhi",(TString("(")+TString(wwselcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && finalstate==3)*weight*puweight")).Data());

  treeZjets->Project("neeLooseInH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0)*weight*puweight")).Data());
  treeZjets->Project("neeLooseOutH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)>15 && finalstate==0)*weight*puweight")).Data());
  treeZjets->Project("nmmLooseInH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1)*weight*puweight")).Data());
  treeZjets->Project("nmmLooseOutH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)>15 && finalstate==1)*weight*puweight")).Data());

  treeData->Project("neeLooseInHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)<15 && finalstate==0)")).Data());
  treeData->Project("neeLooseOutHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)>15 && finalstate==0)")).Data());
  treeData->Project("nmmLooseInHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)<15 && finalstate==1)")).Data());
  treeData->Project("nmmLooseOutHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)>15 && finalstate==1)")).Data());
  treeData->Project("nemLooseInHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)<15 && (finalstate==2 || finalstate==3))")).Data());
  treeData->Project("nemLooseOutHData","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>20 && ")+TString(njcut)+TString( " && bTagTrackCount<2.1 && eleInvMass>12 && abs(eleInvMass-91.1876)>15 && (finalstate==2 || finalstate==3))")).Data());

  // contribution under the peak
  treeOthers->Project("OthersEEH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>40 &&")+TString(njcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0)*weight*puweight")).Data());
  treeOthers->Project("OthersMMH","deltaPhi",(TString("(finalLeptons && pfMet>20 && projMet>40 &&")+TString(njcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1)*weight*puweight")).Data());
  
  treeData->Project("neeInH","deltaPhi",(TString("eleInvMass>12 && finalLeptons && pfMet>20 && projMet>40 && ")+TString(njcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0")).Data()); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nmmInH","deltaPhi",(TString("eleInvMass>12 && finalLeptons && pfMet>20 && projMet>40 && ")+TString(njcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1")).Data()); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nemInH","deltaPhi",(TString("eleInvMass>12 && finalLeptons && pfMet>20 && projMet>40 && ")+TString(njcut)+TString(" && ((leadingJetPt>15 && abs(deltaPhi_LL_JET)<165) || leadingJetPt<=15) && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && (finalstate==2 || finalstate==3)")).Data()); // missing softmu... not available in red trees... hopefully small contrib here

  // R and k estimations
  float R[2], R_err[2], k[2], k_err[2];
  R[ee] = neeLooseOutH->Integral() / neeLooseInH->Integral();
  R_err[ee] = R[ee] * quadrSum( yieldErrPoisson(neeLooseOutH->Integral(),neeLooseOutH->GetEntries()) / neeLooseOutH->Integral(),
                                yieldErrPoisson(neeLooseInH->Integral(),neeLooseInH->GetEntries()) / neeLooseInH->Integral() );

  R[mm] = nmmLooseOutH->Integral() / nmmLooseInH->Integral();
  R_err[mm] = R[mm] * quadrSum( yieldErrPoisson(nmmLooseOutH->Integral(),nmmLooseOutH->GetEntries()) / nmmLooseOutH->Integral(),
                                yieldErrPoisson(nmmLooseInH->Integral(),nmmLooseInH->GetEntries()) / nmmLooseInH->Integral() );

  k[ee] = sqrt(neeLooseInH->Integral() / nmmLooseInH->Integral());
  k_err[ee] = 0.01; // from other studies
  k[mm] = sqrt(nmmLooseInH->Integral() / neeLooseInH->Integral());
  k_err[mm] = 0.01; // from other studies

  // R and k estimations (data)
  float R_data[2], R_data_err[2], k_data[2], k_data_err[2];
  R_data[ee] = (neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  R_data_err[ee] = R_data[ee] * quadrSum( yieldErrPoisson(neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral(),neeLooseOutHData->GetEntries()+0.5*nemLooseOutHData->GetEntries()) / (neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()),
                                          yieldErrPoisson(neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral(),neeLooseInHData->GetEntries()+0.5*nemLooseInHData->GetEntries()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));

  R_data[mm] = (nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  R_data_err[mm] = R_data[mm] * quadrSum( yieldErrPoisson(nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral(),nmmLooseOutHData->GetEntries()+0.5*nemLooseOutHData->GetEntries()) / (nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()),
                                          yieldErrPoisson(nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral(),nmmLooseInHData->GetEntries()+0.5*nemLooseInHData->GetEntries()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));

  k_data[ee] = sqrt((neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) );
  k_data_err[ee] = 0.01; // from other studies
  k_data[mm] = sqrt((nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));
  k_data_err[mm] = 0.01; // from other studies

  std::cout << "My estimation of R and k (MC): " << std::endl;
  std::cout << "Number of events in MC to estimate R (not reweighted): " << std::endl;
  std::cout << "Zee (out/in) = " <<  neeLooseOutH->GetEntries() << " / " << neeLooseInH->GetEntries() << std::endl;
  std::cout << "Zmm (out/in) = " <<  nmmLooseOutH->GetEntries() << " / " << nmmLooseInH->GetEntries() << std::endl;
  std::cout << "Ree = " << R[ee] << " +/- " <<  R_err[ee] << std::endl;
  std::cout << "Rmm = " << R[mm] << " +/- " <<  R_err[mm] << std::endl;
  std::cout << "ke = " << k[ee] << " +/- " <<  k_err[ee] << std::endl;
  std::cout << "km = " << k[mm] << " +/- " <<  k_err[mm] << std::endl;
  std::cout << "--------------------------" << std::endl;
  std::cout << "My estimation of R and k (data): " << std::endl;
  std::cout << "Number of events in data to estimate R (not reweighted): " << std::endl;
  std::cout << "Zee (out/in) = " <<  neeLooseOutHData->GetEntries() - 0.5*nemLooseOutHData->GetEntries() << " / " << neeLooseInHData->GetEntries()-0.5*nemLooseInHData->GetEntries() << std::endl;
  std::cout << "Zmm (out/in) = " <<  nmmLooseOutHData->GetEntries() - 0.5*nemLooseOutHData->GetEntries() << " / " << nmmLooseInHData->GetEntries()-0.5*nemLooseInHData->GetEntries() << std::endl;
  std::cout << "Zem (out/in) to subtract bkg = " <<  nemLooseOutHData->GetEntries() << " / " << nemLooseInHData->GetEntries() << std::endl;
  std::cout << "Ree = " << R_data[ee] << " +/- " <<  R_data_err[ee] << std::endl;
  std::cout << "Rmm = " << R_data[mm] << " +/- " <<  R_data_err[mm] << std::endl;
  std::cout << "ke = " << k_data[ee] << " +/- " <<  k_data_err[ee] << std::endl;
  std::cout << "km = " << k_data[mm] << " +/- " <<  k_data_err[mm] << std::endl;
  std::cout << "--------------------------" << std::endl;

  // DY estimation /////
  std::cout << "DY ESTIMATION..." << std::endl;
  float nZeejetsMC = ZeejetsH->Integral();
  float nZeejetsMC_err = yieldErrPoisson(nZeejetsMC,ZeejetsH->GetEntries());
  float nZmmjetsMC = ZmmjetsH->Integral();
  float nZmmjetsMC_err = yieldErrPoisson(nZmmjetsMC,ZmmjetsH->GetEntries());
  float nZemjetsMC = ZemjetsH->Integral();
  float nZemjetsMC_err = yieldErrPoisson(nZemjetsMC,ZemjetsH->GetEntries());
  float nZmejetsMC = ZmejetsH->Integral();
  float nZmejetsMC_err = yieldErrPoisson(nZmejetsMC,ZmejetsH->GetEntries());
  
  float nDYMC[4], nDYMC_err[4];
  nDYMC[ee] = ZeejetsH->Integral();
  nDYMC_err[ee] = yieldErrPoisson(nDYMC[ee],ZeejetsH->GetEntries());
  nDYMC[mm] = ZmmjetsH->Integral();
  nDYMC_err[mm] = yieldErrPoisson(nDYMC[mm],ZmmjetsH->GetEntries());
  nDYMC[em] = ZemjetsH->Integral();
  nDYMC_err[em] = yieldErrPoisson(nDYMC[em],ZemjetsH->GetEntries());
  nDYMC[me] = ZmejetsH->Integral();
  nDYMC_err[me] = yieldErrPoisson(nDYMC[me],ZmejetsH->GetEntries());

  // WW,ZZ
  float nOthers[4], nOthers_err[4];
  nOthers[ee] = OthersEEH->Integral();
  nOthers_err[ee] = yieldErrPoisson(nOthers[ee],OthersEEH->Integral());
  nOthers[mm] = OthersMMH->Integral();
  nOthers_err[mm] = yieldErrPoisson(nOthers[mm],OthersMMH->Integral());

  float neeIn = neeInH->Integral();
  float nmmIn = nmmInH->Integral();
  float nemIn = (nemInH->Integral());

  std::cout << "Events under the peak in data (no subtraction): " << std::endl;
  std::cout << "ee = " << neeIn << std::endl;
  std::cout << "mm = " << nmmIn << std::endl;
  std::cout << "em (for bkg estimation) = " << nemIn << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Events of WW,ZZ under the peak to subtract:" << std::endl;
  std::cout << "ee = " << nOthers[ee] << " +/- " << nOthers_err[ee] << " (stat)" << std::endl;
  std::cout << "mm = " << nOthers[mm] << " +/- " << nOthers_err[mm] << " (stat)" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  // ee and mm from data
  float neeExp, neeExp_err, nmmExp, nmmExp_err;
  if(useDataRk) {
    neeExp = (nDYout(neeIn, nemIn, R_data[ee], R_data_err[ee], k_data[ee], k_data_err[ee], nOthers[ee])).first;
    neeExp_err = (nDYout(neeIn, nemIn, R_data[ee], R_data_err[ee], k_data[ee], k_data_err[ee], nOthers[ee])).second;
    nmmExp = (nDYout(nmmIn, nemIn, R_data[mm], R_data_err[mm], k_data[mm], k_data_err[mm], nOthers[mm])).first;
    nmmExp_err = (nDYout(nmmIn, nemIn, R_data[mm], R_data_err[mm], k_data[mm], k_data_err[mm], nOthers[mm])).second;
  } else {
    neeExp = (nDYout(neeIn, nemIn, R[ee], R_err[ee], k[ee], k_err[ee], nOthers[ee])).first;
    neeExp_err = (nDYout(neeIn, nemIn, R[ee], R_err[ee], k[ee], k_err[ee], nOthers[ee])).second;
    nmmExp = (nDYout(nmmIn, nemIn, R[mm], R_err[mm], k[mm], k_err[mm], nOthers[mm])).first;
    nmmExp_err = (nDYout(nmmIn, nemIn, R[mm], R_err[mm], k[mm], k_err[mm], nOthers[mm])).second;
  }

  // and em from MC
  float nemmeExp = nZemjetsMC+nZmejetsMC;
  float nemmeExp_err = quadrSum(nZemjetsMC_err,nZmejetsMC_err);

  std::cout << "Number of DY->ll events at W+W- level: " << std::endl;
  std::cout << "neeIn = "  << neeIn << "\tnmmIn = "  << nmmIn << "\tnemIn = " << nemIn << std::endl;
  std::cout << "nEE MC = " << nZeejetsMC  << " +/- " << nZeejetsMC_err
            << "\tData = " << neeExp      << " +/- " << neeExp_err << std::endl;
  std::cout << "nMM MC = " << nZmmjetsMC  << " +/- " << nZmmjetsMC_err
            << "\tData = " << nmmExp      << " +/- " << nmmExp_err << std::endl;
  std::cout << "nEM+nME MC = " << nemmeExp << " +/- " << nemmeExp_err << std::endl; 

  std::cout << "END DY ESTIMATION." << std::endl;
  ////////// END DY ///////////

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK, float nZV) {
  float val = R * (nDYin - 0.5 * nemu * K - nZV);
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

