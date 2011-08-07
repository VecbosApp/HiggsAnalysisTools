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
#include "massDependentCuts.cc"

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK, float nZV);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void estimateDY(int mass, int njets, bool useDataRk, TString addCutIn="1", TString addCutOut="1");

void estimateDYMassDependent(int njets, bool useDataRk) {

  std::cout << "===> ESTIMATION AT WW LEVEL <===" << std::endl;
  estimateDY(0,njets,useDataRk);

  // these are for limits
  ofstream cardfile[2][2]; //[cha][jetbin]
  for(int icha=0; icha<2; icha++) {
    for(int j=0; j<2; j++) {
      char fileName[50];
      if(icha==ee) sprintf(fileName,"DYCard_ee_%dj.txt",j);
      if(icha==mm) sprintf(fileName,"DYCard_mm_%dj.txt",j);
      cardfile[icha][j].open(fileName, ios_base::trunc);
      cardfile[icha][j].precision(3);
    }
  }

  int mH[18] = {115,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  for(int i=0; i<18;i++) {
    std::cout << "mH = " << mH[i] << std::endl;
    TString addCutIn = higgsCuts(mH[i],false);
    TString addCutOut = higgsCuts(mH[i],true);
    estimateDY(mH[i],njets,useDataRk,addCutIn,addCutOut);
  }
  std::cout << "DONE." << std::endl;
  
}

void estimateDY(int mass, int njets, bool useDataRk, TString addCutIn, TString addCutOut) {

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

  TFile *fileData  = TFile::Open("results_data/datasets_trees/dataset_ll.root");
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

  treeZjets->Project("ZeejetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==1)*baseW*puW")).Data());
  treeZjets->Project("ZmmjetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==0)*baseW*puW")).Data());
  treeZjets->Project("ZemjetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==2)*baseW*puW")).Data());
  treeZjets->Project("ZmejetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==3)*baseW*puW")).Data());

  treeZjets->Project("neeLooseInH","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && !zveto && channel==1)*baseW*puW")).Data());
  treeZjets->Project("neeLooseOutH","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && zveto && channel==1)*baseW*puW")).Data());
  treeZjets->Project("nmmLooseInH","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && !zveto && channel==0)*baseW*puW")).Data());
  treeZjets->Project("nmmLooseOutH","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && zveto && channel==0)*baseW*puW")).Data());

  treeData->Project("neeLooseInHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && !zveto && channel==1)")).Data());
  treeData->Project("neeLooseOutHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && zveto && channel==1)")).Data());
  treeData->Project("nmmLooseInHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && !zveto && channel==0)")).Data());
  treeData->Project("nmmLooseOutHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && zveto && channel==0)")).Data());
  treeData->Project("nemLooseInHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && !zveto && (channel==2 || channel==3))")).Data());
  treeData->Project("nemLooseOutHData","dphill",(TString("(finalLeptons && met>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && mll>12 && zveto && (channel==2 || channel==3))")).Data());
  
  // contribution under the peak
  treeOthers->Project("OthersEEH","dphill",(TString("(finalLeptons && met>20 && mpmet>40 &&")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && bveto && !zveto && channel==1)*baseW*puW")).Data());
  treeOthers->Project("OthersMMH","dphill",(TString("(finalLeptons && met>20 && mpmet>40 &&")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && bveto && !zveto && channel==0)*baseW*puW")).Data());
  
  treeData->Project("neeInH","dphill",(TString("mll>12 && finalLeptons && met>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && bveto && !zveto && channel==1")).Data()); 
  treeData->Project("nmmInH","dphill",(TString("mll>12 && finalLeptons && met>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && bveto && !zveto && channel==0")).Data());
  treeData->Project("nemInH","dphill",(TString("mll>12 && finalLeptons && met>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && bveto && !zveto && (channel==2 || channel==3)")).Data());

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

  if(addCutIn==TString("1")) {
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
  }
  // DY estimation /////
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

  if(addCutIn==TString("1")) {
    std::cout << "DY ESTIMATION..." << std::endl;
    std::cout << "Events under the peak in data (no subtraction): " << std::endl;
    std::cout << "ee = " << neeIn << std::endl;
    std::cout << "mm = " << nmmIn << std::endl;
    std::cout << "em (for bkg estimation) = " << nemIn << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Events of WW,ZZ under the peak to subtract:" << std::endl;
    std::cout << "ee = " << nOthers[ee] << " +/- " << nOthers_err[ee] << " (stat)" << std::endl;
    std::cout << "mm = " << nOthers[mm] << " +/- " << nOthers_err[mm] << " (stat)" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
  }

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

  if(addCutIn==TString("1")) {
    std::cout << "Number of DY->ll events at W+W- level: " << std::endl;
    std::cout << "neeIn = "  << neeIn << "\tnmmIn = "  << nmmIn << "\tnemIn = " << nemIn << std::endl;
    std::cout << "nEE MC = " << nZeejetsMC  << " +/- " << nZeejetsMC_err
              << "\tData = " << neeExp      << " +/- " << neeExp_err << std::endl;
    std::cout << "nMM MC = " << nZmmjetsMC  << " +/- " << nZmmjetsMC_err
              << "\tData = " << nmmExp      << " +/- " << nmmExp_err << std::endl;
    std::cout << "nEM+nME MC = " << nemmeExp << " +/- " << nemmeExp_err << std::endl; 
    
    std::cout << "END DY ESTIMATION." << std::endl;
    ////////// END DY ///////////
  } else {
    // short printout
    std::cout << "mode\t" << "nIn\t" << "n(OF)\t" << "R\t" << "n(DY) meas." << std::endl;
    std::cout << "ee:\t" << neeIn << "\t" << nOthers[ee] << "\t" << R_data[ee] << " +/- " << R_data_err[ee] << "\t" <<  neeExp << " +/- " << neeExp_err << std::endl;
    std::cout << "mm:\t" << nmmIn << "\t" << nOthers[mm] << "\t" << R_data[mm] << " +/- " << R_data_err[mm] << "\t" <<  nmmExp << " +/- " << nmmExp_err << std::endl;

    // for the datacards
    float alpha[2], alpha_err[2];
    alpha[ee] = neeExp / neeIn;
    alpha_err[ee] = alpha[ee] * neeExp_err/neeExp;

    alpha[mm] = nmmExp / nmmIn;
    alpha_err[mm] = alpha[mm] * nmmExp_err/nmmExp;

    ofstream cardfile;
    char fileName[50];

    sprintf(fileName,"DYCard_ee_%dj.txt",njets);
    cardfile.open(fileName, ios_base::app);
    cardfile << mass
             << "\t" << neeIn << "\t" << alpha[ee]
             << "\t" << alpha_err[ee]
             << std::endl;
    cardfile.close();

    sprintf(fileName,"DYCard_mm_%dj.txt",njets);
    cardfile.open(fileName, ios_base::app);
    cardfile << mass
             << "\t" << nmmIn << "\t" << alpha[mm]
             << "\t" << alpha_err[mm]
             << std::endl;
    cardfile.close();
  }

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

