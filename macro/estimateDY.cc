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

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK, float nZV, float nZVerr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void estimateDY(float lumiInInvFb, int mass, int njets, bool useDataRk, TString addCutIn="1", TString addCutOut="1");

void estimateDYMassDependent(float lumiInInvFb, int njets, bool useDataRk) {

  std::cout << "===> ESTIMATION AT WW LEVEL <===" << std::endl;
  estimateDY(lumiInInvFb, 0,njets,useDataRk);

  // these are for limits
  ofstream cardfile[2][2]; //[cha][jetbin]
  for(int icha=0; icha<2; icha++) {
    char fileName[50];
    if(icha==ee) sprintf(fileName,"DYCard_ee_%dj.txt",njets);
    if(icha==mm) sprintf(fileName,"DYCard_mm_%dj.txt",njets);
    cardfile[icha][njets].open(fileName, ios_base::trunc);
    cardfile[icha][njets].precision(3);
  }

  ofstream cardfilell[2]; // [jetbin]
  char fileName[50];
  sprintf(fileName,"DYCard_ll_%dj.txt",njets);
  cardfilell[njets].open(fileName, ios_base::trunc);
  cardfilell[njets].precision(3);

  char nameFileTableEE[100];
  sprintf(nameFileTableEE, "DYYieldsData_ForTable_EE_%dj.txt",njets);
  ofstream tablefileEE;
  tablefileEE.open(nameFileTableEE, ios_base::trunc);
  
  tablefileEE << "\\begin{table}[p]" << endl;
  tablefileEE << "\\begin{small}" << endl;
  tablefileEE << "\\begin{center}" << endl;
  tablefileEE << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << endl;
  tablefileEE << "\\hline" << endl;
  tablefileEE << "$m_{H}$ [GeV] \t\t \t & $n^{data}_{in}$ \t\t & \t\t $n^{data e\\mu}_{in}$ \t\t & \t $n^{data}_{in}(sub)$ \t & \t $n^{MC}_{VZ}$ \t & \t $\\sigma_k$ \t & \t $R_{data}$ \t & \t $n^{data}_{DY}$ \t & \t $n^{MC}_{DY}$ \\\\" << endl;
  tablefileEE << "\\hline" << endl;

  char nameFileTableMM[100];
  sprintf(nameFileTableMM, "DYYieldsData_ForTable_MM_%dj.txt",njets);
  ofstream tablefileMM;
  tablefileMM.open(nameFileTableMM, ios_base::trunc);

  tablefileMM << "\\begin{table}[p]" << endl;
  tablefileMM << "\\begin{small}" << endl;
  tablefileMM << "\\begin{center}" << endl;
  tablefileMM << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  tablefileMM << "\\hline" << endl;
  tablefileMM << "$m_{H}$ [GeV] \t & $n^{data}_{in}$ \t & \t $n^{MC}_{VZ}$ \t & \t $R_{data}$ \t & \t $n^{data}_{DY}$ \t & \t $n^{MC}_{DY}$ \\\\" << endl;
  tablefileMM << "\\hline" << endl;


  char nameFileTableLL[100];
  sprintf(nameFileTableLL, "DYYieldsData_ForTable_LL_%dj.txt",njets);
  ofstream tablefileLL;
  tablefileLL.open(nameFileTableLL, ios_base::trunc);

  tablefileLL << "\\begin{table}[p]" << endl;
  tablefileLL << "\\begin{small}" << endl;
  tablefileLL << "\\begin{center}" << endl;
  tablefileLL << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  tablefileLL << "\\hline" << endl;
  tablefileLL << "$m_{H}$ [GeV] \t & $n^{data}_{in}$ \t & \t $n^{MC}_{VZ}$ \t & \t $R_{data}$ \t & \t $n^{data}_{DY}$ \t & \t $n^{MC}_{DY}$ \\\\" << endl;
  tablefileLL << "\\hline" << endl;

  int mH[18] = {115,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  for(int i=0; i<18;i++) {
    std::cout << "mH = " << mH[i] << std::endl;
    TString addCutIn = higgsCuts(mH[i],false);
    TString addCutOut = higgsCuts(mH[i],true);
    estimateDY(lumiInInvFb, mH[i],njets,useDataRk,addCutIn,addCutOut);
  }

  ofstream tablefileEEEnd;
  tablefileEEEnd.open(nameFileTableEE, ios_base::app);

  tablefileEEEnd << "\\hline" << endl;
  tablefileEEEnd << "\\end{tabular}" << endl;
  tablefileEEEnd << "\\caption{DY $\\to e^+e^-$ yields at Higgs selection level from data control sample in the " << njets << " jet bin.}" << std::endl;
  tablefileEEEnd << "\\end{center}" << endl;
  tablefileEEEnd << "\\end{small}" << endl;
  tablefileEEEnd << "\\end{table}" << endl;

  ofstream tablefileMMEnd;
  tablefileMMEnd.open(nameFileTableMM, ios_base::app);    
  tablefileMMEnd << "\\hline" << endl;
  tablefileMMEnd << "\\end{tabular}" << endl;
  tablefileMMEnd << "\\caption{DY $\\to \\mu^+\\mu^-$ yields at Higgs selection level from data control sample in the " << njets << " jet bin.}" << std::endl;
  tablefileMMEnd << "\\end{center}" << endl;
  tablefileMMEnd << "\\end{small}" << endl;
  tablefileMMEnd << "\\end{table}" << endl;

  ofstream tablefileLLEnd;
  tablefileLLEnd.open(nameFileTableLL, ios_base::app);    
  tablefileLLEnd << "\\hline" << endl;
  tablefileLLEnd << "\\end{tabular}" << endl;
  tablefileLLEnd << "\\caption{DY $\\to \\ell^+\\ell^-$ yields at Higgs selection level from data control sample in the " << njets << " jet bin.}" << std::endl;
  tablefileLLEnd << "\\end{center}" << endl;
  tablefileLLEnd << "\\end{small}" << endl;
  tablefileLLEnd << "\\end{table}" << endl;

  std::cout << "DONE." << std::endl;
}

void estimateDY(float lumiInInvFb, int mass, int njets, bool useDataRk, TString addCutIn, TString addCutOut) {

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

  // input files and trees
  TFile *fileData   = TFile::Open("results_data/datasets_trees/dataset_ll.root");
  TFile *fileZjets  = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileOthers = TFile::Open("results/datasets_trees/others_ll.root");

  TTree *treeData   = (TTree*)fileData->Get("latino");
  TTree *treeZjets  = (TTree*)fileZjets->Get("latino");
  TTree *treeOthers = (TTree*)fileOthers->Get("latino");

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

  // to decide the cut for R
  TH1F *neeCheckInH = new TH1F("neeCheckInH","",50,0,180);
  TH1F *nmmCheckInH = new TH1F("nmmCheckInH","",50,0,180);

  // to estimate R and k from MC
  TH1F *neeLooseInH  = new TH1F("neeLooseInH","",50,0,180);
  TH1F *neeLooseOutH = new TH1F("neeLooseOutH","",50,0,180);
  TH1F *nmmLooseInH  = new TH1F("nmmLooseInH","",50,0,180);
  TH1F *nmmLooseOutH = new TH1F("nmmLooseOutH","",50,0,180);

  // to estimate R and k from data
  TH1F *neeLooseInHData  = new TH1F("neeLooseInHData","",50,0,180);
  TH1F *neeLooseOutHData = new TH1F("neeLooseOutHData","",50,0,180);
  TH1F *nmmLooseInHData  = new TH1F("nmmLooseInHData","",50,0,180);
  TH1F *nmmLooseOutHData = new TH1F("nmmLooseOutHData","",50,0,180);  
  TH1F *nemLooseInHData  = new TH1F("nemLooseInHData","",50,0,180);
  TH1F *nemLooseOutHData = new TH1F("nemLooseOutHData","",50,0,180);  

  treeZjets->Project("ZeejetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("ZmmjetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && channel==0)*baseW*puW*effW")).Data());
  treeZjets->Project("ZemjetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && channel==2)*baseW*puW*effW")).Data());
  treeZjets->Project("ZmejetsH","dphill",(TString("(")+TString(wwselcut)+TString(" && ")+addCutOut+TString(" && channel==3)*baseW*puW*effW")).Data());

  // to choose between the last and the last-1 bin for R computation
  treeZjets->Project("neeCheckInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("nmmCheckInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==0)*baseW*puW*effW")).Data());
  int forTheCheck = neeCheckInH->GetEntries() + nmmCheckInH->GetEntries();
  char mpmetcut[30];
  if(forTheCheck<50)  { std::cout << "less than 50 events: cutting at MET>30" << std::endl; sprintf(mpmetcut,"mpmet>30"); }
  if(forTheCheck>=50) { std::cout << "more than 50 events: cutting at MET>40" << std::endl; sprintf(mpmetcut,"mpmet>40"); }

  // nominal R computed with cut at 40 or 30 (according to the statistics at denominator) - MC
  treeZjets->Project("neeLooseInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("neeLooseOutH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && zveto && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("nmmLooseInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==0)*baseW*puW*effW")).Data());
  treeZjets->Project("nmmLooseOutH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && zveto && channel==0)*baseW*puW*effW")).Data());

  // nominal R computed with cut at 40 or 30 (according to the statistics at denominator) - data
  treeData->Project("neeLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && !zveto && channel==1)")).Data());
  treeData->Project("neeLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && zveto && channel==1)")).Data());
  treeData->Project("nmmLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && !zveto && channel==0)")).Data());
  treeData->Project("nmmLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && zveto && channel==0)")).Data());
  treeData->Project("nemLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && !zveto && (channel==2 || channel==3))")).Data());
  treeData->Project("nemLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && ")+TString(mpmetcut)+TString(" && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && mll>12 && zveto && (channel==2 || channel==3))")).Data());
  
  // contribution of dibosons under the peak (from MC) with nominal selection
  treeOthers->Project("OthersEEH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>40 &&")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && !zveto && channel==1)*baseW*puW*effW")).Data());
  treeOthers->Project("OthersMMH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>40 &&")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && !zveto && channel==0)*baseW*puW*effW")).Data());
  
  // number of events measured in data under the peak (control region with nominal selection)
  treeData->Project("neeInH","dphill",(TString("mll>12 && finalLeptons && pfmet>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && !zveto && channel==1")).Data()); 
  treeData->Project("nmmInH","dphill",(TString("mll>12 && finalLeptons && pfmet>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && !zveto && channel==0")).Data());
  treeData->Project("nemInH","dphill",(TString("mll>12 && finalLeptons && pfmet>20 && mpmet>40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString(" && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0  && bveto && nSoftMu==0 && !zveto && (channel==2 || channel==3)")).Data());

  // R and k estimations (from MC)
  float R[3], R_err[3], k[3], k_err[3];
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

  cout << "from loose analysis (cut 40 or 30) => num EE = " << neeLooseOutH->Integral()  << ", denom EE = " << neeLooseInH->Integral() << endl;
  cout << "from loose analysis (cut 40 or 30) => num MM = " << nmmLooseOutH->Integral()  << ", denom MM = " << nmmLooseInH->Integral() << endl;
  cout << "from loose analysis (cut 40 or 30) => Ree = "    << R[ee]  << ", Rmm = " << R[mm] << endl;

  // combine ee+mm because of low stat (MC) =>
  // R combined
  float nllLooseOut     = neeLooseOutH->Integral() + nmmLooseOutH->Integral();
  float nllLooseOut_err = yieldErrPoisson(nllLooseOut, neeLooseOutH->GetEntries()+nmmLooseOutH->GetEntries());
  float nllLooseIn      = neeLooseInH->Integral() + nmmLooseInH->Integral();
  float nllLooseIn_err  = yieldErrPoisson(nllLooseIn, neeLooseInH->GetEntries()+nmmLooseInH->GetEntries());
  float Rll             = nllLooseOut / nllLooseIn;
  float Rll_err         = Rll * quadrSum(nllLooseOut_err/nllLooseOut, nllLooseIn_err/nllLooseIn);

  // k combined 
  float kll     = 0.5 * (k[mm]+1.0/k[mm]);
  float kll_err = k_err[mm] * (1. - 1./k[mm]/k[mm]);
  
  

  // R and k estimations (data)
  float R_data[2], R_data_err[2], k_data[2], k_data_err[2];
  R_data[ee]     = (neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  R_data_err[ee] = R_data[ee] * quadrSum( yieldErrPoisson(neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral(),neeLooseOutHData->GetEntries()+0.5*nemLooseOutHData->GetEntries()) / (neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()),
                                          yieldErrPoisson(neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral(),neeLooseInHData->GetEntries()+0.5*nemLooseInHData->GetEntries()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));

  R_data[mm]     = (nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  R_data_err[mm] = R_data[mm] * quadrSum( yieldErrPoisson(nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral(),nmmLooseOutHData->GetEntries()+0.5*nemLooseOutHData->GetEntries()) / (nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()),
                                          yieldErrPoisson(nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral(),nmmLooseInHData->GetEntries()+0.5*nemLooseInHData->GetEntries()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));

  k_data[ee]     = sqrt((neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) );
  k_data_err[ee] = 0.01; // from other studies
  k_data[mm]     = sqrt((nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral()));
  k_data_err[mm] = 0.01; // from other studies


  // combine ee+mm because of low stat (from data)
  // R
  float nllLooseOutData     = neeLooseOutHData->Integral() + nmmLooseOutHData->Integral();
  float nllLooseOutData_err = yieldErrPoisson(nllLooseOutData, neeLooseOutHData->GetEntries()+nmmLooseOutHData->GetEntries());
  float numR     = nllLooseOutData - nemLooseOutHData->Integral();
  float numR_err = quadrSum(nllLooseOutData_err,nemLooseOutHData->Integral());

  float nllLooseInData     = neeLooseInHData->Integral() + nmmLooseInHData->Integral();
  float nllLooseInData_err = yieldErrPoisson(nllLooseInData, neeLooseInHData->GetEntries()+nmmLooseInHData->GetEntries());
  float denomR     = nllLooseInData - nemLooseInHData->Integral();
  float denomR_err = quadrSum(nllLooseInData_err,nemLooseInHData->Integral());

  float Rll_data     = numR / denomR;
  float Rll_data_err = Rll_data * quadrSum(numR_err/numR, denomR_err/denomR);

  // k combined is 
  //   float kll_data = 0.5*(k_data[mm]+1.0/k_data[mm]);
  //   float kll_data_err = k_data_err[mm] * (1. - 1./k_data[mm]/k_data[mm]);


  // for the systematics: mpMET>20 in MC
  treeZjets->Project("neeLooseInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("nmmLooseInH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && !zveto && channel==0)*baseW*puW*effW")).Data());
  treeZjets->Project("neeLooseOutH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && zveto && channel==1)*baseW*puW*effW")).Data());
  treeZjets->Project("nmmLooseOutH","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && nextra==0 && bveto && nSoftMu==0 && zveto && channel==0)*baseW*puW*effW")).Data());

  float R_interm[2];
  R_interm[ee] = (neeLooseOutH->Integral()) / (neeLooseInH->Integral());
  R_interm[mm] = (nmmLooseOutH->Integral()) / (nmmLooseInH->Integral());
  
  cout << "from interm: cut 20 => num EE = " << neeLooseOutH->Integral()  << ", denom EE = " << neeLooseInH->Integral() << endl;
  cout << "from interm: cut 20 => num MM = " << nmmLooseOutH->Integral()  << ", denom MM = " << nmmLooseInH->Integral() << endl;
  cout << "from interm: cut 20 => Ree = "    << R_interm[ee]  << ", Rmm = " << R_interm[mm] << endl;
  
  float Rll_interm = (neeLooseOutH->Integral() + nmmLooseOutH->Integral()) / (neeLooseInH->Integral() + nmmLooseInH->Integral());

  float R_systerr[2], R_toterr[2], Rll_systerr, Rll_toterr;
  for(int icha=0; icha<2; icha++) {
    R_systerr[icha] = fabs(R[icha] - R_interm[icha]);
    R_toterr[icha]  = quadrSum(R_err[icha],R_systerr[icha]);
  }
  Rll_systerr = fabs(Rll - Rll_interm);
  Rll_toterr  = quadrSum(Rll_err,Rll_systerr);

  // for the systematics: mpMET>20 in data
  treeData->Project("neeLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && !zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==1)")).Data());
  treeData->Project("neeLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==1)")).Data());
  treeData->Project("nmmLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && !zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==0)")).Data());
  treeData->Project("nmmLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && channel==0)")).Data());
  treeData->Project("nemLooseInHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutIn+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && !zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && (channel==2 || channel==3))")).Data());
  treeData->Project("nemLooseOutHData","dphill",(TString("(mll>12 && finalLeptons && pfmet>20 && mpmet>20 && mpmet<40 && ")+addCutOut+TString(" && ")+TString(njcut)+TString( " && bveto && nSoftMu==0 && mll>12 && zveto && ((jetpt1>15 && abs(dphilljet)<165) || jetpt1<=15) && (channel==2 || channel==3))")).Data());

  float R_data_interm[2];
  R_data_interm[ee] = (neeLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (neeLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  R_data_interm[mm] = (nmmLooseOutHData->Integral()-0.5*nemLooseOutHData->Integral()) / (nmmLooseInHData->Integral()-0.5*nemLooseInHData->Integral());
  float Rll_data_interm = (neeLooseOutHData->Integral() + nmmLooseOutHData->Integral() - nemLooseOutHData->Integral()) /
    (neeLooseInHData->Integral() + nmmLooseInHData->Integral() - nemLooseInHData->Integral());

  float R_data_systerr[2], R_data_toterr[2], Rll_data_systerr, Rll_data_toterr;
  for(int icha=0; icha<2; icha++) {
    R_data_systerr[icha] = fabs(R_data[icha] - R_data_interm[icha]);
    R_data_toterr[icha]  = quadrSum(R_data_err[icha],R_data_systerr[icha]);
  }
  Rll_data_systerr = fabs(Rll_data - Rll_data_interm);
  Rll_data_toterr  = quadrSum(Rll_data_err,Rll_data_systerr);


  // printouts
  if(addCutIn==TString("1")) {
    std::cout << "My estimation of R and k (MC): " << std::endl;
    std::cout << "Number of events in MC to estimate R (not reweighted): " << std::endl;
    std::cout << "Zee (out/in) = " <<  neeLooseOutH->GetEntries() << " / " << neeLooseInH->GetEntries() << std::endl;
    std::cout << "Zmm (out/in) = " <<  nmmLooseOutH->GetEntries() << " / " << nmmLooseInH->GetEntries() << std::endl;
    std::cout << "Ree = " << R[ee] << " +/- " <<  R_err[ee] << " +/- " << R_systerr[ee] << std::endl;
    std::cout << "Rmm = " << R[mm] << " +/- " <<  R_err[mm] << " +/- " << R_systerr[mm] << std::endl;
    std::cout << "ke = "  << k[ee] << " +/- " <<  k_err[ee] << std::endl;
    std::cout << "km = "  << k[mm] << " +/- " <<  k_err[mm] << std::endl;
    std::cout << "Combined ee+mm R estimation: " << std::endl;
    std::cout << "Rll = " << Rll << " +/- " << Rll_err << " (stat.) +/- " << Rll_systerr << " (syst.) " << std::endl;
    std::cout << "--------------------------" << std::endl;
    std::cout << "My estimation of R and k (data): " << std::endl;
    std::cout << "Number of events in data to estimate R (not reweighted): " << std::endl;
    std::cout << "Zee (out/in) = " <<  neeLooseOutHData->GetEntries() - 0.5*nemLooseOutHData->GetEntries() << " / " << neeLooseInHData->GetEntries()-0.5*nemLooseInHData->GetEntries() << std::endl;
    std::cout << "Zmm (out/in) = " <<  nmmLooseOutHData->GetEntries() - 0.5*nemLooseOutHData->GetEntries() << " / " << nmmLooseInHData->GetEntries()-0.5*nemLooseInHData->GetEntries() << std::endl;
    std::cout << "Zem (out/in) to subtract bkg = " <<  nemLooseOutHData->GetEntries() << " / " << nemLooseInHData->GetEntries() << std::endl;
    std::cout << "Ree = " << R_data[ee] << " +/- " <<  R_data_err[ee] << " (stat.) +/- " << R_data_systerr[ee] << " (syst.) " << std::endl;
    std::cout << "Rmm = " << R_data[mm] << " +/- " <<  R_data_err[mm] << " (stat.) +/- " << R_data_systerr[ee] << " (syst.) " << std::endl;
    std::cout << "ke = "  << k_data[ee] << " +/- " <<  k_data_err[ee] << std::endl;
    std::cout << "km = "  << k_data[mm] << " +/- " <<  k_data_err[mm] << std::endl;
    std::cout << "Combined ee+mm R estimation: " << std::endl;
    std::cout << "Rll = " << Rll_data << " +/- " << Rll_data_err << " (stat.) +/- " << Rll_data_systerr << " (syst.) " << std::endl;
    std::cout << "--------------------------" << std::endl;
  }


  // --------------------------------------------------------------------------------------
  // DY estimation 
  float nZeejetsMC     = lumiInInvFb * ZeejetsH->Integral();
  float nZeejetsMC_err = lumiInInvFb * yieldErrPoisson(nZeejetsMC,ZeejetsH->GetEntries());
  float nZmmjetsMC     = lumiInInvFb * ZmmjetsH->Integral();
  float nZmmjetsMC_err = lumiInInvFb * yieldErrPoisson(nZmmjetsMC,ZmmjetsH->GetEntries());
  float nZemjetsMC     = lumiInInvFb * ZemjetsH->Integral();
  float nZemjetsMC_err = lumiInInvFb * yieldErrPoisson(nZemjetsMC,ZemjetsH->GetEntries());
  float nZmejetsMC     = lumiInInvFb * ZmejetsH->Integral();
  float nZmejetsMC_err = lumiInInvFb * yieldErrPoisson(nZmejetsMC,ZmejetsH->GetEntries());
  float nZlljetsMC     = lumiInInvFb * (ZeejetsH->Integral() + ZmmjetsH->Integral());
  float nZlljetsMC_err = quadrSum(nZeejetsMC_err, nZmmjetsMC_err);

  // WW,ZZ (assuming 10% error on xsec)
  float nOthers[4], nOthers_err[4], nOthers_systerr[4], nOthers_toterr[4];
  nOthers[ee]     = lumiInInvFb * OthersEEH->Integral();
  nOthers_err[ee] = lumiInInvFb * yieldErrPoisson(nOthers[ee],OthersEEH->GetEntries());
  nOthers_systerr[ee] = 0.10 * nOthers[ee];
  nOthers_toterr[ee]  = quadrSum(nOthers_err[ee], nOthers_systerr[ee]);
  nOthers[mm]     = lumiInInvFb * OthersMMH->Integral();
  nOthers_err[mm] = lumiInInvFb * yieldErrPoisson(nOthers[mm],OthersMMH->GetEntries());
  nOthers_systerr[mm] = 0.10 * nOthers[mm];
  nOthers_toterr[mm]  = quadrSum(nOthers_err[mm], nOthers_systerr[mm]);
  float nOthersll     = nOthers[ee] + nOthers[mm];
  float nOthersll_err = quadrSum(nOthers_err[ee],nOthers_err[mm]);
  float nOthersll_systerr = 0.10 * nOthersll;
  float nOthersll_toterr  = quadrSum(nOthersll_err,nOthersll_systerr);

  // under the peak in data
  float neeIn = neeInH->Integral();
  float nmmIn = nmmInH->Integral();
  float nemIn = nemInH->Integral();
  float nllIn = neeIn + nmmIn;

  if(addCutIn==TString("1")) {
    std::cout << "DY ESTIMATION..." << std::endl;
    std::cout << "Events under the peak in data (no subtraction): " << std::endl;
    std::cout << "ee = " << neeIn << std::endl;
    std::cout << "mm = " << nmmIn << std::endl;
    std::cout << "em (for bkg estimation) = " << nemIn << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Events of WW,ZZ under the peak to subtract:" << std::endl;
    std::cout << "ee = " << nOthers[ee] << " +/- " << nOthers_err[ee] << " (stat)" << " +/- " << nOthers_systerr[ee] << std::endl;
    std::cout << "mm = " << nOthers[mm] << " +/- " << nOthers_err[mm] << " (stat)" << " +/- " << nOthers_systerr[mm] << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
  }

  // ee and mm from data or from MC
  float neeExp, neeExp_err, nmmExp, nmmExp_err;
  float nllExp, nllExp_err;
  if(useDataRk) {  
    neeExp     = (nDYout(neeIn, nemIn, R_data[ee], R_data_toterr[ee], k_data[ee], k_data_err[ee], nOthers[ee], nOthers_toterr[ee])).first;
    neeExp_err = (nDYout(neeIn, nemIn, R_data[ee], R_data_toterr[ee], k_data[ee], k_data_err[ee], nOthers[ee], nOthers_toterr[ee])).second;
    nmmExp     = (nDYout(nmmIn, nemIn, R_data[mm], R_data_toterr[mm], k_data[mm], k_data_err[mm], nOthers[mm], nOthers_toterr[mm])).first;
    nmmExp_err = (nDYout(nmmIn, nemIn, R_data[mm], R_data_toterr[mm], k_data[mm], k_data_err[mm], nOthers[mm], nOthers_toterr[mm])).second;
  } else {
    neeExp     = (nDYout(neeIn, nemIn, R[ee], R_err[ee], k[ee], k_err[ee], nOthers[ee], nOthers_toterr[ee])).first;
    neeExp_err = (nDYout(neeIn, nemIn, R[ee], R_err[ee], k[ee], k_err[ee], nOthers[ee], nOthers_toterr[ee])).second;
    nmmExp     = (nDYout(nmmIn, nemIn, R[mm], R_err[mm], k[mm], k_err[mm], nOthers[mm], nOthers_toterr[mm])).first;
    nmmExp_err = (nDYout(nmmIn, nemIn, R[mm], R_err[mm], k[mm], k_err[mm], nOthers[mm], nOthers_toterr[mm])).second;
  }
  // in any case for the moment using R from MC, as smurfs => to be (maybe) fixed with more data
  nllExp       = (nDYout(nllIn, 2*nemIn, Rll, Rll_toterr, kll, kll_err, nOthersll, nOthersll_toterr)).first;
  nllExp_err   = (nDYout(nllIn, 2*nemIn, Rll, Rll_toterr, kll, kll_err, nOthersll, nOthersll_toterr)).second;

  // and em from MC
  float nemmeExp     = nZemjetsMC+nZmejetsMC;
  float nemmeExp_err = quadrSum(nZemjetsMC_err,nZmejetsMC_err);


  if(addCutIn==TString("1")) {
    std::cout << "Number of DY->ll events at W+W- level: " << std::endl;
    std::cout << "neeIn = "  << neeIn       << "\tnmmIn = "  << nmmIn << "\tnemIn = " << nemIn << std::endl;
    std::cout << "nEE MC = " << nZeejetsMC  << " +/- " << nZeejetsMC_err
              << "\tData = " << neeExp      << " +/- " << neeExp_err << std::endl;
    std::cout << "nMM MC = " << nZmmjetsMC  << " +/- " << nZmmjetsMC_err
              << "\tData = " << nmmExp      << " +/- " << nmmExp_err << std::endl;
    std::cout << "nEM+nME MC = " << nemmeExp    << " +/- " << nemmeExp_err   << std::endl; 
    std::cout << "nLL MC = "     << nZlljetsMC  << " +/- " << nZlljetsMC_err
              << "\tData = " << nllExp      << " +/- " << nllExp_err << std::endl;
    std::cout << "END DY ESTIMATION." << std::endl;
    ////////// END DY ///////////
  } else {
    
    char nameFileTableEE[100];
    sprintf(nameFileTableEE, "DYYieldsData_ForTable_EE_%dj.txt",njets);
    ofstream tablefileEE;
    tablefileEE.open(nameFileTableEE, ios_base::app);
    tablefileEE.precision(2);

    // short printout
    tablefileEE << mass << "\t&\t" << neeIn << "\t&\t" << nemIn << "\t&\t" << nOthers[ee] 
                << "\t&\t" << neeIn - 0.5*nemIn * k[ee] - nOthers[ee]
                << "\t&\t" << k[ee] << " $\\pm$ " << k_err[ee] << " $\\pm$ "
                << "\t&\t" << R[ee] << " $\\pm$ " << R_err[ee] << " $\\pm$ " << R_data_systerr[ee]
                << "\t&\t" <<  neeExp << " $\\pm$ " << neeExp_err 
                << "\t&\t" << nZeejetsMC << " $\\pm$ " << nZeejetsMC_err 
                << "\t \\\\" << std::endl;

    char nameFileTableMM[100];
    sprintf(nameFileTableMM, "DYYieldsData_ForTable_MM_%dj.txt",njets);
    ofstream tablefileMM;
    tablefileMM.open(nameFileTableMM, ios_base::app);
    tablefileMM.precision(2);

    tablefileMM << mass << "\t&\t" << nmmIn << "\t&\t" << nemIn << "\t&\t" << nOthers[mm] 
                << "\t&\t" << nmmIn - 0.5*nemIn * k[mm] - nOthers[mm]
                << "\t&\t" << k[mm] << " $\\pm$ " << k_err[mm]
                << "\t&\t" << R_data[mm] << " $\\pm$ " << R_data_err[mm] << " $\\pm$ " << R_data_systerr[mm] 
                << "\t&\t" <<  nmmExp << " $\\pm$ " << nmmExp_err 
                << "\t&\t" << nZmmjetsMC << " $\\pm$ " << nZmmjetsMC_err 
                << "\t \\\\" << std::endl;


    char nameFileTableLL[100];
    sprintf(nameFileTableLL, "DYYieldsData_ForTable_LL_%dj.txt",njets);
    ofstream tablefileLL;
    tablefileLL.open(nameFileTableLL, ios_base::app);
    tablefileLL.precision(2);

    tablefileLL << mass << "\t&\t" << nllIn << "\t&\t" << nemIn << "\t&\t" << nOthersll
                << "\t&\t" << nllIn - nemIn * kll - nOthersll 
                << "\t&\t" << kll << " $\\pm$ " << kll_err
                << "\t&\t" << Rll << " $\\pm$ " << Rll_err << " $\\pm$ " << Rll_systerr 
                << "\t&\t" <<  nllExp << " $\\pm$ " << nllExp_err 
                << "\t&\t" << nZlljetsMC << " $\\pm$ " << nZlljetsMC_err 
                << "\t \\\\" << std::endl;

    std::cout << "mode\t" << "nIn\t" << "n(OF)\t" << "n(VZ)\t" << "k\t\t" << "R\t" << "n(DY) meas." << std::endl;
    std::cout << "ll:\t" << nllIn << "\t" << nemIn << "\t" << nOthersll 
              << "\t" << kll << " +/- " << kll_err
              << "\t" << Rll << " +/- " << Rll_err << " +/- " << Rll_systerr 
              << "\t" <<  nllExp << " +/- " << nllExp_err 
              << std::endl;



    // for the datacards
    float alpha[2], alpha_err[2];
    float alphall, alphall_err;
    alpha[ee]     = neeExp / neeIn;
    alpha_err[ee] = alpha[ee] * neeExp_err/neeExp;
    
    alpha[mm]     = nmmExp / nmmIn;
    alpha_err[mm] = alpha[mm] * nmmExp_err/nmmExp;

    alphall     = nllExp / nllIn;
    alphall_err = alphall * nllExp_err/nllExp;

    ofstream cardfile;
    char fileName[50];

    sprintf(fileName,"DYCard_ee_%dj.txt",njets);
    cardfile.open(fileName, ios_base::app);
    cardfile << mass
             << "\t" << neeIn << "\t" << alpha[ee] << "\t" << alpha_err[ee]
             << std::endl;
    cardfile.close();

    sprintf(fileName,"DYCard_mm_%dj.txt",njets);
    cardfile.open(fileName, ios_base::app);
    cardfile << mass
             << "\t" << nmmIn << "\t" << alpha[mm] << "\t" << alpha_err[mm]
             << std::endl;
    cardfile.close();

    sprintf(fileName,"DYCard_ll_%dj.txt",njets);
    cardfile.open(fileName, ios_base::app);
    cardfile << mass
             << "\t" << nllIn << "\t" << alphall << "\t" << alphall_err
             << std::endl;
    cardfile.close();
  }
}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK, float nZV, float nZVerr) {
  float val = R * (nDYin - 0.5 * nemu * K - nZV);
  float err = quadrSum(sigmaR*(nDYin - 0.5 * nemu * K -nZV),
                       R * sqrt(nDYin),
                       0.5 * R * K * sqrt(nemu),
                       0.5 * R * nemu * sigmaK,
                       R * nZVerr); 
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

