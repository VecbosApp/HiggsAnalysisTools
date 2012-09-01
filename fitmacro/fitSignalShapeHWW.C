#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "Math/MinimizerOptions.h"

using namespace RooFit;



#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

enum channels { of0j, of1j, sf0j, sf1j };
  
string getStringChannel(int channel) {
  if(channel==of0j) return string("channel>=2 && njet==0");
  if(channel==of1j) return string("channel>=2 && njet==1");
  if(channel==sf0j) return string("channel<2 && njet==0");
  if(channel==sf1j) return string("channel<2 && njet==1");
  return string("ERROR! Unclassified channel!");
}

void fitSignalShapeMR(int massBin, int channel, 
                      double rangeLow, double rangeHigh,
                      double bwSigma,
                      double fitValues[5], double fitErrors[5]);
                        
void all(int channel=0) {
  double bwSigma[30];
  int mass[30]; double xLow[30]; double xHigh[30];  
  int maxMassBin;

  mass[0] = 115; xLow[0] = 50; xHigh[0] = 200; bwSigma[0] = 3.1/1000.;
  mass[1] = 120; xLow[1] = 50; xHigh[1] = 200; bwSigma[1] = 3.5/1000.;
  mass[2] = 125; xLow[1] = 50; xHigh[1] = 200; bwSigma[1] = 3.5/1000.;
  mass[3] = 130; xLow[2] = 50; xHigh[2] = 200; bwSigma[2] = 4.9/1000.;
  mass[4] = 135; xLow[2] = 50; xHigh[2] = 200; bwSigma[2] = 4.9/1000.;
  mass[5] = 140; xLow[3] = 50; xHigh[3] = 200; bwSigma[3] = 8.1/1000.;
  mass[6] = 150; xLow[4] = 50; xHigh[4] = 250; bwSigma[4] = 1.7/100.;
  mass[7] = 160; xLow[5] = 50; xHigh[5] = 250; bwSigma[5] = 8.3/100.;
  mass[8] = 170; xLow[6] = 50; xHigh[6] = 250; bwSigma[6] = 3.8/10.;
  mass[9] = 180; xLow[7] = 60; xHigh[7] = 300; bwSigma[7] = 6.3/10.;
  mass[10] = 190; xLow[8] = 70; xHigh[8] = 300; bwSigma[8] = 1.04;
  mass[11] = 200; xLow[9] = 80; xHigh[9] = 300; bwSigma[9] = 1.43;
  mass[12] = 210; xLow[10] = 80; xHigh[10] = 350; bwSigma[10] = 1.85;
  mass[13] = 220; xLow[11] = 80; xHigh[11] = 350; bwSigma[11] = 2.31;
  maxMassBin = 14;

  double massV[30],massE[30];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double aVal[30],aErr[30];
  double nVal[30],nErr[30];
  double meanCBVal[30],meanCBErr[30];
  double sigmaCBVal[30],sigmaCBErr[30];
  double meanBWVal[30],meanBWErr[30];

  double fitValues[5];
  double fitErrors[5];

//   double extendL(1),extendH(1);

//   if(channels==0) {extendL=1.0;extendH=1.0;}
//   if(channels==1) {extendL=0.90;extendH=1.05;}
//   if(channels==2) {extendL=0.95;extendH=1.04;}

//  for(int i=0; i<maxMassBin;++i){
  
  for(int i=2; i<3;++i) {
    fitSignalShapeMR(mass[i],channel,xLow[i],xHigh[i],bwSigma[i],
		    fitValues,fitErrors);  
  
    cout << "a value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    aVal[i]=fitValues[0]; aErr[i]=fitErrors[0];

    cout << "n value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    nVal[i]=fitValues[3]; nErr[i]=fitErrors[3];

    cout << "meanCB value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
    meanCBVal[i]=fitValues[1]; meanCBErr[i]=fitErrors[1];
    
    cout << "sigmaCB value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    sigmaCBVal[i]=fitValues[4]; sigmaCBErr[i]=fitErrors[4];

    cout << "meanBW value,error: " << fitValues[2] << " , " << fitErrors[2] << endl; 
    meanBWVal[i]=fitValues[2]; meanBWErr[i]=fitErrors[2];

    //Wait();
  }

}

void fitSignalShapeMR(int massBin, int channel, 
                     double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[5], double fitErrors[5]){
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);

  stringstream hFileName;
  hFileName << "results/datasets_trees/H" << massBin << "_ll.root";

  cout << "Opening ROOT file: " << hFileName.str() << endl;

  TFile* hFile = TFile::Open(hFileName.str().c_str());
  TTree* hTree = (TTree*) hFile->Get("latinoFitSkim");
  
  float mr, dphillr;
  float puW,effW,baseW;

  hTree->SetBranchAddress("mr",&mr);
  hTree->SetBranchAddress("dphillr",&dphillr);
  hTree->SetBranchAddress("puW",&puW);
  hTree->SetBranchAddress("baseW",&baseW);
  hTree->SetBranchAddress("effW",&effW);

  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;

  TCut cut = getStringChannel(channel).c_str();
  RooRealVar x("mr","M_{R}",xInit,xMin,xMax,"GeV");
  RooRealVar w("baseW","baseW",1.0,0.,1000.);
  RooRealVar cha("channel","channel",0,-0.5,3.5);
  RooRealVar njet("njet","njet",0,-0.5,10.);

  RooArgSet varset(x,w,cha,njet);
  RooDataSet dataset("mass","mass",varset,Import(*hTree),WeightVar("baseW"),Cut(cut));

  //--- simple CrystalBall
  RooRealVar mean("mean","mean of gaussian",0,-5.0,5.0) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,30.); 
  RooRealVar a("a","a",1.46,0.,4.);
  RooRealVar n("n","n",1.92,0.,25.);   
  RooCBShape CBall("CBall","Crystal ball",x, mean,sigma, a,n);
  

  //--- simple Gaussian
  RooRealVar mean2("mean2","mean2 of gaussian",xInit,xInit*0.80,xInit*1.2) ;
  RooRealVar sigma2("sigma2","width2 of gaussian",10.,5.,500.); 
  RooGaussian tailCatcher("tailCatcher","tailCatcher",x,mean2,sigma2);
  RooRealVar fsig("fsig","signal fraction",0.95,0.7,1.);


  //--- Breit-Wigner
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooBreitWigner bw("bw","bw",x,mean3,sigma3);

  //RooAddPdf model("model","model",RooArgList(CBall,tailCatcher),fsig);
  //RooCBShape model("model","model",x, mean,sigma, a,n);
  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,CBall);

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channel==of0j){frameTitle << "e#mu,0-j, m_{H} = ";}
  if(channel==of1j){frameTitle << "e#mu,1-j, m_{H} = ";}
  if(channel==sf0j){frameTitle << "#ell#ell,0-j, m_{H} = ";}
  if(channel==sf1j){frameTitle << "#ell#ell,1-j, m_{H} = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  model.plotOn(xframe);
  model.paramOn(xframe);

  stringstream nameFile;
  nameFile << "fitM" << massBin << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a.getVal();
    fitValues[1] = mean.getVal();
    fitValues[2] = mean3.getVal();
    fitValues[3] = n.getVal();
    fitValues[4] = sigma.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = a.getError();
    fitErrors[1] = mean.getError();
    fitErrors[2] = mean3.getError();
    fitErrors[3] = n.getError();
    fitErrors[4] = sigma.getError();
  }

  return;

}

