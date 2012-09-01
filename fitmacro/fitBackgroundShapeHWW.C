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
#include "RooLandau.h"
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

int Wait() {
     cout << " Continue [<RET>|q]?  ";
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

enum channels { of0j, of1j, sf0j, sf1j };

  
string getChannelSuffix(int channel) {
  if(channel==of0j) return string("of_0j");
  if(channel==of1j) return string("of_1j");
  if(channel==sf0j) return string("sf_0j");
  if(channel==sf1j) return string("sf_1j");
  return string("ERROR! Unclassified channel!");
}

string getStringChannel(int channel) {
  if(channel==of0j) return string("channel>=2 && njet==0");
  if(channel==of1j) return string("channel>=2 && njet==1");
  if(channel==sf0j) return string("channel<2 && njet==0");
  if(channel==sf1j) return string("channel<2 && njet==1");
  return string("ERROR! Unclassified channel!");
}

void fitWWShapeMR(int channel, 
		  double rangeLow, double rangeHigh,
		  double fitValues[2], double fitErrors[2]);
void all(int channel=0);
                        
void doAllChannels() {
  for(int i=0; i<4; ++i) all(i);
}

void all(int channel) {

  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[2];
  double fitErrors[2];

  fitWWShapeMR(channel,xLow,xHigh,fitValues,fitErrors);
  cout << "mean value,error = " << fitValues[0] << " , " << fitErrors[0] << endl;
  cout << "sigma value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
  

}

void fitWWShapeMR(int channel, 
                  double rangeLow, double rangeHigh,
                  double fitValues[2], double fitErrors[2]){
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
  hFileName << "results/datasets_trees/WW_ll.root";

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
  xInit = 140;
  xMin = rangeLow;
  xMax = rangeHigh ;
  

  TCut cut1 = getStringChannel(channel).c_str();
  stringstream fitrangecut;
  fitrangecut << "mr > " << xMin << " && mr < " << xMax;
  TCut cut2 = fitrangecut.str().c_str();
  TCut cut = cut1 && cut2;

  RooRealVar x("mr","M_{R}",xInit,xMin,xMax,"GeV");
  RooRealVar w("baseW","baseW",1.0,0.,1000.);
  RooRealVar cha("channel","channel",0,-0.5,3.5);
  RooRealVar njet("njet","njet",0,-0.5,10.);

  RooArgSet varset(x,w,cha,njet);
  RooDataSet dataset("mass","mass",varset,Import(*hTree),WeightVar("baseW"),Cut(cut));


  //--- Landau
  RooRealVar mean("mean","mean",140,100,200) ;
  RooRealVar sigma("#sigma","width",50,10,100); 
  RooLandau landau("landau","landau",x,mean,sigma);

  x.setBins(10000,"fft");
  landau.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channel==of0j){frameTitle << "e#mu,0-j";}
  if(channel==of1j){frameTitle << "e#mu,1-j";}
  if(channel==sf0j){frameTitle << "ee+#mu#mu,0-j";}
  if(channel==sf1j){frameTitle << "ee+#mu#mu,1-j";}

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  landau.plotOn(xframe);
  landau.paramOn(xframe);

  stringstream nameFile;
  nameFile << "fitWW_" << getChannelSuffix(channel) << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = mean.getVal();
    fitValues[1] = sigma.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = mean.getError();
    fitErrors[1] = sigma.getError();
  }

  return;

}

