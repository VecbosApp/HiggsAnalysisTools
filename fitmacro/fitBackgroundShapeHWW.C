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
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"

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

void fitLandauShapeMR(int channel, string sample, bool looseloose,
		      double rangeLow, double rangeHigh,
		      double fitValues[2], double fitErrors[2]);

void fitGaussianShapeMR(int channel, string sample,
                        double rangeLow, double rangeHigh,
                        double fitValues[2], double fitErrors[2]);

void fitOthersSFShapeMR(int channel, string sample, bool looseloose,
			double rangeLow, double rangeHigh,
			double fitValues[14], double fitErrors[14]);

void allWW(int channel=0);
void allTop(int channel=0);
void allDY(int channel=0);
void allOthers(int channel=0);
void allWjets(int channel=0);
                        
void doAllChannels() {
  cout << "==> Fitting WW sample..." << endl;
  for(int i=0; i<4; ++i) allWW(i);
  cout << "### Done WW sample ###" << endl;

  cout << "==> Fitting top sample..." << endl;
  for(int i=0; i<4; ++i) allTop(i);
  cout << "### Done top sample ###" << endl;

  cout << "==> Fitting DY sample..." << endl;
  for(int i=0; i<4; ++i) allDY(i);
  cout << "### Done DY sample ###" << endl;

  cout << "==> Fitting others sample..." << endl;
  for(int i=0; i<4; ++i) allOthers(i);
  cout << "### Done Others sample ###" << endl;

  cout << "==> Fitting wjets sample..." << endl;
  for(int i=0; i<4; ++i) allWjets(i);
  cout << "### Done Wjets sample ###" << endl;
}

void allWW(int channel) {
  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[2];
  double fitErrors[2];

  fitLandauShapeMR(channel,"WW",false,xLow,xHigh,fitValues,fitErrors);
  cout << "mean value,error = " << fitValues[0] << " , " << fitErrors[0] << endl;
  cout << "sigma value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
}

void allTop(int channel) {
  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[2];
  double fitErrors[2];

  fitLandauShapeMR(channel,"top",false,xLow,xHigh,fitValues,fitErrors);
  cout << "mean value,error = " << fitValues[0] << " , " << fitErrors[0] << endl;
  cout << "sigma value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
}

void allOthers(int channel) {
  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[14];
  double fitErrors[14];

  if(channel<2) {
    fitLandauShapeMR(channel,"others",false,xLow,xHigh,fitValues,fitErrors);
    cout << "mean value,error = " << fitValues[0] << " , " << fitErrors[0] << endl;
    cout << "sigma value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
  } else {
    fitOthersSFShapeMR(channel,"others",false,xLow,xHigh,fitValues,fitErrors);    
    for(int i=0;i<14;++i) cout << "a" << i << " value,error = " << fitValues[i] << " , " << fitErrors[i] << endl;
  }
}

void allDY(int channel) {
  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[14];
  double fitErrors[14];

  if(channel<2) {
    fitGaussianShapeMR(channel,"dataset_embeddedtt",xLow,xHigh,fitValues,fitErrors);
    cout << "mean1 value,error =  " << fitValues[0] << " , " << fitErrors[0] << endl;
    cout << "sigmq1 value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
    cout << "mean2 value,error =  " << fitValues[2] << " , " << fitErrors[2] << endl;
    cout << "sigma2 value,error = " << fitValues[3] << " , " << fitErrors[3] << endl;
    cout << "fsig value,error   = " << fitValues[4] << " , " << fitErrors[4] << endl;
  } else {
    fitOthersSFShapeMR(channel,"Zjets",false,xLow,xHigh,fitValues,fitErrors);    
    for(int i=0;i<14;++i) cout << "a" << i << " value,error = " << fitValues[i] << " , " << fitErrors[i] << endl;
  }
}

void allWjets(int channel) {
  double xLow, xHigh;
  xLow = 50; xHigh = 500;

  double fitValues[14];
  double fitErrors[14];

  if(channel<2) {
    fitLandauShapeMR(channel,"wjets",true,xLow,xHigh,fitValues,fitErrors);
    cout << "mean value,error = " << fitValues[0] << " , " << fitErrors[0] << endl;
    cout << "sigma value,error = " << fitValues[1] << " , " << fitErrors[1] << endl;
  } else {
    fitOthersSFShapeMR(channel,"wjets",true,xLow,xHigh,fitValues,fitErrors);    
    for(int i=0;i<14;++i) cout << "a" << i << " value,error = " << fitValues[i] << " , " << fitErrors[i] << endl;
  }
}

void fitLandauShapeMR(int channel, string sample, bool looseloose,
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
  if(!looseloose) hFileName << "results/datasets_trees/" << sample << "_ll.root";
  else hFileName << "results_data/datasets_trees/dataset_looseloose_wwbits.root";

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

  stringstream weight;
  if(!looseloose) weight << "baseW";
  else weight << "fake2W";
  RooRealVar x("mr","M_{R}",xInit,xMin,xMax,"GeV");
  RooRealVar w(weight.str().c_str(),weight.str().c_str(),1.0,-1000.,1000.);
  RooRealVar cha("channel","channel",0,-0.5,3.5);
  RooRealVar njet("njet","njet",0,-0.5,10.);

  RooArgSet varset(x,w,cha,njet);
  RooDataSet dataset("mass","mass",varset,Import(*hTree),WeightVar(weight.str().c_str()),Cut(cut));


  //--- Landau
  RooRealVar mean("mean","mean",140,70,200) ;
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
  nameFile << "fit" << sample << "_" << getChannelSuffix(channel) << ".pdf";
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

void fitGaussianShapeMR(int channel, string sample,
                        double rangeLow, double rangeHigh,
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
  hFileName << "results/datasets_trees/" << sample << "_ll.root";

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


  //--- simple Gaussian 1
  RooRealVar mean("mean_{1}","mean",140,50,200) ;
  RooRealVar sigma("#sigma_{1}","width",20,10,30); 
  RooGaussian gauss("gauss","gauss",x,mean,sigma);

  RooRealVar mean2("mean_{2}","mean2",140,50,300) ;
  RooRealVar sigma2("#sigma_{2}","width2",35,30,200); 
  RooGaussian gauss2("gauss2","gauss2",x,mean2,sigma2);
  RooRealVar fsig("f_{1}","signal fraction",0.95,0.7,1.);

  RooAddPdf dgauss("dgauss","model",gauss,gauss2,fsig);

  dgauss.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));


  stringstream frameTitle;
  if(channel==of0j){frameTitle << "e#mu,0-j";}
  if(channel==of1j){frameTitle << "e#mu,1-j";}
  if(channel==sf0j){frameTitle << "ee+#mu#mu,0-j";}
  if(channel==sf1j){frameTitle << "ee+#mu#mu,1-j";}

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  dgauss.plotOn(xframe);
  dgauss.paramOn(xframe);

  stringstream nameFile;
  nameFile << "fit" << sample << "_" << getChannelSuffix(channel) << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = mean.getVal();
    fitValues[1] = sigma.getVal();
    fitValues[2] = mean2.getVal();
    fitValues[3] = sigma2.getVal();
    fitValues[4] = fsig.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = mean.getError();
    fitErrors[1] = sigma.getError();
    fitValues[2] = mean2.getError();
    fitValues[3] = sigma2.getError();
    fitValues[4] = fsig.getError();
  }

  return;

}


void fitOthersSFShapeMR(int channel, string sample, bool looseloose,
			double rangeLow, double rangeHigh,
			double fitValues[14], double fitErrors[14]){
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
  if(!looseloose) hFileName << "results/datasets_trees/" << sample << "_ll.root";
  else hFileName << "results_data/datasets_trees/dataset_looseloose_wwbits.root";

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

  stringstream weight;
  if(!looseloose) weight << "baseW";
  else weight << "fake2W";

  RooRealVar x("mr","M_{R}",xInit,xMin,xMax,"GeV");
  RooRealVar w(weight.str().c_str(),weight.str().c_str(),1.0,-1000.,1000.); 
  RooRealVar cha("channel","channel",0,-0.5,3.5);
  RooRealVar njet("njet","njet",0,-0.5,10.);

  RooArgSet varset(x,w,cha,njet);
  RooDataSet dataset("mass","mass",varset,Import(*hTree),WeightVar(weight.str().c_str()),Cut(cut));


  //--- RooqqZZPdf
  RooRealVar a0("a0","a0",100.,0.,200.);
  RooRealVar a1("a1","a1",15.,0.,50.);
  RooRealVar a2("a2","a2",120.,20.,200.);
  RooRealVar a3("a3","a3",0.04,0.,1.);
  RooRealVar a4("a4","a4",185.,100.,400.);
  RooRealVar a5("a5","a5",10.,0.,150.);
  RooRealVar a6("a6","a6",36.,0.,100.);
  RooRealVar a7("a7","a7",0.11,0.,1.);
  RooRealVar a8("a8","a8",60.,0.,150.);
  RooRealVar a9("a9","a9",0.06,0.,1.);
  RooRealVar a10("a10","a10",95.,20.,200.);
  RooRealVar a11("a11","a11",-6.,-20.,20.);
  RooRealVar a12("a12","a12",1000.,0.,10000.);
  RooRealVar a13("a13","a13",0.1,0.,1.);
  RooqqZZPdf_v2 othersSFPdf("othersSFPdf","othersSFPdf",x,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13);

  x.setBins(10000,"fft");
  othersSFPdf.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channel==of0j){frameTitle << "e#mu,0-j";}
  if(channel==of1j){frameTitle << "e#mu,1-j";}
  if(channel==sf0j){frameTitle << "ee+#mu#mu,0-j";}
  if(channel==sf1j){frameTitle << "ee+#mu#mu,1-j";}

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  othersSFPdf.plotOn(xframe);
  //  othersSFPdf.paramOn(xframe);

  stringstream nameFile;
  nameFile << "fit" << sample << "_" << getChannelSuffix(channel) << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a0.getVal();
    fitValues[1] = a1.getVal();
    fitValues[2] = a2.getVal();
    fitValues[3] = a3.getVal();
    fitValues[4] = a4.getVal();
    fitValues[5] = a5.getVal();
    fitValues[6] = a6.getVal();
    fitValues[7] = a7.getVal();
    fitValues[8] = a8.getVal();
    fitValues[9] = a9.getVal();
    fitValues[10] = a10.getVal();
    fitValues[11] = a11.getVal();
    fitValues[12] = a12.getVal();
    fitValues[13] = a13.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = a0.getError();
    fitErrors[1] = a1.getError();
    fitErrors[2] = a2.getError();
    fitErrors[3] = a3.getError();
    fitErrors[4] = a4.getError();
    fitErrors[5] = a5.getError();
    fitErrors[6] = a6.getError();
    fitErrors[7] = a7.getError();
    fitErrors[8] = a8.getError();
    fitErrors[9] = a9.getError();
    fitErrors[10] = a10.getError();
    fitErrors[11] = a11.getError();
    fitErrors[12] = a12.getError();
    fitErrors[13] = a13.getError();
  }

  return;

}

