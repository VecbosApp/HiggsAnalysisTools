#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TPad.h>

#include "YieldMaker.h"
#include "FitSelection.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

enum prodmode {gg, vbf};

int Wait() {
     cout << " Continue [<RET>|q]?  ";
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

string getChannelSuffix(int channel) {
  if(channel==of0j) return string("of_0j");
  if(channel==of1j) return string("of_1j");
  if(channel==sf0j) return string("sf_0j");
  if(channel==sf1j) return string("sf_1j");
  return string("ERROR! Unclassified channel!");
}

float getYield(int mH, int cha, int prod, float lumi, bool barecount);
float getYieldStatError(float yield, float n) { return yield/sqrt(n); }
void allmasses(int prod, int cha, float lumi);

void all(int lumi) {
  for(int i=0;i<4;++i) {
    for(int j=0;j<2;++j) {
      allmasses(j,i,lumi);
    }
  }
}

void allmasses(int prod, int cha, float lumi) {
  int mass[30];
  float yieldV[30], yieldE[30];
  float massV[30], massE[30];
  int maxMassBin;

  mass[0] = 115;
  mass[1] = 120;
  mass[2] = 125;
  mass[3] = 130;
  mass[4] = 135;
  mass[5] = 140;
  mass[6] = 150;
  mass[7] = 160;
  mass[8] = 170;
  mass[9] = 180;
  maxMassBin = 10;

  for(int i=0;i<maxMassBin;++i) {
    yieldV[i] = getYield(mass[i],cha,prod,lumi,false);
    // stat error is negligible
    // float staterr = getYieldStatError(yieldV[i],getYield(mass[i],cha,prod,lumi,true));
    float systerr = 0.15 *  yieldV[i];
    yieldE[i] = systerr;
    massV[i] = mass[i];
    massE[i]=0;
  }

  TGraphErrors* gY = new TGraphErrors(maxMassBin,massV,yieldV,massE,yieldE);
  gY->SetMarkerStyle(20);   gY->SetMarkerSize(1);
  gY->SetTitle("");
  gY->GetXaxis()->SetTitle("mass (GeV)");
  gY->GetYaxis()->SetTitle("signal yield");

  gStyle->SetOptFit(111111);
  stringstream nameFile;
  nameFile << "higgsYield_" <<  getChannelSuffix(cha) << (prod==gg ? "_ggH" : "_vbfH") << "_lumi" << lumi << "invfb";
  gY->Fit("pol3"); gY->Draw("Ap"); gPad->Update(); gPad->Print((nameFile.str()+string(".pdf")).c_str()); Wait();

}

float getYield(int mH, int cha, int prod, float lumi, bool barecount) {
 stringstream hFileName;
 if(prod==gg) hFileName << "latinos_tree_skim_of/nominals/latino_1" << mH << "_ggToH" << mH << "toWWTo2LAndTau2Nu.root";
 else if(prod==vbf) hFileName << "latinos_tree_skim_of/nominals/latino_2" << mH << "_vbfToH" << mH << "toWWTo2LAndTau2Nu.root";
 cout << "Opening ROOT file: " << hFileName.str() << endl;

 FitSelection sel;

 YieldMaker ymaker_hi;
 ymaker_hi.fill(hFileName.str().c_str());

 float yield = 0.0;
 if(barecount) yield = ymaker_hi.getCount(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax); 
 else yield = ymaker_hi.getYield(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax) * lumi;

 return yield;
}


