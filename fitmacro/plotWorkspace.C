#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>

using namespace RooFit;

void plotWsp2D(char *inputfile, char *figfile) {

  TFile *hww2l2nu = TFile::Open(inputfile);
  RooWorkspace *w = (RooWorkspace*)hww2l2nu->Get("w");

  RooRealVar *mr = w->var("CMS_ww2l_mr_1D");
  RooPlot *mrplot = mr->frame();

  // take the pdfs from the workspace
  RooAbsPdf *qqww = (RooAbsPdf*)w->pdf("bkg_qqww");
  RooAbsPdf *ggww = (RooAbsPdf*)w->pdf("bkg_ggww");
  RooAbsPdf *top = (RooAbsPdf*)w->pdf("bkg_top");
  RooAbsPdf *dy = (RooAbsPdf*)w->pdf("bkg_dy");
  RooAbsPdf *wj = (RooAbsPdf*)w->pdf("bkg_wj");
  RooAbsPdf *others = (RooAbsPdf*)w->pdf("bkg_others");
  RooAbsPdf *ggH = (RooAbsPdf*)w->pdf("ggH");
  RooAbsPdf *qqH = (RooAbsPdf*)w->pdf("qqH");

  qqww->plotOn(mrplot,LineColor(kAzure+1));
  ggww->plotOn(mrplot,LineColor(kAzure+2));
  top->plotOn(mrplot,LineColor(kYellow+2));
  dy->plotOn(mrplot,LineColor(kGreen+2));
  others->plotOn(mrplot,LineColor(kOrange+2));
  wj->plotOn(mrplot,LineColor(kGray+2));
  ggH->plotOn(mrplot,LineColor(kRed));
  qqH->plotOn(mrplot,LineColor(kRed));
  
  // for the legend
  TH1F *hqqww = new TH1F("hqqww","",1,0,1);
  TH1F *hggww = new TH1F("hggww","",1,0,1);
  TH1F *htop = new TH1F("htop","",1,0,1);
  TH1F *hdy = new TH1F("hdy","",1,0,1);
  TH1F *hwj = new TH1F("hwj","",1,0,1);
  TH1F *hothers = new TH1F("hothers","",1,0,1);
  TH1F *hggH = new TH1F("hggH","",1,0,1);

  hqqww->SetLineColor(kAzure+1);
  hggww->SetLineColor(kAzure+2);
  htop->SetLineColor(kYellow+2);
  hdy->SetLineColor(kGreen+2);
  hothers->SetLineColor(kOrange+2);
  hwj->SetLineColor(kGray+2);
  hggH->SetLineColor(kRed);

  TLegend* legend = new TLegend(0.64, 0.64, 0.87, 0.90);
    
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  ( 0.05);

  legend->AddEntry(hqqww,"qq #rightarrow WW","l");
  legend->AddEntry(hggww,"gg #rightarrow WW","l");
  legend->AddEntry(htop,"top","l");
  legend->AddEntry(hdy,"W+jets","l");
  legend->AddEntry(hothers,"WW,ZZ,W#gamma","l");
  legend->AddEntry(hggH,"qq,gg #rightarrow H","l");

  TCanvas *c1 = new TCanvas("c1","c1");
  mrplot->Draw();
  legend->Draw();

  c1->SaveAs(figfile);

}
