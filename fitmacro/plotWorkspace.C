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
#include <sstream>

using namespace RooFit;

void plotWsp2D(const char *inputfile, const char *figfile) {

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
  RooAbsPdf *wgstar = (RooAbsPdf*)w->pdf("bkg_wgstar");
  RooAbsPdf *others = (RooAbsPdf*)w->pdf("bkg_others");
  RooAbsPdf *ggH = (RooAbsPdf*)w->pdf("ggH");
  RooAbsPdf *qqH = (RooAbsPdf*)w->pdf("qqH");

  qqww->plotOn(mrplot,LineColor(kAzure+1));
  ggww->plotOn(mrplot,LineColor(kAzure+2));
  top->plotOn(mrplot,LineColor(kYellow+2));
  dy->plotOn(mrplot,LineColor(kGreen+2));
  wj->plotOn(mrplot,LineColor(kGray+2));
  wgstar->plotOn(mrplot,LineColor(kOrange+7));
  others->plotOn(mrplot,LineColor(kOrange+2));
  ggH->plotOn(mrplot,LineColor(kRed));
  qqH->plotOn(mrplot,LineColor(kRed));
  
  // for the legend
  TH1F *hqqww = new TH1F("hqqww","",1,0,1);
  TH1F *hggww = new TH1F("hggww","",1,0,1);
  TH1F *htop = new TH1F("htop","",1,0,1);
  TH1F *hdy = new TH1F("hdy","",1,0,1);
  TH1F *hwj = new TH1F("hwj","",1,0,1);
  TH1F *hwgstar = new TH1F("hwgstar","",1,0,1);
  TH1F *hothers = new TH1F("hothers","",1,0,1);
  TH1F *hggH = new TH1F("hggH","",1,0,1);

  hqqww->SetLineColor(kAzure+1);
  hggww->SetLineColor(kAzure+2);
  htop->SetLineColor(kYellow+2);
  hdy->SetLineColor(kGreen+2);
  hwj->SetLineColor(kGray+2);
  hwgstar->SetLineColor(kOrange+7);
  hothers->SetLineColor(kOrange+2);
  hggH->SetLineColor(kRed);

  TLegend* legend = new TLegend(0.64, 0.24, 0.87, 0.90);
    
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  ( 0.05);

  legend->AddEntry(hqqww,"qq #rightarrow WW","l");
  legend->AddEntry(hggww,"gg #rightarrow WW","l");
  legend->AddEntry(htop,"top","l");
  legend->AddEntry(hdy,"DY","l");
  legend->AddEntry(hwj,"W+jets","l");
  legend->AddEntry(hothers,"WW,ZZ,W#gamma","l");
  legend->AddEntry(hwgstar,"W#gamma^{(*)}","l");
  legend->AddEntry(hothers,"WW,ZZ,W#gamma","l");
  legend->AddEntry(hggH,"qq,gg #rightarrow H","l");

  TCanvas *c1 = new TCanvas("c1","c1");
  mrplot->Draw();
  legend->Draw();

  c1->SaveAs(figfile);

}

void plotMRAllSignals() {

  stringstream fss0;
  fss0 << "datacards/card_1D_m114_8TeV_of_1j_workspace.root";
  
  TFile *hww2l2nu0 = TFile::Open(fss0.str().c_str());
  RooWorkspace *w = (RooWorkspace*)hww2l2nu0->Get("w");

  TCanvas *c1 = new TCanvas("c1","c1");
  
  RooRealVar *mr = w->var("CMS_ww2l_mr_1D");
  RooPlot *mrplot = mr->frame();
  mrplot->Draw();

  TLegend* legend = new TLegend(0.64, 0.64, 0.87, 0.90);
    
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  ( 0.05);

  int j=0;
  for (float i = 115.; i <= 180.; i += 10.) {
    j++;

    stringstream fss;
    fss << "datacards/card_1D_m" << i << "_8TeV_of_1j_workspace.root";

    cout << "Opening file " << fss.str() << endl;

    TFile *hww2l2nu = TFile::Open(fss.str().c_str());
    w = (RooWorkspace*)hww2l2nu->Get("w");
    
    // take the pdfs from the workspace
    RooAbsPdf *ggH = (RooAbsPdf*)w->pdf("ggH");
    ggH->plotOn(mrplot,LineColor(kAzure+j));  

    TH1F *h = new TH1F("h","h",0,0,1);
    h->SetLineColor(kAzure+j);

    stringstream lab;
    lab << "m_{H}=" << i << " GeV";
    legend->AddEntry(h,lab.str().c_str(),"l");

  }

  mrplot->Draw();
  legend->Draw();
  c1->SaveAs("severalHiggses.png");
}

void plotAll() {
  
  plotWsp2D("datacards/card_1D_m125_8TeV_of_0j_workspace.root","pdfs_of_0j.pdf");
  plotWsp2D("datacards/card_1D_m125_8TeV_of_1j_workspace.root","pdfs_of_1j.pdf");
  plotMRAllSignals();

}
