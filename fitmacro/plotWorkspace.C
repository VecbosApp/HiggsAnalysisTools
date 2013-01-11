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

#include "YieldMaker.h"
#include "DatacardParser.hh"

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
  fss0 << "datacards/hww-12.1fb.mH114.of_0j_shape_workspace.root";
  
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
    fss << "datacards/hww-12.1fb.mH" << i << ".of_0j_shape_workspace.root";

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

void plotOneShapeSyst(string process, string syst, int ch) {

  std::string chstr;
  if (ch == of0j) chstr = "of_0j";
  if (ch == of1j) chstr = "of_1j";
  if (ch == sf0j) chstr = "sf_0j";
  if (ch == sf1j) chstr = "sf_1j";

  stringstream fss;
  fss << "datacards/hww-12.1fb.mH125." << chstr << "_shape_workspace.root";

  TFile *hww2l2nu = TFile::Open(fss.str().c_str());

  RooWorkspace *w = (RooWorkspace*)hww2l2nu->Get("w");

  RooRealVar *mr = w->var("CMS_ww2l_mr_1D");
  RooPlot *mrplot = mr->frame();

  // take the pdfs from the workspace
  RooAbsPdf *pdfnom = (RooAbsPdf*)w->pdf(process.c_str());
  if(pdfnom==0) {
    cout << "PDF for process: " << process << " not found." << endl;
    return;
  }
  pdfnom->plotOn(mrplot,LineColor(kBlack));

  string tevstr = "_8TeV";

  // only special case where name in DC!=name in the workspace. A bit error prone... to be improved 
  string systvar=syst;
  if(syst.compare("scaleup_qcd")!=string::npos && syst.length()==11) systvar="scaleup";
  if(syst.compare("scaledn_qcd")!=string::npos && syst.length()==11) systvar="scaledn";
  
  RooRealVar *WW_mean_err = (RooRealVar*)w->var((process+"_"+chstr+tevstr+"_mean_err_"+systvar).c_str());
  RooRealVar *WW_sigma_err = (RooRealVar*)w->var((process+"_"+chstr+tevstr+"_sigma_err_"+systvar).c_str());

  if(WW_mean_err==0 || WW_sigma_err==0) {
    cout << WW_mean_err->GetName() << " not found " << endl;
    cout << WW_sigma_err->GetName() << " not found " << endl;
    return;
  }

  stringstream fdc;
  fdc << "datacards/hww-12.1fb.mH125." << chstr << "_shape.txt";
  DatacardParser dcp(fdc.str());

  string dcm(WW_mean_err->GetName());
  string dcs(WW_sigma_err->GetName());
  if(syst.compare("qcd")==0) { dcm+="-qcd"; dcs+="-qcd"; }

  float meanShift = dcp.getRelUncertainty(dcm);
  float sigmaShift = dcp.getRelUncertainty(dcs);
  cout << "Landau error on mean/sigma =  " << meanShift << " / " << sigmaShift << endl;
  if(meanShift<-100 || sigmaShift<-100) {
    cout << dcm << " or " << dcs << "  not found in the datacard txt file" << endl;
    return;
  }

  WW_mean_err->setVal(meanShift);
  WW_sigma_err->setVal(sigmaShift);

  pdfnom->plotOn(mrplot,LineColor(kRed+1));

  TCanvas *c1 = new TCanvas("c1","c1");
  mrplot->Draw();
  c1->SaveAs("syst.png");

}

void plotAll() {
  
  plotWsp2D("datacards/card_1D_m125_8TeV_of_0j_workspace.root","pdfs_of_0j.pdf");
  plotWsp2D("datacards/card_1D_m125_8TeV_of_1j_workspace.root","pdfs_of_1j.pdf");
  plotMRAllSignals();

}
