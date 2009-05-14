#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>

using namespace std;

void createDatasetsEE() {

  RooRealVar *nJetsVar = new RooRealVar("jetCat","jetCat",-2,2);
  RooRealVar *METVar = new RooRealVar("MET","MET",0,200,"GeV");
  RooRealVar *deltaPhiVar = new RooRealVar("deltaPhi","deltaPhi",0,180,"#deg");
  RooRealVar *maxPtVar = new RooRealVar("MaxPt","maxPt",0,200,"GeV");
  RooRealVar *minPtVar = new RooRealVar("minPt","minPt",0,200,"GeV");
  RooRealVar *invMassVar = new RooRealVar("invMass","invMass",0,200,"GeV");
  RooRealVar *dxyEVTVar = new RooRealVar("dxyEVT","dxyEVT",0,5000,"#mum");
  RooRealVar *dszEVTVar = new RooRealVar("dszEVT","dszEVT",0,5000,"#mum");
  RooRealVar *weightVar = new RooRealVar("weight","weight",0,100);
   
  RooArgSet setHiggs(*nJetsVar,*METVar,*deltaPhiVar,*maxPtVar,*minPtVar,*invMassVar);
  setHiggs.add(*dxyEVTVar);
  setHiggs.add(*dszEVTVar);
  setHiggs.add(*weightVar);

  // signal
  TFile *file = TFile::Open("hww_2e_Tree_Higgs_21X.root");
  TTree *tree = (TTree*)file->Get("data");
  RooDataSet *dataset_Higgs = new RooDataSet("Higgs","Higgs",tree,setHiggs);

  char outname[200];
  sprintf(outname,"hww_2e_Datasets_21X.root");

  TFile *fileout = TFile::Open(outname,"recreate");
  dataset_Higgs->Write();
  fileout->Close();

  delete dataset_Higgs;

  // WW bkg
  file = TFile::Open("hww_2e_Tree_WW_2l_21X.root");
  tree = (TTree*)file->Get("data");
  RooDataSet *dataset_WW = new RooDataSet("WW","WW",tree,setHiggs);

  fileout = TFile::Open(outname,"update");
  dataset_WW->Write();
  fileout->Close();

  delete dataset_WW;


  // ttbar bkg
  file = TFile::Open("hww_2e_Tree_ttjetsMADGRAPH_21X.root");
  tree = (TTree*)file->Get("data");
  RooDataSet *dataset_ttbar = new RooDataSet("ttbar","ttbar",tree,setHiggs);

  fileout = TFile::Open(outname,"update");
  dataset_ttbar->Write();
  fileout->Close();

  delete dataset_ttbar;
  

  // other bkg
  file = TFile::Open("hww_2e_Tree_other_21X.root");
  tree = (TTree*)file->Get("data");
  RooDataSet *dataset_other = new RooDataSet("other","other",tree,setHiggs);

  fileout = TFile::Open(outname,"update");
  dataset_other->Write();
  fileout->Close();

  delete dataset_other;

}
