#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TSystem.h>

#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>

#include <iostream>

using namespace std;

void createSingleDataset(const char *treefile, const char *roofitfile);

void createAll() {

  cout << "CREATING ROODATSETS FOR EE CHANNEL..." << endl;
  createSingleDataset("results/datasets_trees/H130_ee.root","results/datasets/H130_ee.root");
  createSingleDataset("results/datasets_trees/H160_ee.root","results/datasets/H160_ee.root");
  createSingleDataset("results/datasets_trees/H190_ee.root","results/datasets/H190_ee.root");
  createSingleDataset("results/datasets_trees/others_ee.root","results/datasets/others_ee.root");
  createSingleDataset("results/datasets_trees/top_ee.root","results/datasets/top_ee.root");
  createSingleDataset("results/datasets_trees/WW_ee.root","results/datasets/WW_ee.root");

  cout << "CREATING ROODATSETS FOR MM CHANNEL..." << endl;
  createSingleDataset("results/datasets_trees/H130_mm.root","results/datasets/H130_mm.root");
  createSingleDataset("results/datasets_trees/H160_mm.root","results/datasets/H160_mm.root");
  createSingleDataset("results/datasets_trees/H190_mm.root","results/datasets/H190_mm.root");
  createSingleDataset("results/datasets_trees/others_mm.root","results/datasets/others_mm.root");
  createSingleDataset("results/datasets_trees/top_mm.root","results/datasets/top_mm.root");
  createSingleDataset("results/datasets_trees/WW_mm.root","results/datasets/WW_mm.root");

  cout << "CREATING ROODATSETS FOR EM CHANNEL..." << endl;
  createSingleDataset("results/datasets_trees/H130_em.root","results/datasets/H130_em.root");
  createSingleDataset("results/datasets_trees/H160_em.root","results/datasets/H160_em.root");
  createSingleDataset("results/datasets_trees/H190_em.root","results/datasets/H190_em.root");
  createSingleDataset("results/datasets_trees/others_em.root","results/datasets/others_em.root");
  createSingleDataset("results/datasets_trees/top_em.root","results/datasets/top_em.root");
  createSingleDataset("results/datasets_trees/WW_em.root","results/datasets/WW_em.root");

}

void createSingleDataset(const char *treefile, const char *roofitfile) {

  cout << "roofitting file " << treefile << " in " << roofitfile << endl;

  RooRealVar *jetcat = new RooRealVar("jetcat","jetcat",-1,1); // cut the -2 ( >1 jet )
  RooRealVar *met    = new RooRealVar("met","E_{T}^{miss}",0,200,"GeV");
  RooRealVar *deltaPhi = new RooRealVar("deltaPhi","#Delta#phi",0,180,"#deg");
  RooRealVar *maxPtEle = new RooRealVar("maxPtEle","p_{T}^{max}",20,200,"GeV");
  RooRealVar *eleInvMass = new RooRealVar("eleInvMass","m(l^{+}l^{-})",12,150,"GeV");
  RooRealVar *bTagImpPar = new RooRealVar("bTagImpPar","b-tag",-1001.,2.0);
  RooRealVar *weight = new RooRealVar("weight","weight",0,10000);
  RooRealVar *event = new RooRealVar("event","progressive event number",0,1e+07);

  RooArgSet setHiggs(*jetcat,*met,*deltaPhi,*maxPtEle,*eleInvMass,*bTagImpPar,*weight,*event);

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setHiggs);

  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();

}
