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
  createSingleDataset("results/datasets_trees/H120_ee.root","results/datasets/H120_ee.root");
  createSingleDataset("results/datasets_trees/H130_ee.root","results/datasets/H130_ee.root");
  createSingleDataset("results/datasets_trees/H140_ee.root","results/datasets/H140_ee.root");
  createSingleDataset("results/datasets_trees/H150_ee.root","results/datasets/H150_ee.root");
  createSingleDataset("results/datasets_trees/H155_ee.root","results/datasets/H155_ee.root");
  createSingleDataset("results/datasets_trees/H160_ee.root","results/datasets/H160_ee.root");
  createSingleDataset("results/datasets_trees/H165_ee.root","results/datasets/H165_ee.root");
  createSingleDataset("results/datasets_trees/H170_ee.root","results/datasets/H170_ee.root");
  createSingleDataset("results/datasets_trees/H175_ee.root","results/datasets/H175_ee.root");
  createSingleDataset("results/datasets_trees/H180_ee.root","results/datasets/H180_ee.root");
  createSingleDataset("results/datasets_trees/H190_ee.root","results/datasets/H190_ee.root");
  createSingleDataset("results/datasets_trees/H200_ee.root","results/datasets/H200_ee.root");
  createSingleDataset("results/datasets_trees/H210_ee.root","results/datasets/H210_ee.root");
  createSingleDataset("results/datasets_trees/H220_ee.root","results/datasets/H220_ee.root");
  createSingleDataset("results/datasets_trees/H230_ee.root","results/datasets/H230_ee.root");
  createSingleDataset("results/datasets_trees/H240_ee.root","results/datasets/H240_ee.root");
  createSingleDataset("results/datasets_trees/H250_ee.root","results/datasets/H250_ee.root");
  createSingleDataset("results/datasets_trees/H275_ee.root","results/datasets/H275_ee.root");
  createSingleDataset("results/datasets_trees/H300_ee.root","results/datasets/H300_ee.root");
  createSingleDataset("results/datasets_trees/H350_ee.root","results/datasets/H350_ee.root");
  createSingleDataset("results/datasets_trees/H400_ee.root","results/datasets/H400_ee.root");
  createSingleDataset("results/datasets_trees/H500_ee.root","results/datasets/H500_ee.root");
  createSingleDataset("results/datasets_trees/H550_ee.root","results/datasets/H550_ee.root");
  createSingleDataset("results/datasets_trees/H600_ee.root","results/datasets/H600_ee.root");
  createSingleDataset("results/datasets_trees/others_ee.root","results/datasets/others_ee.root");
  createSingleDataset("results/datasets_trees/TTbar_ee.root","results/datasets/TTbar_ee.root");
  createSingleDataset("results/datasets_trees/TTbarJetsMadgraph_ee.root","results/datasets/TTbarJetsMadgraph_ee.root");
  createSingleDataset("results/datasets_trees/WW_ee.root","results/datasets/WW_ee.root");
  createSingleDataset("results/datasets_trees/Zspecies_ee.root","results/datasets/Zspecies_ee.root");

  cout << "CREATING ROODATSETS FOR MM CHANNEL..." << endl;
  createSingleDataset("results/datasets_trees/H120_mm.root","results/datasets/H120_mm.root");
  createSingleDataset("results/datasets_trees/H130_mm.root","results/datasets/H130_mm.root");
  createSingleDataset("results/datasets_trees/H140_mm.root","results/datasets/H140_mm.root");
  createSingleDataset("results/datasets_trees/H150_mm.root","results/datasets/H150_mm.root");
  createSingleDataset("results/datasets_trees/H155_mm.root","results/datasets/H155_mm.root");
  createSingleDataset("results/datasets_trees/H160_mm.root","results/datasets/H160_mm.root");
  createSingleDataset("results/datasets_trees/H165_mm.root","results/datasets/H165_mm.root");
  createSingleDataset("results/datasets_trees/H170_mm.root","results/datasets/H170_mm.root");
  createSingleDataset("results/datasets_trees/H175_mm.root","results/datasets/H175_mm.root");
  createSingleDataset("results/datasets_trees/H180_mm.root","results/datasets/H180_mm.root");
  createSingleDataset("results/datasets_trees/H190_mm.root","results/datasets/H190_mm.root");
  createSingleDataset("results/datasets_trees/H200_mm.root","results/datasets/H200_mm.root");
  createSingleDataset("results/datasets_trees/H210_mm.root","results/datasets/H210_mm.root");
  createSingleDataset("results/datasets_trees/H220_mm.root","results/datasets/H220_mm.root");
  createSingleDataset("results/datasets_trees/H230_mm.root","results/datasets/H230_mm.root");
  createSingleDataset("results/datasets_trees/H240_mm.root","results/datasets/H240_mm.root");
  createSingleDataset("results/datasets_trees/H250_mm.root","results/datasets/H250_mm.root");
  createSingleDataset("results/datasets_trees/H275_mm.root","results/datasets/H275_mm.root");
  createSingleDataset("results/datasets_trees/H300_mm.root","results/datasets/H300_mm.root");
  createSingleDataset("results/datasets_trees/H350_mm.root","results/datasets/H350_mm.root");
  createSingleDataset("results/datasets_trees/H400_mm.root","results/datasets/H400_mm.root");
  createSingleDataset("results/datasets_trees/H500_mm.root","results/datasets/H500_mm.root");
  createSingleDataset("results/datasets_trees/H550_mm.root","results/datasets/H550_mm.root");
  createSingleDataset("results/datasets_trees/H600_mm.root","results/datasets/H600_mm.root");
  createSingleDataset("results/datasets_trees/others_mm.root","results/datasets/others_mm.root");
  createSingleDataset("results/datasets_trees/TTbar_mm.root","results/datasets/TTbar_mm.root");
  createSingleDataset("results/datasets_trees/TTbarJetsMadgraph_mm.root","results/datasets/TTbarJetsMadgraph_mm.root");
  createSingleDataset("results/datasets_trees/WW_mm.root","results/datasets/WW_mm.root");
  createSingleDataset("results/datasets_trees/Zspecies_mm.root","results/datasets/Zspecies_mm.root");

  cout << "CREATING ROODATSETS FOR EM CHANNEL..." << endl;
  createSingleDataset("results/datasets_trees/H120_em.root","results/datasets/H120_em.root");
  createSingleDataset("results/datasets_trees/H130_em.root","results/datasets/H130_em.root");
  createSingleDataset("results/datasets_trees/H140_em.root","results/datasets/H140_em.root");
  createSingleDataset("results/datasets_trees/H150_em.root","results/datasets/H150_em.root");
  createSingleDataset("results/datasets_trees/H155_em.root","results/datasets/H155_em.root");
  createSingleDataset("results/datasets_trees/H160_em.root","results/datasets/H160_em.root");
  createSingleDataset("results/datasets_trees/H165_em.root","results/datasets/H165_em.root");
  createSingleDataset("results/datasets_trees/H170_em.root","results/datasets/H170_em.root");
  createSingleDataset("results/datasets_trees/H175_em.root","results/datasets/H175_em.root");
  createSingleDataset("results/datasets_trees/H180_em.root","results/datasets/H180_em.root");
  createSingleDataset("results/datasets_trees/H190_em.root","results/datasets/H190_em.root");
  createSingleDataset("results/datasets_trees/H200_em.root","results/datasets/H200_em.root");
  createSingleDataset("results/datasets_trees/H210_em.root","results/datasets/H210_em.root");
  createSingleDataset("results/datasets_trees/H220_em.root","results/datasets/H220_em.root");
  createSingleDataset("results/datasets_trees/H230_em.root","results/datasets/H230_em.root");
  createSingleDataset("results/datasets_trees/H240_em.root","results/datasets/H240_em.root");
  createSingleDataset("results/datasets_trees/H250_em.root","results/datasets/H250_em.root");
  createSingleDataset("results/datasets_trees/H275_em.root","results/datasets/H275_em.root");
  createSingleDataset("results/datasets_trees/H300_em.root","results/datasets/H300_em.root");
  createSingleDataset("results/datasets_trees/H350_em.root","results/datasets/H350_em.root");
  createSingleDataset("results/datasets_trees/H400_em.root","results/datasets/H400_em.root");
  createSingleDataset("results/datasets_trees/H500_em.root","results/datasets/H500_em.root");
  createSingleDataset("results/datasets_trees/H550_em.root","results/datasets/H550_em.root");
  createSingleDataset("results/datasets_trees/H600_em.root","results/datasets/H600_em.root");
  createSingleDataset("results/datasets_trees/others_em.root","results/datasets/others_em.root");
  createSingleDataset("results/datasets_trees/TTbar_em.root","results/datasets/TTbar_em.root");
  createSingleDataset("results/datasets_trees/TTbarJetsMadgraph_em.root","results/datasets/TTbarJetsMadgraph_em.root");
  createSingleDataset("results/datasets_trees/WW_em.root","results/datasets/WW_em.root");

}

void createSingleDataset(const char *treefile, const char *roofitfile) {

  cout << "roofitting file " << treefile << " in " << roofitfile << endl;

  RooRealVar *jetcat = new RooRealVar("jetcat","jetcat",-1,1); // cut the -2 ( >1 jet )
  RooRealVar *met    = new RooRealVar("met","E_{T}^{miss}",0,200,"GeV");
  RooRealVar *deltaPhi = new RooRealVar("deltaPhi","#Delta#phi",0,180,"#deg");
  RooRealVar *minPtEle = new RooRealVar("minPtEle","p_{T}^{min}",0,200,"GeV");
  RooRealVar *eleInvMass = new RooRealVar("eleInvMass","m(l^{+}l^{-})",0,150,"GeV");
  RooRealVar *bTagImpPar = new RooRealVar("bTagImpPar","b-tag",-1001.,2.0);
  RooRealVar *weight = new RooRealVar("weight","weight",0,10000);

  RooArgSet setHiggs(*jetcat,*met,*deltaPhi,*minPtEle,*eleInvMass,*bTagImpPar,*weight);

  TFile *file = TFile::Open(treefile);
  TTree *tree = (TTree*)file->Get("T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setHiggs);

  TFile *roofitFile = TFile::Open(roofitfile,"recreate");
  data->Write();
  roofitFile->Close();

}
