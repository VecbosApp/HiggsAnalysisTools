#include "HiggsAnalysisTools/include/RedHiggsTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedHiggsTree::RedHiggsTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  // GENERAL block
  myTree->Branch("met",                 &myMet,                 "met/F");  
  myTree->Branch("deltaPhi",            &myDeltaPhi,            "deltaPhi/F");  
  myTree->Branch("transvMass",          &myTransvMass,          "transvMass/F");  
  myTree->Branch("eleInvMass",          &myEleInvMass,          "eleInvMass/F");  
  myTree->Branch("maxPtEle",            &maxPtEle,              "maxPtEle/F");  
  myTree->Branch("minPtEle",            &minPtEle,              "minPtEle/F");  
  myTree->Branch("detaLeptons",         &myDetaLeptons,         "detaLeptons/F");  
  myTree->Branch("finalLeptons",        &myFinalLeptons,        "finalLeptons/B");
  myTree->Branch("jetVeto",             &myJetVeto,             "jetVeto/B");
  myTree->Branch("preDeltaPhi",         &myPreDeltaPhi,         "preDeltaPhi/B");
  myTree->Branch("finalSelection",      &myFinalSelection,      "finalSelection/B");

}

RedHiggsTree::~RedHiggsTree() 
{
  delete myFile;
}

void RedHiggsTree::addMLVars() {
  myTree->Branch("maxPtLh",             &myMaxPtLh,             "maxPtLh/F");
  myTree->Branch("minPtLh",             &myMinPtLh,             "minPtLh/F");
  myTree->Branch("njets",               &myNjets,               "njets/I");
  myTree->Branch("dxyEVT",              &myDxyEVT,              "dxyEVT/F");
  myTree->Branch("dszEVT",              &myDszEVT,              "dszEVT/F");
}

void RedHiggsTree::addCSA07Infos() {

  myTree->Branch("CSA07weight",              &myWeight,              "CSA07weight/D");
  myTree->Branch("CSA07processId",           &myProcesId,            "CSA07processId/D");
  myTree->Branch("CSA07lumi",                &myLumi,                "CSA07lumi/F");

}

void RedHiggsTree::addKFactor() {
  
  myTree->Branch("KFactor",  &myKFactor,   "KFactor/F");

}

void RedHiggsTree::addMcTruthInfos() {

  myTree->Branch("promptDecay",         &myPromptDecay,         "promptDecay/B");

}

void RedHiggsTree::addHLTElectronsInfos() {

  myTree->Branch("HLTSingleElectron",        &myHLTSingleElectron,        "HLTSingleElectron/B");
  myTree->Branch("HLTSingleElectronRelaxed", &myHLTSingleElectronRelaxed, "HLTSingleElectronRelaxed/B");
  myTree->Branch("HLTSingleElectronOR",      &myHLTSingleElectronOR,      "HLTSingleElectronOR/B");

}

void RedHiggsTree::addHLTMuonsInfos() {

  myTree->Branch("HLTSingleMuon",        &myHLTSingleMuon,        "HLTSingleMuon/B");
  myTree->Branch("HLTSingleMuonRelaxed", &myHLTSingleMuonRelaxed, "HLTSingleMuonRelaxed/B");
  myTree->Branch("HLTSingleMuonOR",      &myHLTSingleMuonOR,      "HLTSingleMuonOR/B");

}

void RedHiggsTree::store()
{
  myTree->Fill();
}


void RedHiggsTree::save() 
{
  myFile->cd();
  myTree->Write();
  myFile->Close();
}


void RedHiggsTree::fillAll(float mt, float dphi, float tmass, float mee, float max, float min, float deta,
			   bool finalLeptons, bool jetVeto, bool preDeltaPhi, bool finalSelection)
{

  myMet         = mt;
  myDeltaPhi    = dphi;
  myTransvMass  = tmass;
  myEleInvMass  = mee;
  maxPtEle      = max;
  minPtEle      = min;
  myDetaLeptons = deta;
  myFinalLeptons = finalLeptons;
  myJetVeto = jetVeto;
  myPreDeltaPhi = preDeltaPhi;
  myFinalSelection = finalSelection;

}

void RedHiggsTree::fillMLVars(float maxlh, float minlh, int njets, float dxyEVT, float dszEVT) {

  myMaxPtLh = maxlh;
  myMinPtLh = minlh;
  myNjets = njets;
  myDxyEVT = dxyEVT;
  myDszEVT = dszEVT;

}

void RedHiggsTree::fillCSA07(double weight, double processId, float lumi) 
{

  myWeight = weight;
  myProcesId = processId;
  myLumi = lumi;

}

void RedHiggsTree::fillKFactor(double kfactor) {

  myKFactor = kfactor;

}

void RedHiggsTree::fillMcTruth(bool prompt) {

  myPromptDecay = prompt;

}

void RedHiggsTree::fillHLTElectrons(bool singleEle, bool singleEleRelaxed, bool singleEleOR) {

  myHLTSingleElectron = singleEle;
  myHLTSingleElectronRelaxed = singleEleRelaxed;
  myHLTSingleElectronOR = singleEleOR;

}

void RedHiggsTree::fillHLTMuons(bool singleMuon, bool singleMuonRelaxed, bool singleMuonOR) {

  myHLTSingleMuon = singleMuon;
  myHLTSingleMuonRelaxed = singleMuonRelaxed;
  myHLTSingleMuonOR = singleMuonOR;

}
