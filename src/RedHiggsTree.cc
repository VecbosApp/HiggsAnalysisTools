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
  myTree->Branch("finalSelection",      &myFinalSelection,      "finalSelection/B");

}

RedHiggsTree::~RedHiggsTree() 
{
  delete myFile;
}

void RedHiggsTree::addCSA07Infos() {

  myTree->Branch("CSA07weight",              &myWeight,              "CSA07weight/D");
  myTree->Branch("CSA07processId",           &myProcesId,            "CSA07processId/D");
  myTree->Branch("CSA07lumi",                &myLumi,                "CSA07lumi/F");

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
			   bool finalLeptons, bool jetVeto, bool finalSelection)
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
  myFinalSelection = finalSelection;

}

void RedHiggsTree::fillCSA07(double weight, double processId, float lumi) 
{

  myWeight = weight;
  myProcesId = processId;
  myLumi = lumi;

}
