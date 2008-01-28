#include "./RedEleIDTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedEleIDTree::RedEleIDTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  // GENERAL block
  myTree->Branch("sampleOk",                 &mySampleOk,                 "sampleOk/I");  
  myTree->Branch("singleElePassedTrg",       &mySingleElePassedTrg,       "singleElePassedTrg/I");  
  myTree->Branch("singleEleRelaxPassedTrg",  &mySingleEleRelaxPassedTrg,  "singleEleRelaxPassedTrg/I");  
  myTree->Branch("doubleElePassedTrg",       &myDoubleElePassedTrg,       "doubleElePassedTrg/I");  
  myTree->Branch("doubleEleRelaxPassedTrg",  &myDoubleEleRelaxPassedTrg,  "doubleEleRelaxPassedTrg/I");  
  myTree->Branch("chargeEle",                &myChargeEle,                "chargeEle/I");  
  myTree->Branch("energyEle",                &myEnergyEle,                "energyEle/F");  
  myTree->Branch("etEle",                    &myEtEle,                    "etEle/F");  
  myTree->Branch("momentumEle",              &myMomentumEle,              "momentumEle/F");
  myTree->Branch("thetaEle",                 &myThetaEle,                 "thetaEle/F");  
  myTree->Branch("etaEle",                   &myEtaEle,                   "etaEle/F");  
  myTree->Branch("phiEle",                   &myPhiEle,                   "phiEle/F");  
  myTree->Branch("latEle",                   &myLatEle,                   "latEle/F");  
  myTree->Branch("a20Ele",                   &myA20Ele,                   "a20Ele/F");  
  myTree->Branch("s9s25Ele",                 &myS9s25Ele,                 "s9s25Ele/F");
  myTree->Branch("covEtaEtaEle",             &myCovEtaEtaEle,             "covEtaEtaEle/F");  
  myTree->Branch("eleClassEle",              &myEleClassEle,              "eleClassEle/I");  
  myTree->Branch("eleHoEEle",                &myEleHoEEle,                "eleHoEEle/F");  
  myTree->Branch("eleCorrEoPEle",            &myEleCorrEoPEle,            "eleCorrEoPEle/F");  
  myTree->Branch("eleCorrEoPoutEle",         &myEleCorrEoPoutEle,         "eleCorrEoPEle/F");  
  myTree->Branch("eleDeltaEtaAtVtxEle",      &myEleDeltaEtaAtVtxEle,      "eleDeltaEtaAtVtxEle/F");  
  myTree->Branch("eleDeltaPhiAtVtxEle",      &myEleDeltaPhiAtVtxEle,      "eleDeltaPhiAtVtxEle/F");  
  myTree->Branch("eleTrackerIso_sumPtEle",   &myEleTrackerIso_SumPtEle,   "eleTrackerIso_sumPtEle/F");  
  myTree->Branch("eleLikelihoodEle",         &myEleLikelihoodEle,         "eleLikelihoodEle/F");  
  myTree->Branch("eleFisherEle",             &myEleFisherEle,             "eleFisherEle/F");  
}


RedEleIDTree::~RedEleIDTree() 
{
  delete myFile;
}

void RedEleIDTree::store()
{
  myTree->Fill();
}


void RedEleIDTree::save() 
{
  myFile->cd();
  myTree->Write();
  myFile->Close();
}


void RedEleIDTree::fillAll(int isok, int seHLT, int serHLT, int deHLT, int derHLT, int charge, float ene, float et, float mom, float theta, float eta, float phi, float lat, float a20, float s9s25, float covEE, int theclass, float hoe, float eop, float eopout, float deta, float dphi, float iso, float like, float fis)
{
  mySampleOk                = isok;
  mySingleElePassedTrg      = seHLT;
  mySingleEleRelaxPassedTrg = serHLT; 
  myDoubleElePassedTrg      = deHLT;       
  myDoubleEleRelaxPassedTrg = derHLT;   
  myChargeEle               = charge;
  myEnergyEle               = ene; 
  myEtEle                   = et;
  myMomentumEle             = mom;
  myThetaEle                = theta; 
  myEtaEle                  = eta; 
  myPhiEle                  = phi; 
  myLatEle                  = lat; 
  myA20Ele                  = a20;
  myS9s25Ele                = s9s25; 
  myCovEtaEtaEle            = covEE; 
  myEleClassEle             = theclass;
  myEleHoEEle               = hoe; 
  myEleCorrEoPEle           = eop; 
  myEleCorrEoPoutEle        = eopout;
  myEleDeltaEtaAtVtxEle     = deta; 
  myEleDeltaPhiAtVtxEle     = dphi;
  myEleTrackerIso_SumPtEle  = iso;
  myEleLikelihoodEle        = like; 
  myEleFisherEle            = fis; 
}
