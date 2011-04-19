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
  myTree->Branch("run",                 &myRun,                 "run/I");
  myTree->Branch("ls",                  &myLS,                  "ls/I");
  myTree->Branch("event",               &myEvent,               "event/I");
  myTree->Branch("met",                 &myMet,                 "met/F");  
  myTree->Branch("pfMet",               &myPFMet,               "pfMet/F");  
  myTree->Branch("caloMet",             &myCaloMet,             "caloMet/F");  
  myTree->Branch("projMet",             &myProjectedMet,        "projMet/F");
  myTree->Branch("deltaPhi",            &myDeltaPhi,            "deltaPhi/F");  
  myTree->Branch("deltaR",              &myDeltaR,              "deltaR/F");  
  myTree->Branch("transvMass",          &myTransvMass,          "transvMass/F");  
  myTree->Branch("eleInvMass",          &myEleInvMass,          "eleInvMass/F");  
  myTree->Branch("maxPtEle",            &maxPtEle,              "maxPtEle/F");  
  myTree->Branch("minPtEle",            &minPtEle,              "minPtEle/F");  
  myTree->Branch("detaLeptons",         &myDetaLeptons,         "detaLeptons/F");  
  myTree->Branch("finalLeptons",        &myFinalLeptons,        "finalLeptons/B");
  myTree->Branch("jetVeto",             &myJetVeto,             "jetVeto/B");
  myTree->Branch("uncorrJetVeto",       &myUncorrJetVeto,       "uncorrJetVeto/B");
  myTree->Branch("preDeltaPhi",         &myPreDeltaPhi,         "preDeltaPhi/B");
  myTree->Branch("finalSelection",      &myFinalSelection,      "finalSelection/B");

}

RedHiggsTree::~RedHiggsTree() 
{
  delete myFile;
}

void RedHiggsTree::addMLVars() {
  myTree->Branch("njets",               &myNjets,               "njets/I");
  myTree->Branch("nuncorrjets",         &myNuncorrjets,         "nuncorrjets/I");
  myTree->Branch("dxyEVT",              &myDxyEVT,              "dxyEVT/F");
  myTree->Branch("dszEVT",              &myDszEVT,              "dszEVT/F");
  myTree->Branch("bTagTrackCount",      &myBTagTrackCount,      "bTagTrackCount/F");
  myTree->Branch("bTagImpPar",          &myBTagImpPar,          "bTagImpPar/F");
  myTree->Branch("bTagSecVertex",       &myBTagSecVertex,       "bTagSecVertex/F");
}

void RedHiggsTree::addElectronInfos() {
  
  myTree->Branch("recoflag", myRecoflag, "recoflag[2]/I");
  myTree->Branch("pt", myPt, "pt[2]/F");
  myTree->Branch("eta", myEta, "eta[2]/F");
  myTree->Branch("phi", myPhi, "phi[2]/F");
  myTree->Branch("classification", myClassification, "classification[2]/I");
  myTree->Branch("nbrems", myNBremClusters, "nbrems[2]/I");
  myTree->Branch("deta", myDeta, "deta[2]/F");
  myTree->Branch("dphi", myDphi, "dphi[2]/F");
  myTree->Branch("hoe", myHoe, "hoe[2]/F");
  myTree->Branch("see", mySee, "see[2]/F");
  myTree->Branch("spp", mySpp, "spp[2]/F");
  myTree->Branch("eop", myEop, "eop[2]/F");
  myTree->Branch("fbrem", myFbrem, "fbrem[2]/F");
  myTree->Branch("trackerIso", myTrackerIso, "trackerIso[2]/F");
  myTree->Branch("hcalIso", myHcalIso, "hcalIso[2]/F");
  myTree->Branch("ecalJIso", myEcalJIso, "ecalJIso[2]/F");
  myTree->Branch("ecalGTIso", myEcalGTIso, "ecalGTIso[2]/F");
  myTree->Branch("combinedIso", myCombinedIso, "combinedIso[2]/F");
  myTree->Branch("charge", myCharge, "charge[2]/I");
  myTree->Branch("missHits", myMissHits, "missHits[2]/I");
  myTree->Branch("dist", myDist, "dist[2]/F");
  myTree->Branch("dcot", myDcot, "dcot[2]/F");
  myTree->Branch("lh", myLh, "lh[2]/F");
  myTree->Branch("matched", myMatched, "matched[2]/I");
}

void RedHiggsTree::addCSA07Infos() {

  myTree->Branch("CSA07weight",              &myWeight,              "CSA07weight/D");
  myTree->Branch("CSA07processId",           &myProcesId,            "CSA07processId/D");
  myTree->Branch("CSA07lumi",                &myLumi,                "CSA07lumi/F");

}

void RedHiggsTree::addLatinos() {

  myTree->Branch("step1",              &myStep1,              "step1/B");
  myTree->Branch("step2",              &myStep2,              "step2/B");
  myTree->Branch("step3",              &myStep3,              "step3/B");
  myTree->Branch("step4",              &myStep4,              "step4/B");
  myTree->Branch("step5",              &myStep5,              "step5/B");
  myTree->Branch("step6",              &myStep6,              "step6/B");
  myTree->Branch("step7",              &myStep7,              "step7/B");
  myTree->Branch("step8",              &myStep8,              "step8/B");
  myTree->Branch("step9",              &myStep9,              "step9/B");
  myTree->Branch("step10",             &myStep10,             "step10/B");
  myTree->Branch("step11",             &myStep11,             "step11/B");
  myTree->Branch("step12",             &myStep12,             "step12/B");
  myTree->Branch("step13",             &myStep13,             "step13/B");
  myTree->Branch("step14",             &myStep14,             "step14/B");
  myTree->Branch("step15",             &myStep15,             "step15/B");
  myTree->Branch("step16",             &myStep16,             "step16/B");
  myTree->Branch("step17",             &myStep17,             "step17/B");

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

void RedHiggsTree::addRunInfos() {

  myTree->Branch("run", &myRun,     "run/I");
  myTree->Branch("lumi", &myLS,     "lumi/I");
  myTree->Branch("event", &myEvent, "event/I");

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


void RedHiggsTree::fillAll(float met, float pfmet, float cmet, float projmet, 
			   float dphi, float derre, float tmass, float mee, float max, float min, float deta,
			   bool finalLeptons, bool jetVeto, bool uncorrjetVeto, bool preDeltaPhi, bool finalSelection)
{

  myMet         = met;
  myPFMet       = pfmet;
  myCaloMet     = cmet;
  myProjectedMet = projmet;
  myDeltaPhi    = dphi;
  myDeltaR      = derre;
  myTransvMass  = tmass;
  myEleInvMass  = mee;
  maxPtEle      = max;
  minPtEle      = min;
  myDetaLeptons = deta;
  myFinalLeptons = finalLeptons;
  myJetVeto       = jetVeto;
  myUncorrJetVeto = uncorrjetVeto;
  myPreDeltaPhi = preDeltaPhi;
  myFinalSelection = finalSelection;

}

void RedHiggsTree::fillMLVars(int njets, int nuncorrjets, float dxyEVT, float dszEVT, 
                              float bTagTrackCount, float bTagImpPar, float bTagSecVertex) {

  myNjets   = njets;
  myNuncorrjets = nuncorrjets;
  myDxyEVT = dxyEVT;
  myDszEVT = dszEVT;
  myBTagTrackCount = bTagTrackCount;
  myBTagImpPar = bTagImpPar;
  myBTagSecVertex = bTagSecVertex;

}

void RedHiggsTree::fillLatinos(bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16, bool s17) {

  myStep1  = s1;
  myStep2  = s2;
  myStep3  = s3;
  myStep4  = s4;
  myStep5  = s5;
  myStep6  = s6;
  myStep7  = s7;
  myStep8  = s8;
  myStep9  = s9;
  myStep10 = s10;
  myStep11 = s11;
  myStep12 = s12;
  myStep13 = s13;
  myStep14 = s14;
  myStep15 = s15;
  myStep16 = s16;
  myStep17 = s17;
}

void RedHiggsTree::fillElectrons(int recoflag[2], float pt[2], float eta[2], float phi[2],
                                 int classification[2], int nbrems[2], float deta[2], float dphi[2], float hoe[2], float see[2], float spp[2], float eop[2], float fbrem[2],
                                 float trackerIso[2], float hcalIso[2], float ecalJIso[2], float ecalGTIso[2], float combinedIso[2], int charge[2],
                                 int missHits[2], float dist[2], float dcot[2], float lh[2], int matched[2]) {

  for(int i=0; i<2; i++) {
    myRecoflag[i] = recoflag[i];
    myPt[i] = pt[i];
    myEta[i] = eta[i];
    myPhi[i] = phi[i];
    myClassification[i] = classification[i];
    myNBremClusters[i] = nbrems[i];
    myDeta[i] = deta[i];
    myDphi[i] = dphi[i];
    myHoe[i] = hoe[i];
    mySee[i] = see[i];
    mySpp[i] = spp[i];
    myEop[i] = eop[i];
    myFbrem[i] = fbrem[i];
    myTrackerIso[i] = trackerIso[i];
    myHcalIso[i] = hcalIso[i];
    myEcalJIso[i] = ecalJIso[i];
    myEcalGTIso[i] = ecalGTIso[i];
    myCombinedIso[i] = combinedIso[i];
    myCharge[i] = charge[i];
    myMissHits[i] = missHits[i];
    myDist[i] = dist[i];
    myDcot[i] = dcot[i];
    myLh[i] = lh[i];
    myMatched[i] = matched[i];
  }
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

void RedHiggsTree::fillRunInfos(int run, int lumi, int event) {

  myRun = run;
  myLS = lumi;
  myEvent = event;

}


