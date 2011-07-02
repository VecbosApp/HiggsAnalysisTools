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
  myTree->Branch("puweight",            &myPUWeight,            "puweight/F");
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
  myTree->Branch("maxEtaEle",           &maxEtaEle,             "maxEtaEle/F");  
  myTree->Branch("minEtaEle",           &minEtaEle,             "minEtaEle/F");  
  myTree->Branch("detaLeptons",         &myDetaLeptons,         "detaLeptons/F");  
  myTree->Branch("nVtx",                &myNVtx,                "nVtx/I");
  myTree->Branch("finalLeptons",        &myFinalLeptons,        "finalLeptons/O");
  myTree->Branch("jetVeto",             &myJetVeto,             "jetVeto/O");
  myTree->Branch("uncorrJetVeto",       &myUncorrJetVeto,       "uncorrJetVeto/O");
  myTree->Branch("preDeltaPhi",         &myPreDeltaPhi,         "preDeltaPhi/O");
  myTree->Branch("finalSelection",      &myFinalSelection,      "finalSelection/O");

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
  myTree->Branch("nSoftMu",             &myNSoftMu,             "nSoftMu/I");
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

  myTree->Branch("step",              mySteps,              "step[25]/O"); 

}

void RedHiggsTree::addKinematics() {
  myTree->Branch("pxTkMet", &myPxTkMet, "pxTkMet/F");
  myTree->Branch("pyTkMet", &myPyTkMet, "pyTkMet/F");
  myTree->Branch("pzTkMet", &myPzTkMet, "pzTkMet/F");
  myTree->Branch("pxLeadJet", &myPxLeadJet, "pxLeadJet/F");
  myTree->Branch("pyLeadJet", &myPyLeadJet, "pyLeadJet/F");
  myTree->Branch("pzLeadJet", &myPzLeadJet, "pzLeadJet/F");
  myTree->Branch("pxL1", &myPxL1, "pxL1/F");
  myTree->Branch("pyL1", &myPyL1, "pyL1/F");
  myTree->Branch("pzL1", &myPzL1, "pzL1/F");
  myTree->Branch("pxL2", &myPxL2, "pxL2/F");
  myTree->Branch("pyL2", &myPyL2, "pyL2/F");
  myTree->Branch("pzL2", &myPzL2, "pzL2/F");
}

void RedHiggsTree::addKFactor() {
  
  myTree->Branch("KFactor",       &myKFactor,      "KFactor/F");
  myTree->Branch("GenHPt",        &myGenHPt,       "GenHPt/F");
  myTree->Branch("leadingJetPt",  &myLeadingJetPt, "leadingJetPt/F");
}

void RedHiggsTree::addMcTruthInfos() {

  myTree->Branch("promptDecay",         &myPromptDecay,         "promptDecay/O");

}

void RedHiggsTree::addHLTElectronsInfos() {

  myTree->Branch("HLTSingleElectron",        &myHLTSingleElectron,        "HLTSingleElectron/O");
  myTree->Branch("HLTSingleElectronRelaxed", &myHLTSingleElectronRelaxed, "HLTSingleElectronRelaxed/O");
  myTree->Branch("HLTSingleElectronOR",      &myHLTSingleElectronOR,      "HLTSingleElectronOR/O");

}

void RedHiggsTree::addHLTMuonsInfos() {

  myTree->Branch("HLTSingleMuon",        &myHLTSingleMuon,        "HLTSingleMuon/O");
  myTree->Branch("HLTSingleMuonRelaxed", &myHLTSingleMuonRelaxed, "HLTSingleMuonRelaxed/O");
  myTree->Branch("HLTSingleMuonOR",      &myHLTSingleMuonOR,      "HLTSingleMuonOR/O");

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
			   float dphi, float derre, float tmass, float mee, float max, float min, float deta, int nvtx,
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
  myNVtx        = nvtx;
  myFinalLeptons = finalLeptons;
  myJetVeto       = jetVeto;
  myUncorrJetVeto = uncorrjetVeto;
  myPreDeltaPhi = preDeltaPhi;
  myFinalSelection = finalSelection;

}

void RedHiggsTree::fillAll(float met, float pfmet, float cmet, float projmet, 
			   float dphi, float derre, float tmass, float mee, 
			   float max, float min, float maxEta, float minEta, float deta, int nvtx,
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
  maxEtaEle     = maxEta;
  minEtaEle     = minEta;
  myDetaLeptons = deta;
  myNVtx        = nvtx;
  myFinalLeptons = finalLeptons;
  myJetVeto       = jetVeto;
  myUncorrJetVeto = uncorrjetVeto;
  myPreDeltaPhi = preDeltaPhi;
  myFinalSelection = finalSelection;

}

void RedHiggsTree::fillMLVars(int njets, int nuncorrjets, float dxyEVT, float dszEVT, 
                              float bTagTrackCount, float bTagImpPar, float bTagSecVertex, int nsoftmu) {

  myNjets   = njets;
  myNuncorrjets = nuncorrjets;
  myDxyEVT = dxyEVT;
  myDszEVT = dszEVT;
  myBTagTrackCount = bTagTrackCount;
  myBTagImpPar = bTagImpPar;
  myBTagSecVertex = bTagSecVertex;
  myNSoftMu = nsoftmu;

}

void RedHiggsTree::fillLatinos(bool s0, bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16, bool s17,
                               bool s18, bool s19, bool s20, bool s21, bool s22, bool s23, bool s24) {
  mySteps[0]  = s0;
  mySteps[1]  = s1;
  mySteps[2]  = s2;
  mySteps[3]  = s3;
  mySteps[4]  = s4;
  mySteps[5]  = s5;
  mySteps[6]  = s6;
  mySteps[7]  = s7;
  mySteps[8]  = s8;
  mySteps[9]  = s9;
  mySteps[10] = s10;
  mySteps[11] = s11;
  mySteps[12] = s12;
  mySteps[13] = s13;
  mySteps[14] = s14;
  mySteps[15] = s15;
  mySteps[16] = s16;
  mySteps[17] = s17;
  mySteps[18] = s18;
  mySteps[19] = s19;
  mySteps[20] = s20;
  mySteps[21] = s21;
  mySteps[22] = s22;
  mySteps[23] = s23;
  mySteps[24] = s24;

}

void RedHiggsTree::fillKinematics(float pxTkMet, float pyTkMet, float pzTkMet,
                                  float pxLeadJet, float pyLeadJet, float pzLeadJet,
                                  float pxL1, float pyL1, float pzL1,
                                  float pxL2, float pyL2, float pzL2) {
  myPxTkMet = pxTkMet;
  myPyTkMet = pyTkMet;
  myPzTkMet = pzTkMet;
  myPxLeadJet = pxLeadJet;
  myPyLeadJet = pyLeadJet;
  myPzLeadJet = pzLeadJet;
  myPxL1 = pxL1;
  myPyL1 = pyL1;
  myPzL1 = pzL1;
  myPxL2 = pxL2;
  myPyL2 = pyL2;
  myPzL2 = pzL2;
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

void RedHiggsTree::fillKFactor(float kfactor, float genh, float ptlj ) {

  myKFactor      = kfactor;
  myGenHPt       = genh;
  myLeadingJetPt = ptlj;
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

void RedHiggsTree::fillRunInfos(int run, int lumi, int event, float puweight) {

  myRun = run;
  myLS = lumi;
  myEvent = event;
  myPUWeight = puweight;

}


