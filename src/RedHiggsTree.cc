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
  myTree->Branch("hlt",                 &myHLT,                 "hlt/O");
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
  myTree->Branch("nSoftMuNoJets",       &myNSoftMuNoJets,       "nSoftMuNoJets/I");
  myTree->Branch("leadingJetBTagTrackCount", &myLeadingJetBTagTrackCount,    "leadingJetBTagTrackCount/F");
  myTree->Branch("subleadingJetBTagTrackCount", &mySubleadingJetBTagTrackCount,    "subleadingJetBTagTrackCount/F");    
  myTree->Branch("subleadingJetsMaxBTagTrackCount", &mySubleadingJetsMaxBTagTrackCount,    "subleadingJetsMaxBTagTrackCount/F");    
  myTree->Branch("numExtraLep", &myNumExtraLep, "numExtraLep/I");   
}

void RedHiggsTree::addSystematics() {
  
  myTree->Branch("scEnergy", myScEnergy, "scEnergy[2]/F");
  myTree->Branch("R9", myR9, "R9[2]/F");
  myTree->Branch("eneL1",  &myEneL1,  "eneL1/F");
  myTree->Branch("eneL2",  &myEneL2, " eneL2/F");
  myTree->Branch("typeL1", &myTypeL1, "typeL1/I");
  myTree->Branch("typeL2", &myTypeL2, "typeL2/I");
  myMetFromJets = 0;
  myPfMetUp = 0;
  myPfMetDown = 0;
  myTree->Branch("metFromJets", "TVector3", &myMetFromJets);
  myTree->Branch("pfMetUp", "TVector3", &myPfMetUp);
  myTree->Branch("pfMetDown", "TVector3", &myPfMetDown);
  myTree->Branch("transvMassUp", &myMtUp, "transvMassUp/F");
  myTree->Branch("transvMassDown", &myMtDown, "transvMassDown/F");

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

  myTree->Branch("step",              mySteps,              "step[29]/O"); 
}

void RedHiggsTree::addRazor() {

  myTree->Branch("mtr",  &myMTR,  "mtr/F");
  myTree->Branch("mr",  &myMR,  "mr/F");
  myTree->Branch("gammamr", &myGammaMR, "gammamr/F");
}

void RedHiggsTree::addFake() {
  
  myTree->Branch("tight",        &myTight,        "tight/I");
  myTree->Branch("weightFP",     &myWeightFP,     "weightFP/F");
  myTree->Branch("weightStatFP", &myWeightStatFP, "weightStatFP/F");

  myTree->Branch("weightFP15",     &myWeightFP15,     "weightFP15/F");
  myTree->Branch("weightStatFP15", &myWeightStatFP15, "weightStatFP15/F");
  myTree->Branch("weightFF15",     &myWeightFF15,     "weightFF15/F");
  myTree->Branch("weightStatFF15", &myWeightStatFF15, "weightStatFF15/F");
  myTree->Branch("weightPP15",     &myWeightPP15,     "weightPP15/F");
  myTree->Branch("weightStatPP15", &myWeightStatPP15, "weightStatPP15/F");

  myTree->Branch("weightFP30",     &myWeightFP30,     "weightFP30/F");
  myTree->Branch("weightStatFP30", &myWeightStatFP30, "weightStatFP30/F");
  myTree->Branch("weightFF30",     &myWeightFF30,     "weightFF30/F");
  myTree->Branch("weightStatFF30", &myWeightStatFF30, "weightStatFF30/F");
  myTree->Branch("weightPP30",     &myWeightPP30,     "weightPP30/F");
  myTree->Branch("weightStatPP30", &myWeightStatPP30, "weightStatPP30/F");

  myTree->Branch("weightFP35",     &myWeightFP35,     "weightFP35/F");
  myTree->Branch("weightStatFP35", &myWeightStatFP35, "weightStatFP35/F");
  myTree->Branch("weightFF35",     &myWeightFF35,     "weightFF35/F");
  myTree->Branch("weightStatFF35", &myWeightStatFF35, "weightStatFF35/F");
  myTree->Branch("weightPP35",     &myWeightPP35,     "weightPP35/F");
  myTree->Branch("weightStatPP35", &myWeightStatPP35, "weightStatPP35/F");

  myTree->Branch("weightFP50",     &myWeightFP50,     "weightFP50/F");
  myTree->Branch("weightStatFP50", &myWeightStatFP50, "weightStatFP50/F");
  myTree->Branch("weightFF50",     &myWeightFF50,     "weightFF50/F");
  myTree->Branch("weightStatFF50", &myWeightStatFF50, "weightStatFF50/F");
  myTree->Branch("weightPP50",     &myWeightPP50,     "weightPP50/F");
  myTree->Branch("weightStatPP50", &myWeightStatPP50, "weightStatPP50/F");

  myTree->Branch("weightFPQCD",     &myWeightFPQCD,     "weightFPQCD/F");
  myTree->Branch("weightStatFPQCD", &myWeightStatFPQCD, "weightStatFPQCD/F");
  myTree->Branch("weightFFQCD",     &myWeightFFQCD,     "weightFFQCD/F");
  myTree->Branch("weightStatFFQCD", &myWeightStatFFQCD, "weightStatFFQCD/F");
  myTree->Branch("weightPPQCD",     &myWeightPPQCD,     "weightPPQCD/F");
  myTree->Branch("weightStatPPQCD", &myWeightStatPPQCD, "weightStatPPQCD/F");
}


void RedHiggsTree::addKinematics() {

  myTree->Branch("pxTkMet", &myPxTkMet, "pxTkMet/F");
  myTree->Branch("pyTkMet", &myPyTkMet, "pyTkMet/F");
  myTree->Branch("pzTkMet", &myPzTkMet, "pzTkMet/F");
  myTree->Branch("pxLeadJet", myPxLeadJet, "pxLeadJet[3]/F"); // 0th is nominal JEC, 1st=+1sigma JES, 2nd=-1sigma JES   
  myTree->Branch("pyLeadJet", myPyLeadJet, "pyLeadJet[3]/F");
  myTree->Branch("pzLeadJet", myPzLeadJet, "pzLeadJet[3]/F");
  myTree->Branch("pxSecondJet", myPxSecondJet, "pxSecondJet[3]/F");
  myTree->Branch("pySecondJet", myPySecondJet, "pySecondJet[3]/F");
  myTree->Branch("pzSecondJet", myPzSecondJet, "pzSecondJet[3]/F");
  myTree->Branch("pxL1", &myPxL1, "pxL1/F");
  myTree->Branch("pyL1", &myPyL1, "pyL1/F");
  myTree->Branch("pzL1", &myPzL1, "pzL1/F");
  myTree->Branch("pxL2", &myPxL2, "pxL2/F");
  myTree->Branch("pyL2", &myPyL2, "pyL2/F");
  myTree->Branch("pzL2", &myPzL2, "pzL2/F");
  myTree->Branch("ch", myLepCharge, "ch[2]/I");
  myTree->Branch("lh", myEleLh, "lh[2]/F");
  myTree->Branch("bdt", myEleBdt, "bdt[2]/F");
  myTree->Branch("iso", myIso, "iso[2]/F");
  myTree->Branch("chmajority", myMajority, "chmajority[2]/I");

  myJetsSum = 0;
  myUncorrJetsSum = 0;
  myPfMet = 0;
  myTree->Branch("sumJetsV4", "TLorentzVector", &myJetsSum);
  myTree->Branch("uncorrSumJetsV4", "TLorentzVector", &myUncorrJetsSum);
  myTree->Branch("pfmetV", "TVector3", &myPfMet);
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

void RedHiggsTree::addMetStudies() {

  myTree->Branch("projPFMet",        &myProjPFMet,        "projPFMet/F");
  myTree->Branch("projPFChargedMet", &myProjPFChargedMet, "projPFChargedMet/F");
  myTree->Branch("signPFMet",        &mySignPFMet,        "signPFMet/F");
  myTree->Branch("signPFChargedMet", &mySignPFChargedMet, "signPFChargedMet/F");
  myTree->Branch("mtrchargedMet",    &myMTRchargedMet,    "mtrchargedMet/F");
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
                              float bTagTrackCount, float bTagImpPar, float bTagSecVertex, int nsoftmu, 
                              float leadJetBTagTrackCount, float subleadJetBTagTrackCount, float subleadJetsMaxBTagTrackCount,
                              int numExtraLep, int nsoftmunojets) {
 
  myNjets   = njets;
  myNuncorrjets = nuncorrjets;
  myDxyEVT = dxyEVT;
  myDszEVT = dszEVT;
  myBTagTrackCount = bTagTrackCount;
  myBTagImpPar = bTagImpPar;
  myBTagSecVertex = bTagSecVertex;
  myNSoftMu = nsoftmu;
  myNSoftMuNoJets = nsoftmunojets;
  myLeadingJetBTagTrackCount = leadJetBTagTrackCount;
  mySubleadingJetBTagTrackCount = subleadJetBTagTrackCount;
  mySubleadingJetsMaxBTagTrackCount = subleadJetsMaxBTagTrackCount;
  myNumExtraLep = numExtraLep;
}

void RedHiggsTree::fillLatinos(bool s0, bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16, bool s17,
                               bool s18, bool s19, bool s20, bool s21, bool s22, bool s23, bool s24, bool s25, bool s26, bool s27, bool s28) {
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
  mySteps[25] = s25;
  mySteps[26] = s26;
  mySteps[27] = s27;
  mySteps[28] = s28;

}

void RedHiggsTree::fillRazor(float MTR, float mR, float gammaMR) {

  myMTR = MTR;
  myMR = mR;
  myGammaMR = gammaMR;
}

void RedHiggsTree::fillFake(int ntigh, float wfp, float wsfp, 
			    float wfp15, float wsfp15, float wff15, float wsff15, float wpp15, float wspp15,
			    float wfp30, float wsfp30, float wff30, float wsff30, float wpp30, float wspp30,
			    float wfp35, float wsfp35, float wff35, float wsff35, float wpp35, float wspp35,
			    float wfp50, float wsfp50, float wff50, float wsff50, float wpp50, float wspp50,
			    float wfpQCD, float wsfpQCD, float wffQCD, float wsffQCD, float wppQCD, float wsppQCD) {
  
  myTight        = ntigh;
  myWeightFP     = wfp;
  myWeightStatFP = wsfp;
  //
  myWeightFP15     = wfp15;
  myWeightStatFP15 = wsfp15;
  myWeightFF15     = wff15;
  myWeightStatFF15 = wsff15;
  myWeightPP15     = wpp15;
  myWeightStatPP15 = wspp15;
  //
  myWeightFP30     = wfp30;
  myWeightStatFP30 = wsfp30;
  myWeightFF30     = wff30;
  myWeightStatFF30 = wsff30;
  myWeightPP30     = wpp30;
  myWeightStatPP30 = wspp30;
  //
  myWeightFP35     = wfp35;
  myWeightStatFP35 = wsfp35;
  myWeightFF35     = wff35;
  myWeightStatFF35 = wsff35;
  myWeightPP35     = wpp35;
  myWeightStatPP35 = wspp35;
  //
  myWeightFP50     = wfp50;
  myWeightStatFP50 = wsfp50;
  myWeightFF50     = wff50;
  myWeightStatFF50 = wsff50;
  myWeightPP50     = wpp50;
  myWeightStatPP50 = wspp50;
  // 
  myWeightFPQCD     = wfpQCD;
  myWeightStatFPQCD = wsfpQCD;
  myWeightFFQCD     = wffQCD;
  myWeightStatFFQCD = wsffQCD;
  myWeightPPQCD     = wppQCD;
  myWeightStatPPQCD = wsppQCD;
}


void RedHiggsTree::fillKinematics(float pxTkMet, float pyTkMet, float pzTkMet,
                                  float pxLeadJet[3], float pyLeadJet[3], float pzLeadJet[3], 
                                  float pxSecJet[3], float pySecJet[3], float pzSecJet[3],    
                                  float pxL1, float pyL1, float pzL1,
                                  float pxL2, float pyL2, float pzL2,
                                  int ch[2], float lh[2], float iso[2], int majority[2], float bdt[2],
                                  TLorentzVector *jetSum, TLorentzVector *uncorrJetSum, TVector3 *pfmet) {

  myPxTkMet = pxTkMet;
  myPyTkMet = pyTkMet;
  myPzTkMet = pzTkMet;

  for(int jes=0; jes<3; jes++) {
    myPxLeadJet[jes] = pxLeadJet[jes];
    myPyLeadJet[jes] = pyLeadJet[jes];
    myPzLeadJet[jes] = pzLeadJet[jes];
    myPxSecondJet[jes] = pxSecJet[jes];
    myPySecondJet[jes] = pySecJet[jes];
    myPzSecondJet[jes] = pzSecJet[jes];
  }

  myPxL1 = pxL1;
  myPyL1 = pyL1;
  myPzL1 = pzL1;
  myPxL2 = pxL2;
  myPyL2 = pyL2;
  myPzL2 = pzL2;

  myJetsSum = jetSum;
  myUncorrJetsSum = uncorrJetSum;
  myPfMet = pfmet;
  for(int i=0; i<2; i++) {
    myLepCharge[i] = ch[i];
    myIso[i] = iso[i];
    myEleLh[i] = lh[i];
    myEleBdt[i] = bdt[i];
    myMajority[i] = majority[i];
  }
}

void RedHiggsTree::fillSystematics(float scE[2], float r9[2], float ene1, float ene2, int ty1, int ty2, 
                                   TVector3 *metFromJets, TVector3 *pfMetUp, TVector3 *pfMetDown, float mtUp, float mtDown) {

  for(int i=0; i<2; i++) {
    myScEnergy[i] = scE[i];
    myR9[i] = r9[i];
  }

  myEneL1  = ene1;
  myEneL2  = ene2;
  myTypeL1 = ty1;
  myTypeL2 = ty2;
  myMetFromJets = metFromJets;
  myPfMetUp = pfMetUp;
  myPfMetDown = pfMetDown;  
  myMtUp = mtUp;
  myMtDown = mtDown;
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

void RedHiggsTree::fillRunInfos(int run, int lumi, int event, float puweight, bool HLT) {

  myRun = run;
  myLS = lumi;
  myEvent = event;
  myPUWeight = puweight;
  myHLT = HLT;
}

void RedHiggsTree::fillMetStudies(float projPF, float projTk, float signPFMet, float signChMet, float m_MTRcha ) {
  
  myProjPFMet        = projPF;  
  myProjPFChargedMet = projTk;  
  mySignPFMet        = signPFMet;
  mySignPFChargedMet = signChMet;
  myMTRchargedMet    = m_MTRcha;
}
