#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

int fullFormat = 1;

void setReducedFormat() { fullFormat = 1; }

void addWeights(const char* filename, float weight, int processId, int finalstate) {

  cout << "Adding weight branch to file " << filename << " with weight " << weight << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  if ( treeOrig ) {
    int nentriesOrig = treeOrig->GetEntries();

    TFile *fileNew = TFile::Open(filename,"recreate");
    TTree *treeNew = new TTree("T1","tree with only selected events");

    TString skimFile(filename);
    skimFile.ReplaceAll("merged","merged_skim");
    TFile *fileNewSkim = TFile::Open(skimFile.Data(),"recreate");
    TTree *treeNewSkim = new TTree("T1","tree with only selected events");

    std::vector<TTree*> trees; 
    trees.push_back(treeNew);
    trees.push_back(treeNewSkim);
    
    // add also a branch with jet category (1 for njets=0, -1 for njets=1: useful for the fit)
    // and a branch with float final selection bool (for roofit)
    Int_t           run;
    // Int_t           lumi;
    Int_t           event;
    Float_t         puweight;
    Float_t         met;
    Float_t         pfMet;
    Float_t         caloMet;
    Float_t         projMet;
    Float_t         deltaPhi;
    Float_t         deltaR;
    Float_t         transvMass;
    Float_t         eleInvMass;
    Float_t         maxPtEle;
    Float_t         minPtEle;
    Float_t         detaLeptons;
    Bool_t          finalLeptons;
    Bool_t          jetVeto;
    Bool_t          uncorrJetVeto;
    Bool_t          preDeltaPhi;
    Bool_t          finalSelection;
    Bool_t          step[25];
    Bool_t          hlt;
    Float_t         KFactor;
    Bool_t          promptDecay;
    // Float_t         maxPtLh;
    // Float_t         minPtLh;
    Int_t           njets;
    Int_t           nuncorrjets;
    Int_t           nVtx;
    Float_t         dxyEVT;
    Float_t         dszEVT;
    Float_t         bTagTrackCount;
    Float_t         bTagImpPar;
    Float_t         bTagSecVertex;
    Float_t         leadingJetBTagTrackCount;
    Float_t         subleadingJetBTagTrackCount;
    Float_t         pt[2];
    Float_t         eta[2];
    Float_t         deta[2];
    Float_t         dphi[2];
    Float_t         hoe[2];
    Float_t         see[2];
    Float_t         scEnergy[2];
    Float_t         R9[2];
    Int_t           matched[2];
    Float_t         pxTkMet;
    Float_t         pyTkMet;
    Float_t         pzTkMet;
    Float_t         pxLeadJet;
    Float_t         pyLeadJet;
    Float_t         pzLeadJet;
    Float_t         pxSecondJet;
    Float_t         pySecondJet;
    Float_t         pzSecondJet;
    Float_t         pxL1;
    Float_t         pyL1;
    Float_t         pzL1;
    Float_t         pxL2;
    Float_t         pyL2;
    Float_t         pzL2;
    Float_t         eneL1;
    Float_t         eneL2;
    Int_t           typeL1;
    Int_t           typeL2;
    Int_t           nSoftMu;
    Float_t         mtr;
    Float_t         mr;
    Float_t         gammamr;

    TLorentzVector *sumJetsV4 = 0;
    TLorentzVector *uncorrSumJetsV4 = 0;
    TVector3        *pfmetV = 0;

    treeOrig->SetBranchAddress("run", &run);
    // treeOrig->SetBranchAddress("lumi", &lumi);
    treeOrig->SetBranchAddress("event", &event);
    treeOrig->SetBranchAddress("puweight", &puweight);
    treeOrig->SetBranchAddress("hlt", &hlt);
    treeOrig->SetBranchAddress("met", &met);  // default MET is tcMET for WW
    treeOrig->SetBranchAddress("pfMet", &pfMet);
    treeOrig->SetBranchAddress("caloMet", &caloMet);
    treeOrig->SetBranchAddress("projMet", &projMet);
    treeOrig->SetBranchAddress("deltaPhi", &deltaPhi);
    treeOrig->SetBranchAddress("deltaR", &deltaR);
    treeOrig->SetBranchAddress("transvMass", &transvMass);
    treeOrig->SetBranchAddress("eleInvMass", &eleInvMass);
    treeOrig->SetBranchAddress("maxPtEle", &maxPtEle);
    treeOrig->SetBranchAddress("minPtEle", &minPtEle);
    treeOrig->SetBranchAddress("detaLeptons", &detaLeptons);
    treeOrig->SetBranchAddress("nVtx", &nVtx);
    treeOrig->SetBranchAddress("finalLeptons", &finalLeptons);
    treeOrig->SetBranchAddress("jetVeto", &jetVeto);
    treeOrig->SetBranchAddress("uncorrJetVeto", &uncorrJetVeto);
    treeOrig->SetBranchAddress("preDeltaPhi", &preDeltaPhi);
    treeOrig->SetBranchAddress("finalSelection", &finalSelection);
    treeOrig->SetBranchAddress("step", step);
    treeOrig->SetBranchAddress("KFactor", &KFactor);
    treeOrig->SetBranchAddress("promptDecay", &promptDecay);
    // treeOrig->SetBranchAddress("maxPtLh", &maxPtLh);
    // treeOrig->SetBranchAddress("minPtLh", &minPtLh);
    treeOrig->SetBranchAddress("njets", &njets);
    treeOrig->SetBranchAddress("nuncorrjets", &nuncorrjets);
    treeOrig->SetBranchAddress("dxyEVT", &dxyEVT);
    treeOrig->SetBranchAddress("dszEVT", &dszEVT);
    treeOrig->SetBranchAddress("bTagTrackCount", &bTagTrackCount);
    treeOrig->SetBranchAddress("bTagImpPar", &bTagImpPar);
    treeOrig->SetBranchAddress("bTagSecVertex", &bTagSecVertex);
    treeOrig->SetBranchAddress("leadingJetBTagTrackCount", &leadingJetBTagTrackCount);
    treeOrig->SetBranchAddress("subleadingJetBTagTrackCount", &subleadingJetBTagTrackCount);
    treeOrig->SetBranchAddress("pt", pt);
    treeOrig->SetBranchAddress("eta", eta);
    treeOrig->SetBranchAddress("deta", deta);
    treeOrig->SetBranchAddress("dphi", dphi);
    treeOrig->SetBranchAddress("hoe", hoe);
    treeOrig->SetBranchAddress("see", see);
    treeOrig->SetBranchAddress("scEnergy", scEnergy);
    treeOrig->SetBranchAddress("R9", R9);
    treeOrig->SetBranchAddress("matched", matched);
    treeOrig->SetBranchAddress("pxTkMet", &pxTkMet);
    treeOrig->SetBranchAddress("pyTkMet", &pyTkMet);
    treeOrig->SetBranchAddress("pzTkMet", &pzTkMet);
    treeOrig->SetBranchAddress("pxLeadJet", &pxLeadJet);
    treeOrig->SetBranchAddress("pyLeadJet", &pyLeadJet);
    treeOrig->SetBranchAddress("pzLeadJet", &pzLeadJet);
    treeOrig->SetBranchAddress("pxSecondJet", &pxSecondJet);
    treeOrig->SetBranchAddress("pySecondJet", &pySecondJet);
    treeOrig->SetBranchAddress("pzSecondJet", &pzSecondJet);
    treeOrig->SetBranchAddress("pxL1", &pxL1);
    treeOrig->SetBranchAddress("pyL1", &pyL1);
    treeOrig->SetBranchAddress("pzL1", &pzL1);
    treeOrig->SetBranchAddress("pxL2", &pxL2);
    treeOrig->SetBranchAddress("pyL2", &pyL2);
    treeOrig->SetBranchAddress("pzL2", &pzL2);
    treeOrig->SetBranchAddress("eneL1", &eneL1);
    treeOrig->SetBranchAddress("eneL2", &eneL2);
    treeOrig->SetBranchAddress("typeL1", &typeL1);
    treeOrig->SetBranchAddress("typeL2", &typeL2);
    treeOrig->SetBranchAddress("nSoftMu", &nSoftMu);
    treeOrig->SetBranchAddress("mtr", &mtr);
    treeOrig->SetBranchAddress("mr", &mr);
    treeOrig->SetBranchAddress("gammamr", &gammamr);
    treeOrig->SetBranchAddress("sumJetsV4", &sumJetsV4);
    treeOrig->SetBranchAddress("uncorrSumJetsV4", &uncorrSumJetsV4);
    treeOrig->SetBranchAddress("pfmetV", &pfmetV);

    // 
    Float_t pt_1,       pt_2;
    Float_t eta_1,      eta_2;
    Float_t deta_1,     deta_2;
    Float_t dphi_1,     dphi_2;
    Float_t hoe_1,      hoe_2;
    Float_t see_1,      see_2;
    Float_t R9_1,       R9_2;
    Float_t scEnergy_1, scEnergy_2;
    Int_t   matched_1,  matched_2;

    // convert the booleans into integers (to insert in RooDataset)
    Int_t         i_finalLeptons;
    Int_t         i_jetVeto;
    Int_t         i_uncorrJetVeto;
    Int_t         i_preDeltaPhi;
    Int_t         i_finalSelection;
    Int_t         i_promptDecay;
    Int_t         i_WWSel;
    Int_t         i_WWSel1j;
    Int_t         i_hlt;
    float deltaPhi_LL;    
    float deltaPhi_LL_MET;
    float deltaPhi_LLJ1_MET;
    float deltaPhi_LL_JET1;
    float deltaPhi_LL_JET2;
    float deltaPhi_MET_JET1;
    float deltaPhi_MET_JET2;
    float deltaPhi_LL_JJ;
    float leadingJetPt;
    float secondJetPt;
    float L1eta, L1phi;
    float L2eta, L2phi;
    float dileptonPt;
    float R;
    float dgammamr;
    float jetcat = 0;
    float consecevent = -1;

    for(int i=0; i<(int)trees.size();i++) {
      TTree *theTreeNew = trees[i];

    // the selected final state: ee=0, mm=1, em=2
    theTreeNew->Branch("finalstate", &finalstate, "finalstate/I");

    // one integer containing the process identifier (for MC, 0 for data)
    theTreeNew->Branch("process", &processId, "process/I");

    // copy branches
    theTreeNew->Branch("run", &run, "run/I");
    // theTreeNew->Branch("lumi", &lumi, "lumi/I");
    theTreeNew->Branch("event", &event, "event/I");
    theTreeNew->Branch("puweight", &puweight, "puweight/F");
    theTreeNew->Branch("met", &met, "met/F");  // default MET is tcMET for WW
    theTreeNew->Branch("hlt", &i_hlt, "hlt/I");
    theTreeNew->Branch("pfMet", &pfMet, "pfMet/F");
    theTreeNew->Branch("caloMet", &caloMet, "caloMet/F");
    theTreeNew->Branch("projMet", &projMet, "projMet/F");
    theTreeNew->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
    theTreeNew->Branch("deltaR", &deltaR, "deltaR/F");
    theTreeNew->Branch("eleInvMass", &eleInvMass, "eleInvMass/F");
    theTreeNew->Branch("transvMass", &transvMass, "transvMass/F");
    theTreeNew->Branch("maxPtEle", &maxPtEle, "maxPtEle/F");
    theTreeNew->Branch("minPtEle", &minPtEle, "minPtEle/F");
    theTreeNew->Branch("detaLeptons", &detaLeptons, "detaLeptons/F");
    theTreeNew->Branch("nVtx", &nVtx, "nVtx/I");
    theTreeNew->Branch("finalLeptons", &i_finalLeptons, "finalLeptons/I");
    theTreeNew->Branch("jetVeto", &i_jetVeto, "jetVeto/I");
    theTreeNew->Branch("uncorrJetVeto", &i_uncorrJetVeto, "uncorrJetVeto/I");
    theTreeNew->Branch("preDeltaPhi", &i_preDeltaPhi, "preDeltaPhi/I");
    theTreeNew->Branch("finalSelection", &i_finalSelection, "finalSelection/I");
    theTreeNew->Branch("WWSel", &i_WWSel, "WWSel/I");
    theTreeNew->Branch("WWSel1j", &i_WWSel1j, "WWSel1j/I");
    theTreeNew->Branch("KFactor", &KFactor, "KFactor/F");
    theTreeNew->Branch("promptDecay", &i_promptDecay, "promptDecay/I");
    // theTreeNew->Branch("maxPtLh", &maxPtLh, "maxPtLh/F");
    // theTreeNew->Branch("minPtLh", &minPtLh, "minPtLh/F");
    theTreeNew->Branch("njets", &njets, "njets/I");
    theTreeNew->Branch("nuncorrjets", &nuncorrjets, "nuncorrjets/I");
    theTreeNew->Branch("dxyEVT", &dxyEVT, "dxyEVT/F");
    theTreeNew->Branch("dszEVT", &dszEVT, "dszEVT/F");
    theTreeNew->Branch("bTagTrackCount", &bTagTrackCount, "bTagTrackCount/F");
    theTreeNew->Branch("bTagImpPar", &bTagImpPar, "bTagImpPar/F");
    theTreeNew->Branch("bTagSecVertex", &bTagSecVertex, "bTagSecVertex/F");
    theTreeNew->Branch("leadingJetBTagTrackCount", &leadingJetBTagTrackCount, "leadingJetBTagTrackCount/F");
    theTreeNew->Branch("subleadingJetBTagTrackCount", &subleadingJetBTagTrackCount, "subleadingJetBTagTrackCount/F");
    theTreeNew->Branch("pt1", &pt_1, "pt1/F");
    theTreeNew->Branch("eta1", &eta_1, "eta1/F");
    theTreeNew->Branch("deta1", &deta_1, "deta1/F");
    theTreeNew->Branch("dphi1", &dphi_1, "dphi1/F");
    theTreeNew->Branch("hoe1", &hoe_1, "hoe1/F");
    theTreeNew->Branch("see1", &see_1, "see1/F");
    theTreeNew->Branch("R91", &R9_1, "R91/F");
    theTreeNew->Branch("scEnergy1", &scEnergy_1, "scEnergy1/F");
    theTreeNew->Branch("matched1", &matched_1, "matched1/I");
    theTreeNew->Branch("pt2", &pt_2, "pt2/F");
    theTreeNew->Branch("eta2", &eta_2, "eta2/F");
    theTreeNew->Branch("deta2", &deta_2, "deta2/F");
    theTreeNew->Branch("dphi2", &dphi_2, "dphi2/F");
    theTreeNew->Branch("hoe2", &hoe_2, "hoe2/F");
    theTreeNew->Branch("see2", &see_2, "see2/F");
    theTreeNew->Branch("R92", &R9_2, "R92/F");
    theTreeNew->Branch("scEnergy2", &scEnergy_2, "scEnergy2/F");
    theTreeNew->Branch("matched2", &matched_2, "matched2/I");
    theTreeNew->Branch("step", step, "step[25]/O");
    theTreeNew->Branch("mtr", &mtr, "mtr/F");
    theTreeNew->Branch("mr", &mr, "mr/F");
    theTreeNew->Branch("dgammamr", &dgammamr, "dgammamr/F");
    theTreeNew->Branch("R", &R, "R/F");
    if(fullFormat) {
      theTreeNew->Branch("pxTkMet", &pxTkMet, "pxTkMet/F");
      theTreeNew->Branch("pyTkMet", &pyTkMet, "pyTkMet/F");
      theTreeNew->Branch("pzTkMet", &pzTkMet, "pzTkMet/F");
      theTreeNew->Branch("pxLeadJet", &pxLeadJet, "pxLeadJet/F");
      theTreeNew->Branch("pyLeadJet", &pyLeadJet, "pyLeadJet/F");
      theTreeNew->Branch("pzLeadJet", &pzLeadJet, "pzLeadJet/F");
      theTreeNew->Branch("pxSecondJet", &pxSecondJet, "pxSecondJet/F");
      theTreeNew->Branch("pySecondJet", &pySecondJet, "pySecondJet/F");
      theTreeNew->Branch("pzSecondJet", &pzSecondJet, "pzSecondJet/F");
      theTreeNew->Branch("pxL1", &pxL1, "pxL1/F");
      theTreeNew->Branch("pyL1", &pyL1, "pyL1/F");
      theTreeNew->Branch("pzL1", &pzL1, "pzL1/F");
      theTreeNew->Branch("eneL1",  &eneL1,  "eneL1/F");
      theTreeNew->Branch("typeL1", &typeL1, "typeL1/I");
      theTreeNew->Branch("pxL2", &pxL2, "pxL2/F");
      theTreeNew->Branch("pyL2", &pyL2, "pyL2/F");
      theTreeNew->Branch("pzL2", &pzL2, "pzL2/F");
      theTreeNew->Branch("eneL2",  &eneL2,  "eneL2/F");
      theTreeNew->Branch("typeL2", &typeL2, "typeL2/I");
      theTreeNew->Branch("deltaPhi_LL", &deltaPhi_LL, "deltaPhi_LL/F");
      theTreeNew->Branch("deltaPhi_LL_MET", &deltaPhi_LL_MET, "deltaPhi_LL_MET/F");
      theTreeNew->Branch("deltaPhi_LLJ1_MET", &deltaPhi_LLJ1_MET, "deltaPhi_LLJ1_MET/F");
      theTreeNew->Branch("deltaPhi_LL_JET1", &deltaPhi_LL_JET1, "deltaPhi_LL_JET1/F");
      theTreeNew->Branch("deltaPhi_LL_JET2", &deltaPhi_LL_JET2, "deltaPhi_LL_JET2/F");
      theTreeNew->Branch("deltaPhi_MET_JET1", &deltaPhi_MET_JET1, "deltaPhi_MET_JET1/F");
      theTreeNew->Branch("deltaPhi_MET_JET2", &deltaPhi_MET_JET2, "deltaPhi_MET_JET2/F");
      theTreeNew->Branch("deltaPhi_LL_JJ", &deltaPhi_LL_JJ, "deltaPhi_LL_JJ/F");
      theTreeNew->Branch("dileptonPt", &dileptonPt, "dileptonPt/F");
      theTreeNew->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
      theTreeNew->Branch("secondJetPt", &secondJetPt, "secondJetPt/F");
      theTreeNew->Branch("L1eta", &L1eta, "L1eta/F");
      theTreeNew->Branch("L1phi", &L1phi, "L1phi/F");
      theTreeNew->Branch("L2eta", &L2eta, "L2eta/F");
      theTreeNew->Branch("L2phi", &L2phi, "L2phi/F");
      theTreeNew->Branch("nSoftMu", &nSoftMu, "nSoftMu/I");
      theTreeNew->Branch("sumJetsV4", "TLorentzVector", &sumJetsV4);
      theTreeNew->Branch("uncorrSumJetsV4", "TLorentzVector", &uncorrSumJetsV4);
      theTreeNew->Branch("pfmetV", "TVector3", &pfmetV);
    }

    theTreeNew->Branch("jetcat", &jetcat,  "jetcat/F");

    theTreeNew->Branch("consecevent", &consecevent, "consecevent/F");
    theTreeNew->Branch("weight", &weight,  "weight/F");
    }

    int j =0;

    for(int i=0; i<nentriesOrig; i++) {
      if (i%10000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl;
      treeOrig->GetEntry(i);

      if (njets==0) jetcat = 1;
      else if(njets==1) jetcat = -1;
      else jetcat = -2;

      pt_1       = pt[0];
      eta_1      = eta[0];
      deta_1     = deta[0];
      dphi_1     = dphi[0];
      hoe_1      = hoe[0];
      see_1      = see[0];
      scEnergy_1 = scEnergy[0];
      R9_1       = R9[0];
      matched_1  = matched[0];
      pt_2       = pt[1];
      eta_2      = eta[1];
      deta_2     = deta[1];
      dphi_2     = dphi[1];
      hoe_2      = hoe[1];
      see_2      = see[1];
      scEnergy_2 = scEnergy[1];
      R9_2       = R9[1];
      matched_2  = matched[1];
      R = mtr/mr;
      dgammamr = 2*gammamr;

      // consider only events with 0 or 1 jet
      // and fit variables within fit range
      //     if (njets<=1 && 
      //         met>=0 && met<=200 &&
      //         deltaPhi>=0 && deltaPhi<=180 &&
      //         maxPtEle>=20 && maxPtEle<=200 &&
      //         eleInvMass>=12 && eleInvMass<=150 &&
      //         bTagImpPar>=-1001 && bTagImpPar<=2 &&
      //         finalSelection) {
      
      i_finalLeptons = (finalLeptons) ? 1 : 0;
      i_jetVeto = (jetVeto) ? 1 : 0;
      i_uncorrJetVeto = (uncorrJetVeto) ? 1 : 0;
      i_preDeltaPhi = (preDeltaPhi) ? 1 : 0;
      i_finalSelection = (finalSelection) ? 1 : 0;
      i_promptDecay = (promptDecay) ? 1 : 0;
      i_WWSel = (step[13]) ? 1 : 0;
      i_WWSel1j = (step[19] && njets==1) ? 1 : 0;
      i_hlt = (hlt) ? 1 : 0;

      if (finalLeptons && fullFormat) {
        TVector3 TV_L1( pxL1, pyL1, pzL1 );
        TVector3 TV_L2( pxL2, pyL2, pzL2 );
        TVector3 TV_L1p2 = TV_L1 + TV_L2;
        TVector3 TV_met( pxTkMet, pyTkMet, pzTkMet );
        TVector3 TV_jet1( pxLeadJet, pyLeadJet, pzLeadJet );
        TVector3 TV_jet2( pxSecondJet, pySecondJet, pzSecondJet );
        TVector3 TV_J1p2  = TV_jet1 + TV_jet2;
        TVector3 TV_L12pJ1 = TV_L1p2 + TV_jet1;
        deltaPhi_LL       = (180./3.14) * TV_L1.DeltaPhi(TV_L2);
        deltaPhi_LL_MET   = (180./3.14) * TV_met.DeltaPhi(TV_L1p2);

        if(TV_jet1.Pt()>15) {
          deltaPhi_LLJ1_MET = (180./TMath::Pi()) * fabs(TV_met.DeltaPhi(TV_L12pJ1));   
          deltaPhi_MET_JET1 = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_met));
          deltaPhi_LL_JET1  = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_L1p2));
        } else {
          deltaPhi_LLJ1_MET = -1.;
          deltaPhi_MET_JET1 = -1.;
          deltaPhi_LL_JET1 = -1.;
        }

        if(TV_jet2.Pt()>15) {
          deltaPhi_MET_JET2 = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_met));
          deltaPhi_LL_JET2  = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_L1p2));
        } else {
          deltaPhi_MET_JET2 = -1.;
          deltaPhi_LL_JET2 = -1.;
        }

        if(TV_jet1.Pt()>15 && TV_jet2.Pt()>15) {
          deltaPhi_LL_JJ    = (180./3.14) * TV_L1p2.DeltaPhi(TV_J1p2);
        } else {
          deltaPhi_LL_JJ = -1.;
        }

        leadingJetPt = sqrt(pxLeadJet*pxLeadJet + pyLeadJet*pyLeadJet);
        secondJetPt  = sqrt(pxSecondJet*pxSecondJet + pySecondJet*pySecondJet);
        dileptonPt   = TV_L1p2.Pt();
        L1eta = TV_L1.Eta();
        L2eta = TV_L2.Eta();
        L1phi = TV_L1.Phi();
        L2phi = TV_L2.Phi();
      } else { 
        L1eta = 100.;
        L2eta = 100.;
        L1phi = 100.;
        L2phi = 100.;
        deltaPhi_LL = -9999.;
        deltaPhi_LL_MET = -9999.;
        deltaPhi_LLJ1_MET = -9999.;
        deltaPhi_LL_JET1 = -9999.;
        deltaPhi_LL_JET2 = -9999.;
        deltaPhi_MET_JET1 = -9999.;
        deltaPhi_MET_JET2 = -9999.;
	deltaPhi_LL_JJ = -9999.;	
        leadingJetPt = -9999.;
        secondJetPt = -9999.;
        dileptonPt = -9999.;
      } 

      consecevent = (float)j;
      if(finalLeptons) {
        if(processId>0) { // MC
          treeNew->Fill();
          if(i_WWSel || i_WWSel1j) treeNewSkim->Fill();
        } else { // data: apply the trigger only for the single lepton trigger datasets
          if((processId==-1) || (processId==-2 && hlt)) {
            treeNew->Fill();
	    if(i_WWSel || i_WWSel1j) treeNewSkim->Fill();
          }
        }
      }
      j++;
      //    }
    }
  
    fileNew->cd();
    treeNew->Write();
    fileNew->Close();

    fileNewSkim->cd();
    treeNewSkim->Write();
    fileNewSkim->Close();

    fileOrig->cd();
    fileOrig->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }
}
