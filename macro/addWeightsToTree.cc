#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

int FRWeights = 0;

float GetProjectedMet(TVector3 met, TVector3 p1, TVector3 p2);
float calcMT(TVector3 met, TVector3 lepton);
void addFRWeights() { FRWeights = 1;};

void addWeights(const char* filename, float baseW, int processId, int finalstate) {

  cout << "Adding weight branch to file " << filename << " with weight " << baseW << endl;

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
    Int_t           lumi;
    Int_t           event;
    Float_t         puweight;
    Float_t         met;
    Float_t         pfMet;
    Float_t         caloMet;
    Float_t         projMet;
    Float_t         deltaPhi;
    Float_t         deltaR;
    Float_t         transvMass, transvMassUp, transvMassDown;
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
    Float_t         pxLeadJet[3];
    Float_t         pyLeadJet[3];
    Float_t         pzLeadJet[3];
    Float_t         pxSecondJet[3];
    Float_t         pySecondJet[3];
    Float_t         pzSecondJet[3];
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
    Int_t           numExtraLep;
    // lepton quantities (independent from lepton flavor)
    // 0 leading, 1 is trailing
    Int_t           ch[2];
    Float_t         lh[2];
    Float_t         iso[2];

    TLorentzVector *sumJetsV4 = 0;
    TLorentzVector *uncorrSumJetsV4 = 0;
    TVector3        *pfmetV = 0;
    float pfmetX, pfmetY;

    // Fake rate weights (only for W+jets selection)
    Float_t         weightFP;
    Float_t         weightStatFP;
    Float_t         weightFF;
    Float_t         weightStatFF;
    Float_t         weightPP;
    Float_t         weightStatPP;

    treeOrig->SetBranchAddress("run", &run);
    treeOrig->SetBranchAddress("lumi", &lumi);
    treeOrig->SetBranchAddress("event", &event);
    treeOrig->SetBranchAddress("puweight", &puweight);
    treeOrig->SetBranchAddress("hlt", &hlt);
    treeOrig->SetBranchAddress("met", &met);
    treeOrig->SetBranchAddress("pfMet", &pfMet);
    treeOrig->SetBranchAddress("caloMet", &caloMet);
    treeOrig->SetBranchAddress("projMet", &projMet);
    treeOrig->SetBranchAddress("deltaPhi", &deltaPhi);
    treeOrig->SetBranchAddress("deltaR", &deltaR);
    treeOrig->SetBranchAddress("transvMass", &transvMass);
    treeOrig->SetBranchAddress("transvMassUp", &transvMassUp);
    treeOrig->SetBranchAddress("transvMassDown", &transvMassDown);
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
    treeOrig->SetBranchAddress("pxLeadJet", pxLeadJet);
    treeOrig->SetBranchAddress("pyLeadJet", pyLeadJet);
    treeOrig->SetBranchAddress("pzLeadJet", pzLeadJet);
    treeOrig->SetBranchAddress("pxSecondJet", pxSecondJet);
    treeOrig->SetBranchAddress("pySecondJet", pySecondJet);
    treeOrig->SetBranchAddress("pzSecondJet", pzSecondJet);
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
    treeOrig->SetBranchAddress("numExtraLep", &numExtraLep);
    treeOrig->SetBranchAddress("mtr", &mtr);
    treeOrig->SetBranchAddress("mr", &mr);
    treeOrig->SetBranchAddress("gammamr", &gammamr);
    treeOrig->SetBranchAddress("sumJetsV4", &sumJetsV4);
    treeOrig->SetBranchAddress("uncorrSumJetsV4", &uncorrSumJetsV4);
    treeOrig->SetBranchAddress("pfmetV", &pfmetV);
    treeOrig->SetBranchAddress("ch", ch);
    treeOrig->SetBranchAddress("lh", lh);
    treeOrig->SetBranchAddress("iso", iso);
    if(FRWeights) {
      treeOrig->SetBranchAddress("weightFP", &weightFP);
      treeOrig->SetBranchAddress("weightStatFP", &weightStatFP);
      treeOrig->SetBranchAddress("weightFF", &weightFF);
      treeOrig->SetBranchAddress("weightStatFF", &weightStatFF);
      treeOrig->SetBranchAddress("weightPP", &weightPP);
      treeOrig->SetBranchAddress("weightStatPP", &weightStatPP);
    }

    // electron ID (only filled for electrons)
    Float_t pt_1,       pt_2;
    Float_t eta_1,      eta_2;
    Float_t deta_1,     deta_2;
    Float_t dphi_1,     dphi_2;
    Float_t hoe_1,      hoe_2;
    Float_t see_1,      see_2;
    Float_t R9_1,       R9_2;
    Float_t scEnergy_1, scEnergy_2;
    Int_t   matched_1,  matched_2;
    
    Float_t ellh;

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

    // other kinematics
    float deltaPhi_LL;    
    float deltaPhi_LL_MET;
    float deltaPhi_LLJ1_MET;
    float deltaPhi_LL_JET1;
    float deltaPhi_LL_JET2;
    float deltaPhi_MET_JET1;
    float deltaPhi_MET_JET2;
    float deltaPhi_LL_JJ;
    float dphilmet1, dphilmet2;
    float jetpt1, jeteta1, jetphi1;
    float jetpt2, jeteta2, jetphi2;
    float tkMET;
    float pmet, pmet2;
    float L1eta, L1phi;
    float L2eta, L2phi;
    float dileptonPt;
    float mtw1, mtw2;
    float R;
    float dgammamr;
    float jetcat = 0;
    float consecevent = -1;

    // vetoes
    int zveto, bveto_ip, bveto_mu, bveto;

    // additional (dummy for the moment)
    Float_t effW = 1.0;
    Float_t triggW = 1.0;

    for(int i=0; i<(int)trees.size();i++) {
      TTree *theTreeNew = trees[i];

      // the selected final state: ee=0, mm=1, em=2
      theTreeNew->Branch("channel", &finalstate, "channel/I");

      // one integer containing the process identifier (for MC, 0 for data)
      theTreeNew->Branch("process", &processId, "process/I");

      // copy branches
      theTreeNew->Branch("run", &run, "run/I");
      theTreeNew->Branch("lumi", &lumi, "lumi/I");
      theTreeNew->Branch("event", &event, "event/I");
      theTreeNew->Branch("puW", &puweight, "puW/F");
      theTreeNew->Branch("effW", &effW, "effW/F");
      theTreeNew->Branch("triggW", &triggW, "triggW/F");
      theTreeNew->Branch("met", &pfMet, "met/F");
      theTreeNew->Branch("met2", &tkMET, "met2/F");
      theTreeNew->Branch("hlt", &i_hlt, "hlt/I");
      theTreeNew->Branch("met", &pfMet, "met/F");
      theTreeNew->Branch("caloMet", &caloMet, "caloMet/F");
      theTreeNew->Branch("pmet", &pmet, "pmet/F");
      theTreeNew->Branch("pmet2", &pmet2, "pmet2/F");
      theTreeNew->Branch("mpmet", &projMet, "mpmet/F");
      theTreeNew->Branch("dphill", &deltaPhi, "dphill/F");
      theTreeNew->Branch("drll", &deltaR, "drll/F");
      theTreeNew->Branch("mll", &eleInvMass, "mll/F");
      theTreeNew->Branch("mth", &transvMass, "mth/F");
      theTreeNew->Branch("mthUp", &transvMassUp, "mthUp/F");
      theTreeNew->Branch("mthDown", &transvMassDown, "mthDown/F");
      theTreeNew->Branch("mtw1", &mtw1, "mtw1/F");
      theTreeNew->Branch("mtw2", &mtw2, "mtw2/F");
      theTreeNew->Branch("pt1", &maxPtEle, "pt1/F");
      theTreeNew->Branch("eta1", &L1eta, "eta1/F");
      theTreeNew->Branch("phi1", &L1phi, "phi1/F");
      theTreeNew->Branch("pt2", &minPtEle, "pt2/F");
      theTreeNew->Branch("eta2", &L2eta, "eta2/F");
      theTreeNew->Branch("phi2", &L2phi, "phi2/F");
      theTreeNew->Branch("detaLeptons", &detaLeptons, "detaLeptons/F");
      theTreeNew->Branch("nvtx", &nVtx, "nvtx/I");
      theTreeNew->Branch("finalLeptons", &i_finalLeptons, "finalLeptons/I");
      theTreeNew->Branch("jetVeto", &i_jetVeto, "jetVeto/I");
      theTreeNew->Branch("uncorrJetVeto", &i_uncorrJetVeto, "uncorrJetVeto/I");
      theTreeNew->Branch("preDeltaPhi", &i_preDeltaPhi, "preDeltaPhi/I");
      theTreeNew->Branch("finalSelection", &i_finalSelection, "finalSelection/I");
      theTreeNew->Branch("WWSel", &i_WWSel, "WWSel/I");
      theTreeNew->Branch("WWSel1j", &i_WWSel1j, "WWSel1j/I");
      theTreeNew->Branch("kfW", &KFactor, "kfW/F");
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
      theTreeNew->Branch("step", step, "step[25]/O");
      theTreeNew->Branch("zveto", &zveto, "zveto/I");
      theTreeNew->Branch("bveto_ip", &bveto_ip, "bveto_ip/I");
      theTreeNew->Branch("bveto_mu", &bveto_mu, "bveto_mu/I");
      theTreeNew->Branch("bveto", &bveto, "bveto/I");
      theTreeNew->Branch("mtr", &mtr, "mtr/F");
      theTreeNew->Branch("mr", &mr, "mr/F");
      theTreeNew->Branch("dgammamr", &dgammamr, "dgammamr/F");
      theTreeNew->Branch("R", &R, "R/F");
      theTreeNew->Branch("pt1_eleid", &pt_1, "pt1_eleid/F");
      theTreeNew->Branch("eta1_eleid", &eta_1, "eta1_eleid/F");
      theTreeNew->Branch("deta1_eleid", &deta_1, "deta1_eleid/F");
      theTreeNew->Branch("dphi1_eleid", &dphi_1, "dphi1_eleid/F");
      theTreeNew->Branch("hoe1_eleid", &hoe_1, "hoe1_eleid/F");
      theTreeNew->Branch("see1_eleid", &see_1, "see1_eleid/F");
      theTreeNew->Branch("R91_eleid", &R9_1, "R91_eleid/F");
      theTreeNew->Branch("scEnergy1_eleid", &scEnergy_1, "scEnergy1_eleid/F");
      theTreeNew->Branch("matched1_eleid", &matched_1, "matched1_eleid/I");
      theTreeNew->Branch("pt2_eleid", &pt_2, "pt2_eleid/F");
      theTreeNew->Branch("eta2_eleid", &eta_2, "eta2_eleid/F");
      theTreeNew->Branch("deta2_eleid", &deta_2, "deta2_eleid/F");
      theTreeNew->Branch("dphi2_eleid", &dphi_2, "dphi2_eleid/F");
      theTreeNew->Branch("hoe2_eleid", &hoe_2, "hoe2_eleid/F");
      theTreeNew->Branch("see2_eleid", &see_2, "see2_eleid/F");
      theTreeNew->Branch("R92_eleid", &R9_2, "R92_eleid/F");
      theTreeNew->Branch("scEnergy2_eleid", &scEnergy_2, "scEnergy2_eleid/F");
      theTreeNew->Branch("matched2_eleid", &matched_2, "matched2_eleid/I");
      theTreeNew->Branch("pxTkMet", &pxTkMet, "pxTkMet/F");
      theTreeNew->Branch("pyTkMet", &pyTkMet, "pyTkMet/F");
      theTreeNew->Branch("pzTkMet", &pzTkMet, "pzTkMet/F");
      theTreeNew->Branch("pxLeadJet", pxLeadJet, "pxLeadJet[3]/F");
      theTreeNew->Branch("pyLeadJet", pyLeadJet, "pyLeadJet[3]/F");
      theTreeNew->Branch("pzLeadJet", pzLeadJet, "pzLeadJet[3]/F");
      theTreeNew->Branch("pxSecondJet", pxSecondJet, "pxSecondJet[3]/F");
      theTreeNew->Branch("pySecondJet", pySecondJet, "pySecondJet[3]/F");
      theTreeNew->Branch("pzSecondJet", pzSecondJet, "pzSecondJet[3]/F");
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
      theTreeNew->Branch("dphilmet1", &dphilmet1, "dphilmet1/F");
      theTreeNew->Branch("dphilmet2", &dphilmet2, "dphilmet2/F");
      theTreeNew->Branch("deltaPhi_LL", &deltaPhi_LL, "deltaPhi_LL/F");
      theTreeNew->Branch("deltaPhi_LL_MET", &deltaPhi_LL_MET, "deltaPhi_LL_MET/F");
      theTreeNew->Branch("deltaPhi_LLJ1_MET", &deltaPhi_LLJ1_MET, "deltaPhi_LLJ1_MET/F");
      theTreeNew->Branch("dphilljet", &deltaPhi_LL_JET1, "dphilljet/F");
      theTreeNew->Branch("dphilljet2", &deltaPhi_LL_JET2, "dphilljet2/F");
      theTreeNew->Branch("deltaPhi_MET_JET1", &deltaPhi_MET_JET1, "deltaPhi_MET_JET1/F");
      theTreeNew->Branch("deltaPhi_MET_JET2", &deltaPhi_MET_JET2, "deltaPhi_MET_JET2/F");
      theTreeNew->Branch("deltaPhi_LL_JJ", &deltaPhi_LL_JJ, "deltaPhi_LL_JJ/F");
      theTreeNew->Branch("ptll", &dileptonPt, "ptll/F");
      theTreeNew->Branch("jetpt1", &jetpt1, "jetpt1/F");
      theTreeNew->Branch("jeteta1", &jeteta1, "jeteta1/F");
      theTreeNew->Branch("jetphi1", &jetphi1, "jetphi1/F");
      theTreeNew->Branch("jetpt2", &jetpt2, "jetpt2/F");
      theTreeNew->Branch("jeteta2", &jeteta2, "jeteta2/F");
      theTreeNew->Branch("jetphi2", &jetphi2, "jetphi2/F");
      theTreeNew->Branch("nSoftMu", &nSoftMu, "nSoftMu/I");
      theTreeNew->Branch("nExtraLep", &numExtraLep, "nExtraLep/I");
      theTreeNew->Branch("sumJetsV4", "TLorentzVector", &sumJetsV4);
      theTreeNew->Branch("uncorrSumJetsV4", "TLorentzVector", &uncorrSumJetsV4);
      theTreeNew->Branch("pfmetV", "TVector3", &pfmetV);
      theTreeNew->Branch("pfmetX", &pfmetX, "pfmetX/F");
      theTreeNew->Branch("pfmetY", &pfmetY, "pfmetY/F");
      theTreeNew->Branch("ch1", &(ch[0]), "ch1/I");
      theTreeNew->Branch("ch2", &(ch[1]), "ch2/I");
      theTreeNew->Branch("iso1", &(iso[0]), "iso1/F");
      theTreeNew->Branch("iso2", &(iso[1]), "iso2/F");

      theTreeNew->Branch("lh1", &(lh[0]), "lh1/F");
      theTreeNew->Branch("lh2", &(lh[1]), "lh2/F");
      // likelihood of the electron (or worst likelihood  if 2 electrons in the elel)
      theTreeNew->Branch("ellh", &ellh, "ellh/F");

      theTreeNew->Branch("jetcat", &jetcat,  "jetcat/F");

      theTreeNew->Branch("consecevent", &consecevent, "consecevent/F");
      theTreeNew->Branch("baseW", &baseW,  "baseW/F");

      if(FRWeights) {
        theTreeNew->Branch("weightFP", &weightFP, "weightFP/F");
        theTreeNew->Branch("weightStatFP", &weightStatFP, "weightStatFP/F");
        theTreeNew->Branch("weightFF", &weightFF, "weightFF/F");
        theTreeNew->Branch("weightStatFF", &weightStatFF, "weightStatFF/F");
        theTreeNew->Branch("weightPP", &weightPP, "weightPP/F");
        theTreeNew->Branch("weightStatPP", &weightStatPP, "weightStatPP/F");
      }

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

      zveto = (fabs(eleInvMass-91.1876)>15) ? 1 : 0;
      bveto_ip = (bTagTrackCount<2.1) ? 1 : 0;
      bveto_mu = (nSoftMu==0) ? 1 : 0;
      bveto = (bveto_ip && bveto_mu) ? 1 : 0;

      pfmetX = pfmetV->Px();
      pfmetY = pfmetV->Py();

      if(lh[0]!=-999 && lh[1]!=-999) ellh = TMath::Min(lh[0],lh[1]); // ee
      if(lh[0]!=-999 && lh[1]==-999) ellh = lh[0]; // emu
      if(lh[0]==-999 && lh[1]!=-999) ellh = lh[1]; // mue
      if(lh[0]==-999 && lh[1]==-999) ellh = -999; // mumu

      if (finalLeptons) {
        TVector3 TV_L1( pxL1, pyL1, pzL1 );
        TVector3 TV_L2( pxL2, pyL2, pzL2 );
        TVector3 TV_L1p2 = TV_L1 + TV_L2;
        TVector3 TV_tkmet( pxTkMet, pyTkMet, pzTkMet );
        TVector3 TV_jet1( pxLeadJet[0], pyLeadJet[0], pzLeadJet[0] );
        TVector3 TV_jet2( pxSecondJet[0], pySecondJet[0], pzSecondJet[0] );
        deltaPhi_LL       = (180./3.14) * TV_L1.DeltaPhi(TV_L2);
        deltaPhi_LL_MET   = (180./3.14) * TV_tkmet.DeltaPhi(TV_L1p2);

        dphilmet1 = fabs(pfmetV->DeltaPhi(TV_L1));
        dphilmet2 = fabs(pfmetV->DeltaPhi(TV_L2));

        tkMET = TV_tkmet.Pt();
        pmet = GetProjectedMet(*pfmetV,TV_L1,TV_L2);
        pmet2 = GetProjectedMet(tkMET,TV_L1,TV_L2);

        mtw1 = calcMT(*pfmetV,TV_L1);
        mtw2 = calcMT(*pfmetV,TV_L2);

        if(TV_jet1.Pt()>15) {
          TVector3 TV_L12pJ1 = TV_L1p2 + TV_jet1;
          deltaPhi_LLJ1_MET = (180./TMath::Pi()) * fabs(TV_tkmet.DeltaPhi(TV_L12pJ1));   
          deltaPhi_MET_JET1 = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_tkmet));
          deltaPhi_LL_JET1  = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_L1p2));
          jetpt1 = TV_jet1.Pt();
          jeteta1 = TV_jet1.Eta();
          jetphi1 = TV_jet1.Phi();
        } else {
          deltaPhi_LLJ1_MET = -1.;
          deltaPhi_MET_JET1 = -1.;
          deltaPhi_LL_JET1 = -1.;
          jetpt1 = -999.;
          jeteta1 = -999.;
          jetphi1 = -999.;
        }

        if(TV_jet2.Pt()>15) {
          deltaPhi_MET_JET2 = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_tkmet));
          deltaPhi_LL_JET2  = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_L1p2));
          jetpt2 = TV_jet2.Pt();
          jeteta2 = TV_jet2.Eta();
          jetphi2 = TV_jet2.Phi();
        } else {
          deltaPhi_MET_JET2 = -1.;
          deltaPhi_LL_JET2 = -1.;
          jetpt2 = -999.;
          jeteta2 = -999.;
          jetphi2 = -999.;
        }

        if(TV_jet1.Pt()>15 && TV_jet2.Pt()>15) {
          TVector3 TV_J1p2  = TV_jet1 + TV_jet2;
          deltaPhi_LL_JJ    = (180./3.14) * TV_L1p2.DeltaPhi(TV_J1p2);
        } else {
          deltaPhi_LL_JJ = -1.;
        }

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
        jetpt1 = -9999.;
        jetpt2 = -9999.;
        dileptonPt = -9999.;
      } 
      consecevent = (float)j;
      if(finalLeptons) {
        if(processId>0) { // MC
          treeNew->Fill();
          if(i_WWSel || i_WWSel1j) treeNewSkim->Fill();
        } else { // data: apply the trigger 
          if(hlt) {
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

float GetProjectedMet(TVector3 met, TVector3 p1, TVector3 p2) {

  float projMET = 0.0;
  float deltaPhi1 = fabs(p1.DeltaPhi(met));
  float deltaPhi2 = fabs(p2.DeltaPhi(met));
  float deltaphi = TMath::Min(deltaPhi1,deltaPhi2);
  if(deltaphi<TMath::Pi()/2.) projMET = met.Mag() * sin(deltaphi);
  else projMET = met.Mag();

  return projMET;
}

float calcMT(TVector3 met, TVector3 lepton) {
  return sqrt( 2.*(lepton.Pt())*(met.Pt())*( 1 - cos(met.DeltaPhi(lepton))) );
}
