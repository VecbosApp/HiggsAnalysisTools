#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

int fullFormat = 1;

void setReducedFormat() { fullFormat = 0; }

void addWeights(const char* filename, float weight, int finalstate) {

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
    Float_t         KFactor;
    Bool_t          promptDecay;
    Float_t         maxPtLh;
    Float_t         minPtLh;
    Int_t           njets;
    Int_t           nuncorrjets;
    Int_t           nVtx;
    Float_t         dxyEVT;
    Float_t         dszEVT;
    Float_t         bTagTrackCount;
    Float_t         bTagImpPar;
    Float_t         bTagSecVertex;
    Float_t         leadingJetBTagTrackCount;
    Float_t         pt[2];
    Float_t         eta[2];
    Float_t         deta[2];
    Float_t         dphi[2];
    Float_t         hoe[2];
    Float_t         see[2];
    Int_t           matched[2];
    Float_t         pxTkMet;
    Float_t         pyTkMet;
    Float_t         pzTkMet;
    Float_t         pxLeadJet;
    Float_t         pyLeadJet;
    Float_t         pzLeadJet;
    Float_t         pxL1;
    Float_t         pyL1;
    Float_t         pzL1;
    Float_t         pxL2;
    Float_t         pyL2;
    Float_t         pzL2;
    Int_t           nSoftMu;

    treeOrig->SetBranchAddress("run", &run);
    treeOrig->SetBranchAddress("lumi", &lumi);
    treeOrig->SetBranchAddress("event", &event);
    treeOrig->SetBranchAddress("puweight", &puweight);
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
    treeOrig->SetBranchAddress("maxPtLh", &maxPtLh);
    treeOrig->SetBranchAddress("minPtLh", &minPtLh);
    treeOrig->SetBranchAddress("njets", &njets);
    treeOrig->SetBranchAddress("nuncorrjets", &nuncorrjets);
    treeOrig->SetBranchAddress("dxyEVT", &dxyEVT);
    treeOrig->SetBranchAddress("dszEVT", &dszEVT);
    treeOrig->SetBranchAddress("bTagTrackCount", &bTagTrackCount);
    treeOrig->SetBranchAddress("bTagImpPar", &bTagImpPar);
    treeOrig->SetBranchAddress("bTagSecVertex", &bTagSecVertex);
    treeOrig->SetBranchAddress("leadingJetBTagTrackCount", &leadingJetBTagTrackCount);
    treeOrig->SetBranchAddress("pt", pt);
    treeOrig->SetBranchAddress("eta", eta);
    treeOrig->SetBranchAddress("deta", deta);
    treeOrig->SetBranchAddress("dphi", dphi);
    treeOrig->SetBranchAddress("hoe", hoe);
    treeOrig->SetBranchAddress("see", see);
    treeOrig->SetBranchAddress("matched", matched);
    treeOrig->SetBranchAddress("pxTkMet", &pxTkMet);
    treeOrig->SetBranchAddress("pyTkMet", &pyTkMet);
    treeOrig->SetBranchAddress("pzTkMet", &pzTkMet);
    treeOrig->SetBranchAddress("pxLeadJet", &pxLeadJet);
    treeOrig->SetBranchAddress("pyLeadJet", &pyLeadJet);
    treeOrig->SetBranchAddress("pzLeadJet", &pzLeadJet);
    treeOrig->SetBranchAddress("pxL1", &pxL1);
    treeOrig->SetBranchAddress("pyL1", &pyL1);
    treeOrig->SetBranchAddress("pzL1", &pzL1);
    treeOrig->SetBranchAddress("pxL2", &pxL2);
    treeOrig->SetBranchAddress("pyL2", &pyL2);
    treeOrig->SetBranchAddress("pzL2", &pzL2);
    treeOrig->SetBranchAddress("nSoftMu", &nSoftMu);

    // 
    Float_t pt_1,   pt_2;
    Float_t eta_1,  eta_2;
    Float_t deta_1, deta_2;
    Float_t dphi_1, dphi_2;
    Float_t hoe_1,  hoe_2;
    Float_t see_1,  see_2;
    Int_t   matched_1,  matched_2;
    Float_t expCosDphi;
    Float_t gammaStMRSt;

    // convert the booleans into integers (to insert in RooDataset)
    Int_t         i_finalLeptons;
    Int_t         i_jetVeto;
    Int_t         i_uncorrJetVeto;
    Int_t         i_preDeltaPhi;
    Int_t         i_finalSelection;
    Int_t         i_promptDecay;
    Int_t         i_WWSel;
    Int_t         i_WWSel1j;
    float deltaPhi_LL;    
    float deltaPhi_LL_MET;
    float deltaPhi_LL_JET;
    float deltaPhi_MET_JET;
    float leadingJetPt;
    float L1eta, L1phi;
    float L2eta, L2phi;
    float dileptonPt;

    // the selected final state: ee=0, mm=1, em=2
    treeNew->Branch("finalstate", &finalstate, "finalstate/I");

    // copy branches
    treeNew->Branch("run", &run, "run/I");
    treeNew->Branch("lumi", &lumi, "lumi/I");
    treeNew->Branch("event", &event, "event/I");
    treeNew->Branch("puweight", &puweight, "puweight/F");
    treeNew->Branch("met", &met, "met/F");  // default MET is tcMET for WW
    treeNew->Branch("pfMet", &pfMet, "pfMet/F");
    treeNew->Branch("caloMet", &caloMet, "caloMet/F");
    treeNew->Branch("projMet", &projMet, "projMet/F");
    treeNew->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
    treeNew->Branch("deltaR", &deltaR, "deltaR/F");
    treeNew->Branch("gammaStMRSt", &gammaStMRSt, "gammaStMRSt/F");
    treeNew->Branch("eleInvMass", &eleInvMass, "eleInvMass/F");
    treeNew->Branch("transvMass", &transvMass, "transvMass/F");
    treeNew->Branch("maxPtEle", &maxPtEle, "maxPtEle/F");
    treeNew->Branch("minPtEle", &minPtEle, "minPtEle/F");
    treeNew->Branch("detaLeptons", &detaLeptons, "detaLeptons/F");
    treeNew->Branch("nVtx", &nVtx, "nVtx/I");
    treeNew->Branch("finalLeptons", &i_finalLeptons, "finalLeptons/I");
    treeNew->Branch("jetVeto", &i_jetVeto, "jetVeto/I");
    treeNew->Branch("uncorrJetVeto", &i_uncorrJetVeto, "uncorrJetVeto/I");
    treeNew->Branch("preDeltaPhi", &i_preDeltaPhi, "preDeltaPhi/I");
    treeNew->Branch("finalSelection", &i_finalSelection, "finalSelection/I");
    treeNew->Branch("WWSel", &i_WWSel, "WWSel/I");
    treeNew->Branch("WWSel1j", &i_WWSel1j, "WWSel1j/I");
    treeNew->Branch("KFactor", &KFactor, "KFactor/F");
    treeNew->Branch("promptDecay", &i_promptDecay, "promptDecay/I");
    treeNew->Branch("maxPtLh", &maxPtLh, "maxPtLh/F");
    treeNew->Branch("minPtLh", &minPtLh, "minPtLh/F");
    treeNew->Branch("njets", &njets, "njets/I");
    treeNew->Branch("nuncorrjets", &nuncorrjets, "nuncorrjets/I");
    treeNew->Branch("dxyEVT", &dxyEVT, "dxyEVT/F");
    treeNew->Branch("dszEVT", &dszEVT, "dszEVT/F");
    treeNew->Branch("bTagTrackCount", &bTagTrackCount, "bTagTrackCount/F");
    treeNew->Branch("bTagImpPar", &bTagImpPar, "bTagImpPar/F");
    treeNew->Branch("bTagSecVertex", &bTagSecVertex, "bTagSecVertex/F");
    treeNew->Branch("leadingJetBTagTrackCount", &leadingJetBTagTrackCount, "leadingJetBTagTrackCount/F");
    treeNew->Branch("pt1", &pt_1, "pt1/F");
    treeNew->Branch("eta1", &eta_1, "eta1/F");
    treeNew->Branch("deta1", &deta_1, "deta1/F");
    treeNew->Branch("dphi1", &dphi_1, "dphi1/F");
    treeNew->Branch("hoe1", &hoe_1, "hoe1/F");
    treeNew->Branch("see1", &see_1, "see1/F");
    treeNew->Branch("matched1", &matched_1, "matched1/I");
    treeNew->Branch("pt2", &pt_2, "pt2/F");
    treeNew->Branch("eta2", &eta_2, "eta2/F");
    treeNew->Branch("deta2", &deta_2, "deta2/F");
    treeNew->Branch("dphi2", &dphi_2, "dphi2/F");
    treeNew->Branch("hoe2", &hoe_2, "hoe2/F");
    treeNew->Branch("see2", &see_2, "see2/F");
    treeNew->Branch("matched2", &matched_2, "matched2/I");
    treeNew->Branch("expCosDphi", &expCosDphi, "expCosDphi/F");
    treeNew->Branch("step", step, "step[25]/O");
    if(fullFormat) {
      treeNew->Branch("pxTkMet", &pxTkMet, "pxTkMet/F");
      treeNew->Branch("pyTkMet", &pyTkMet, "pyTkMet/F");
      treeNew->Branch("pzTkMet", &pzTkMet, "pzTkMet/F");
      treeNew->Branch("pxLeadJet", &pxLeadJet, "pxLeadJet/F");
      treeNew->Branch("pyLeadJet", &pyLeadJet, "pyLeadJet/F");
      treeNew->Branch("pzLeadJet", &pzLeadJet, "pzLeadJet/F");
      treeNew->Branch("pxL1", &pxL1, "pxL1/F");
      treeNew->Branch("pyL1", &pyL1, "pyL1/F");
      treeNew->Branch("pzL1", &pzL1, "pzL1/F");
      treeNew->Branch("pxL2", &pxL2, "pxL2/F");
      treeNew->Branch("pyL2", &pyL2, "pyL2/F");
      treeNew->Branch("pzL2", &pzL2, "pzL2/F");
      treeNew->Branch("deltaPhi_LL", &deltaPhi_LL, "deltaPhi_LL/F");
      treeNew->Branch("deltaPhi_LL_MET", &deltaPhi_LL_MET, "deltaPhi_LL_MET/F");
      treeNew->Branch("deltaPhi_LL_JET", &deltaPhi_LL_JET, "deltaPhi_LL_JET/F");
      treeNew->Branch("deltaPhi_MET_JET", &deltaPhi_MET_JET, "deltaPhi_MET_JET/F");
      treeNew->Branch("dileptonPt", &dileptonPt, "dileptonPt/F");
      treeNew->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/F");
      treeNew->Branch("L1eta", &L1eta, "L1eta/F");
      treeNew->Branch("L1phi", &L1phi, "L1phi/F");
      treeNew->Branch("L2eta", &L2eta, "L2eta/F");
      treeNew->Branch("L2phi", &L2phi, "L2phi/F");
      treeNew->Branch("nSoftMu", &nSoftMu, "nSoftMu/I");
    }

    float jetcat = 0;
    treeNew->Branch("jetcat", &jetcat,  "jetcat/F");
    float consecevent = -1;
    treeNew->Branch("consecevent", &consecevent, "consecevent/F");
    treeNew->Branch("weight", &weight,  "weight/F");

    int j =0;

    for(int i=0; i<nentriesOrig; i++) {
      if (i%10000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl;
      treeOrig->GetEntry(i);

      if (njets==0) jetcat = 1;
      else if(njets==1) jetcat = -1;
      else jetcat = -2;

      pt_1      = pt[0];
      eta_1     = eta[0];
      deta_1    = deta[0];
      dphi_1    = dphi[0];
      hoe_1     = hoe[0];
      see_1     = see[0];
      matched_1 = matched[0];
      pt_2      = pt[1];
      eta_2     = eta[1];
      deta_2    = deta[1];
      dphi_2    = dphi[1];
      hoe_2     = hoe[1];
      see_2     = see[1];
      matched_2 = matched[1];
      expCosDphi = exp(cos(TMath::Pi()*deltaPhi/180.));
      gammaStMRSt = 2*transvMass;

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

      if (finalLeptons && fullFormat) {
        TVector3 TV_L1( pxL1, pyL1, pzL1 );
        TVector3 TV_L2( pxL2, pyL2, pzL2 );
        TVector3 TV_L1p2 = TV_L1 + TV_L2;
        TVector3 TV_met( pxTkMet, pyTkMet, pzTkMet );
        TVector3 TV_jet( pxLeadJet, pyLeadJet, pzLeadJet );
        deltaPhi_LL     = (180./3.14) * TV_L1.DeltaPhi(TV_L2);
        deltaPhi_LL_MET = (180./3.14) * TV_met.DeltaPhi(TV_L1p2);
        deltaPhi_LL_JET = (180./3.14) * TV_jet.DeltaPhi(TV_L1p2);
        deltaPhi_MET_JET = (180./3.14) * TV_jet.DeltaPhi(TV_met);
        leadingJetPt = sqrt(pxLeadJet*pxLeadJet + pyLeadJet*pyLeadJet);
        dileptonPt = TV_L1p2.Pt();
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
        deltaPhi_LL_JET = -9999.;
        deltaPhi_MET_JET = -9999.;
        leadingJetPt = -9999.;
        dileptonPt = -9999.;
      } 

      consecevent = (float)j;
      treeNew->Fill();
      j++;
      //    }
    }
  
    fileNew->cd();
    treeNew->Write();
    fileNew->Close();

    fileOrig->cd();
    fileOrig->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }
}
