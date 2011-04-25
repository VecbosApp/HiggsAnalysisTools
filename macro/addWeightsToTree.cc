#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <iostream>

using namespace std;

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
    Float_t         KFactor;
    Bool_t          promptDecay;
    Float_t         maxPtLh;
    Float_t         minPtLh;
    Int_t           njets;
    Int_t           nuncorrjets;
    Float_t         dxyEVT;
    Float_t         dszEVT;
    Float_t         bTagTrackCount;
    Float_t         bTagImpPar;
    Float_t         bTagSecVertex;
    Float_t         pt[2];
    Float_t         eta[2];
    Float_t         deta[2];
    Float_t         dphi[2];
    Float_t         hoe[2];
    Float_t         see[2];
    Int_t           matched[2];

    treeOrig->SetBranchAddress("run", &run);
    treeOrig->SetBranchAddress("lumi", &lumi);
    treeOrig->SetBranchAddress("event", &event);
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
    treeOrig->SetBranchAddress("finalLeptons", &finalLeptons);
    treeOrig->SetBranchAddress("jetVeto", &jetVeto);
    treeOrig->SetBranchAddress("uncorrJetVeto", &uncorrJetVeto);
    treeOrig->SetBranchAddress("preDeltaPhi", &preDeltaPhi);
    treeOrig->SetBranchAddress("finalSelection", &finalSelection);
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
    treeOrig->SetBranchAddress("pt", pt);
    treeOrig->SetBranchAddress("eta", eta);
    treeOrig->SetBranchAddress("deta", deta);
    treeOrig->SetBranchAddress("dphi", dphi);
    treeOrig->SetBranchAddress("hoe", hoe);
    treeOrig->SetBranchAddress("see", see);
    treeOrig->SetBranchAddress("matched", matched);

    // 
    Float_t pt_1,   pt_2;
    Float_t eta_1,  eta_2;
    Float_t deta_1, deta_2;
    Float_t dphi_1, dphi_2;
    Float_t hoe_1,  hoe_2;
    Float_t see_1,  see_2;
    Int_t   matched_1,  matched_2;
    Float_t expCosDphi;

    // convert the booleans into integers (to insert in RooDataset)
    Int_t         i_finalLeptons;
    Int_t         i_jetVeto;
    Int_t         i_uncorrJetVeto;
    Int_t         i_preDeltaPhi;
    Int_t         i_finalSelection;
    Int_t         i_promptDecay;

    // the selected final state: ee=0, mm=1, em=2
    treeNew->Branch("finalstate", &finalstate, "finalstate/I");

    // copy branches
    treeNew->Branch("run", &run, "run/I");
    treeNew->Branch("lumi", &lumi, "lumi/I");
    treeNew->Branch("event", &event, "event/I");
    treeNew->Branch("met", &met, "met/F");  // default MET is tcMET for WW
    treeNew->Branch("pfMet", &pfMet, "pfMet/F");
    treeNew->Branch("caloMet", &caloMet, "caloMet/F");
    treeNew->Branch("projMet", &projMet, "projMet/F");
    treeNew->Branch("deltaPhi", &deltaPhi, "deltaPhi/F");
    treeNew->Branch("deltaR", &deltaR, "deltaR/F");
    treeNew->Branch("transvMass", &transvMass, "transvMass/F");
    treeNew->Branch("eleInvMass", &eleInvMass, "eleInvMass/F");
    treeNew->Branch("maxPtEle", &maxPtEle, "maxPtEle/F");
    treeNew->Branch("minPtEle", &minPtEle, "minPtEle/F");
    treeNew->Branch("detaLeptons", &detaLeptons, "detaLeptons/F");
    treeNew->Branch("finalLeptons", &i_finalLeptons, "finalLeptons/I");
    treeNew->Branch("jetVeto", &i_jetVeto, "jetVeto/I");
    treeNew->Branch("uncorrJetVeto", &i_uncorrJetVeto, "uncorrJetVeto/I");
    treeNew->Branch("preDeltaPhi", &i_preDeltaPhi, "preDeltaPhi/I");
    treeNew->Branch("finalSelection", &i_finalSelection, "finalSelection/I");
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
