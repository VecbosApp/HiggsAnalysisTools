#define addWeightsToTreeLooseLoose_cxx
#include "addWeightsToTreeLooseLoose.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

void addWeightsToTreeLooseLoose::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L addWeightsToTreeLooseLoose.C
  //      Root > addWeightsToTreeLooseLoose t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  float WWSel, WWSel1j;

  TFile *fileNew = TFile::Open("results_data/datasets_trees_looseloose/looseloose.root","recreate");
  TTree *treeNew = new TTree("latino","tree with only selected events");

  TFile *fileNewSkim = TFile::Open("results_data/datasets_trees_looseloose_skim/looseloose.root","recreate");
  TTree *treeNewSkim = new TTree("latino","tree with only selected events");

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);
  trees.push_back(treeNewSkim);

  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];

    theTreeNew->Branch("baseW", &baseW, "baseW/F");
    theTreeNew->Branch("bdt1", &bdt1, "bdt1/F");
    theTreeNew->Branch("bdt2", &bdt2, "bdt2/F");
    theTreeNew->Branch("ch1", &ch1, "ch1/F");
    theTreeNew->Branch("ch2", &ch2, "ch2/F");
    theTreeNew->Branch("channel", &channel, "channel/F");
    theTreeNew->Branch("chmet", &chmet, "chmet/F");
    theTreeNew->Branch("chmetphi", &chmetphi, "chmetphi/F");
    theTreeNew->Branch("dataset", &dataset, "dataset/F");
    theTreeNew->Branch("detajj", &detajj, "detajj/F");
    theTreeNew->Branch("dphill", &dphill, "dphill/F");
    theTreeNew->Branch("dphilljet", &dphilljet, "dphilljet/F");
    theTreeNew->Branch("dphilljetjet", &dphilljetjet, "dphilljetjet/F");
    theTreeNew->Branch("dphillmet", &dphillmet, "dphillmet/F");
    theTreeNew->Branch("dphilmet", &dphilmet, "dphilmet/F");
    theTreeNew->Branch("dphilmet1", &dphilmet1, "dphilmet1/F");
    theTreeNew->Branch("dphilmet2", &dphilmet2, "dphilmet2/F");
    theTreeNew->Branch("drll", &drll, "drll/F");
    theTreeNew->Branch("effAW", &effAW, "effAW/F");
    theTreeNew->Branch("effBW", &effBW, "effBW/F");
    theTreeNew->Branch("effW", &effW, "effW/F");
    theTreeNew->Branch("eta1", &eta1, "eta1/F");
    theTreeNew->Branch("eta2", &eta2, "eta2/F");
    theTreeNew->Branch("fakeAW", &fakeAW, "fakeAW/F");
    theTreeNew->Branch("fakeBW", &fakeBW, "fakeBW/F");
    theTreeNew->Branch("fakeW", &fakeW, "fakeW/F");
    theTreeNew->Branch("fermiW", &fermiW, "fermiW/F");
    theTreeNew->Branch("fourW", &fourW, "fourW/F");
    theTreeNew->Branch("gammaMRStar", &gammaMRStar, "gammaMRStar/F");
    theTreeNew->Branch("hardbdisc", &hardbdisc, "hardbdisc/F");
    theTreeNew->Branch("hypo", &hypo, "hypo/F");
    theTreeNew->Branch("imet", &imet, "imet/F");
    theTreeNew->Branch("iso1", &iso1, "iso1/F");
    theTreeNew->Branch("iso2", &iso2, "iso2/F");
    theTreeNew->Branch("jeteta1", &jeteta1, "jeteta1/F");
    theTreeNew->Branch("jeteta2", &jeteta2, "jeteta2/F");
    theTreeNew->Branch("jetphi1", &jetphi1, "jetphi1/F");
    theTreeNew->Branch("jetphi2", &jetphi2, "jetphi2/F");
    theTreeNew->Branch("jetpt1", &jetpt1, "jetpt1/F");
    theTreeNew->Branch("jetpt2", &jetpt2, "jetpt2/F");
    theTreeNew->Branch("jettche1", &jettche1, "jettche1/F");
    theTreeNew->Branch("jettche2", &jettche2, "jettche2/F");
    theTreeNew->Branch("jettchp1", &jettchp1, "jettchp1/F");
    theTreeNew->Branch("jettchp2", &jettchp2, "jettchp2/F");
    theTreeNew->Branch("kfW", &kfW, "kfW/F");
    theTreeNew->Branch("lh1", &lh1, "lh1/F");
    theTreeNew->Branch("lh2", &lh2, "lh2/F");
    theTreeNew->Branch("mjj", &mjj, "mjj/F");
    theTreeNew->Branch("mll", &mll, "mll/F");
    theTreeNew->Branch("mpmet", &mpmet, "mpmet/F");
    theTreeNew->Branch("mth", &mth, "mth/F");
    theTreeNew->Branch("mtw1", &mtw1, "mtw1/F");
    theTreeNew->Branch("mtw2", &mtw2, "mtw2/F");
    theTreeNew->Branch("nbjet", &nbjet, "nbjet/F");
    theTreeNew->Branch("nbrem1", &nbrem1, "nbrem1/F");
    theTreeNew->Branch("nbrem2", &nbrem2, "nbrem2/F");
    theTreeNew->Branch("nextra", &nextra, "nextra/F");
    theTreeNew->Branch("njet", &njet, "njet/F");
    theTreeNew->Branch("njetid", &njetid, "njetid/F");
    theTreeNew->Branch("njetvbf", &njetvbf, "njetvbf/F");
    theTreeNew->Branch("pchmet", &pchmet, "pchmet/F");
    theTreeNew->Branch("peaking", &peaking, "peaking/F");
    theTreeNew->Branch("pfmet", &pfmet, "pfmet/F");
    theTreeNew->Branch("pfmetphi", &pfmetphi, "pfmetphi/F");
    theTreeNew->Branch("phi1", &phi1, "phi1/F");
    theTreeNew->Branch("phi2", &phi2, "phi2/F");
    theTreeNew->Branch("ppfmet", &ppfmet, "ppfmet/F");
    theTreeNew->Branch("predmet", &predmet, "predmet/F");
    theTreeNew->Branch("pt1", &pt1, "pt1/F");
    theTreeNew->Branch("pt2", &pt2, "pt2/F");
    theTreeNew->Branch("ptcmet", &ptcmet, "ptcmet/F");
    theTreeNew->Branch("ptll", &ptll, "ptll/F");
    theTreeNew->Branch("puAW", &puAW, "puAW/F");
    theTreeNew->Branch("puBW", &puBW, "puBW/F");
    theTreeNew->Branch("puW", &puW, "puW/F");
    theTreeNew->Branch("redmet", &redmet, "redmet/F");
    theTreeNew->Branch("sceta1", &sceta1, "sceta1/F");
    theTreeNew->Branch("sceta2", &sceta2, "sceta2/F");
    theTreeNew->Branch("softbdisc", &softbdisc, "softbdisc/F");
    theTreeNew->Branch("tcmet", &tcmet, "tcmet/F");
    theTreeNew->Branch("tcmetphi", &tcmetphi, "tcmetphi/F");
    theTreeNew->Branch("tightmu", &tightmu, "tightmu/F");
    theTreeNew->Branch("triggAW", &triggAW, "triggAW/F");
    theTreeNew->Branch("triggBW", &triggBW, "triggBW/F");
    theTreeNew->Branch("triggW", &triggW, "triggW/F");
    theTreeNew->Branch("trigger", &trigger, "trigger/F");
    theTreeNew->Branch("worstJetLepPt", &worstJetLepPt, "worstJetLepPt/F");
    theTreeNew->Branch("yll", &yll, "yll/F");
    theTreeNew->Branch("nvtx", &nvtx, "nvtx/F");
    theTreeNew->Branch("bveto", &bveto, "bveto/F");
    theTreeNew->Branch("bveto_ip", &bveto_ip, "bveto_ip/F");
    theTreeNew->Branch("bveto_mu", &bveto_mu, "bveto_mu/F");
    theTreeNew->Branch("bveto_munj", &bveto_munj, "bveto_munj/F");
    theTreeNew->Branch("bveto_nj", &bveto_nj, "bveto_nj/F");
    theTreeNew->Branch("dphiveto", &dphiveto, "dphiveto/F");
    theTreeNew->Branch("passBDT1", &passBDT1, "passBDT1/F");
    theTreeNew->Branch("passBDT2", &passBDT2, "passBDT2/F");
    theTreeNew->Branch("passCB1", &passCB1, "passCB1/F");
    theTreeNew->Branch("passCB2", &passCB2, "passCB2/F");
    theTreeNew->Branch("passCBOld1", &passCBOld1, "passCBOld1/F");
    theTreeNew->Branch("passCBOld2", &passCBOld2, "passCBOld2/F");
    theTreeNew->Branch("passLH1", &passLH1, "passLH1/F");
    theTreeNew->Branch("passLH2", &passLH2, "passLH2/F");
    theTreeNew->Branch("sameflav", &sameflav, "sameflav/F");
    theTreeNew->Branch("zveto", &zveto, "zveto/F");
    theTreeNew->Branch("run", &run, "run/F");
    theTreeNew->Branch("lumi", &lumi, "lumi/F");
    theTreeNew->Branch("event", &event, "event/F");
    theTreeNew->Branch("fakeW", &fakeW, "fakeW/F");
    theTreeNew->Branch("WWSel", &WWSel, "WWSel/F");
    theTreeNew->Branch("WWSel1j", &WWSel1j, "WWSel1j/F");
  }

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000 == 0) std::cout << ">>> Weighting event # " << jentry << " / " << nentries << " entries" << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    WWSel = trigger && pfmet > 20 && mll > (12 + 8*sameflav) && zveto && mpmet > (20+(17+nvtx/2.)*sameflav)&& njet==0 && (dphiveto || ! sameflav) && bveto_mu && nextra == 0 && bveto_ip && (pt2 > 15||!sameflav) && ptll > 45;
    WWSel1j = trigger && pfmet > 20 && mll > (12 + 8*sameflav) && zveto && mpmet > (20+(17+nvtx/2.)*sameflav)&& njet==1 && (dphiveto || ! sameflav) && bveto_mu && nextra == 0 && bveto_ip && nbjet==0 && (pt2 > 15||!sameflav) && ptll > 45 ;
    gammaMRStar = 2*gammaMRStar;
    treeNew->Fill();
    if(WWSel || WWSel1j) treeNewSkim->Fill();

  }

  fileNew->cd();
  treeNew->Write();
  fileNew->Close();

  fileNewSkim->cd();
  treeNewSkim->Write();
  fileNewSkim->Close();

}
