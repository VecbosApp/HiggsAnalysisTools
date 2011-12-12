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

   treeNew->Branch("baseW", &baseW, "baseW/F");
   treeNew->Branch("bdt1", &bdt1, "bdt1/F");
   treeNew->Branch("bdt2", &bdt2, "bdt2/F");
   treeNew->Branch("ch1", &ch1, "ch1/F");
   treeNew->Branch("ch2", &ch2, "ch2/F");
   treeNew->Branch("channel", &channel, "channel/F");
   treeNew->Branch("chmet", &chmet, "chmet/F");
   treeNew->Branch("chmetphi", &chmetphi, "chmetphi/F");
   treeNew->Branch("dataset", &dataset, "dataset/F");
   treeNew->Branch("detajj", &detajj, "detajj/F");
   treeNew->Branch("dphill", &dphill, "dphill/F");
   treeNew->Branch("dphilljet", &dphilljet, "dphilljet/F");
   treeNew->Branch("dphilljetjet", &dphilljetjet, "dphilljetjet/F");
   treeNew->Branch("dphillmet", &dphillmet, "dphillmet/F");
   treeNew->Branch("dphilmet", &dphilmet, "dphilmet/F");
   treeNew->Branch("dphilmet1", &dphilmet1, "dphilmet1/F");
   treeNew->Branch("dphilmet2", &dphilmet2, "dphilmet2/F");
   treeNew->Branch("drll", &drll, "drll/F");
   treeNew->Branch("effAW", &effAW, "effAW/F");
   treeNew->Branch("effBW", &effBW, "effBW/F");
   treeNew->Branch("effW", &effW, "effW/F");
   treeNew->Branch("eta1", &eta1, "eta1/F");
   treeNew->Branch("eta2", &eta2, "eta2/F");
   treeNew->Branch("fakeAW", &fakeAW, "fakeAW/F");
   treeNew->Branch("fakeBW", &fakeBW, "fakeBW/F");
   treeNew->Branch("fakeW", &fakeW, "fakeW/F");
   treeNew->Branch("fermiW", &fermiW, "fermiW/F");
   treeNew->Branch("fourW", &fourW, "fourW/F");
   treeNew->Branch("gammaMRStar", &gammaMRStar, "gammaMRStar/F");
   treeNew->Branch("hardbdisc", &hardbdisc, "hardbdisc/F");
   treeNew->Branch("hypo", &hypo, "hypo/F");
   treeNew->Branch("imet", &imet, "imet/F");
   treeNew->Branch("iso1", &iso1, "iso1/F");
   treeNew->Branch("iso2", &iso2, "iso2/F");
   treeNew->Branch("jeteta1", &jeteta1, "jeteta1/F");
   treeNew->Branch("jeteta2", &jeteta2, "jeteta2/F");
   treeNew->Branch("jetphi1", &jetphi1, "jetphi1/F");
   treeNew->Branch("jetphi2", &jetphi2, "jetphi2/F");
   treeNew->Branch("jetpt1", &jetpt1, "jetpt1/F");
   treeNew->Branch("jetpt2", &jetpt2, "jetpt2/F");
   treeNew->Branch("jettche1", &jettche1, "jettche1/F");
   treeNew->Branch("jettche2", &jettche2, "jettche2/F");
   treeNew->Branch("jettchp1", &jettchp1, "jettchp1/F");
   treeNew->Branch("jettchp2", &jettchp2, "jettchp2/F");
   treeNew->Branch("kfW", &kfW, "kfW/F");
   treeNew->Branch("lh1", &lh1, "lh1/F");
   treeNew->Branch("lh2", &lh2, "lh2/F");
   treeNew->Branch("mjj", &mjj, "mjj/F");
   treeNew->Branch("mll", &mll, "mll/F");
   treeNew->Branch("mpmet", &mpmet, "mpmet/F");
   treeNew->Branch("mth", &mth, "mth/F");
   treeNew->Branch("mtw1", &mtw1, "mtw1/F");
   treeNew->Branch("mtw2", &mtw2, "mtw2/F");
   treeNew->Branch("nbjet", &nbjet, "nbjet/F");
   treeNew->Branch("nbrem1", &nbrem1, "nbrem1/F");
   treeNew->Branch("nbrem2", &nbrem2, "nbrem2/F");
   treeNew->Branch("nextra", &nextra, "nextra/F");
   treeNew->Branch("njet", &njet, "njet/F");
   treeNew->Branch("njetid", &njetid, "njetid/F");
   treeNew->Branch("njetvbf", &njetvbf, "njetvbf/F");
   treeNew->Branch("pchmet", &pchmet, "pchmet/F");
   treeNew->Branch("peaking", &peaking, "peaking/F");
   treeNew->Branch("pfmet", &pfmet, "pfmet/F");
   treeNew->Branch("pfmetphi", &pfmetphi, "pfmetphi/F");
   treeNew->Branch("phi1", &phi1, "phi1/F");
   treeNew->Branch("phi2", &phi2, "phi2/F");
   treeNew->Branch("ppfmet", &ppfmet, "ppfmet/F");
   treeNew->Branch("predmet", &predmet, "predmet/F");
   treeNew->Branch("pt1", &pt1, "pt1/F");
   treeNew->Branch("pt2", &pt2, "pt2/F");
   treeNew->Branch("ptcmet", &ptcmet, "ptcmet/F");
   treeNew->Branch("ptll", &ptll, "ptll/F");
   treeNew->Branch("puAW", &puAW, "puAW/F");
   treeNew->Branch("puBW", &puBW, "puBW/F");
   treeNew->Branch("puW", &puW, "puW/F");
   treeNew->Branch("redmet", &redmet, "redmet/F");
   treeNew->Branch("sceta1", &sceta1, "sceta1/F");
   treeNew->Branch("sceta2", &sceta2, "sceta2/F");
   treeNew->Branch("softbdisc", &softbdisc, "softbdisc/F");
   treeNew->Branch("tcmet", &tcmet, "tcmet/F");
   treeNew->Branch("tcmetphi", &tcmetphi, "tcmetphi/F");
   treeNew->Branch("tightmu", &tightmu, "tightmu/F");
   treeNew->Branch("triggAW", &triggAW, "triggAW/F");
   treeNew->Branch("triggBW", &triggBW, "triggBW/F");
   treeNew->Branch("triggW", &triggW, "triggW/F");
   treeNew->Branch("trigger", &trigger, "trigger/F");
   treeNew->Branch("worstJetLepPt", &worstJetLepPt, "worstJetLepPt/F");
   treeNew->Branch("yll", &yll, "yll/F");
   treeNew->Branch("nvtx", &nvtx, "nvtx/F");
   treeNew->Branch("bveto", &bveto, "bveto/F");
   treeNew->Branch("bveto_ip", &bveto_ip, "bveto_ip/F");
   treeNew->Branch("bveto_mu", &bveto_mu, "bveto_mu/F");
   treeNew->Branch("bveto_munj", &bveto_munj, "bveto_munj/F");
   treeNew->Branch("bveto_nj", &bveto_nj, "bveto_nj/F");
   treeNew->Branch("dphiveto", &dphiveto, "dphiveto/F");
   treeNew->Branch("passBDT1", &passBDT1, "passBDT1/F");
   treeNew->Branch("passBDT2", &passBDT2, "passBDT2/F");
   treeNew->Branch("passCB1", &passCB1, "passCB1/F");
   treeNew->Branch("passCB2", &passCB2, "passCB2/F");
   treeNew->Branch("passCBOld1", &passCBOld1, "passCBOld1/F");
   treeNew->Branch("passCBOld2", &passCBOld2, "passCBOld2/F");
   treeNew->Branch("passLH1", &passLH1, "passLH1/F");
   treeNew->Branch("passLH2", &passLH2, "passLH2/F");
   treeNew->Branch("sameflav", &sameflav, "sameflav/F");
   treeNew->Branch("zveto", &zveto, "zveto/F");
   treeNew->Branch("run", &run, "run/F");
   treeNew->Branch("lumi", &lumi, "lumi/F");
   treeNew->Branch("event", &event, "event/F");
   treeNew->Branch("fakeW", &fakeW, "fakeW/F");
   treeNew->Branch("WWSel", &WWSel, "WWSel/F");
   treeNew->Branch("WWSel1j", &WWSel1j, "WWSel1j/F");

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
