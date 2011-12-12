//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 11 16:51:28 2011 by ROOT version 5.27/06b
// from TTree latino/latino
// found on file: results_data/datasets_trees_looseloose/Full2011.root
//////////////////////////////////////////////////////////

#ifndef addWeightsToTreeLooseLoose_h
#define addWeightsToTreeLooseLoose_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class addWeightsToTreeLooseLoose {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         baseW;
   Float_t         bdt1;
   Float_t         bdt2;
   Float_t         ch1;
   Float_t         ch2;
   Float_t         channel;
   Float_t         chmet;
   Float_t         chmetphi;
   Float_t         dataset;
   Float_t         detajj;
   Float_t         dphill;
   Float_t         dphilljet;
   Float_t         dphilljetjet;
   Float_t         dphillmet;
   Float_t         dphilmet;
   Float_t         dphilmet1;
   Float_t         dphilmet2;
   Float_t         drll;
   Float_t         effAW;
   Float_t         effBW;
   Float_t         effW;
   Float_t         eta1;
   Float_t         eta2;
   Float_t         fakeAW;
   Float_t         fakeBW;
   Float_t         fakeW;
   Float_t         fake2W;
   Float_t         fermiW;
   Float_t         fourW;
   Float_t         gammaMRStar;
   Float_t         hardbdisc;
   Float_t         hypo;
   Float_t         imet;
   Float_t         iso1;
   Float_t         iso2;
   Float_t         jeteta1;
   Float_t         jeteta2;
   Float_t         jetphi1;
   Float_t         jetphi2;
   Float_t         jetpt1;
   Float_t         jetpt2;
   Float_t         jettche1;
   Float_t         jettche2;
   Float_t         jettchp1;
   Float_t         jettchp2;
   Float_t         kfW;
   Float_t         lh1;
   Float_t         lh2;
   Float_t         mjj;
   Float_t         mll;
   Float_t         mpmet;
   Float_t         mth;
   Float_t         mtw1;
   Float_t         mtw2;
   Float_t         nbjet;
   Float_t         nbrem1;
   Float_t         nbrem2;
   Float_t         nextra;
   Float_t         njet;
   Float_t         njetid;
   Float_t         njetvbf;
   Float_t         pchmet;
   Float_t         peaking;
   Float_t         pfmet;
   Float_t         pfmetphi;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         ppfmet;
   Float_t         predmet;
   Float_t         pt1;
   Float_t         pt2;
   Float_t         ptcmet;
   Float_t         ptll;
   Float_t         puAW;
   Float_t         puBW;
   Float_t         puW;
   Float_t         redmet;
   Float_t         sceta1;
   Float_t         sceta2;
   Float_t         softbdisc;
   Float_t         tcmet;
   Float_t         tcmetphi;
   Float_t         tightmu;
   Float_t         triggAW;
   Float_t         triggBW;
   Float_t         triggW;
   Float_t         trigger;
   Float_t         worstJetLepPt;
   Float_t         yll;
   Float_t         nvtx;
   Int_t           bveto;
   Int_t           bveto_ip;
   Int_t           bveto_mu;
   Int_t           bveto_munj;
   Int_t           bveto_nj;
   Int_t           dphiveto;
   Int_t           passBDT1;
   Int_t           passBDT2;
   Int_t           passCB1;
   Int_t           passCB2;
   Int_t           passCBOld1;
   Int_t           passCBOld2;
   Int_t           passLH1;
   Int_t           passLH2;
   Int_t           sameflav;
   Int_t           zveto;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;

   // List of branches
   TBranch        *b_baseW;   //!
   TBranch        *b_bdt1;   //!
   TBranch        *b_bdt2;   //!
   TBranch        *b_ch1;   //!
   TBranch        *b_ch2;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_chmet;   //!
   TBranch        *b_chmetphi;   //!
   TBranch        *b_dataset;   //!
   TBranch        *b_detajj;   //!
   TBranch        *b_dphill;   //!
   TBranch        *b_dphilljet;   //!
   TBranch        *b_dphilljetjet;   //!
   TBranch        *b_dphillmet;   //!
   TBranch        *b_dphilmet;   //!
   TBranch        *b_dphilmet1;   //!
   TBranch        *b_dphilmet2;   //!
   TBranch        *b_drll;   //!
   TBranch        *b_effAW;   //!
   TBranch        *b_effBW;   //!
   TBranch        *b_effW;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_fakeAW;   //!
   TBranch        *b_fakeBW;   //!
   TBranch        *b_fakeW;   //!
   TBranch        *b_fermiW;   //!
   TBranch        *b_fourW;   //!
   TBranch        *b_gammaMRStar;   //!
   TBranch        *b_hardbdisc;   //!
   TBranch        *b_hypo;   //!
   TBranch        *b_imet;   //!
   TBranch        *b_iso1;   //!
   TBranch        *b_iso2;   //!
   TBranch        *b_jeteta1;   //!
   TBranch        *b_jeteta2;   //!
   TBranch        *b_jetphi1;   //!
   TBranch        *b_jetphi2;   //!
   TBranch        *b_jetpt1;   //!
   TBranch        *b_jetpt2;   //!
   TBranch        *b_jettche1;   //!
   TBranch        *b_jettche2;   //!
   TBranch        *b_jettchp1;   //!
   TBranch        *b_jettchp2;   //!
   TBranch        *b_kfW;   //!
   TBranch        *b_lh1;   //!
   TBranch        *b_lh2;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_mll;   //!
   TBranch        *b_mpmet;   //!
   TBranch        *b_mth;   //!
   TBranch        *b_mtw1;   //!
   TBranch        *b_mtw2;   //!
   TBranch        *b_nbjet;   //!
   TBranch        *b_nbrem1;   //!
   TBranch        *b_nbrem2;   //!
   TBranch        *b_nextra;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_njetid;   //!
   TBranch        *b_njetvbf;   //!
   TBranch        *b_pchmet;   //!
   TBranch        *b_peaking;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pfmetphi;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_ppfmet;   //!
   TBranch        *b_predmet;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_ptcmet;   //!
   TBranch        *b_ptll;   //!
   TBranch        *b_puAW;   //!
   TBranch        *b_puBW;   //!
   TBranch        *b_puW;   //!
   TBranch        *b_redmet;   //!
   TBranch        *b_sceta1;   //!
   TBranch        *b_sceta2;   //!
   TBranch        *b_softbdisc;   //!
   TBranch        *b_tcmet;   //!
   TBranch        *b_tcmetphi;   //!
   TBranch        *b_tightmu;   //!
   TBranch        *b_triggAW;   //!
   TBranch        *b_triggBW;   //!
   TBranch        *b_triggW;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_worstJetLepPt;   //!
   TBranch        *b_yll;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_bveto;   //!
   TBranch        *b_bveto_ip;   //!
   TBranch        *b_bveto_mu;   //!
   TBranch        *b_bveto_munj;   //!
   TBranch        *b_bveto_nj;   //!
   TBranch        *b_dphiveto;   //!
   TBranch        *b_passBDT1;   //!
   TBranch        *b_passBDT2;   //!
   TBranch        *b_passCB1;   //!
   TBranch        *b_passCB2;   //!
   TBranch        *b_passCBOld1;   //!
   TBranch        *b_passCBOld2;   //!
   TBranch        *b_passLH1;   //!
   TBranch        *b_passLH2;   //!
   TBranch        *b_sameflav;   //!
   TBranch        *b_zveto;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!

   addWeightsToTreeLooseLoose(TTree *tree=0);
   virtual ~addWeightsToTreeLooseLoose();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef addWeightsToTreeLooseLoose_cxx
addWeightsToTreeLooseLoose::addWeightsToTreeLooseLoose(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmsrm/pc23_2/crovelli/data/Higgs4.2.X/thirdTestFull2011_latinosTree/Full2011.root");
      if (!f) {
         f = new TFile("/cmsrm/pc23_2/crovelli/data/Higgs4.2.X/thirdTestFull2011_latinosTree/Full2011.root");
      }
      tree = (TTree*)gDirectory->Get("latino");

   }
   Init(tree);
}

addWeightsToTreeLooseLoose::~addWeightsToTreeLooseLoose()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t addWeightsToTreeLooseLoose::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t addWeightsToTreeLooseLoose::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void addWeightsToTreeLooseLoose::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("baseW", &baseW, &b_baseW);
   fChain->SetBranchAddress("bdt1", &bdt1, &b_bdt1);
   fChain->SetBranchAddress("bdt2", &bdt2, &b_bdt2);
   fChain->SetBranchAddress("ch1", &ch1, &b_ch1);
   fChain->SetBranchAddress("ch2", &ch2, &b_ch2);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("chmet", &chmet, &b_chmet);
   fChain->SetBranchAddress("chmetphi", &chmetphi, &b_chmetphi);
   fChain->SetBranchAddress("dataset", &dataset, &b_dataset);
   fChain->SetBranchAddress("detajj", &detajj, &b_detajj);
   fChain->SetBranchAddress("dphill", &dphill, &b_dphill);
   fChain->SetBranchAddress("dphilljet", &dphilljet, &b_dphilljet);
   fChain->SetBranchAddress("dphilljetjet", &dphilljetjet, &b_dphilljetjet);
   fChain->SetBranchAddress("dphillmet", &dphillmet, &b_dphillmet);
   fChain->SetBranchAddress("dphilmet", &dphilmet, &b_dphilmet);
   fChain->SetBranchAddress("dphilmet1", &dphilmet1, &b_dphilmet1);
   fChain->SetBranchAddress("dphilmet2", &dphilmet2, &b_dphilmet2);
   fChain->SetBranchAddress("drll", &drll, &b_drll);
   fChain->SetBranchAddress("effAW", &effAW, &b_effAW);
   fChain->SetBranchAddress("effBW", &effBW, &b_effBW);
   fChain->SetBranchAddress("effW", &effW, &b_effW);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("fakeAW", &fakeAW, &b_fakeAW);
   fChain->SetBranchAddress("fakeBW", &fakeBW, &b_fakeBW);
   fChain->SetBranchAddress("fermiW", &fermiW, &b_fermiW);
   fChain->SetBranchAddress("fourW", &fourW, &b_fourW);
   fChain->SetBranchAddress("gammaMRStar", &gammaMRStar, &b_gammaMRStar);
   fChain->SetBranchAddress("hardbdisc", &hardbdisc, &b_hardbdisc);
   fChain->SetBranchAddress("hypo", &hypo, &b_hypo);
   fChain->SetBranchAddress("imet", &imet, &b_imet);
   fChain->SetBranchAddress("iso1", &iso1, &b_iso1);
   fChain->SetBranchAddress("iso2", &iso2, &b_iso2);
   fChain->SetBranchAddress("jeteta1", &jeteta1, &b_jeteta1);
   fChain->SetBranchAddress("jeteta2", &jeteta2, &b_jeteta2);
   fChain->SetBranchAddress("jetphi1", &jetphi1, &b_jetphi1);
   fChain->SetBranchAddress("jetphi2", &jetphi2, &b_jetphi2);
   fChain->SetBranchAddress("jetpt1", &jetpt1, &b_jetpt1);
   fChain->SetBranchAddress("jetpt2", &jetpt2, &b_jetpt2);
   fChain->SetBranchAddress("jettche1", &jettche1, &b_jettche1);
   fChain->SetBranchAddress("jettche2", &jettche2, &b_jettche2);
   fChain->SetBranchAddress("jettchp1", &jettchp1, &b_jettchp1);
   fChain->SetBranchAddress("jettchp2", &jettchp2, &b_jettchp2);
   fChain->SetBranchAddress("kfW", &kfW, &b_kfW);
   fChain->SetBranchAddress("lh1", &lh1, &b_lh1);
   fChain->SetBranchAddress("lh2", &lh2, &b_lh2);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
   fChain->SetBranchAddress("mpmet", &mpmet, &b_mpmet);
   fChain->SetBranchAddress("mth", &mth, &b_mth);
   fChain->SetBranchAddress("mtw1", &mtw1, &b_mtw1);
   fChain->SetBranchAddress("mtw2", &mtw2, &b_mtw2);
   fChain->SetBranchAddress("nbjet", &nbjet, &b_nbjet);
   fChain->SetBranchAddress("nbrem1", &nbrem1, &b_nbrem1);
   fChain->SetBranchAddress("nbrem2", &nbrem2, &b_nbrem2);
   fChain->SetBranchAddress("nextra", &nextra, &b_nextra);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("njetid", &njetid, &b_njetid);
   fChain->SetBranchAddress("njetvbf", &njetvbf, &b_njetvbf);
   fChain->SetBranchAddress("pchmet", &pchmet, &b_pchmet);
   fChain->SetBranchAddress("peaking", &peaking, &b_peaking);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pfmetphi", &pfmetphi, &b_pfmetphi);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("ppfmet", &ppfmet, &b_ppfmet);
   fChain->SetBranchAddress("predmet", &predmet, &b_predmet);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("ptcmet", &ptcmet, &b_ptcmet);
   fChain->SetBranchAddress("ptll", &ptll, &b_ptll);
   fChain->SetBranchAddress("puAW", &puAW, &b_puAW);
   fChain->SetBranchAddress("puBW", &puBW, &b_puBW);
   fChain->SetBranchAddress("puW", &puW, &b_puW);
   fChain->SetBranchAddress("redmet", &redmet, &b_redmet);
   fChain->SetBranchAddress("sceta1", &sceta1, &b_sceta1);
   fChain->SetBranchAddress("sceta2", &sceta2, &b_sceta2);
   fChain->SetBranchAddress("softbdisc", &softbdisc, &b_softbdisc);
   fChain->SetBranchAddress("tcmet", &tcmet, &b_tcmet);
   fChain->SetBranchAddress("tcmetphi", &tcmetphi, &b_tcmetphi);
   fChain->SetBranchAddress("tightmu", &tightmu, &b_tightmu);
   fChain->SetBranchAddress("triggAW", &triggAW, &b_triggAW);
   fChain->SetBranchAddress("triggBW", &triggBW, &b_triggBW);
   fChain->SetBranchAddress("triggW", &triggW, &b_triggW);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("worstJetLepPt", &worstJetLepPt, &b_worstJetLepPt);
   fChain->SetBranchAddress("yll", &yll, &b_yll);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("bveto", &bveto, &b_bveto);
   fChain->SetBranchAddress("bveto_ip", &bveto_ip, &b_bveto_ip);
   fChain->SetBranchAddress("bveto_mu", &bveto_mu, &b_bveto_mu);
   fChain->SetBranchAddress("bveto_munj", &bveto_munj, &b_bveto_munj);
   fChain->SetBranchAddress("bveto_nj", &bveto_nj, &b_bveto_nj);
   fChain->SetBranchAddress("dphiveto", &dphiveto, &b_dphiveto);
   fChain->SetBranchAddress("passBDT1", &passBDT1, &b_passBDT1);
   fChain->SetBranchAddress("passBDT2", &passBDT2, &b_passBDT2);
   fChain->SetBranchAddress("passCB1", &passCB1, &b_passCB1);
   fChain->SetBranchAddress("passCB2", &passCB2, &b_passCB2);
   fChain->SetBranchAddress("passCBOld1", &passCBOld1, &b_passCBOld1);
   fChain->SetBranchAddress("passCBOld2", &passCBOld2, &b_passCBOld2);
   fChain->SetBranchAddress("passLH1", &passLH1, &b_passLH1);
   fChain->SetBranchAddress("passLH2", &passLH2, &b_passLH2);
   fChain->SetBranchAddress("sameflav", &sameflav, &b_sameflav);
   fChain->SetBranchAddress("zveto", &zveto, &b_zveto);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("fake2W", &fakeW, &b_fakeW);
   Notify();
}

Bool_t addWeightsToTreeLooseLoose::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void addWeightsToTreeLooseLoose::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t addWeightsToTreeLooseLoose::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef addWeightsToTreeLooseLoose_cxx
