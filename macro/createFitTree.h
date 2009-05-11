//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  6 14:32:09 2009 by ROOT version 5.18/00a
// from TTree T1/vecbos analysis tree for Z studies
// found on file: WjetsMADGRAPH_Fall08/root/WjetsMADGRAPH_Fall08_102_VecBosOutput-out-Wenu.root
//////////////////////////////////////////////////////////

#ifndef createFitTree_h
#define createFitTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

class createFitTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         met;
   Float_t         deltaPhi;
   Float_t         transvMass;
   Float_t         eleInvMass;
   Float_t         maxPtEle;
   Float_t         minPtEle;
   Float_t         detaLeptons;
   Char_t          finalLeptons;
   Char_t          jetVeto;
   Char_t          preDeltaPhi;
   Char_t          finalSelection;
   Char_t          HLTSingleElectron;
   Char_t          HLTSingleElectronRelaxed;
   Char_t          HLTSingleElectronOR;
   Float_t         KFactor;
   Char_t          promptDecay;
   Float_t         maxPtLh;
   Float_t         minPtLh;
   Int_t           njets;
   Float_t         dxyEVT;
   Float_t         dszEVT;

   // List of branches
   TBranch        *b_met;   //!
   TBranch        *b_deltaPhi;   //!
   TBranch        *b_transvMass;   //!
   TBranch        *b_eleInvMass;   //!
   TBranch        *b_maxPtEle;   //!
   TBranch        *b_minPtEle;   //!
   TBranch        *b_detaLeptons;   //!
   TBranch        *b_finalLeptons;   //!
   TBranch        *b_jetVeto;   //!
   TBranch        *b_preDeltaPhi;   //!
   TBranch        *b_finalSelection;   //!
   TBranch        *b_HLTSingleElectron;   //!
   TBranch        *b_HLTSingleElectronRelaxed;   //!
   TBranch        *b_HLTSingleElectronOR;   //!
   TBranch        *b_KFactor;   //!
   TBranch        *b_promptDecay;   //!
   TBranch        *b_maxPtLh;   //!
   TBranch        *b_minPtLh;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_dxyEVT;   //!
   TBranch        *b_dszEVT;   //!

   createFitTree(const char *option="recreate");
   virtual ~createFitTree();
   virtual void     AddFiles(const char *wildcard);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SetDatasetName(const char *datasetname);
   virtual void     SaveRootTree(bool save) { _saveTree = save; }
   
   TChain *_chain;
   char _datasetname[200];
   char _option[100];
   bool _saveTree;

};

#endif

#ifdef createFitTree_cxx
createFitTree::createFitTree(const char *option) {

  _chain = new TChain("T1","");
  sprintf(_option, "%s", option);
  _saveTree = false;

}
void createFitTree::AddFiles(const char *wildcard) {

  _chain->Add(wildcard);

}
createFitTree::~createFitTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t createFitTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t createFitTree::LoadTree(Long64_t entry)
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

void createFitTree::Init()
{
  
  TTree *tree = _chain;

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

   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("deltaPhi", &deltaPhi, &b_deltaPhi);
   fChain->SetBranchAddress("transvMass", &transvMass, &b_transvMass);
   fChain->SetBranchAddress("eleInvMass", &eleInvMass, &b_eleInvMass);
   fChain->SetBranchAddress("maxPtEle", &maxPtEle, &b_maxPtEle);
   fChain->SetBranchAddress("minPtEle", &minPtEle, &b_minPtEle);
   fChain->SetBranchAddress("detaLeptons", &detaLeptons, &b_detaLeptons);
   fChain->SetBranchAddress("finalLeptons", &finalLeptons, &b_finalLeptons);
   fChain->SetBranchAddress("jetVeto", &jetVeto, &b_jetVeto);
   fChain->SetBranchAddress("preDeltaPhi", &preDeltaPhi, &b_preDeltaPhi);
   fChain->SetBranchAddress("finalSelection", &finalSelection, &b_finalSelection);
   fChain->SetBranchAddress("HLTSingleElectron", &HLTSingleElectron, &b_HLTSingleElectron);
   fChain->SetBranchAddress("HLTSingleElectronRelaxed", &HLTSingleElectronRelaxed, &b_HLTSingleElectronRelaxed);
   fChain->SetBranchAddress("HLTSingleElectronOR", &HLTSingleElectronOR, &b_HLTSingleElectronOR);
   fChain->SetBranchAddress("KFactor", &KFactor, &b_KFactor);
   fChain->SetBranchAddress("promptDecay", &promptDecay, &b_promptDecay);
   fChain->SetBranchAddress("maxPtLh", &maxPtLh, &b_maxPtLh);
   fChain->SetBranchAddress("minPtLh", &minPtLh, &b_minPtLh);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("dxyEVT", &dxyEVT, &b_dxyEVT);
   fChain->SetBranchAddress("dszEVT", &dszEVT, &b_dszEVT);

   Notify();
}

Bool_t createFitTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void createFitTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t createFitTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void createFitTree::SetDatasetName(const char *datasetname)
{
  sprintf(_datasetname,"%s",datasetname);
}
#endif // #ifdef createFitTree_cxx
