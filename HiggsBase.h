//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 13 10:17:06 2007 by ROOT version 5.14/00b
// from TTree ntp1/ntp1
// found on file: H160_WW_2l/default_1.root
//////////////////////////////////////////////////////////

#ifndef HiggsBase_h
#define HiggsBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class HiggsBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           nMc;
   Float_t         pMc[101];   //[nMc]
   Float_t         massMc[101];   //[nMc]
   Float_t         thetaMc[101];   //[nMc]
   Float_t         etaMc[101];   //[nMc]
   Float_t         phiMc[101];   //[nMc]
   Int_t           idMc[101];   //[nMc]
   Int_t           mothMc[101];   //[nMc]
   Int_t           nDauMc[101];   //[nMc]
   Float_t         xMc[101];   //[nMc]
   Float_t         yMc[101];   //[nMc]
   Float_t         zMc[101];   //[nMc]
   UChar_t         singleElePassedTrg;
   UChar_t         singleEleRelaxPassedTrg;
   UChar_t         doubleElePassedTrg;
   UChar_t         doubleEleRelaxPassedTrg;
   UChar_t         singlePhoPassedTrg;
   UChar_t         singlePhoRelaxPassedTrg;
   UChar_t         doublePhoPassedTrg;
   UChar_t         doublePhoRelaxPassedTrg;
   UChar_t         highEMPassedTrg;
   UChar_t         veryHighEMPassedTrg;
   Int_t           nEle;
   Int_t           chargeEle[50];   //[nEle]
   Float_t         energyEle[50];   //[nEle]
   Float_t         etEle[50];   //[nEle]
   Float_t         momentumEle[50];   //[nEle]
   Float_t         thetaEle[50];   //[nEle]
   Float_t         etaEle[50];   //[nEle]
   Float_t         phiEle[50];   //[nEle]
   Float_t         pxEle[50];   //[nEle]
   Float_t         pyEle[50];   //[nEle]
   Float_t         pzEle[50];   //[nEle]
   Float_t         vertexXEle[50];   //[nEle]
   Float_t         vertexYEle[50];   //[nEle]
   Float_t         vertexZEle[50];   //[nEle]
   Float_t         massEle[50];   //[nEle]
   Float_t         mtEle[50];   //[nEle]
   Int_t           pdgIdEle[50];   //[nEle]
   Float_t         ecalEle[50];   //[nEle]
   Int_t           nCluEle[50];   //[nEle]
   Int_t           nCryEle[50];   //[nEle]
   Float_t         e3x3Ele[50];   //[nEle]
   Float_t         e5x5Ele[50];   //[nEle]
   Float_t         eMaxEle[50];   //[nEle]
   Float_t         latEle[50];   //[nEle]
   Float_t         erawEle[50];   //[nEle]
   Float_t         caloEtaEle[50];   //[nEle]
   Float_t         caloPhiEle[50];   //[nEle]
   Float_t         e2x2Ele[50];   //[nEle]
   Float_t         e2ndEle[50];   //[nEle]
   Float_t         s1s9Ele[50];   //[nEle]
   Float_t         s9s25Ele[50];   //[nEle]
   Float_t         covEtaEtaEle[50];   //[nEle]
   Float_t         covEtaPhiEle[50];   //[nEle]
   Float_t         covPhiPhiEle[50];   //[nEle]
   Float_t         a00Ele[50];   //[nEle]
   Float_t         a11Ele[50];   //[nEle]
   Float_t         a20Ele[50];   //[nEle]
   Float_t         a22Ele[50];   //[nEle]
   Float_t         a42Ele[50];   //[nEle]
   Float_t         pxAtOuterEle[50];   //[nEle]
   Float_t         pyAtOuterEle[50];   //[nEle]
   Float_t         pzAtOuterEle[50];   //[nEle]
   Float_t         xAtOuterEle[50];   //[nEle]
   Float_t         yAtOuterEle[50];   //[nEle]
   Float_t         zAtOuterEle[50];   //[nEle]
   Float_t         pxAtInnerEle[50];   //[nEle]
   Float_t         pyAtInnerEle[50];   //[nEle]
   Float_t         pzAtInnerEle[50];   //[nEle]
   Float_t         xAtInnerEle[50];   //[nEle]
   Float_t         yAtInnerEle[50];   //[nEle]
   Float_t         zAtInnerEle[50];   //[nEle]
   Float_t         eleFullCorrEEle[50];   //[nEle]
   Float_t         eleCaloCorrEEle[50];   //[nEle]
   Float_t         eleNxtalCorrEEle[50];   //[nEle]
   Float_t         eleRawEEle[50];   //[nEle]
   Float_t         eleTrackerPEle[50];   //[nEle]
   Int_t           eleClassEle[50];   //[nEle]
   Float_t         eleHoEEle[50];   //[nEle]
   Float_t         eleCorrEoPEle[50];   //[nEle]
   Float_t         eleNotCorrEoPEle[50];   //[nEle]
   Float_t         eleCorrEoPoutEle[50];   //[nEle]
   Float_t         eleNotCorrEoPoutEle[50];   //[nEle]
   Float_t         eleDeltaEtaAtVtxEle[50];   //[nEle]
   Float_t         eleDeltaPhiAtVtxEle[50];   //[nEle]
   Float_t         eleDeltaEtaAtCaloEle[50];   //[nEle]
   Float_t         eleDeltaPhiAtCaloEle[50];   //[nEle]
   Float_t         eleTrackerIso_minDREle[50];   //[nEle]
   Float_t         eleTrackerIso_minDR_vetoEle[50];   //[nEle]
   Float_t         eleTrackerIso_sumPtEle[50];   //[nEle]
   Float_t         eleCaloIso_minDREle[50];   //[nEle]
   Float_t         eleCaloIso_sumPtEle[50];   //[nEle]
   Float_t         eleLikelihoodEle[50];   //[nEle]
   Float_t         eleTipEle[50];   //[nEle]
   Int_t           nMet;
   Int_t           chargeMet[1];   //[nMet]
   Float_t         energyMet[1];   //[nMet]
   Float_t         etMet[1];   //[nMet]
   Float_t         momentumMet[1];   //[nMet]
   Float_t         thetaMet[1];   //[nMet]
   Float_t         etaMet[1];   //[nMet]
   Float_t         phiMet[1];   //[nMet]
   Float_t         pxMet[1];   //[nMet]
   Float_t         pyMet[1];   //[nMet]
   Float_t         pzMet[1];   //[nMet]
   Float_t         vertexXMet[1];   //[nMet]
   Float_t         vertexYMet[1];   //[nMet]
   Float_t         vertexZMet[1];   //[nMet]
   Float_t         massMet[1];   //[nMet]
   Float_t         mtMet[1];   //[nMet]
   Int_t           pdgIdMet[1];   //[nMet]
   Int_t           nGenMet;
   Int_t           chargeGenMet[1];   //[nGenMet]
   Float_t         energyGenMet[1];   //[nGenMet]
   Float_t         etGenMet[1];   //[nGenMet]
   Float_t         momentumGenMet[1];   //[nGenMet]
   Float_t         thetaGenMet[1];   //[nGenMet]
   Float_t         etaGenMet[1];   //[nGenMet]
   Float_t         phiGenMet[1];   //[nGenMet]
   Float_t         pxGenMet[1];   //[nGenMet]
   Float_t         pyGenMet[1];   //[nGenMet]
   Float_t         pzGenMet[1];   //[nGenMet]
   Float_t         vertexXGenMet[1];   //[nGenMet]
   Float_t         vertexYGenMet[1];   //[nGenMet]
   Float_t         vertexZGenMet[1];   //[nGenMet]
   Float_t         massGenMet[1];   //[nGenMet]
   Float_t         mtGenMet[1];   //[nGenMet]
   Int_t           pdgIdGenMet[1];   //[nGenMet]
   Int_t           nJet;
   Int_t           chargeJet[500];   //[nJet]
   Float_t         energyJet[500];   //[nJet]
   Float_t         etJet[500];   //[nJet]
   Float_t         momentumJet[500];   //[nJet]
   Float_t         thetaJet[500];   //[nJet]
   Float_t         etaJet[500];   //[nJet]
   Float_t         phiJet[500];   //[nJet]
   Float_t         pxJet[500];   //[nJet]
   Float_t         pyJet[500];   //[nJet]
   Float_t         pzJet[500];   //[nJet]
   Float_t         vertexXJet[500];   //[nJet]
   Float_t         vertexYJet[500];   //[nJet]
   Float_t         vertexZJet[500];   //[nJet]
   Float_t         massJet[500];   //[nJet]
   Float_t         mtJet[500];   //[nJet]
   Int_t           pdgIdJet[500];   //[nJet]
   Float_t         alphaJet[500];   //[nJet]
   Float_t         emFracJet[500];   //[nJet]
   Float_t         hadFracJet[500];   //[nJet]
   Int_t           nGenJet;
   Int_t           chargeGenJet[500];   //[nGenJet]
   Float_t         energyGenJet[500];   //[nGenJet]
   Float_t         etGenJet[500];   //[nGenJet]
   Float_t         momentumGenJet[500];   //[nGenJet]
   Float_t         thetaGenJet[500];   //[nGenJet]
   Float_t         etaGenJet[500];   //[nGenJet]
   Float_t         phiGenJet[500];   //[nGenJet]
   Float_t         pxGenJet[500];   //[nGenJet]
   Float_t         pyGenJet[500];   //[nGenJet]
   Float_t         pzGenJet[500];   //[nGenJet]
   Float_t         vertexXGenJet[500];   //[nGenJet]
   Float_t         vertexYGenJet[500];   //[nGenJet]
   Float_t         vertexZGenJet[500];   //[nGenJet]
   Float_t         massGenJet[500];   //[nGenJet]
   Float_t         mtGenJet[500];   //[nGenJet]
   Int_t           pdgIdGenJet[500];   //[nGenJet]
   Double_t        genProcessId;
   Double_t        genPtHat;
   Double_t        genFilterEff;
   Double_t        genXsec;
   Double_t        genWeight;

   // List of branches
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_massMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_nDauMc;   //!
   TBranch        *b_xMc;   //!
   TBranch        *b_yMc;   //!
   TBranch        *b_zMc;   //!
   TBranch        *b_singleElePassedTrg;   //!
   TBranch        *b_singleEleRelaxPassedTrg;   //!
   TBranch        *b_doubleElePassedTrg;   //!
   TBranch        *b_doubleEleRelaxPassedTrg;   //!
   TBranch        *b_singlePhoPassedTrg;   //!
   TBranch        *b_singlePhoRelaxPassedTrg;   //!
   TBranch        *b_doublePhoPassedTrg;   //!
   TBranch        *b_doublePhoRelaxPassedTrg;   //!
   TBranch        *b_highEMPassedTrg;   //!
   TBranch        *b_veryHighEMPassedTrg;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_etEle;   //!
   TBranch        *b_momentumEle;   //!
   TBranch        *b_thetaEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_pxEle;   //!
   TBranch        *b_pyEle;   //!
   TBranch        *b_pzEle;   //!
   TBranch        *b_vertexXEle;   //!
   TBranch        *b_vertexYEle;   //!
   TBranch        *b_vertexZEle;   //!
   TBranch        *b_massEle;   //!
   TBranch        *b_mtEle;   //!
   TBranch        *b_pdgIdEle;   //!
   TBranch        *b_ecalEle;   //!
   TBranch        *b_nCluEle;   //!
   TBranch        *b_nCryEle;   //!
   TBranch        *b_e3x3Ele;   //!
   TBranch        *b_e5x5Ele;   //!
   TBranch        *b_eMaxEle;   //!
   TBranch        *b_latEle;   //!
   TBranch        *b_erawEle;   //!
   TBranch        *b_caloEtaEle;   //!
   TBranch        *b_caloPhiEle;   //!
   TBranch        *b_e2x2Ele;   //!
   TBranch        *b_e2ndEle;   //!
   TBranch        *b_s1s9Ele;   //!
   TBranch        *b_s9s25Ele;   //!
   TBranch        *b_covEtaEtaEle;   //!
   TBranch        *b_covEtaPhiEle;   //!
   TBranch        *b_covPhiPhiEle;   //!
   TBranch        *b_a00Ele;   //!
   TBranch        *b_a11Ele;   //!
   TBranch        *b_a20Ele;   //!
   TBranch        *b_a22Ele;   //!
   TBranch        *b_a42Ele;   //!
   TBranch        *b_pxAtOuterEle;   //!
   TBranch        *b_pyAtOuterEle;   //!
   TBranch        *b_pzAtOuterEle;   //!
   TBranch        *b_xAtOuterEle;   //!
   TBranch        *b_yAtOuterEle;   //!
   TBranch        *b_zAtOuterEle;   //!
   TBranch        *b_pxAtInnerEle;   //!
   TBranch        *b_pyAtInnerEle;   //!
   TBranch        *b_pzAtInnerEle;   //!
   TBranch        *b_xAtInnerEle;   //!
   TBranch        *b_yAtInnerEle;   //!
   TBranch        *b_zAtInnerEle;   //!
   TBranch        *b_eleFullCorrEEle;   //!
   TBranch        *b_eleCaloCorrEEle;   //!
   TBranch        *b_eleNxtalCorrEEle;   //!
   TBranch        *b_eleRawEEle;   //!
   TBranch        *b_eleTrackerPEle;   //!
   TBranch        *b_eleClassEle;   //!
   TBranch        *b_eleHoEEle;   //!
   TBranch        *b_eleCorrEoPEle;   //!
   TBranch        *b_eleNotCorrEoPEle;   //!
   TBranch        *b_eleCorrEoPoutEle;   //!
   TBranch        *b_eleNotCorrEoPoutEle;   //!
   TBranch        *b_eleDeltaEtaAtVtxEle;   //!
   TBranch        *b_eleDeltaPhiAtVtxEle;   //!
   TBranch        *b_eleDeltaEtaAtCaloEle;   //!
   TBranch        *b_eleDeltaPhiAtCaloEle;   //!
   TBranch        *b_eleTrackerIso_minDREle;   //!
   TBranch        *b_eleTrackerIso_minDR_vetoEle;   //!
   TBranch        *b_eleTrackerIso_sumPtEle;   //!
   TBranch        *b_eleCaloIso_minDREle;   //!
   TBranch        *b_eleCaloIso_sumPtEle;   //!
   TBranch        *b_eleLikelihoodEle;   //!
   TBranch        *b_eleTipEle;   //!
   TBranch        *b_nMet;   //!
   TBranch        *b_chargeMet;   //!
   TBranch        *b_energyMet;   //!
   TBranch        *b_etMet;   //!
   TBranch        *b_momentumMet;   //!
   TBranch        *b_thetaMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_vertexXMet;   //!
   TBranch        *b_vertexYMet;   //!
   TBranch        *b_vertexZMet;   //!
   TBranch        *b_massMet;   //!
   TBranch        *b_mtMet;   //!
   TBranch        *b_pdgIdMet;   //!
   TBranch        *b_nGenMet;   //!
   TBranch        *b_chargeGenMet;   //!
   TBranch        *b_energyGenMet;   //!
   TBranch        *b_etGenMet;   //!
   TBranch        *b_momentumGenMet;   //!
   TBranch        *b_thetaGenMet;   //!
   TBranch        *b_etaGenMet;   //!
   TBranch        *b_phiGenMet;   //!
   TBranch        *b_pxGenMet;   //!
   TBranch        *b_pyGenMet;   //!
   TBranch        *b_pzGenMet;   //!
   TBranch        *b_vertexXGenMet;   //!
   TBranch        *b_vertexYGenMet;   //!
   TBranch        *b_vertexZGenMet;   //!
   TBranch        *b_massGenMet;   //!
   TBranch        *b_mtGenMet;   //!
   TBranch        *b_pdgIdGenMet;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_chargeJet;   //!
   TBranch        *b_energyJet;   //!
   TBranch        *b_etJet;   //!
   TBranch        *b_momentumJet;   //!
   TBranch        *b_thetaJet;   //!
   TBranch        *b_etaJet;   //!
   TBranch        *b_phiJet;   //!
   TBranch        *b_pxJet;   //!
   TBranch        *b_pyJet;   //!
   TBranch        *b_pzJet;   //!
   TBranch        *b_vertexXJet;   //!
   TBranch        *b_vertexYJet;   //!
   TBranch        *b_vertexZJet;   //!
   TBranch        *b_massJet;   //!
   TBranch        *b_mtJet;   //!
   TBranch        *b_pdgIdJet;   //!
   TBranch        *b_alphaJet;   //!
   TBranch        *b_emFracJet;   //!
   TBranch        *b_hadFracJet;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_chargeGenJet;   //!
   TBranch        *b_energyGenJet;   //!
   TBranch        *b_etGenJet;   //!
   TBranch        *b_momentumGenJet;   //!
   TBranch        *b_thetaGenJet;   //!
   TBranch        *b_etaGenJet;   //!
   TBranch        *b_phiGenJet;   //!
   TBranch        *b_pxGenJet;   //!
   TBranch        *b_pyGenJet;   //!
   TBranch        *b_pzGenJet;   //!
   TBranch        *b_vertexXGenJet;   //!
   TBranch        *b_vertexYGenJet;   //!
   TBranch        *b_vertexZGenJet;   //!
   TBranch        *b_massGenJet;   //!
   TBranch        *b_mtGenJet;   //!
   TBranch        *b_pdgIdGenJet;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genPtHat;   //!
   TBranch        *b_genFilterEff;   //!
   TBranch        *b_genXsec;   //!
   TBranch        *b_genWeight;   //!


   HiggsBase(TTree *tree=0);
   virtual ~HiggsBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HiggsBase_cxx
HiggsBase::HiggsBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("H160_WW_2l/default_1.root");
      if (!f) {
         f = new TFile("H160_WW_2l/default_1.root");
      }
      tree = (TTree*)gDirectory->Get("ntp1");

   }
   Init(tree);
}

HiggsBase::~HiggsBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HiggsBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HiggsBase::LoadTree(Long64_t entry)
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

void HiggsBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("massMc", massMc, &b_massMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("nDauMc", nDauMc, &b_nDauMc);
   fChain->SetBranchAddress("xMc", xMc, &b_xMc);
   fChain->SetBranchAddress("yMc", yMc, &b_yMc);
   fChain->SetBranchAddress("zMc", zMc, &b_zMc);
   fChain->SetBranchAddress("singleElePassedTrg", &singleElePassedTrg, &b_singleElePassedTrg);
   fChain->SetBranchAddress("singleEleRelaxPassedTrg", &singleEleRelaxPassedTrg, &b_singleEleRelaxPassedTrg);
   fChain->SetBranchAddress("doubleElePassedTrg", &doubleElePassedTrg, &b_doubleElePassedTrg);
   fChain->SetBranchAddress("doubleEleRelaxPassedTrg", &doubleEleRelaxPassedTrg, &b_doubleEleRelaxPassedTrg);
   fChain->SetBranchAddress("singlePhoPassedTrg", &singlePhoPassedTrg, &b_singlePhoPassedTrg);
   fChain->SetBranchAddress("singlePhoRelaxPassedTrg", &singlePhoRelaxPassedTrg, &b_singlePhoRelaxPassedTrg);
   fChain->SetBranchAddress("doublePhoPassedTrg", &doublePhoPassedTrg, &b_doublePhoPassedTrg);
   fChain->SetBranchAddress("doublePhoRelaxPassedTrg", &doublePhoRelaxPassedTrg, &b_doublePhoRelaxPassedTrg);
   fChain->SetBranchAddress("highEMPassedTrg", &highEMPassedTrg, &b_highEMPassedTrg);
   fChain->SetBranchAddress("veryHighEMPassedTrg", &veryHighEMPassedTrg, &b_veryHighEMPassedTrg);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("etEle", etEle, &b_etEle);
   fChain->SetBranchAddress("momentumEle", momentumEle, &b_momentumEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("massEle", massEle, &b_massEle);
   fChain->SetBranchAddress("mtEle", mtEle, &b_mtEle);
   fChain->SetBranchAddress("pdgIdEle", pdgIdEle, &b_pdgIdEle);
   fChain->SetBranchAddress("ecalEle", ecalEle, &b_ecalEle);
   fChain->SetBranchAddress("nCluEle", nCluEle, &b_nCluEle);
   fChain->SetBranchAddress("nCryEle", nCryEle, &b_nCryEle);
   fChain->SetBranchAddress("e3x3Ele", e3x3Ele, &b_e3x3Ele);
   fChain->SetBranchAddress("e5x5Ele", e5x5Ele, &b_e5x5Ele);
   fChain->SetBranchAddress("eMaxEle", eMaxEle, &b_eMaxEle);
   fChain->SetBranchAddress("latEle", latEle, &b_latEle);
   fChain->SetBranchAddress("erawEle", erawEle, &b_erawEle);
   fChain->SetBranchAddress("caloEtaEle", caloEtaEle, &b_caloEtaEle);
   fChain->SetBranchAddress("caloPhiEle", caloPhiEle, &b_caloPhiEle);
   fChain->SetBranchAddress("e2x2Ele", e2x2Ele, &b_e2x2Ele);
   fChain->SetBranchAddress("e2ndEle", e2ndEle, &b_e2ndEle);
   fChain->SetBranchAddress("s1s9Ele", s1s9Ele, &b_s1s9Ele);
   fChain->SetBranchAddress("s9s25Ele", s9s25Ele, &b_s9s25Ele);
   fChain->SetBranchAddress("covEtaEtaEle", covEtaEtaEle, &b_covEtaEtaEle);
   fChain->SetBranchAddress("covEtaPhiEle", covEtaPhiEle, &b_covEtaPhiEle);
   fChain->SetBranchAddress("covPhiPhiEle", covPhiPhiEle, &b_covPhiPhiEle);
   fChain->SetBranchAddress("a00Ele", a00Ele, &b_a00Ele);
   fChain->SetBranchAddress("a11Ele", a11Ele, &b_a11Ele);
   fChain->SetBranchAddress("a20Ele", a20Ele, &b_a20Ele);
   fChain->SetBranchAddress("a22Ele", a22Ele, &b_a22Ele);
   fChain->SetBranchAddress("a42Ele", a42Ele, &b_a42Ele);
   fChain->SetBranchAddress("pxAtOuterEle", pxAtOuterEle, &b_pxAtOuterEle);
   fChain->SetBranchAddress("pyAtOuterEle", pyAtOuterEle, &b_pyAtOuterEle);
   fChain->SetBranchAddress("pzAtOuterEle", pzAtOuterEle, &b_pzAtOuterEle);
   fChain->SetBranchAddress("xAtOuterEle", xAtOuterEle, &b_xAtOuterEle);
   fChain->SetBranchAddress("yAtOuterEle", yAtOuterEle, &b_yAtOuterEle);
   fChain->SetBranchAddress("zAtOuterEle", zAtOuterEle, &b_zAtOuterEle);
   fChain->SetBranchAddress("pxAtInnerEle", pxAtInnerEle, &b_pxAtInnerEle);
   fChain->SetBranchAddress("pyAtInnerEle", pyAtInnerEle, &b_pyAtInnerEle);
   fChain->SetBranchAddress("pzAtInnerEle", pzAtInnerEle, &b_pzAtInnerEle);
   fChain->SetBranchAddress("xAtInnerEle", xAtInnerEle, &b_xAtInnerEle);
   fChain->SetBranchAddress("yAtInnerEle", yAtInnerEle, &b_yAtInnerEle);
   fChain->SetBranchAddress("zAtInnerEle", zAtInnerEle, &b_zAtInnerEle);
   fChain->SetBranchAddress("eleFullCorrEEle", eleFullCorrEEle, &b_eleFullCorrEEle);
   fChain->SetBranchAddress("eleCaloCorrEEle", eleCaloCorrEEle, &b_eleCaloCorrEEle);
   fChain->SetBranchAddress("eleNxtalCorrEEle", eleNxtalCorrEEle, &b_eleNxtalCorrEEle);
   fChain->SetBranchAddress("eleRawEEle", eleRawEEle, &b_eleRawEEle);
   fChain->SetBranchAddress("eleTrackerPEle", eleTrackerPEle, &b_eleTrackerPEle);
   fChain->SetBranchAddress("eleClassEle", eleClassEle, &b_eleClassEle);
   fChain->SetBranchAddress("eleHoEEle", eleHoEEle, &b_eleHoEEle);
   fChain->SetBranchAddress("eleCorrEoPEle", eleCorrEoPEle, &b_eleCorrEoPEle);
   fChain->SetBranchAddress("eleNotCorrEoPEle", eleNotCorrEoPEle, &b_eleNotCorrEoPEle);
   fChain->SetBranchAddress("eleCorrEoPoutEle", eleCorrEoPoutEle, &b_eleCorrEoPoutEle);
   fChain->SetBranchAddress("eleNotCorrEoPoutEle", eleNotCorrEoPoutEle, &b_eleNotCorrEoPoutEle);
   fChain->SetBranchAddress("eleDeltaEtaAtVtxEle", eleDeltaEtaAtVtxEle, &b_eleDeltaEtaAtVtxEle);
   fChain->SetBranchAddress("eleDeltaPhiAtVtxEle", eleDeltaPhiAtVtxEle, &b_eleDeltaPhiAtVtxEle);
   fChain->SetBranchAddress("eleDeltaEtaAtCaloEle", eleDeltaEtaAtCaloEle, &b_eleDeltaEtaAtCaloEle);
   fChain->SetBranchAddress("eleDeltaPhiAtCaloEle", eleDeltaPhiAtCaloEle, &b_eleDeltaPhiAtCaloEle);
   fChain->SetBranchAddress("eleTrackerIso_minDREle", eleTrackerIso_minDREle, &b_eleTrackerIso_minDREle);
   fChain->SetBranchAddress("eleTrackerIso_minDR_vetoEle", eleTrackerIso_minDR_vetoEle, &b_eleTrackerIso_minDR_vetoEle);
   fChain->SetBranchAddress("eleTrackerIso_sumPtEle", eleTrackerIso_sumPtEle, &b_eleTrackerIso_sumPtEle);
   fChain->SetBranchAddress("eleCaloIso_minDREle", eleCaloIso_minDREle, &b_eleCaloIso_minDREle);
   fChain->SetBranchAddress("eleCaloIso_sumPtEle", eleCaloIso_sumPtEle, &b_eleCaloIso_sumPtEle);
   fChain->SetBranchAddress("eleLikelihoodEle", eleLikelihoodEle, &b_eleLikelihoodEle);
   fChain->SetBranchAddress("eleTipEle", eleTipEle, &b_eleTipEle);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("etMet", etMet, &b_etMet);
   fChain->SetBranchAddress("momentumMet", momentumMet, &b_momentumMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("massMet", massMet, &b_massMet);
   fChain->SetBranchAddress("mtMet", mtMet, &b_mtMet);
   fChain->SetBranchAddress("pdgIdMet", pdgIdMet, &b_pdgIdMet);
   fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
   fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
   fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
   fChain->SetBranchAddress("etGenMet", etGenMet, &b_etGenMet);
   fChain->SetBranchAddress("momentumGenMet", momentumGenMet, &b_momentumGenMet);
   fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
   fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
   fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
   fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
   fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
   fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
   fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
   fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
   fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   fChain->SetBranchAddress("massGenMet", massGenMet, &b_massGenMet);
   fChain->SetBranchAddress("mtGenMet", mtGenMet, &b_mtGenMet);
   fChain->SetBranchAddress("pdgIdGenMet", pdgIdGenMet, &b_pdgIdGenMet);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("chargeJet", chargeJet, &b_chargeJet);
   fChain->SetBranchAddress("energyJet", energyJet, &b_energyJet);
   fChain->SetBranchAddress("etJet", etJet, &b_etJet);
   fChain->SetBranchAddress("momentumJet", momentumJet, &b_momentumJet);
   fChain->SetBranchAddress("thetaJet", thetaJet, &b_thetaJet);
   fChain->SetBranchAddress("etaJet", etaJet, &b_etaJet);
   fChain->SetBranchAddress("phiJet", phiJet, &b_phiJet);
   fChain->SetBranchAddress("pxJet", pxJet, &b_pxJet);
   fChain->SetBranchAddress("pyJet", pyJet, &b_pyJet);
   fChain->SetBranchAddress("pzJet", pzJet, &b_pzJet);
   fChain->SetBranchAddress("vertexXJet", vertexXJet, &b_vertexXJet);
   fChain->SetBranchAddress("vertexYJet", vertexYJet, &b_vertexYJet);
   fChain->SetBranchAddress("vertexZJet", vertexZJet, &b_vertexZJet);
   fChain->SetBranchAddress("massJet", massJet, &b_massJet);
   fChain->SetBranchAddress("mtJet", mtJet, &b_mtJet);
   fChain->SetBranchAddress("pdgIdJet", pdgIdJet, &b_pdgIdJet);
   fChain->SetBranchAddress("alphaJet", alphaJet, &b_alphaJet);
   fChain->SetBranchAddress("emFracJet", emFracJet, &b_emFracJet);
   fChain->SetBranchAddress("hadFracJet", hadFracJet, &b_hadFracJet);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("chargeGenJet", chargeGenJet, &b_chargeGenJet);
   fChain->SetBranchAddress("energyGenJet", energyGenJet, &b_energyGenJet);
   fChain->SetBranchAddress("etGenJet", etGenJet, &b_etGenJet);
   fChain->SetBranchAddress("momentumGenJet", momentumGenJet, &b_momentumGenJet);
   fChain->SetBranchAddress("thetaGenJet", thetaGenJet, &b_thetaGenJet);
   fChain->SetBranchAddress("etaGenJet", etaGenJet, &b_etaGenJet);
   fChain->SetBranchAddress("phiGenJet", phiGenJet, &b_phiGenJet);
   fChain->SetBranchAddress("pxGenJet", pxGenJet, &b_pxGenJet);
   fChain->SetBranchAddress("pyGenJet", pyGenJet, &b_pyGenJet);
   fChain->SetBranchAddress("pzGenJet", pzGenJet, &b_pzGenJet);
   fChain->SetBranchAddress("vertexXGenJet", vertexXGenJet, &b_vertexXGenJet);
   fChain->SetBranchAddress("vertexYGenJet", vertexYGenJet, &b_vertexYGenJet);
   fChain->SetBranchAddress("vertexZGenJet", vertexZGenJet, &b_vertexZGenJet);
   fChain->SetBranchAddress("massGenJet", massGenJet, &b_massGenJet);
   fChain->SetBranchAddress("mtGenJet", mtGenJet, &b_mtGenJet);
   fChain->SetBranchAddress("pdgIdGenJet", pdgIdGenJet, &b_pdgIdGenJet);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genPtHat", &genPtHat, &b_genPtHat);
   fChain->SetBranchAddress("genFilterEff", &genFilterEff, &b_genFilterEff);
   fChain->SetBranchAddress("genXsec", &genXsec, &b_genXsec);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);

   Notify();
}

Bool_t HiggsBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HiggsBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HiggsBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HiggsBase_cxx
