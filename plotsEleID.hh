#ifndef plotsEleID_h
#define plotsEleID_h

#include <vector>
#include "Monitor.hh"
#include "HiggsBase.h"

#include <TVector3.h>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"

class plotsEleID : public HiggsBase {
 public:
    plotsEleID(TTree *tree=0);
    virtual ~plotsEleID();
    void Loop();

private:

  // private members
  bool _verbose;

  // sample
  char* category;

  // output files
  TFile *tfilePdfs, *tfileIso, *tfileEff;

  // histos
  // ECALsubdet: 0 = EB; 1 = EE
  // ptBin: 0 = < 15GeV; 1 = > 15 GeV 
  // class: 0 = non-showering; 1 = showering
  // fullclass: 0 = golden; 1 = bigbrem; 2 = narrow; 3 = showering+cracks

  // Hadrons: not splitted
  TH1F *dPhiCaloHad[2][2];
  TH1F *dPhiVtxHad[2][2];
  TH1F *dEtaHad[2][2];
  TH1F *EoPoutHad[2][2];
  TH1F *HoEHad[2][2];  
  TH1F *shapeFisherHad[2][2];  

  // Electrons: not splitted
  // histo[ecalsubdet][ptbin]
  TH1F *dPhiCaloUnsplitEle[2][2];
  TH1F *dPhiVtxUnsplitEle[2][2];
  TH1F *dEtaUnsplitEle[2][2];
  TH1F *EoPoutUnsplitEle[2][2];
  TH1F *HoEUnsplitEle[2][2];  
  TH1F *shapeFisherUnsplitEle[2][2];  
  
  // Electrons class-splitted
  // histo[ecalsubdet][ptbin][class]
  TH1F *dPhiCaloClassEle[2][2][2];
  TH1F *dPhiVtxClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPoutClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *shapeFisherClassEle[2][2][2];

  // Electrons fullclass-splitted
  // histo[ecalsubdet][ptbin][fullclass]
  TH1F *dPhiCaloFullclassEle[2][2][4];
  TH1F *dPhiVtxFullclassEle[2][2][4];
  TH1F *dEtaFullclassEle[2][2][4];
  TH1F *EoPoutFullclassEle[2][2][4];
  TH1F *HoEFullclassEle[2][2][4];
  TH1F *shapeFisherFullclassEle[2][2][4];

  TH1F *H_dRmin_tracker_withVeto,      *E_dRmin_tracker_withVeto;
  TH1F *H_dRmin_tracker_withVeto_zoom, *E_dRmin_tracker_withVeto_zoom;
  TH1F *H_dRmin_tracker_noVeto,        *E_dRmin_tracker_noVeto;
  TH1F *H_dRmin_tracker_noVeto_zoom,   *E_dRmin_tracker_noVeto_zoom;
  TH1F *H_ptFrac_tracker,              *E_ptFrac_tracker;
  TH1F *H_ptFrac_tracker_zoom,         *E_ptFrac_tracker_zoom;
  TH1F *H_dRmin_calo_noVeto,           *E_dRmin_calo_noVeto;
  TH1F *H_ptFrac_calo,                 *E_ptFrac_calo;

  TH1F *H_Reco_eta_wcm,    *E_Reco_eta_wcm;
  TH1F *H_Reco_eta_wgm,    *E_Reco_eta_wgm;
  TH1F *H_Reco_eta_wogm,   *E_Reco_eta_wogm;
  TH1F *H_Gene_eta,        *E_Gene_eta;

  // counters
  int ntotEve;
  int ntot, ntotEB, ntotEE;
  int nGsfClass0[2], nGsfClass1[2], nGsfClass2[2], nGsfClass3[2];
  int trackIso, caloIso;

};
#endif
