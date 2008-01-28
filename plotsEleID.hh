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
  bool splitelectrons;
  bool isGt15GeV;

  // output files
  TFile *tfileEB, *tfileEE, *tfileIso, *tfileEff;

  // histos
  TH1F *dPhiCalo[2];
  TH1F *dPhiVtx[2];
  TH1F *dEta[2];
  TH1F *EoPout[2];
  TH1F *HoE[2];  
  TH1F *shapeFisher[2];  
  TH1F *dPhiCaloUnsplit[2];
  TH1F *dPhiVtxUnsplit[2];
  TH1F *dEtaUnsplit[2];
  TH1F *EoPoutUnsplit[2];
  TH1F *HoEUnsplit[2];  
  TH1F *shapeFisherUnsplit[2];  
  TH1F *dPhiCalo0[2],    *dPhiCalo1[2],     *dPhiCalo2[2],    *dPhiCalo3[2];
  TH1F *dPhiVtx0[2],     *dPhiVtx1[2],      *dPhiVtx2[2],     *dPhiVtx3[2];
  TH1F *dEta0[2],        *dEta1[2],         *dEta2[2],        *dEta3[2];
  TH1F *EoPout0[2],      *EoPout1[2],       *EoPout2[2],      *EoPout3[2];
  TH1F *HoE0[2],         *HoE1[2],          *HoE2[2],         *HoE3[2];
  TH1F *shapeFisher0[2], *shapeFisher1[2],  *shapeFisher2[2], *shapeFisher3[2];

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
