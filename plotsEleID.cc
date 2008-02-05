#include <iostream>
#include <string>

#include <TTree.h>
#include <TMath.h>

#include "Counters.hh"
#include "Selection.hh"
#include "SprDataFiller.hh"
#include "plotsEleID.hh"

plotsEleID::plotsEleID(TTree *tree) 
  : HiggsBase(tree) {

  int nbins = 100;

  float dPhiCaloMin = -0.3;
  float dPhiCaloMax =  0.3;
  float dPhiVtxMin  = -0.1; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiVtxMax  =  0.1;
  float dEtaMin     = -0.05;
  float dEtaMax     =  0.05;
  float EoPoutMin   =  0.0;
  float EoPoutMax   =  50.0;
  float HoEMin      = -0.2; // pixelMatchGsfElectron pre-selection: |H/E| < 0.2
  float HoEMax      =  0.2;
  float fisherMin   = -15.0;
  float fisherMax   =  15.0;
  float sigmaEtaEtaMin = 0.0;
  float sigmaEtaEtaMax = 0.1;
  float sigmaEtaPhiMin = 0.0;
  float sigmaEtaPhiMax = 0.1;
  float sigmaPhiPhiMin = 0.0;
  float sigmaPhiPhiMax = 0.1;
  float s1s9Min  = 0.0;
  float s1s9Max  = 1.0;
  float s9s25Min = 0.0;
  float s9s25Max = 1.0;
  float LATMin   = 0.0;
  float LATMax   = 1.0;
  float etaLATMin = 0.0;
  float etaLATMax = 1.0;
  float phiLATMin = 0.0;
  float phiLATMax = 1.0;
  float a20Min    = 0.0;
  float a20Max    = 1.0;
  float a42Min    = 0.0;
  float a42Max    = 0.5;

  // sample region
  //category = "hadrons";  
  category = "electrons";  

  // max number of shape variables to be inserted in the Fisher training
  nStatPatternRecognitionVars = 10;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: < 15 GeV
    // iptbin = 1: > 15 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];
      
      sprintf(histo,"dPhiCalo_hadrons_%d_%d",iecal,iptbin);
      dPhiCaloHad[iecal][iptbin]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
      sprintf(histo,"dPhiVtx_hadrons_%d_%d",iecal,iptbin);
      dPhiVtxHad[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
      sprintf(histo,"dEta_hadrons_%d_%d",iecal,iptbin);
      dEtaHad[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPout_hadrons_%d_%d",iecal,iptbin);
      EoPoutHad[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
      sprintf(histo,"HoE_hadrons_%d_%d",iecal,iptbin);
      HoEHad[iecal][iptbin]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      sprintf(histo,"shapeFisher_hadrons_%d_%d",iecal,iptbin);
      shapeFisherHad[iecal][iptbin] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
      sprintf(histo,"sigmaEtaEta_hadrons_%d_%d",iecal,iptbin);
      sigmaEtaEtaHad[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
      sprintf(histo,"sigmaEtaPhi_hadrons_%d_%d",iecal,iptbin);
      sigmaEtaPhiHad[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
      sprintf(histo,"sigmaPhiPhi_hadrons_%d_%d",iecal,iptbin);
      sigmaPhiPhiHad[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
      sprintf(histo,"s1s9_hadrons_%d_%d",iecal,iptbin);
      s1s9Had[iecal][iptbin] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
      sprintf(histo,"s9s25_hadrons_%d_%d",iecal,iptbin);
      s9s25Had[iecal][iptbin] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
      sprintf(histo,"LAT_hadrons_%d_%d",iecal,iptbin);
      LATHad[iecal][iptbin] = new TH1F(histo, histo, nbins, LATMin, LATMax);
      sprintf(histo,"etaLAT_hadrons_%d_%d",iecal,iptbin);
      etaLATHad[iecal][iptbin] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
      sprintf(histo,"phiLAT_hadrons_%d_%d",iecal,iptbin);
      phiLATHad[iecal][iptbin] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
      sprintf(histo,"a20_hadrons_%d_%d",iecal,iptbin);
      a20Had[iecal][iptbin] = new TH1F(histo, histo, nbins, a20Min, a20Max);
      sprintf(histo,"a42_hadrons_%d_%d",iecal,iptbin);
      a42Had[iecal][iptbin] = new TH1F(histo, histo, nbins, a42Min, a42Max);

      sprintf(histo,"dPhiCaloUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiCaloUnsplitEle[iecal][iptbin]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
      sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiVtxUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
      sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
      EoPoutUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
      sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      sprintf(histo,"shapeFisherUnsplit_electrons_%d_%d",iecal,iptbin);
      shapeFisherUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
      sprintf(histo,"sigmaEtaEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
      sprintf(histo,"sigmaEtaPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
      sprintf(histo,"sigmaPhiPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaPhiPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
      sprintf(histo,"s1s9Unsplit_electrons_%d_%d",iecal,iptbin);
      s1s9UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
      sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
      s9s25UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
      sprintf(histo,"LATUnsplit_electrons_%d_%d",iecal,iptbin);
      LATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, LATMin, LATMax);
      sprintf(histo,"etaLATUnsplit_electrons_%d_%d",iecal,iptbin);
      etaLATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
      sprintf(histo,"phiLATUnsplit_electrons_%d_%d",iecal,iptbin);
      phiLATUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
      sprintf(histo,"a20Unsplit_electrons_%d_%d",iecal,iptbin);
      a20UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, a20Min, a20Max);
      sprintf(histo,"a42Unsplit_electrons_%d_%d",iecal,iptbin);
      a42UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, a42Min, a42Max);

      // iclass = 0: non-showering
      // iclass = 1: showering
      for(int iclass=0; iclass<2; iclass++) {
      
	sprintf(histo,"dPhiCaloClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiCaloClassEle[iecal][iptbin][iclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiVtxClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPoutClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"shapeFisherClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	shapeFisherClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
	sprintf(histo,"sigmaEtaEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaEtaEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
	sprintf(histo,"sigmaEtaPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaEtaPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
	sprintf(histo,"sigmaPhiPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaPhiPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
	sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s1s9ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s9s25ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
	sprintf(histo,"LATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	LATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, LATMin, LATMax);
	sprintf(histo,"etaLATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	etaLATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
	sprintf(histo,"phiLATClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	phiLATClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
	sprintf(histo,"a20Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	a20ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, a20Min, a20Max);
	sprintf(histo,"a42Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	a42ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, a42Min, a42Max);

      }

      // iclass = 0: golden
      // iclass = 1: bigbrem
      // iclass = 2: narrow
      // iclass = 3: showering + cracks

      for(int ifullclass=0; ifullclass<4; ifullclass++) {
      
	sprintf(histo,"dPhiCaloFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dEtaFullclassEle[iecal][iptbin][ifullclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	EoPoutFullclassEle[iecal][iptbin][ifullclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	HoEFullclassEle[iecal][iptbin][ifullclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"shapeFisherFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	shapeFisherFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, fisherMin, fisherMax);
	sprintf(histo,"sigmaEtaEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaEtaEtaFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaEtaEtaMin, sigmaEtaEtaMax);
	sprintf(histo,"sigmaEtaPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaEtaPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaEtaPhiMin, sigmaEtaPhiMax);
	sprintf(histo,"sigmaPhiPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaPhiPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaPhiPhiMin, sigmaPhiPhiMax);
	sprintf(histo,"s1s9Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s1s9FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s9s25FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
	sprintf(histo,"LATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	LATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, LATMin, LATMax);
	sprintf(histo,"etaLATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	etaLATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, etaLATMin, etaLATMax);
	sprintf(histo,"phiLATFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	phiLATFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, phiLATMin, phiLATMax);
	sprintf(histo,"a20Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	a20FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, a20Min, a20Max);
	sprintf(histo,"a42Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	a42FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, a42Min, a42Max);

      }

    }

  }

  // booking histos isolation
  H_dRmin_tracker_withVeto      = new TH1F("H_dRmin_tracker_withVeto",      "min #Delta R", 100, 0., 0.5);
  H_dRmin_tracker_noVeto        = new TH1F("H_dRmin_tracker_noVeto",        "min #Delta R", 100, 0., 0.5);
  H_dRmin_tracker_withVeto_zoom = new TH1F("H_dRmin_tracker_withVeto_zoom", "min #Delta R", 100, 0., 0.1);
  H_dRmin_tracker_noVeto_zoom   = new TH1F("H_dRmin_tracker_noVeto_zoom",   "min #Delta R", 100, 0., 0.1);
  H_ptFrac_tracker              = new TH1F("H_ptFrac_tracker",              "Sum tracks p_{T} / ele p_{T}", 100, 0., 0.5);
  H_ptFrac_tracker_zoom         = new TH1F("H_ptFrac_tracker_zoom",         "Sum tracks p_{T} / ele p_{T}", 100, 0., 0.05);
  H_dRmin_calo_noVeto           = new TH1F("H_dRmin_calo_noVeto",           "min #Delta R", 100, 0., 0.5);
  H_ptFrac_calo                 = new TH1F("H_ptFrac_calo",                 "Sum hcal rechits E_{T} / ele E_{T}", 100, 0., 0.5);
  E_dRmin_tracker_withVeto      = new TH1F("E_dRmin_tracker_withVeto",      "min #Delta R", 100, 0., 0.5);
  E_dRmin_tracker_noVeto        = new TH1F("E_dRmin_tracker_noVeto",        "min #Delta R", 100, 0., 0.5);
  E_dRmin_tracker_withVeto_zoom = new TH1F("E_dRmin_tracker_withVeto_zoom", "min #Delta R", 100, 0., 0.1);
  E_dRmin_tracker_noVeto_zoom   = new TH1F("E_dRmin_tracker_noVeto_zoom",   "min #Delta R", 100, 0., 0.1);
  E_ptFrac_tracker              = new TH1F("E_ptFrac_tracker",              "Sum tracks p_{T} / ele p_{T}", 100, 0., 0.5);
  E_ptFrac_tracker_zoom         = new TH1F("E_ptFrac_tracker_zoom",         "Sum tracks p_{T} / ele p_{T}", 100, 0., 0.05);
  E_dRmin_calo_noVeto           = new TH1F("E_dRmin_calo_noVeto",           "min #Delta R", 100, 0., 0.5);
  E_ptFrac_calo                 = new TH1F("E_ptFrac_calo",                 "Sum hcal rechits E_{T} / ele E_{T}", 100, 0., 0.5);

  // booking histos efficiency
  H_Reco_eta_wcm  = new TH1F("H_Reco_eta_wcm",  "reconstructed #eta - with charge matching",    25, -2.5,2.5);
  H_Reco_eta_wgm  = new TH1F("H_Reco_eta_wgm",  "reconstructed #eta - with geom matching only", 25, -2.5,2.5);
  H_Reco_eta_wogm = new TH1F("H_Reco_eta_wogm", "reconstructed #eta - without matching",        25, -2.5,2.5);
  //  H_Reco_eta_tiso = new TH1F("H_Reco_eta_tiso", "reconstructed #eta - without matching, with tracker isolation", 25, -2.5,2.5);
  H_Gene_eta      = new TH1F("H_Gene_eta",      "generated #eta",                               25, -2.5,2.5);
  E_Reco_eta_wcm  = new TH1F("E_Reco_eta_wcm",  "reconstructed #eta - with charge matching",    25, -2.5,2.5);
  E_Reco_eta_wgm  = new TH1F("E_Reco_eta_wgm",  "reconstructed #eta - with geom matching only", 25, -2.5,2.5);
  E_Reco_eta_wogm = new TH1F("E_Reco_eta_wogm", "reconstructed #eta - without matching",        25, -2.5,2.5);
  //  E_Reco_eta_tiso = new TH1F("H_Reco_eta_tiso", "reconstructed #eta - without matching, with tracker isolation", 25, -2.5,2.5);
  E_Gene_eta      = new TH1F("E_Gene_eta",      "generated #eta",                               25, -2.5,2.5);

  // counters
  ntotEve  = 0;
  ntot     = 0;
  ntotEB   = 0;
  ntotEE   = 0;
  caloIso  = 0;
  trackIso = 0;
  for(int ii=0; ii<2; ii++){
    nGsfClass0[ii] = 0;
    nGsfClass1[ii] = 0;
    nGsfClass2[ii] = 0;
    nGsfClass3[ii] = 0;
  }


  // output files
  tfilePdfs  = new TFile("pdfs.root","RECREATE");
  tfileIso = new TFile("isolation.root", "RECREATE");
  tfileEff = new TFile("efficiency.root","RECREATE");
}
 
plotsEleID::~plotsEleID() {

  // fraction of electons x class
  int ntot = ntotEB + ntotEE;
  std::cout << "Processed ntot = " << ntot << " electrons" << std::endl;
  std::cout << std::endl;
  std::cout << "In barrel: "  << ntotEB << " electrons" << std::endl;
  std::cout << "GsfClass0_EB = " << float(nGsfClass0[0])/float(ntotEB)*100 << " %" << std::endl;
  std::cout << "GsfClass1_EB = " << float(nGsfClass1[0])/float(ntotEB)*100 << " %" << std::endl;
  std::cout << "GsfClass2_EB = " << float(nGsfClass2[0])/float(ntotEB)*100 << " %" << std::endl;
  std::cout << "GsfClass3_EB = " << float(nGsfClass3[0])/float(ntotEB)*100 << " %" << std::endl;
  std::cout << std::endl;
  std::cout << "In barrel: "  << ntotEE << " electrons" << std::endl;
  std::cout << "GsfClass0_EE = " << float(nGsfClass0[1])/float(ntotEE)*100 << " %" << std::endl;
  std::cout << "GsfClass1_EE = " << float(nGsfClass1[1])/float(ntotEE)*100 << " %" << std::endl;
  std::cout << "GsfClass2_EE = " << float(nGsfClass2[1])/float(ntotEE)*100 << " %" << std::endl;
  std::cout << "GsfClass3_EE = " << float(nGsfClass3[1])/float(ntotEE)*100 << " %" << std::endl;
  std::cout << std::endl;
  std::cout << "Iso cuts efficiency: " << std::endl;
  std::cout << "tracker: " << trackIso << " events passing the cut = " << float(trackIso)/float(ntot)*100 << " %" << std::endl;
  std::cout << "calo: "    << caloIso  << " events passing the cut = " << float(caloIso)/float(ntot)*100  << " %" << std::endl;
  std::cout << std::endl;

  // saving histos

  // pdf electron ID
  tfilePdfs->cd();
    
  for (int iecal=0; iecal<2; iecal++) {
    for(int iptbin=0; iptbin<2; iptbin++) {
      if(strcmp(category,"hadrons")==0) {
	dPhiCaloHad[iecal][iptbin]    -> Write();
	dPhiVtxHad[iecal][iptbin]     -> Write();
	dEtaHad[iecal][iptbin]        -> Write();
	EoPoutHad[iecal][iptbin]      -> Write();
	HoEHad[iecal][iptbin]         -> Write();
	shapeFisherHad[iecal][iptbin] -> Write();
	sigmaEtaEtaHad[iecal][iptbin] -> Write();
	sigmaEtaPhiHad[iecal][iptbin] -> Write();
	sigmaPhiPhiHad[iecal][iptbin] -> Write();
	s1s9Had[iecal][iptbin]        -> Write();
	s9s25Had[iecal][iptbin]       -> Write();
	LATHad[iecal][iptbin]         -> Write();
	etaLATHad[iecal][iptbin]      -> Write();
	phiLATHad[iecal][iptbin]      -> Write();
	a20Had[iecal][iptbin]         -> Write();
	a42Had[iecal][iptbin]         -> Write();
      }
      else if(strcmp(category,"electrons")==0) {

	dPhiCaloUnsplitEle[iecal][iptbin]    -> Write();
	dPhiVtxUnsplitEle[iecal][iptbin]     -> Write();
	dEtaUnsplitEle[iecal][iptbin]        -> Write();
	EoPoutUnsplitEle[iecal][iptbin]      -> Write();
	HoEUnsplitEle[iecal][iptbin]         -> Write();
	shapeFisherUnsplitEle[iecal][iptbin] -> Write();
	sigmaEtaEtaUnsplitEle[iecal][iptbin] -> Write();
	sigmaEtaPhiUnsplitEle[iecal][iptbin] -> Write();
	sigmaPhiPhiUnsplitEle[iecal][iptbin] -> Write();
	s1s9UnsplitEle[iecal][iptbin]        -> Write();
	s9s25UnsplitEle[iecal][iptbin]       -> Write();
	LATUnsplitEle[iecal][iptbin]         -> Write();
	etaLATUnsplitEle[iecal][iptbin]      -> Write();
	phiLATUnsplitEle[iecal][iptbin]      -> Write();
	a20UnsplitEle[iecal][iptbin]         -> Write();
	a42UnsplitEle[iecal][iptbin]         -> Write();

	for(int iclass=0; iclass<2; iclass++) {

	  dPhiCaloClassEle[iecal][iptbin][iclass]    -> Write();
	  dPhiVtxClassEle[iecal][iptbin][iclass]     -> Write();
	  dEtaClassEle[iecal][iptbin][iclass]        -> Write();
	  EoPoutClassEle[iecal][iptbin][iclass]      -> Write();
	  HoEClassEle[iecal][iptbin][iclass]         -> Write();
	  shapeFisherClassEle[iecal][iptbin][iclass] -> Write();
	  sigmaEtaEtaClassEle[iecal][iptbin][iclass] -> Write();
	  sigmaEtaPhiClassEle[iecal][iptbin][iclass] -> Write();
	  sigmaPhiPhiClassEle[iecal][iptbin][iclass] -> Write();
	  s1s9ClassEle[iecal][iptbin][iclass]        -> Write();
	  s9s25ClassEle[iecal][iptbin][iclass]       -> Write();
	  LATClassEle[iecal][iptbin][iclass]         -> Write();
	  etaLATClassEle[iecal][iptbin][iclass]      -> Write();
	  phiLATClassEle[iecal][iptbin][iclass]      -> Write();
	  a20ClassEle[iecal][iptbin][iclass]         -> Write();
	  a42ClassEle[iecal][iptbin][iclass]         -> Write();

	}

	for(int ifullclass=0; ifullclass<4; ifullclass++) {

	  dPhiCaloFullclassEle[iecal][iptbin][ifullclass]    -> Write();
	  dPhiVtxFullclassEle[iecal][iptbin][ifullclass]     -> Write();
	  dEtaFullclassEle[iecal][iptbin][ifullclass]        -> Write();
	  EoPoutFullclassEle[iecal][iptbin][ifullclass]      -> Write();
	  HoEFullclassEle[iecal][iptbin][ifullclass]         -> Write();
	  shapeFisherFullclassEle[iecal][iptbin][ifullclass] -> Write();
	  sigmaEtaEtaFullclassEle[iecal][iptbin][ifullclass] -> Write();
	  sigmaEtaPhiFullclassEle[iecal][iptbin][ifullclass] -> Write();
	  sigmaPhiPhiFullclassEle[iecal][iptbin][ifullclass] -> Write();
	  s1s9FullclassEle[iecal][iptbin][ifullclass]        -> Write();
	  s9s25FullclassEle[iecal][iptbin][ifullclass]       -> Write();
	  LATFullclassEle[iecal][iptbin][ifullclass]         -> Write();
	  etaLATFullclassEle[iecal][iptbin][ifullclass]      -> Write();
	  phiLATFullclassEle[iecal][iptbin][ifullclass]      -> Write();
	  a20FullclassEle[iecal][iptbin][ifullclass]         -> Write();
	  a42FullclassEle[iecal][iptbin][ifullclass]         -> Write();

	}

      }
      else {
	std::cout << "ERROR! category " << category << " not recognized!" << std::endl;
      }
	
    }

  }
	

  // isolation
  tfileIso->cd();
  H_dRmin_tracker_withVeto      -> Write();
  H_dRmin_tracker_noVeto        -> Write();
  H_ptFrac_tracker              -> Write();
  H_dRmin_tracker_withVeto_zoom -> Write();
  H_dRmin_tracker_noVeto_zoom   -> Write();
  H_ptFrac_tracker_zoom         -> Write();
  H_dRmin_calo_noVeto           -> Write();
  H_ptFrac_calo                 -> Write();

  
  // efficiency
  tfileEff->cd();
  H_Reco_eta_wcm     ->Write();
  H_Reco_eta_wgm     ->Write();
  H_Reco_eta_wogm    ->Write();
  //    H_Reco_eta_tiso    ->Write();
  H_Gene_eta         ->Write();

  tfilePdfs ->Close();
  tfileIso->Close();
  
  // deleting
  for (int iecal=0; iecal<2; iecal++) {
    for(int iptbin=0; iptbin<2; iptbin++) {

      if(strcmp(category,"hadrons")==0) {

	delete dPhiCaloHad[iecal][iptbin];
	delete dPhiVtxHad[iecal][iptbin];
	delete dEtaHad[iecal][iptbin];
	delete EoPoutHad[iecal][iptbin];
	delete HoEHad[iecal][iptbin];
	delete shapeFisherHad[iecal][iptbin];
	delete sigmaEtaEtaHad[iecal][iptbin];
	delete sigmaEtaPhiHad[iecal][iptbin];
	delete sigmaPhiPhiHad[iecal][iptbin];
	delete s1s9Had[iecal][iptbin];
	delete s9s25Had[iecal][iptbin];
	delete LATHad[iecal][iptbin];
	delete etaLATHad[iecal][iptbin];
	delete phiLATHad[iecal][iptbin];
	delete a20Had[iecal][iptbin];
	delete a42Had[iecal][iptbin];

      }
      else if(strcmp(category,"electrons")==0) {

	delete dPhiCaloUnsplitEle[iecal][iptbin];
	delete dPhiVtxUnsplitEle[iecal][iptbin];
	delete dEtaUnsplitEle[iecal][iptbin];
	delete EoPoutUnsplitEle[iecal][iptbin];
	delete HoEUnsplitEle[iecal][iptbin];
	delete shapeFisherUnsplitEle[iecal][iptbin];
	delete sigmaEtaEtaUnsplitEle[iecal][iptbin];
	delete sigmaEtaPhiUnsplitEle[iecal][iptbin];
	delete sigmaPhiPhiUnsplitEle[iecal][iptbin];
	delete s1s9UnsplitEle[iecal][iptbin];
	delete s9s25UnsplitEle[iecal][iptbin];
	delete LATUnsplitEle[iecal][iptbin];
	delete etaLATUnsplitEle[iecal][iptbin];
	delete phiLATUnsplitEle[iecal][iptbin];
	delete a20UnsplitEle[iecal][iptbin];
	delete a42UnsplitEle[iecal][iptbin];
	
	for(int iclass=0; iclass<2; iclass++) {

	  delete dPhiCaloClassEle[iecal][iptbin][iclass];
	  delete dPhiVtxClassEle[iecal][iptbin][iclass];
	  delete dEtaClassEle[iecal][iptbin][iclass];
	  delete EoPoutClassEle[iecal][iptbin][iclass];
	  delete HoEClassEle[iecal][iptbin][iclass];
	  delete shapeFisherClassEle[iecal][iptbin][iclass];
	  delete sigmaEtaEtaClassEle[iecal][iptbin][iclass];
	  delete sigmaEtaPhiClassEle[iecal][iptbin][iclass];
	  delete sigmaPhiPhiClassEle[iecal][iptbin][iclass];
	  delete s1s9ClassEle[iecal][iptbin][iclass];
	  delete s9s25ClassEle[iecal][iptbin][iclass];
	  delete LATClassEle[iecal][iptbin][iclass];
	  delete etaLATClassEle[iecal][iptbin][iclass];
	  delete phiLATClassEle[iecal][iptbin][iclass];
	  delete a20ClassEle[iecal][iptbin][iclass];
	  delete a42ClassEle[iecal][iptbin][iclass];

	}

	for(int ifullclass=0; ifullclass<4; ifullclass++) {

	  delete dPhiCaloFullclassEle[iecal][iptbin][ifullclass];
	  delete dPhiVtxFullclassEle[iecal][iptbin][ifullclass];
	  delete dEtaFullclassEle[iecal][iptbin][ifullclass];
	  delete EoPoutFullclassEle[iecal][iptbin][ifullclass];
	  delete HoEFullclassEle[iecal][iptbin][ifullclass];
	  delete shapeFisherFullclassEle[iecal][iptbin][ifullclass];
	  delete sigmaEtaEtaFullclassEle[iecal][iptbin][ifullclass];
	  delete sigmaEtaPhiFullclassEle[iecal][iptbin][ifullclass];
	  delete sigmaPhiPhiFullclassEle[iecal][iptbin][ifullclass];
	  delete s1s9FullclassEle[iecal][iptbin][ifullclass];
	  delete s9s25FullclassEle[iecal][iptbin][ifullclass];
	  delete LATFullclassEle[iecal][iptbin][ifullclass];
	  delete etaLATFullclassEle[iecal][iptbin][ifullclass];
	  delete phiLATFullclassEle[iecal][iptbin][ifullclass];
	  delete a20FullclassEle[iecal][iptbin][ifullclass];
	  delete a42FullclassEle[iecal][iptbin][ifullclass];

	}
	
      }

    }

  }

  delete E_dRmin_tracker_withVeto,      H_dRmin_tracker_withVeto;
  delete E_dRmin_tracker_withVeto_zoom, H_dRmin_tracker_withVeto_zoom;
  delete E_dRmin_tracker_noVeto,        H_dRmin_tracker_noVeto;
  delete E_dRmin_tracker_noVeto_zoom,   H_dRmin_tracker_noVeto_zoom;
  delete E_ptFrac_tracker,              H_ptFrac_tracker;
  delete E_ptFrac_tracker_zoom,         H_ptFrac_tracker_zoom;
  delete E_dRmin_calo_noVeto,           H_dRmin_calo_noVeto;
  delete E_ptFrac_calo,                 H_ptFrac_calo;

  delete E_Reco_eta_wcm,      H_Reco_eta_wcm;
  delete E_Reco_eta_wgm,      H_Reco_eta_wgm;
  delete E_Reco_eta_wogm,     H_Reco_eta_wogm;  
  //  delete E_Reco_eta_tiso,     H_Reco_eta_tiso;  
  delete E_Gene_eta,          H_Gene_eta;
}


void plotsEleID::Loop() {

  _verbose=true;
  if(fChain == 0) return;
    
  Long64_t nentries = fChain->GetEntries();

  // set up the filler for StatPatternRecognition
  std::string sprName;
  if ( strcmp(category,"hadrons")==0 ) sprName = "clusterShapeVars-hadrons";
  else if ( strcmp(category,"electrons")==0 ) sprName = "clusterShapeVars-electrons";
  float sigmaEtaEta, sigmaEtaPhi, sigmaPhiPhi, s1s9, s9s25, LAT, etaLAT, phiLAT, a20, a42;
  SprDataFiller clusterShapeFiller;
  clusterShapeFiller.setName( sprName.c_str() );
  clusterShapeFiller.setEntries(nentries);
  clusterShapeFiller.add("sigmaEtaEta",&sigmaEtaEta);
  clusterShapeFiller.add("sigmaEtaPhi",&sigmaEtaPhi);
  clusterShapeFiller.add("sigmaPhiPhi",&sigmaPhiPhi);
  clusterShapeFiller.add("s1s9",&s1s9);
  clusterShapeFiller.add("s9s25",&s9s25);
  clusterShapeFiller.add("a20",&a20);
  clusterShapeFiller.add("a42",&a42);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%1000==0) std::cout << "Processing " << jentry << std::endl;

    for(int iele=0;iele<nEle;iele++) {

      int iecal = -1;
      float fisher;
      if ( fabs(etaEle[iele]) <= 1.476 ) { // barrel
	iecal = 0; 
	fisher = 42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele];
      }

      if ( fabs(etaEle[iele]) > 1.476  ) { // endcap
	iecal = 1; 
	fisher = 27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele];
      }

      int iptbin = -1;
      if ( etEle[iele] < 15.0 ) iptbin = 0;
      else iptbin = 1;

      int iclass = -1;
      int ifullclass = -1;
      if ( eleClassEle[iele] == 0 || eleClassEle[iele] == 100 ) { // golden
	iclass = 0;
	ifullclass = 0;
      }
      else if ( eleClassEle[iele] == 10 || eleClassEle[iele] == 110 ) { // bigbrem
	iclass = 0;
	ifullclass = 1;
      }
      else if ( eleClassEle[iele] == 20 || eleClassEle[iele] == 120 ) { // narrow
	iclass = 0;
	ifullclass = 2;
      }
      else if ( (eleClassEle[iele] >= 30 && eleClassEle[iele] <= 40) ||
		(eleClassEle[iele] >= 130 && eleClassEle[iele] <= 140) ) { // showering + cracks
	iclass = 1;
	ifullclass = 3;
      }

      // efficiency
      int recCharge = chargeEle[iele];
      double recEta = etaEle[iele];
      double recPhi = phiEle[iele];
      double recEne = energyEle[iele];
      bool isolated = true;
      if (eleTrackerIso_sumPtEle[iele]<0.05) { isolated = true; }
      if (eleTrackerIso_sumPtEle[iele]>=0.05){ isolated = false; }

      float dr_min_wc       = 10000.;
      float dr_min_woc      = 10000.;
      float dr_min_wiso     = 10000.;
      float EoPm1_best_wc   = 10000.;
      float EoPm1_best_woc  = 10000.;
      float EoPm1_best_wiso = 10000.;
      float thisMcEta_wc    = 0.;
      float thisMcEta_woc   = 0.;
      float thisMcEta_wiso  = 0.;
      for(int imc=0;imc<nMc;imc++) {
	
	double trueEta;
	if (idMc[imc] == 11 || idMc[imc] == -11){
	  
	  trueEta        = etaMc[imc];
	  double truePhi = phiMc[imc];
	  double trueEne = pMc[imc];
	  int trueCharge;
	  if (idMc[imc] == -11){ trueCharge = -1; }
	  if (idMc[imc] == +11){ trueCharge = +1; }
	  
	  double pi = 3.1415927;
	  double deltaEta = recEta - trueEta;
	  double deltaPhi = recPhi - truePhi;
	  if(deltaPhi > pi)  deltaPhi -= 2.*pi;
	  if(deltaPhi < -pi) deltaPhi += 2.*pi;
	  double deltaR  = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
	  double EoEtrue = recEne/trueEne;
	  
	  // without charge matching
	  if (deltaR < 0.5 ){
	    if (fabs(EoEtrue-1) < EoPm1_best_woc){
	      dr_min_woc     = deltaR;
	      EoPm1_best_woc = fabs(EoEtrue-1);
	      thisMcEta_woc  = trueEta;
	    }}

	  // with charge matching
	  if (deltaR < 0.5){
	    if (trueCharge*recCharge > 0){
	      if (fabs(EoEtrue-1) < EoPm1_best_wc){
		dr_min_wc     = deltaR;
		EoPm1_best_wc = fabs(EoEtrue-1);
		thisMcEta_wc  = trueEta;
	      }}}

	  // isolated (in tracker) ele only
	  if (deltaR < 0.5){
	    if (isolated){
	      if (fabs(EoEtrue-1) < EoPm1_best_wc){
		//		dr_min_iso     = deltaR;
		//		EoPm1_best_iso = fabs(EoEtrue-1);
		//		thisMcEta_iso  = trueEta;
	      }}}
	} // e+ or e-
      } // mc particles
      
      
      if(strcmp(category,"hadrons")==0) {

	ntot++;

	// eleID
	dPhiCaloHad    [iecal][iptbin] -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxHad     [iecal][iptbin] -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaHad        [iecal][iptbin] -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutHad      [iecal][iptbin] -> Fill ( eleCorrEoPoutEle[iele] );
	HoEHad         [iecal][iptbin] -> Fill ( eleHoEEle[iele] );
	shapeFisherHad [iecal][iptbin] -> Fill ( fisher );
	sigmaEtaEtaHad [iecal][iptbin] -> Fill ( sqrt(covEtaEtaEle[iele]) );
	sigmaEtaPhiHad [iecal][iptbin] -> Fill ( sqrt(fabs(covEtaPhiEle[iele])) );
	sigmaPhiPhiHad [iecal][iptbin] -> Fill ( sqrt(covPhiPhiEle[iele]) );
	s1s9Had        [iecal][iptbin] -> Fill ( s1s9Ele[iele] );
	s9s25Had       [iecal][iptbin] -> Fill ( s9s25Ele[iele] );
	LATHad         [iecal][iptbin] -> Fill ( latEle[iele] );
	etaLATHad      [iecal][iptbin] -> Fill ( etaLatEle[iele] );
	phiLATHad      [iecal][iptbin] -> Fill ( phiLatEle[iele] );
	a20Had         [iecal][iptbin] -> Fill ( a20Ele[iele] );
	a42Had         [iecal][iptbin] -> Fill ( a20Ele[iele] );

	// isolation
	if(eleTrackerIso_sumPtEle[iele]>0. && eleTrackerIso_sumPtEle[iele]<0.5)           { H_ptFrac_tracker         -> Fill(eleTrackerIso_sumPtEle[iele]);}
	if(eleTrackerIso_sumPtEle[iele]>0. && eleTrackerIso_sumPtEle[iele]<0.05)          { H_ptFrac_tracker_zoom    -> Fill(eleTrackerIso_sumPtEle[iele]);}
	if(eleTrackerIso_minDR_vetoEle[iele]> 0. && eleTrackerIso_minDR_vetoEle[iele]<0.5){ H_dRmin_tracker_withVeto -> Fill(eleTrackerIso_minDR_vetoEle[iele]);}
	if(eleTrackerIso_minDR_vetoEle[iele]> 0. && eleTrackerIso_minDR_vetoEle[iele]<0.1){ H_dRmin_tracker_withVeto_zoom -> Fill(eleTrackerIso_minDR_vetoEle[iele]);}
	if(eleTrackerIso_minDREle[iele]> 0. && eleTrackerIso_minDREle[iele]<0.5)          { H_dRmin_tracker_noVeto   -> Fill(eleTrackerIso_minDREle[iele]);}
	if(eleTrackerIso_minDREle[iele]> 0. && eleTrackerIso_minDREle[iele]<0.1)          { H_dRmin_tracker_noVeto_zoom -> Fill(eleTrackerIso_minDREle[iele]); }
	if(eleCaloIso_minDREle[iele]>0 && eleCaloIso_minDREle[iele]<0.5)                  { H_dRmin_calo_noVeto         -> Fill(eleCaloIso_minDREle[iele]);}
	if(eleCaloIso_sumPtEle[iele]>0 && eleCaloIso_sumPtEle[iele]<0.5)                  { H_ptFrac_calo               -> Fill(eleCaloIso_sumPtEle[iele]);}
	
	// computing the efficiency of iso cuts:
	if (eleTrackerIso_sumPtEle[iele]<0.05){ trackIso++; }
	if (eleCaloIso_sumPtEle[iele]<0.05)   { caloIso++; }

	// efficiency: all reco
	H_Reco_eta_wogm ->Fill(etaEle[iele]);
	
	// efficiency: all gene
	for(int imc=0;imc<nMc;imc++) {
	  if (idMc[imc] == 11 || idMc[imc] == -11){H_Gene_eta->Fill(etaMc[imc]); }} 
	
	// efficiency: matching
	if (dr_min_wc  < 0.05){ H_Reco_eta_wcm -> Fill(recEta); }
	if (dr_min_woc < 0.05){ H_Reco_eta_wgm -> Fill(recEta); }
      }
      
      else if(strcmp(category,"electrons")==0) {
	
	// isolation
	if(eleTrackerIso_sumPtEle[iele]>0. && eleTrackerIso_sumPtEle[iele]<0.5)           { E_ptFrac_tracker         -> Fill(eleTrackerIso_sumPtEle[iele]);}
	if(eleTrackerIso_sumPtEle[iele]>0. && eleTrackerIso_sumPtEle[iele]<0.05)          { E_ptFrac_tracker_zoom    -> Fill(eleTrackerIso_sumPtEle[iele]);}
	if(eleTrackerIso_minDR_vetoEle[iele]> 0. && eleTrackerIso_minDR_vetoEle[iele]<0.5){ E_dRmin_tracker_withVeto -> Fill(eleTrackerIso_minDR_vetoEle[iele]);}
	if(eleTrackerIso_minDR_vetoEle[iele]> 0. && eleTrackerIso_minDR_vetoEle[iele]<0.1){ E_dRmin_tracker_withVeto_zoom -> Fill(eleTrackerIso_minDR_vetoEle[iele]);}
	if(eleTrackerIso_minDREle[iele]> 0. && eleTrackerIso_minDREle[iele]<0.5)          { E_dRmin_tracker_noVeto   -> Fill(eleTrackerIso_minDREle[iele]);}
	if(eleTrackerIso_minDREle[iele]> 0. && eleTrackerIso_minDREle[iele]<0.1)          { E_dRmin_tracker_noVeto_zoom -> Fill(eleTrackerIso_minDREle[iele]); }
	if(eleCaloIso_minDREle[iele]>0 && eleCaloIso_minDREle[iele]<0.5)                  { E_dRmin_calo_noVeto         -> Fill(eleCaloIso_minDREle[iele]);}
	if(eleCaloIso_sumPtEle[iele]>0 && eleCaloIso_sumPtEle[iele]<0.5)                  { E_ptFrac_calo               -> Fill(eleCaloIso_sumPtEle[iele]);}

	// computing the efficiency of iso cuts:
	if (eleTrackerIso_sumPtEle[iele]<0.05){ trackIso++; }
	if (eleCaloIso_sumPtEle[iele]<0.05)   { caloIso++; }

	// computing the efficiency of iso cuts:
	if (eleTrackerIso_sumPtEle[iele]<0.05){ trackIso++; }
	if (eleCaloIso_sumPtEle[iele]<0.05)   { caloIso++; }
	
	// efficiency: all reco
	E_Reco_eta_wogm ->Fill(etaEle[iele]);
	
	// efficiency: all gene
	for(int imc=0;imc<nMc;imc++) {
	  if (idMc[imc] == 11 || idMc[imc] == -11){E_Gene_eta->Fill(etaMc[imc]); }} 
	
	// efficiency: matching
	if (dr_min_wc  < 0.05){ E_Reco_eta_wcm -> Fill(recEta); }
	if (dr_min_woc < 0.05){ E_Reco_eta_wgm -> Fill(recEta); }

	if(eleClassEle[iele]<99) ntotEB++;
	else ntotEE++;

	int GsfClass0[2], GsfClass1[2], GsfClass2[2];
	bool GsfClass3[2];
	if ( ifullclass == 0 ) nGsfClass0[iecal]++;
	else if ( ifullclass == 1 ) nGsfClass1[iecal]++;
	else if ( ifullclass == 2 ) nGsfClass2[iecal]++;
	else if ( ifullclass == 3 ) nGsfClass3[iecal]++;
	
	dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaUnsplitEle        [iecal][iptbin] -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( eleCorrEoPoutEle[iele] );
	HoEUnsplitEle         [iecal][iptbin] -> Fill ( eleHoEEle[iele] );
	shapeFisherUnsplitEle [iecal][iptbin] -> Fill ( fisher ); 
	sigmaEtaEtaUnsplitEle [iecal][iptbin] -> Fill ( sqrt(covEtaEtaEle[iele]) );
	sigmaEtaPhiUnsplitEle [iecal][iptbin] -> Fill ( sqrt(fabs(covEtaPhiEle[iele])) );
	sigmaPhiPhiUnsplitEle [iecal][iptbin] -> Fill ( sqrt(covPhiPhiEle[iele]) );
	s1s9UnsplitEle        [iecal][iptbin] -> Fill ( s1s9Ele[iele] );
	s9s25UnsplitEle       [iecal][iptbin] -> Fill ( s9s25Ele[iele] );
	LATUnsplitEle         [iecal][iptbin] -> Fill ( latEle[iele] );
	etaLATUnsplitEle      [iecal][iptbin] -> Fill ( etaLatEle[iele] );
	phiLATUnsplitEle      [iecal][iptbin] -> Fill ( phiLatEle[iele] );
	a20UnsplitEle         [iecal][iptbin] -> Fill ( a20Ele[iele] );
	a42UnsplitEle         [iecal][iptbin] -> Fill ( a20Ele[iele] );
	
	dPhiCaloClassEle    [iecal][iptbin][iclass] -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxClassEle     [iecal][iptbin][iclass] -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaClassEle        [iecal][iptbin][iclass] -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutClassEle      [iecal][iptbin][iclass] -> Fill ( eleCorrEoPoutEle[iele] );
	HoEClassEle         [iecal][iptbin][iclass] -> Fill ( eleHoEEle[iele] );
	shapeFisherClassEle [iecal][iptbin][iclass] -> Fill ( fisher ); 
	sigmaEtaEtaClassEle [iecal][iptbin][iclass] -> Fill ( sqrt(covEtaEtaEle[iele]) );
	sigmaEtaPhiClassEle [iecal][iptbin][iclass] -> Fill ( sqrt(fabs(covEtaPhiEle[iele])) );
	sigmaPhiPhiClassEle [iecal][iptbin][iclass] -> Fill ( sqrt(covPhiPhiEle[iele]) );
	s1s9ClassEle        [iecal][iptbin][iclass] -> Fill ( s1s9Ele[iele] );
	s9s25ClassEle       [iecal][iptbin][iclass] -> Fill ( s9s25Ele[iele] );
	LATClassEle         [iecal][iptbin][iclass] -> Fill ( latEle[iele] );
	etaLATClassEle      [iecal][iptbin][iclass] -> Fill ( etaLatEle[iele] );
	phiLATClassEle      [iecal][iptbin][iclass] -> Fill ( phiLatEle[iele] );
	a20ClassEle         [iecal][iptbin][iclass] -> Fill ( a20Ele[iele] );
	a42ClassEle         [iecal][iptbin][iclass] -> Fill ( a20Ele[iele] );

	dPhiCaloFullclassEle    [iecal][iptbin][ifullclass] -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxFullclassEle     [iecal][iptbin][ifullclass] -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( eleCorrEoPoutEle[iele] );
	HoEFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( eleHoEEle[iele] );
	shapeFisherFullclassEle [iecal][iptbin][ifullclass] -> Fill ( fisher ); 
	sigmaEtaEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sqrt(covEtaEtaEle[iele]) );
	sigmaEtaPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sqrt(fabs(covEtaPhiEle[iele])) );
	sigmaPhiPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sqrt(covPhiPhiEle[iele]) );
	s1s9FullclassEle        [iecal][iptbin][ifullclass] -> Fill ( s1s9Ele[iele] );
	s9s25FullclassEle       [iecal][iptbin][ifullclass] -> Fill ( s9s25Ele[iele] );
	LATFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( latEle[iele] );
	etaLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( etaLatEle[iele] );
	phiLATFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( phiLatEle[iele] );
	a20FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a20Ele[iele] );
	a42FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( a20Ele[iele] );
	
      }

      // fill the StatPatternRecognition file
      sigmaEtaEta = sqrt(covEtaEtaEle[iele]);
      sigmaEtaPhi = sqrt(fabs(covEtaPhiEle[iele]));
      sigmaPhiPhi = sqrt(covPhiPhiEle[iele]);
      s1s9 = s1s9Ele[iele];
      s9s25 = s9s25Ele[iele];
      LAT = latEle[iele];
      etaLAT = etaLatEle[iele];
      phiLAT = phiLatEle[iele];
      a20 = a20Ele[iele];
      a42 = a42Ele[iele];

      int signal = -1;
      if ( strcmp(category,"hadrons")==0 ) signal = 0;
      else if ( strcmp(category,"electrons")==0 ) signal = 1;
      clusterShapeFiller.fill(jentry,signal);

    }

  }

}



