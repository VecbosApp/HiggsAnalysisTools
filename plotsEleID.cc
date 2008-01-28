#include <iostream>
#include <string>

#include <TTree.h>
#include <TMath.h>

#include "Counters.hh"
#include "Selection.hh"
#include "plotsEleID.hh"

plotsEleID::plotsEleID(TTree *tree) 
  : HiggsBase(tree) {

  int nbins = 50;

  float dPhiCaloMin = -0.1;
  float dPhiCaloMax =  0.1;
  float dPhiVtxMin  = -0.1;
  float dPhiVtxMax  =  0.1;
  float dEtaMin     = -0.02;
  float dEtaMax     =  0.02;
  float EoPoutMin   =  0.0;
  float EoPoutMax   =  10.0;
  float HoEMin      = -0.1;
  float HoEMax      =  0.1;
  float fisherMin   = -15.0;
  float fisherMax   =  15.0;

  // sample region
  //category = "hadrons";  
  category = "electrons";  

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
      }
      else if(strcmp(category,"electrons")==0) {

	dPhiCaloUnsplitEle[iecal][iptbin]    -> Write();
	dPhiVtxUnsplitEle[iecal][iptbin]     -> Write();
	dEtaUnsplitEle[iecal][iptbin]        -> Write();
	EoPoutUnsplitEle[iecal][iptbin]      -> Write();
	HoEUnsplitEle[iecal][iptbin]         -> Write();
	shapeFisherUnsplitEle[iecal][iptbin] -> Write();

	for(int iclass=0; iclass<2; iclass++) {

	  dPhiCaloClassEle[iecal][iptbin][iclass]    -> Write();
	  dPhiVtxClassEle[iecal][iptbin][iclass]     -> Write();
	  dEtaClassEle[iecal][iptbin][iclass]        -> Write();
	  EoPoutClassEle[iecal][iptbin][iclass]      -> Write();
	  HoEClassEle[iecal][iptbin][iclass]         -> Write();
	  shapeFisherClassEle[iecal][iptbin][iclass] -> Write();

	}

	for(int ifullclass=0; ifullclass<4; ifullclass++) {

	  dPhiCaloFullclassEle[iecal][iptbin][ifullclass]    -> Write();
	  dPhiVtxFullclassEle[iecal][iptbin][ifullclass]     -> Write();
	  dEtaFullclassEle[iecal][iptbin][ifullclass]        -> Write();
	  EoPoutFullclassEle[iecal][iptbin][ifullclass]      -> Write();
	  HoEFullclassEle[iecal][iptbin][ifullclass]         -> Write();
	  shapeFisherFullclassEle[iecal][iptbin][ifullclass] -> Write();

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

      }
      else if(strcmp(category,"electrons")==0) {

	delete dPhiCaloUnsplitEle[iecal][iptbin];
	delete dPhiVtxUnsplitEle[iecal][iptbin];
	delete dEtaUnsplitEle[iecal][iptbin];
	delete EoPoutUnsplitEle[iecal][iptbin];
	delete HoEUnsplitEle[iecal][iptbin];
	delete shapeFisherUnsplitEle[iecal][iptbin];
	
	for(int iclass=0; iclass<2; iclass++) {

	  delete dPhiCaloClassEle[iecal][iptbin][iclass];
	  delete dPhiVtxClassEle[iecal][iptbin][iclass];
	  delete dEtaClassEle[iecal][iptbin][iclass];
	  delete EoPoutClassEle[iecal][iptbin][iclass];
	  delete HoEClassEle[iecal][iptbin][iclass];
	  delete shapeFisherClassEle[iecal][iptbin][iclass];

	}

	for(int ifullclass=0; ifullclass<4; ifullclass++) {

	  delete dPhiCaloFullclassEle[iecal][iptbin][ifullclass];
	  delete dPhiVtxFullclassEle[iecal][iptbin][ifullclass];
	  delete dEtaFullclassEle[iecal][iptbin][ifullclass];
	  delete EoPoutFullclassEle[iecal][iptbin][ifullclass];
	  delete HoEFullclassEle[iecal][iptbin][ifullclass];
	  delete shapeFisherFullclassEle[iecal][iptbin][ifullclass];

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
    
  Long64_t nentries = fChain->GetEntriesFast();

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
	dPhiCaloHad    [iecal][iptbin] -> Fill(eleDeltaPhiAtCaloEle[iele]);
	dPhiVtxHad     [iecal][iptbin] -> Fill(eleDeltaPhiAtVtxEle[iele]);
	dEtaHad        [iecal][iptbin] -> Fill(eleDeltaEtaAtCaloEle[iele]);
	EoPoutHad      [iecal][iptbin] -> Fill(eleCorrEoPoutEle[iele]);
	HoEHad         [iecal][iptbin] -> Fill(eleHoEEle[iele]);
	shapeFisherHad [iecal][iptbin] -> Fill(fisher);

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
	
	dPhiCaloUnsplitEle[iecal][iptbin]   -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxUnsplitEle[iecal][iptbin]    -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaUnsplitEle[iecal][iptbin]       -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutUnsplitEle[iecal][iptbin]     -> Fill ( eleCorrEoPoutEle[iele] );
	HoEUnsplitEle[iecal][iptbin]        -> Fill ( eleHoEEle[iele] );
	shapeFisherUnsplitEle[iecal][iptbin] -> Fill ( fisher ); 
	
	dPhiCaloClassEle[iecal][iptbin][iclass]   -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxClassEle[iecal][iptbin][iclass]    -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaClassEle[iecal][iptbin][iclass]       -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutClassEle[iecal][iptbin][iclass]     -> Fill ( eleCorrEoPoutEle[iele] );
	HoEClassEle[iecal][iptbin][iclass]        -> Fill ( eleHoEEle[iele] );
	shapeFisherClassEle[iecal][iptbin][iclass] -> Fill ( fisher ); 
	
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]   -> Fill ( eleDeltaPhiAtCaloEle[iele] );
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]    -> Fill ( eleDeltaPhiAtVtxEle[iele] );
	dEtaFullclassEle[iecal][iptbin][ifullclass]       -> Fill ( eleDeltaEtaAtCaloEle[iele] );
	EoPoutFullclassEle[iecal][iptbin][ifullclass]     -> Fill ( eleCorrEoPoutEle[iele] );
	HoEFullclassEle[iecal][iptbin][ifullclass]        -> Fill ( eleHoEEle[iele] );
	shapeFisherFullclassEle[iecal][iptbin][ifullclass] -> Fill ( fisher ); 
	
      }

    }

  }

}



