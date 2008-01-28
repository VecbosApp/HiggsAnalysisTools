#include <iostream>
#include <string>

#include <TTree.h>
#include <TMath.h>

#include "Counters.hh"
#include "Selection.hh"
#include "plotsEleID.hh"


plotsEleID::plotsEleID(TTree *tree) 
  : HiggsBase(tree) {

  // sample region
  //category = "hadrons";  
  category = "electrons";  
  splitelectrons = false;
  isGt15GeV = false;

  // booking histos eleID
  // ii = 0 --> barrel
  // ii = 1 --> endcap
  for (int ii=0; ii<2; ii++){

    dPhiCalo[ii]    = new TH1F("dPhiCalo_hadrons",   "dPhiCalo_hadrons",   100,-0.1,0.1);
    dPhiVtx[ii]     = new TH1F("dPhiVtx_hadrons",    "dPhiVtx_hadrons",    100,-0.1,0.1);
    dEta[ii]        = new TH1F("dEta_hadrons",       "dEta_hadrons",       100,-0.02,0.02);
    EoPout[ii]      = new TH1F("EoPout_hadrons",     "EoPout_hadrons"     ,100,0.0,4.0);
    HoE[ii]         = new TH1F("HoE_hadrons",        "HoE_hadrons",         50,-0.1,0.1);
    shapeFisher[ii] = new TH1F("shapeFisher_hadrons","shapeFisher_hadrons",100,-15,15);

    dPhiCaloUnsplit[ii]    = new TH1F("dPhiCalo_electrons",   "dPhiCalo_electrons",   100,-0.1,0.1);
    dPhiVtxUnsplit[ii]     = new TH1F("dPhiVtx_electrons",    "dPhiVtx_electrons",    100,-0.1,0.1);
    dEtaUnsplit[ii]        = new TH1F("dEtaCalo_electrons",   "dEtaCalo_electrons",       100,-0.02,0.02);
    EoPoutUnsplit[ii]      = new TH1F("EoPout_electrons",     "EoPout_electrons"     ,100,0.0,4.0);
    HoEUnsplit[ii]         = new TH1F("HoE_electrons",        "HoE_electrons",         50,-0.1,0.1);
    shapeFisherUnsplit[ii] = new TH1F("shapeFisher_electrons","shapeFisher_electrons",100,-15,15);

    dPhiCalo0[ii] = new TH1F("dPhiCalo_electrons_GsfClass0","dPhiCalo_electrons_GsfClass0",100,-0.05,0.05);
    dPhiCalo1[ii] = new TH1F("dPhiCalo_electrons_GsfClass1","dPhiCalo_electrons_GsfClass1",100,-0.05,0.05);
    dPhiCalo2[ii] = new TH1F("dPhiCalo_electrons_GsfClass2","dPhiCalo_electrons_GsfClass2",100,-0.05,0.05);
    dPhiCalo3[ii] = new TH1F("dPhiCalo_electrons_GsfClass3","dPhiCalo_electrons_GsfClass3",100,-0.05,0.05);
    
    dPhiVtx0[ii]  = new TH1F("dPhiVtx_electrons_GsfClass0","dPhiVtx_electrons_GsfClass0",100,-0.05,0.05);
    dPhiVtx1[ii]  = new TH1F("dPhiVtx_electrons_GsfClass1","dPhiVtx_electrons_GsfClass1",100,-0.05,0.05);
    dPhiVtx2[ii]  = new TH1F("dPhiVtx_electrons_GsfClass2","dPhiVtx_electrons_GsfClass2",100,-0.05,0.05);
    dPhiVtx3[ii]  = new TH1F("dPhiVtx_electrons_GsfClass3","dPhiVtx_electrons_GsfClass3",100,-0.05,0.05);
    
    dEta0[ii]   = new TH1F("dEtaCalo_electrons_GsfClass0","dEtaCalo_electrons_GsfClass0",100,-0.02,0.02);
    dEta1[ii]   = new TH1F("dEtaCalo_electrons_GsfClass1","dEtaCalo_electrons_GsfClass1",100,-0.02,0.02);
    dEta2[ii]   = new TH1F("dEtaCalo_electrons_GsfClass2","dEtaCalo_electrons_GsfClass2",100,-0.02,0.02);
    dEta3[ii]   = new TH1F("dEtaCalo_electrons_GsfClass3","dEtaCalo_electrons_GsfClass3",100,-0.02,0.02);
    
    EoPout0[ii] = new TH1F("EoPout_electrons_GsfClass0","EoPout_electrons_GsfClass0",100,0.0,4.0);
    EoPout1[ii] = new TH1F("EoPout_electrons_GsfClass1","EoPout_electrons_GsfClass1",100,0.0,4.0);
    EoPout2[ii] = new TH1F("EoPout_electrons_GsfClass2","EoPout_electrons_GsfClass2",100,0.0,4.0);
    EoPout3[ii] = new TH1F("EoPout_electrons_GsfClass3","EoPout_electrons_GsfClass3",100,0.0,4.0);
    
    HoE0[ii] = new TH1F("HoE_electrons_GsfClass0","HoE_electrons_GsfClass0",50,-0.1,0.1);
    HoE1[ii] = new TH1F("HoE_electrons_GsfClass1","HoE_electrons_GsfClass1",50,-0.1,0.1);
    HoE2[ii] = new TH1F("HoE_electrons_GsfClass2","HoE_electrons_GsfClass2",50,-0.1,0.1);
    HoE3[ii] = new TH1F("HoE_electrons_GsfClass3","HoE_electrons_GsfClass3",50,-0.1,0.1);
    
    shapeFisher0[ii] = new TH1F("shapeFisher_electrons_GsfClass0","shapeFisher_electrons_GsfClass0",100,-15,15);
    shapeFisher1[ii] = new TH1F("shapeFisher_electrons_GsfClass1","shapeFisher_electrons_GsfClass1",100,-15,15);
    shapeFisher2[ii] = new TH1F("shapeFisher_electrons_GsfClass2","shapeFisher_electrons_GsfClass2",100,-15,15);
    shapeFisher3[ii] = new TH1F("shapeFisher_electrons_GsfClass3","shapeFisher_electrons_GsfClass3",100,-15,15);
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
  tfileEB  = new TFile("pdfEB.root","RECREATE");
  tfileEE  = new TFile("pdfEE.root","RECREATE");
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
  if(strcmp(category,"hadrons")==0){

    // barrel
    tfileEB->cd();
    dPhiCalo[0]    -> Write();
    dPhiVtx[0]     -> Write();
    dEta[0]        -> Write();
    EoPout[0]      -> Write();
    HoE[0]         -> Write();
    shapeFisher[0] -> Write();

    // endcap
    tfileEB->cd();
    dPhiCalo[1]    -> Write();
    dPhiVtx[1]     -> Write();
    dEta[1]        -> Write();
    EoPout[1]      -> Write();
    HoE[1]         -> Write();
    shapeFisher[1] -> Write();

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
  }


  if(strcmp(category,"electrons")==0){

    if(splitelectrons) {
      // barrel
      tfileEB->cd();
      dPhiCalo0[0] -> Write(); dPhiVtx0[0] -> Write();  dEta0[0] -> Write();  EoPout0[0] -> Write(); HoE0[0] -> Write(); shapeFisher0[0] -> Write();
      dPhiCalo1[0] -> Write(); dPhiVtx1[0] -> Write();  dEta1[0] -> Write();  EoPout1[0] -> Write(); HoE1[0] -> Write(); shapeFisher1[0] -> Write();
      dPhiCalo2[0] -> Write(); dPhiVtx2[0] -> Write();  dEta2[0] -> Write();  EoPout2[0] -> Write(); HoE2[0] -> Write(); shapeFisher2[0] -> Write();
      dPhiCalo3[0] -> Write(); dPhiVtx3[0] -> Write();  dEta3[0] -> Write();  EoPout3[0] -> Write(); HoE3[0] -> Write(); shapeFisher3[0] -> Write();
      
      // endcap
      tfileEE->cd();
      dPhiCalo0[1] -> Write(); dPhiVtx0[1] -> Write();  dEta0[1] -> Write();  EoPout0[1] -> Write(); HoE0[1] -> Write(); shapeFisher0[1] -> Write();
      dPhiCalo1[1] -> Write(); dPhiVtx1[1] -> Write();  dEta1[1] -> Write();  EoPout1[1] -> Write(); HoE1[1] -> Write(); shapeFisher1[1] -> Write();
      dPhiCalo2[1] -> Write(); dPhiVtx2[1] -> Write();  dEta2[1] -> Write();  EoPout2[1] -> Write(); HoE2[1] -> Write(); shapeFisher2[1] -> Write();
      dPhiCalo3[1] -> Write(); dPhiVtx3[1] -> Write();  dEta3[1] -> Write();  EoPout3[1] -> Write(); HoE3[1] -> Write(); shapeFisher3[1] -> Write();
    }
    else {
      // barrel
      tfileEB->cd();
      dPhiCaloUnsplit[0] -> Write(); dPhiVtxUnsplit[0] -> Write();  dEtaUnsplit[0] -> Write();  EoPoutUnsplit[0] -> Write(); HoEUnsplit[0] -> Write(); shapeFisherUnsplit[0] -> Write();

      // endcap
      tfileEE->cd();
      dPhiCaloUnsplit[1] -> Write(); dPhiVtxUnsplit[1] -> Write();  dEtaUnsplit[1] -> Write();  EoPoutUnsplit[1] -> Write(); HoEUnsplit[1] -> Write(); shapeFisherUnsplit[1] -> Write();
    }

    // isolation
    tfileIso->cd();
    E_dRmin_tracker_withVeto      -> Write();
    E_dRmin_tracker_noVeto        -> Write();
    E_ptFrac_tracker              -> Write();
    E_dRmin_tracker_withVeto_zoom -> Write();
    E_dRmin_tracker_noVeto_zoom   -> Write();
    E_ptFrac_tracker_zoom         -> Write();
    E_dRmin_calo_noVeto           -> Write();
    E_ptFrac_calo                 -> Write();

    tfileEff->cd();
    E_Reco_eta_wcm     ->Write();
    E_Reco_eta_wgm     ->Write();
    E_Reco_eta_wogm    ->Write();
    //    E_Reco_eta_tiso    ->Write();
    E_Gene_eta         ->Write();
  }

  tfileEB ->Close();
  tfileEE ->Close();
  tfileIso->Close();

  // deleting
  for (int ii=0; ii<2; ii++){
    delete dPhiCalo[ii];    
    delete dPhiVtx[ii];    
    delete dEta[ii];        
    delete EoPout[ii];      
    delete HoE[ii];         
    delete shapeFisher[ii];
    delete dPhiCaloUnsplit[ii],    dPhiCaloUnsplit[ii],    dPhiCaloUnsplit[ii],    dPhiCaloUnsplit[ii]; 
    delete dPhiCalo0[ii],    dPhiCalo1[ii],    dPhiCalo2[ii],    dPhiCalo3[ii]; 
    delete dPhiVtx0[ii],     dPhiVtx1[ii],     dPhiVtx2[ii],     dPhiVtx3[ii]; 
    delete dEta0[ii],        dEta1[ii],        dEta2[ii],        dEta3[ii]; 
    delete EoPout0[ii],      EoPout1[ii],      EoPout2[ii],      EoPout3[ii]; 
    delete HoE0[ii],         HoE1[ii],         HoE2[ii],         HoE3[ii]; 
    delete shapeFisher0[ii], shapeFisher1[ii], shapeFisher2[ii], shapeFisher3[ii]; 
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

      int isEndcap = -1;
      float thisEta = etaEle[iele];  
      if (fabs(thisEta) <= 1.476){ isEndcap = 0; }   // barrel
      if (fabs(thisEta) > 1.476) { isEndcap = 1; }   // endcap 

      bool dumpPtBin = false;
      if(isGt15GeV && etEle[iele]>15.0) dumpPtBin = true;
      else if((!isGt15GeV) && etEle[iele]<15.0) dumpPtBin = true;

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
	dPhiCalo[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	dPhiVtx [isEndcap]->Fill(eleDeltaPhiAtVtxEle[iele]);
	dEta    [isEndcap]->Fill(eleDeltaEtaAtCaloEle[iele]);
	EoPout  [isEndcap]->Fill(eleCorrEoPoutEle[iele]);
	HoE     [isEndcap]->Fill(eleHoEEle[iele]);

	if(isEndcap==0)
	  shapeFisher[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	if(isEndcap==1)
	  shapeFisher[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);

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
      
      if(strcmp(category,"electrons")==0) {
	
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
	GsfClass0[0] = 0;
	GsfClass1[0] = 10;
	GsfClass2[0] = 20;
	GsfClass3[0] = eleClassEle[iele]==30 || eleClassEle[iele]==31 || eleClassEle[iele]==32 || eleClassEle[iele]==33 || eleClassEle[iele]==34;
	GsfClass0[1] = 100;
	GsfClass1[1] = 110;
	GsfClass2[1] = 120;
	GsfClass3[1] = eleClassEle[iele]==130 || eleClassEle[iele]==131 || eleClassEle[iele]==132 || eleClassEle[iele]==133 || eleClassEle[iele]==134;

	if(dumpPtBin) {
	  dPhiCaloUnsplit[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	  dPhiVtxUnsplit[isEndcap] ->Fill(eleDeltaPhiAtVtxEle[iele]);
	  dEtaUnsplit[isEndcap]    ->Fill(eleDeltaEtaAtCaloEle[iele]);
	  EoPoutUnsplit[isEndcap]  ->Fill(eleCorrEoPoutEle[iele]);
	  HoEUnsplit[isEndcap]     ->Fill(eleHoEEle[iele]);

	  if(isEndcap==0)
	    shapeFisherUnsplit[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	  if(isEndcap==1)
	    shapeFisherUnsplit[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);

	  if(eleClassEle[iele]==GsfClass0[isEndcap]) {
	    nGsfClass0[isEndcap]++;
	    dPhiCalo0[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	    dPhiVtx0[isEndcap] ->Fill(eleDeltaPhiAtVtxEle[iele]);
	    dEta0[isEndcap]    ->Fill(eleDeltaEtaAtCaloEle[iele]);
	    EoPout0[isEndcap]  ->Fill(eleCorrEoPoutEle[iele]);
	    HoE0[isEndcap]     ->Fill(eleHoEEle[iele]);
	    if(isEndcap==0)
	      shapeFisher0[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	    if(isEndcap==1)
	      shapeFisher0[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);
	  }
	  else if(eleClassEle[iele]==GsfClass1[isEndcap]) {
	    nGsfClass1[isEndcap]++;
	    dPhiCalo1[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	    dPhiVtx1[isEndcap] ->Fill(eleDeltaPhiAtVtxEle[iele]);
	    dEta1[isEndcap]    ->Fill(eleDeltaEtaAtCaloEle[iele]);
	    EoPout1[isEndcap]  ->Fill(eleCorrEoPoutEle[iele]);
	    HoE1[isEndcap]     ->Fill(eleHoEEle[iele]);
	    if(isEndcap==0)
	      shapeFisher1[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	    if(isEndcap==1)
	      shapeFisher1[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);
	  }
	  else if(eleClassEle[iele]==GsfClass2[isEndcap]) {
	    nGsfClass2[isEndcap]++;
	    dPhiCalo2[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	    dPhiVtx2[isEndcap] ->Fill(eleDeltaPhiAtVtxEle[iele]);
	    dEta2[isEndcap]    ->Fill(eleDeltaEtaAtCaloEle[iele]);
	    EoPout2[isEndcap]  ->Fill(eleCorrEoPoutEle[iele]);
	    HoE2[isEndcap]     ->Fill(eleHoEEle[iele]);
	    if(isEndcap==0)
	      shapeFisher2[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	    if(isEndcap==1)
	      shapeFisher2[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);
	  }
	  else if(GsfClass3[isEndcap]) {
	    nGsfClass3[isEndcap]++;
	    dPhiCalo3[isEndcap]->Fill(eleDeltaPhiAtCaloEle[iele]);
	    dPhiVtx3[isEndcap] ->Fill(eleDeltaPhiAtVtxEle[iele]);
	    dEta3[isEndcap]    ->Fill(eleDeltaEtaAtCaloEle[iele]);
	    EoPout3[isEndcap]  ->Fill(eleCorrEoPoutEle[iele]);
	    HoE3[isEndcap]     ->Fill(eleHoEEle[iele]);
	    if(isEndcap==0)
	      shapeFisher3[isEndcap]->Fill(42.0238-3.38943*s9s25Ele[iele]-794.092*sqrt(covEtaEtaEle[iele])-15.3449*latEle[iele]-31.1032*a20Ele[iele]);
	    if(isEndcap==1)
	      shapeFisher3[isEndcap]->Fill(27.2967+2.97453*s9s25Ele[iele]-169.219*sqrt(covEtaEtaEle[iele])-17.0445*latEle[iele]-24.8542*a20Ele[iele]);
	  }
	} 
      }
    }
  }
}



