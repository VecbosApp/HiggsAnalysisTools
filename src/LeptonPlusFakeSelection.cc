#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"
#include "HiggsAnalysisTools/include/LeptonPlusFakeSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"

#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>

LeptonPlusFakeSelection::LeptonPlusFakeSelection(TTree *tree) 
  : HiggsBase(tree) {

  // choose the Higgs Mass
  std::string higgsConfigDir;
  std::string higgsConfigDirMass;
  std::ifstream setfile("config/higgs/higgsMass.txt");
  higgsConfigDir="config/higgs/";
  if(!setfile.good()) {
    std::cout << "Cannot read the higgsMass file to choose the selection: config/higgs/higgsMass.txt" << std::endl
	      << "using default (160 GeV)" << std::endl;
    higgsConfigDirMass="config/higgs/h160/";
  }
  else {
    std::string var, massVal; 
    bool found=false;
    while (1) {
      setfile >> var >> massVal;
      _massVal = atoi(massVal.c_str());
      if(!setfile.good()) break;
      if(var.compare("HiggsMass")==0) { 
	found=true;
	higgsConfigDirMass="config/higgs/h" + massVal + "/";
	std::cout << "Reading configuration for Higgs mass = " << massVal << " GeV/c^2" << std::endl;
	break;
      }
    }
  }

  // selection efficiencies
  std::string fileCutsEE     = higgsConfigDirMass + "2e2nuCuts.txt";
  std::string fileSwitchesEE = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedLeptonPlusFakeSelectionEE.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str()); 
  _selectionEE = CutBasedLeptonPlusFakeSelectionEE.GetSelection();  

  // single electron efficiency
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/looseEleId/"); 
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/tightEleId/"); 
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/hwwAnEleId/");
  EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/newOptimEleId_looseOthers_m160/");
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/newOptimEleId_tightOthers_m160/");
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/newOptimEleId_looseOthers_m190/");
  // EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/newOptimEleId_tightOthers_m190/");

}

LeptonPlusFakeSelection::~LeptonPlusFakeSelection(){

  delete _selectionEE;

}

void LeptonPlusFakeSelection::initialiseFakeRate() {

  // binning
  m_minFakePt[0] = 10.;  m_maxFakePt[0] = 15.;
  m_minFakePt[1] = 15.;  m_maxFakePt[1] = 20.;
  m_minFakePt[2] = 20.;  m_maxFakePt[2] = 25.;
  m_minFakePt[3] = 25.;  m_maxFakePt[3] = 30.;
  m_minFakePt[4] = 30.;  m_maxFakePt[4] = 35.;
  m_minFakePt[5] = 35.;  m_maxFakePt[5] = 40.;
  m_minFakePt[6] = 40.;  m_maxFakePt[6] = 14000.;

  m_fakeRateFromQCD[0] = 0.00174309;  m_fakeRateFromQCD_err[0] = 6.45959e-05;
  m_fakeRateFromQCD[1] = 0.00119603;  m_fakeRateFromQCD_err[1] = 4.88512e-05;
  m_fakeRateFromQCD[2] = 0.00102552;  m_fakeRateFromQCD_err[2] = 4.47325e-05;
  m_fakeRateFromQCD[3] = 0.000842069; m_fakeRateFromQCD_err[3] = 3.85378e-05;
  m_fakeRateFromQCD[4] = 0.000775075; m_fakeRateFromQCD_err[4] = 3.55649e-05;
  m_fakeRateFromQCD[5] = 0.000815457; m_fakeRateFromQCD_err[5] = 4.31591e-05;
  m_fakeRateFromQCD[6] = 0.00284392;  m_fakeRateFromQCD_err[6] = 4.91138e-05;

}

void LeptonPlusFakeSelection::Loop() {

  if(fChain == 0) return;
  
  initialiseFakeRate();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    resetKinematics();

    float weight = 1;
    
    // get the alpgen weight to normalize correctly jet bins
    weight = genWeight; 

    // trigger
    Utils anaUtils;
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

    // get the best L1 + Fake combination
    getBestLplusFakePair();

    if (m_theL1 == -1 || m_theFake == -1) continue;

    // set the kinematic variables for the rest of the selection
    setKinematics();

    // weight with the Fake -> L2 probability
    float fakerate = getFakeRate( m_theFake4Momentum.Pt(), QCD );
    weight = weight * fakerate;
    if (fakerate<0) cout << "BIG ERROR: m_theFake4Momentum.Pt() = " << m_theFake4Momentum.Pt() << endl;


    // ----------------------- selection ----------------------------
    // electron ID (true by default - studied only if ee or emu channel)
    bool theElectronID = true;
    // custom electron ID
    if ( m_theL1 > -1 ) theElectronID = isEleID(m_theL1);

    // tracker isolation for electrons
    float theEleTrackerPtSum = 0.;
    if ( m_theL1 > -1 ) theEleTrackerPtSum = eleSumPtPreselectionEle[m_theL1];

    // hcal isolation for electrons
    float theEleHcalPtSum = 0.;
    if ( m_theL1 > -1 ) theEleHcalPtSum = eleSumHadEt04Ele[m_theL1];

    // ecal isolation for electrons
    float theEleEcalPtSum = 0.;
    if ( m_theL1 > -1 ) theEleEcalPtSum = eleSumEmEt04Ele[m_theL1];

    // jet veto: method gives true if good jet is found (L1 and Fake are excluded)
    bool passedJetVeto = !goodJetFound();

    // selections
    CutBasedLeptonPlusFakeSelectionEE.SetProcessID((int)genAlpgenID/1000);
    CutBasedLeptonPlusFakeSelectionEE.SetWeight( weight );
    CutBasedLeptonPlusFakeSelectionEE.SetHighElePt( m_theL14Momentum.Pt() );
    CutBasedLeptonPlusFakeSelectionEE.SetLowElePt( m_theFake4Momentum.Pt() );
    CutBasedLeptonPlusFakeSelectionEE.SetElectronId(theElectronID);
    CutBasedLeptonPlusFakeSelectionEE.SetPositronId( true );
    CutBasedLeptonPlusFakeSelectionEE.SetEleTrackerPtSum( theEleTrackerPtSum );
    CutBasedLeptonPlusFakeSelectionEE.SetPosTrackerPtSum( 0. );
    CutBasedLeptonPlusFakeSelectionEE.SetEleHcalPtSum( theEleHcalPtSum );
    CutBasedLeptonPlusFakeSelectionEE.SetPosHcalPtSum( 0. );
    CutBasedLeptonPlusFakeSelectionEE.SetEleEcalPtSum( theEleEcalPtSum );
    CutBasedLeptonPlusFakeSelectionEE.SetPosEcalPtSum( 0. );
    CutBasedLeptonPlusFakeSelectionEE.SetJetVeto(passedJetVeto);
    CutBasedLeptonPlusFakeSelectionEE.SetMet( etMet[0] );					
    CutBasedLeptonPlusFakeSelectionEE.SetDeltaPhi( m_deltaPhi );
    CutBasedLeptonPlusFakeSelectionEE.SetInvMass( m_mll );
    CutBasedLeptonPlusFakeSelectionEE.SetDetaLeptons( 0. );
    bool selOutput = CutBasedLeptonPlusFakeSelectionEE.output();

  }

  CutBasedLeptonPlusFakeSelectionEE.diplayEfficiencies();

}

void LeptonPlusFakeSelection::getBestLplusFakePair() {
  int theLep1=-1;
  int theFake=-1;
  float maxPtLep1=-1000.;

  for(int i=0;i<nEle;i++) {
    if( fabs(etaEle[i]) > 2.5 ) continue;
    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();
    if( thisPt < 20 ) continue;
    if (thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
  }
  m_theL1 = theLep1;
  if ( theLep1 > -1 ) 
    m_theL14Momentum.SetXYZT(pxEle[theLep1],pyEle[theLep1],pzEle[theLep1],energyEle[theLep1]);
  else 
    m_theL14Momentum.SetXYZT(0,0,0,0);
  
  float maxPtJet=-1000.;
  for(int j=0;j<nJet;j++) {

    // check if the electron or the positron falls into the jet
    // common cleaning class
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    if ( m_theL14Momentum.Vect().Mag() != 0 ) {
      float deltaR =  p3Jet.DeltaR( m_theL14Momentum.Vect() );
      // taking from ee config file, but jets veto is the same for all the channels
      if(_selectionEE->getSwitch("jetConeWidth") && 
	 _selectionEE->passCut("jetConeWidth",deltaR)
	// && (m_HoEElectronMinus < 0.2)  
	// && (m_CaloEneElectronMinus/energyJet[j] > 0.9)
	 ) continue;
    }

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc",etaJet[j])) continue;
    if(_selectionEE->getSwitch("etJetLowAcc") && !_selectionEE->passCut("etJetLowAcc",etJet[j]) ) continue;

    if(etJet[j]>maxPtJet) { maxPtJet=etJet[j]; theFake=j; }

    break;
  }

  m_theFake = theFake;
  if ( theFake > -1 )
    m_theFake4Momentum.SetXYZT(pxJet[theFake],pyJet[theFake],pzJet[theFake],energyJet[theFake]);
  else
    m_theFake4Momentum.SetXYZT(0,0,0,0);
  
}


bool LeptonPlusFakeSelection::isEleID(int eleIndex) {

  TVector3 pTrkAtOuter(pxAtOuterEle[eleIndex],pyAtOuterEle[eleIndex],pzAtOuterEle[eleIndex]);

  EgammaCutBasedID.SetHOverE( eleHoEEle[eleIndex] );
  EgammaCutBasedID.SetS9S25( s9s25Ele[eleIndex] );
  EgammaCutBasedID.SetDEta( eleDeltaEtaAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiIn( eleDeltaPhiAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiOut( eleDeltaPhiAtCaloEle[eleIndex] );
  EgammaCutBasedID.SetInvEminusInvP( 1./eleCaloCorrEEle[eleIndex]-1./eleTrackerPEle[eleIndex] );
  EgammaCutBasedID.SetBremFraction( fabs(eleTrackerPEle[eleIndex]-pTrkAtOuter.Mag())/eleTrackerPEle[eleIndex] );
  EgammaCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[eleIndex]) );
  EgammaCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[eleIndex]) );
  EgammaCutBasedID.SetEOverPout( eleCorrEoPoutEle[eleIndex] );
  EgammaCutBasedID.SetEOverPin( eleCorrEoPEle[eleIndex] );
  EgammaCutBasedID.SetElectronClass ( eleClassEle[eleIndex] );
  EgammaCutBasedID.SetEgammaCutBasedID ( eleIdCutBasedEle[eleIndex] );
  EgammaCutBasedID.SetLikelihood( eleLikelihoodEle[eleIndex] );

  bool isIdentified = EgammaCutBasedID.output();

  return isIdentified;
}



void LeptonPlusFakeSelection::setKinematics() {

  m_mll = (m_theL14Momentum + m_theFake4Momentum).M();

  // MET
  m_MET = etMet[0]; 

  // compute delta Phi in degrees
  TVector3 dilepPt;
  m_deltaPhi = fabs( 180./TMath::Pi() * m_theL14Momentum.Vect().DeltaPhi(m_theFake4Momentum.Vect()) );

}


bool LeptonPlusFakeSelection::goodJetFound() {

  // first check that kinematics has been set
  bool foundJet=false;
  for(int j=0;j<nJet;j++) {

    if ( j == m_theFake ) continue;

    // check if the electron or the positron falls into the jet
    // common cleaning class
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    if ( m_theL14Momentum.Vect().Mag() != 0 ) {
      float deltaR =  p3Jet.DeltaR( m_theL14Momentum.Vect() );
      // taking from ee config file, but jets veto is the same for all the channels
      if(_selectionEE->getSwitch("jetConeWidth") && 
	 _selectionEE->passCut("jetConeWidth",deltaR)
	// && (m_HoEElectronMinus < 0.2)  
	// && (m_CaloEneElectronMinus/energyJet[j] > 0.9)
	 ) continue;
    }

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc",etaJet[j])) continue;
    if(_selectionEE->getSwitch("etJetLowAcc") && !_selectionEE->passCut("etJetLowAcc",etJet[j]) ) continue;

    if( (_selectionEE->getSwitch("etJetHighAcc") && _selectionEE->passCut("etJetHighAcc",etJet[j])) &&
 	(_selectionEE->getSwitch("alphaJet") && !_selectionEE->passCut("alphaJet",alphaJet[j])) 
	) continue;
    foundJet=true;
    break;
  }

  return foundJet;

}

void LeptonPlusFakeSelection::resetKinematics() {

  m_theL1 = -1;
  m_theFake = -1;
  m_theL14Momentum.SetXYZT(0,0,0,0);
  m_theFake4Momentum.SetXYZT(0,0,0,0);
  m_deltaPhi = 0;
  m_mll = 0;

}


float LeptonPlusFakeSelection::getFakeRate(float fakePt, int source) {
  for (int i=0; i<7; i++) {
    if( fakePt >= m_minFakePt[i] && fakePt < m_maxFakePt[i] ) {
      if( source == QCD ) {
	return m_fakeRateFromQCD[i];
      } else {
	return m_fakeRateFromWjets[i];
      }
    }
  }
  return -1.;
}

float LeptonPlusFakeSelection::getFakeRateError(float fakePt, int source) {
  for (int i=0; i<7; i++) {
    if( fakePt >= m_minFakePt[i] && fakePt < m_maxFakePt[i] ) {
      if( source == QCD ) {
	return m_fakeRateFromQCD_err[i];
      } else {
	return m_fakeRateFromWjets_err[i];
      }
    }
  }
  return -1.;
}
