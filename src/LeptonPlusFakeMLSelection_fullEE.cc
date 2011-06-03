#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/LeptonPlusFakeMLSelection_fullEE.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <TTree.h>

using namespace bits;

LeptonPlusFakeMLSelection_fullEE::LeptonPlusFakeMLSelection_fullEE(TTree *tree) 
  : Higgs(tree) {
  
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

  CutBasedHiggsSelectionEE1.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE");  
  CutBasedHiggsErrorsSelectionEE1.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION ERRORS EE");   

  _selectionEE1    = CutBasedHiggsSelectionEE1.GetSelection();  
  _selectionErrEE1 = CutBasedHiggsErrorsSelectionEE1.GetSelection();  

  //  extra selection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  _selectionEE1->addCut("etaElectronAcc");    
  _selectionEE1->addCut("ptElectronAcc");
  _selectionEE1->addCut("etaMuonAcc");
  _selectionEE1->addCut("ptMuonAcc");
  _selectionEE1->addCut("etUncorrJetAcc");  
  _selectionEE1->addSwitch("apply_kFactor");   
  _selectionEE1->addSwitch("isData");
  _selectionEE1->addSwitch("goodRunLS");
  _selectionEE1->addSwitch("asymmetricID");
  _selectionEE1->addStringParameter("electronIDType");
  _selectionEE1->addStringParameter("electronIDTypeLow");  

  // cut based electron id or likelihood 
  TString selectionString(_selectionEE1->getStringParameter("electronIDType"));
  if (!_selectionEE1->getSwitch("asymmetricID")) 
    cout << "=== CONFIGURING " << selectionString << " CUT BASED SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/higgs/electronId/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/higgs/electronId/");
  
  if (_selectionEE1->getSwitch("asymmetricID")) {
    TString selectionStringLow (_selectionEE1->getStringParameter("electronIDTypeLow"));
    cout << "=== CONFIGURING "  << selectionStringLow << " and " 
	 << selectionString << " for CUT BASED ASYMMETRIC ELECTRON ID ===" << endl;
    EgammaCutBasedIDLow.ConfigureNoClass("config/higgs/electronId/"+selectionStringLow);
    EgammaCutBasedIDLow.ConfigureEcalCleaner("config/higgs/electronId/");
  }
  
  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useOneOverEMinusOneOverP = true;
  LH = new ElectronLikelihood(&(*EB0lt15dir), &(*EB1lt15dir), &(*EElt15dir), &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir), defaultSwitches, std::string("class"),std::string("class"),true,true);
  
  // Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _selectionEE1->getSwitch("goodRunLS") << " isData is " <<  _selectionEE1->getSwitch("isData") << std::endl;

  // To read good run list!
  if (_selectionEE1->getSwitch("goodRunLS") && _selectionEE1->getSwitch("isData")) {
    std::string goodRunJsonFile = "config/json/goodCollisions2011.json";         // chiara
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics
  for(int theChannel=0; theChannel<1; theChannel++) {
    m_p4LeptonPlus[theChannel]  = new TLorentzVector(0.,0.,0.,0.);
    m_p4LeptonMinus[theChannel] = new TLorentzVector(0.,0.,0.,0.);
    m_p3PFMET = new TVector3(0.,0.,0.);
    m_p3TKMET = new TVector3(0.,0.,0.);
  }    

  // b-veto event variables
  m_maxDxyEvt = 0.0;
  m_maxDszEvt = 0.0;
}

LeptonPlusFakeMLSelection_fullEE::~LeptonPlusFakeMLSelection_fullEE(){

  for(int theChannel=0; theChannel<1; theChannel++) {  
    delete m_p4LeptonPlus[theChannel];
    delete m_p4LeptonMinus[theChannel];
  }
  delete m_p3PFMET;
  delete m_p3TKMET;

  delete _selectionEE1;
  delete _selectionErrEE1;
  
  myOutTreeEE1 -> save();
}

void LeptonPlusFakeMLSelection_fullEE::initialiseFakeRate() {

  // binning                                      
  m_minFakePt[0] = 10.;   m_maxFakePt[0] = 15.;
  m_minFakePt[1] = 15.;   m_maxFakePt[1] = 20.;
  m_minFakePt[2] = 20.;   m_maxFakePt[2] = 25.;
  m_minFakePt[3] = 25.;   m_maxFakePt[3] = 50.;
  m_minFakePt[4] = 50.;   m_maxFakePt[4] = 10000.;

  /*
  // fake in the barrel from QCD MC                                  
  m_fakeRateEB[0] = 0.31852;
  m_fakeRateEB[1] = 0.0815463;
  m_fakeRateEB[2] = 0.0707214;
  m_fakeRateEB[3] = 0.0588124;
  m_fakeRateEB[4] = 0.0458273;

  m_fakeRateEB_err[0] = 0.00598496;
  m_fakeRateEB_err[1] = 0.00467256;
  m_fakeRateEB_err[2] = 0.00381354;
  m_fakeRateEB_err[3] = 0.00335443;
  m_fakeRateEB_err[4] = 0.0120015;

  // fake in the endcap from QCD MC                                                                     
  m_fakeRateEE[0] = 0.208365;
  m_fakeRateEE[1] = 0.0621581;
  m_fakeRateEE[2] = 0.0546968; 
  m_fakeRateEE[3] = 0.0522816; 
  m_fakeRateEE[4] = 0.0515749;

  m_fakeRateEE_err[0] = 0.00714442;
  m_fakeRateEE_err[1] = 0.00411465;
  m_fakeRateEE_err[2] = 0.00314477;
  m_fakeRateEE_err[3] = 0.00299393;
  m_fakeRateEE_err[4] = 0.0142672;
  */

  // fake in the barrel from data (jet:30)) - with EWK subtraction (apart from bin0, where no EWK removal)
  m_fakeRateEB[0] = 0.170909;
  m_fakeRateEB[1] = 0.126939;
  m_fakeRateEB[2] = 0.0829306;
  m_fakeRateEB[3] = 0.0804465;
  m_fakeRateEB[4] = 0.0993018;

  m_fakeRateEB_err[0] = 0.00926705;
  m_fakeRateEB_err[1] = 0.00916318;
  m_fakeRateEB_err[2] = 0.00703008;
  m_fakeRateEB_err[3] = 0.00687927;
  m_fakeRateEB_err[4] = 0.0285547;

  // fake in the endcap from data (jet:30) - with EWK subtraction (apart from bin0, where no EWK removal)
  m_fakeRateEE[0] = 0.136997;
  m_fakeRateEE[1] = 0.0628376;
  m_fakeRateEE[2] = 0.0559572;
  m_fakeRateEE[3] = 0.0441833;
  m_fakeRateEE[4] = 0.0255683;

  m_fakeRateEE_err[0] = 0.009566;
  m_fakeRateEE_err[1] = 0.00691408;
  m_fakeRateEE_err[2] = 0.00498304;
  m_fakeRateEE_err[3] = 0.00418404;
  m_fakeRateEE_err[4] = 0.01372;

  /*
  // fake in the barrel from data (jet:15)) - no EWK subtraction
  m_fakeRateEB[0] = 0.176334; 
  m_fakeRateEB[1] = 0.13883;
  m_fakeRateEB[2] = 0.105588;
  m_fakeRateEB[3] = 0.0916634;
  m_fakeRateEB[4] = 0.114754;
  
  m_fakeRateEB_err[0] = 0.00382773;
  m_fakeRateEB_err[1] = 0.00532204;
  m_fakeRateEB_err[2] = 0.00503719;
  m_fakeRateEB_err[3] = 0.00573555;
  m_fakeRateEB_err[4] = 0.028856;

  // fake in the endcap from data (jet:15)
  m_fakeRateEE[0] = 0.138302;
  m_fakeRateEE[1] = 0.0695556;
  m_fakeRateEE[2] = 0.0601607;
  m_fakeRateEE[3] = 0.0486044;
  m_fakeRateEE[4] = 0.0357143;

  m_fakeRateEE_err[0] = 0.00379404;
  m_fakeRateEE_err[1] = 0.00379232;
  m_fakeRateEE_err[2] = 0.00310862;
  m_fakeRateEE_err[3] = 0.00333565;
  m_fakeRateEE_err[4] = 0.0156841;
  */

  /*
  // fake in the barrel from data (jet:50))
  m_fakeRateEB[0] = 0.130682;
  m_fakeRateEB[1] = 0.101523;
  m_fakeRateEB[2] = 0.0920502;
  m_fakeRateEB[3] = 0.0910973;
  m_fakeRateEB[4] = 0.0804598;
  
  m_fakeRateEB_err[0] = 0.0254063;
  m_fakeRateEB_err[1] = 0.021518;
  m_fakeRateEB_err[2] = 0.0187001;
  m_fakeRateEB_err[3] = 0.013093;
  m_fakeRateEB_err[4] = 0.0291619;

  // fake in the endcap from data (jet:50)
  m_fakeRateEE[0] = 0.104575;
  m_fakeRateEE[1] = 0.0753425;
  m_fakeRateEE[2] = 0.0536913;
  m_fakeRateEE[3] = 0.0351438;
  m_fakeRateEE[4] = 0.0222222;

  m_fakeRateEE_err[0] = 0.0247391;
  m_fakeRateEE_err[1] = 0.0218441;
  m_fakeRateEE_err[2] = 0.0130575;
  m_fakeRateEE_err[3] = 0.00735984;
  m_fakeRateEE_err[4] = 0.0155379;
  */
}

void LeptonPlusFakeMLSelection_fullEE::initialisePromptRate() {

  // binning                                      
  m_minPromptPt[0] = 10.;   m_maxPromptPt[0] = 15.;
  m_minPromptPt[1] = 15.;   m_maxPromptPt[1] = 20.;
  m_minPromptPt[2] = 20.;   m_maxPromptPt[2] = 25.;
  m_minPromptPt[3] = 25.;   m_maxPromptPt[3] = 50.;
  m_minPromptPt[4] = 50.;   m_maxPromptPt[4] = 10000.;

  // prompt in the barrel
  m_promptRateEB[0] = 0.917;
  m_promptRateEB[1] = 0.887;
  m_promptRateEB[2] = 0.911;
  m_promptRateEB[3] = 0.94380;
  m_promptRateEB[4] = 0.962;
  
  m_promptRateEB_err[0] = 0.018;
  m_promptRateEB_err[1] = 0.009;
  m_promptRateEB_err[2] = 0.005;
  m_promptRateEB_err[3] = 0.00003;
  m_promptRateEB_err[4] = 0.002;

  // prompt in the endcap
  m_promptRateEE[0] = 0.851;
  m_promptRateEE[1] = 0.817;
  m_promptRateEE[2] = 0.8546;
  m_promptRateEE[3] = 0.9114;
  m_promptRateEE[4] = 0.9482;

  m_promptRateEE_err[0] = 0.021;
  m_promptRateEE_err[1] = 0.011;
  m_promptRateEE_err[2] = 0.0077;
  m_promptRateEE_err[3] = 0.0009;
  m_promptRateEE_err[4] = 0.0041;
}

void LeptonPlusFakeMLSelection_fullEE::Loop() {

  _verbose=false;
  if(fChain == 0) return;

  // initializations
  initialiseFakeRate();
  initialisePromptRate();
  
  // kinematics reduced tree
  std::string reducedTreeNameEE1 = _datasetName+"-datasetEE.root";    
  myOutTreeEE1 = new RedHiggsTree(reducedTreeNameEE1.c_str());

  if ( _selectionEE1->getSwitch("isData")) myOutTreeEE1->addRunInfos();
  myOutTreeEE1->addMLVars();
  myOutTreeEE1->addLatinos();

  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  PUWeight* fPUWeight = new PUWeight();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    resetKinematicsStart();

    // weight for the PU observed in 2011 data
    float tmpWeight = 1.;
    if ( !_selectionEE1->getSwitch("isData") ) tmpWeight *= fPUWeight->GetWeight(nPU);  

    // dummy MC truth
    bool promptEE = true;

    // Good Run selection
    if (_selectionEE1->getSwitch("isData") && _selectionEE1->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_selectionEE1->getSwitch("isData") && _selectionEE1->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    
    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask();
    bool passedHLT[1];
    passedHLT[ee] = hasPassedHLT();

    // ----------------------------------------------------------------------------
    // get the best electrons and best muons ==> tu be used to select ALL the possible channels at the beginning only
    std::pair<int,int> thePreElectrons = getBestElectronPair_acceptance();
    thePreElectron  = thePreElectrons.second;
    thePrePositron  = thePreElectrons.first;

    // reconstructed channel
    m_channel[ee] = false;     

    // at this level the SELECTED channel should have pT > 10 and > 20. So far, at least 2 leptons with pT >20 and 10 in the event
    if (thePreElectron > -1 && thePrePositron > -1) {
      float thisMaxPt = TMath::Max(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      float thisMinPt = TMath::Min(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      if (thisMaxPt>20 && thisMinPt>15) m_channel[ee] = true;    // fixme: hardcoded
    }

    if (_verbose) {
      std::cout << "nEle = "   << nEle << "\tnMuon = " << nMuon << std::endl;
      std::cout << "indices: " << thePreElectron << " " << thePrePositron << std::endl;
      std::cout << "chargeEle = " << chargeEle[thePreElectron] << "\tchargePos = " << chargeEle[thePrePositron] << std::endl;
      std::cout << "ee = " << m_channel[ee] << std::endl;
    }
    


    // -------------------------------------------------------------
    // EE candidates: preparing vectors of candidates and selecting the two highest pT ele- and ele+ after each step - to check the 20 GeV cut after 

    // eleID, for electrons in acceptance
    std::pair<int,int> theBestIdEle = getBestElectronPair_id(_acceptEleAll);   

    // isolation, for identified electrons
    std::pair<int,int> theBestIsolEle = getBestElectronPair_isol(_idEleAll); 

    // conversion rejection, for isolated electrons
    std::pair<int,int> theBestConvEle = getBestElectronPair_conv(_isolEleAll);     

    // transverse impact parameter, for electrons passing conversion rejection
    std::pair<int,int> theBestIpEle = getBestElectronPair_ip(_convEleAll);     

    // electrons passing the denominator definition
    std::pair<int,int> theBestDenomEle = getBestElectronPair_denominator();

    // the two highest pT electrons at this point are those I use for my analysis since they passed the full lepton selection
    int thePositron = theBestIpEle.first;    
    int theElectron = theBestIpEle.second;

    // 0 - 1 - 2 tight candidates
    bool is0tight = false;
    bool is1tight = false;
    bool is2tight = false;

    // here I have two candidates with opposite sign passing the tight selection: is N2T
    if (thePositron>-1 && theElectron>-1) {  
      float ptPositron = GetPt(pxEle[thePositron],pyEle[thePositron]);
      float ptElectron = GetPt(pxEle[theElectron],pyEle[theElectron]);
      if (ptPositron>ptElectron) { 
	theReal = thePositron;    
	theFake = theElectron;
      }
      else { 
	theReal = theElectron;
	theFake = thePositron;
      }
      is2tight = true;
    } 

    // here I have only 1 candidate passing the tight selection: is N1T
    if (thePositron>-1 && theElectron<0) { 
      theReal = thePositron;
      theFake = getBestDenominator(theReal);      
      is1tight = true;
    }
    if (theElectron>-1 && thePositron<0) { 
      theReal = theElectron;
      theFake = getBestDenominator(theReal);      
      is1tight = true;
    }

    // here I have zero candidates passing the tight selection: is N0T
    if (thePositron<0 && theElectron<0) { 
      int theDenomPlus  = theBestDenomEle.first;    
      int theDenomMinus = theBestDenomEle.second;
      float ptPlus  = GetPt(pxEle[theDenomPlus], pyEle[theDenomPlus]);
      float ptMinus = GetPt(pxEle[theDenomMinus],pyEle[theDenomMinus]);
      if (ptPlus>ptMinus) { 
	theReal = theDenomPlus;
	theFake = theDenomMinus;
      }
      else { 
	theReal = theDenomMinus;
	theFake = theDenomPlus;
      }
      is0tight = true;
    }    
    
    // sanity check
    if ( (is0tight && is1tight) || (is0tight && is2tight) || (is1tight && is2tight) ) cout << "questo non puo' succedere mai" << endl;

    // set of kinematics: : now I've all the final leptons 
    resetKinematics();

    // MET is an event variable. Independent o the channel        
    m_p3PFMET->SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);
    m_p3TKMET->SetXYZ(pxChMetPV[0],pyChMetPV[0],pzChMetPV[0]); // the one associated to the 0th vertex     
    m_theMET = m_p3PFMET->Pt();

    setKinematicsEE(theReal, theFake);
        

    // weight with the Fake / Prompt -> L2 probability	
    // I apply both fake rate and prompt rate to my "fake" candidate
    float theFakePt    = GetPt(pxEle[theFake],pyEle[theFake]);
    bool  isFakeBarrel = false;
    if ( fabs(etaEle[theFake])<1.476 ) isFakeBarrel   = true;

    // do both F-P and F-F analysis
    float weightFP      = 1.;
    float weightErrorFP = 1.;
    float weightFF      = 1.;
    float weightErrorFF = 1.;

    if ( theFake>-1 && theReal>-1) {
      float fakerate      = getFakeRate( theFakePt, isFakeBarrel );
      float fakerateErr   = getFakeRateError( theFakePt, isFakeBarrel );
      float promptrate    = getPromptRate( theFakePt, isFakeBarrel );
      float promptrateErr = getPromptRateError( theFakePt, isFakeBarrel );

      float thisPartWeightFP    = 1.;
      float thisPartWeightErrFP = 1.;
      float thisPartWeightFF    = 1.;
      float thisPartWeightErrFF = 1.;

      if (is0tight) {
	// for F-P
	thisPartWeightFP    = (fakerate*promptrate) * ( ( -2.*fakerate*promptrate) / ( (promptrate - fakerate)*(promptrate - fakerate) ) );
	float d1F = (4*fakerate*pow(promptrate,3)) / ( pow((fakerate-promptrate),3) );
	float d1P = (4*promptrate*pow(fakerate,3)) / ( pow((fakerate-promptrate),3) );
	thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);

	// for F-F
	thisPartWeightFF    = (fakerate*fakerate) * ( (promptrate*promptrate) / ( (promptrate - fakerate)*(promptrate - fakerate) ) );
	float d2F = -(2*fakerate*pow(promptrate,3)) / ( pow((fakerate-promptrate),3) );
	float d2P =  (2*promptrate*pow(fakerate,3)) / ( pow((fakerate-promptrate),3) );
	thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
      }

      if (is1tight) {
	// for F-P
	thisPartWeightFP    = (fakerate*promptrate) * (( fakerate*(1.-promptrate) + promptrate*(1.-fakerate) ) / ( (promptrate - fakerate)*(promptrate - fakerate) ));
	float d1F = ((promptrate*promptrate)*(-promptrate + fakerate*(-3+4*promptrate))) / ( pow((fakerate-promptrate),3) );
	float d1P = ((fakerate*fakerate)*(fakerate + promptrate*(3-4*fakerate))) / ( pow((fakerate-promptrate),3) );
	thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);

	// for F-F 
	thisPartWeightFF    = (fakerate*fakerate) * ( (-promptrate*(1.-promptrate)) / ( (promptrate - fakerate)*(promptrate - fakerate) ) );
	float d2F = (2*fakerate*(promptrate-1)*pow(promptrate,2))/( pow((fakerate-promptrate),3) );
	float d2P = (fakerate*fakerate*(fakerate+promptrate-2*fakerate*promptrate)) / ( pow((fakerate-promptrate),3) );
	thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
      }
      
      if (is2tight) {
	// for F-P
	thisPartWeightFP    = (fakerate*promptrate) * (( -2*(1.-promptrate)*(1.-fakerate) ) / ( (promptrate - fakerate)*(promptrate - fakerate) ));
	float d1P = -( 2*(fakerate-1)*fakerate*(-promptrate + fakerate*(-1 + 2*promptrate))) / ( pow((fakerate-promptrate),3) );
	float d1F =  ( 2*(promptrate-1)*promptrate*(-promptrate + fakerate*(-1 + 2*promptrate))) / ( pow((fakerate-promptrate),3) );
	thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);

	// for F-F
	thisPartWeightFF    = (fakerate*fakerate) * ( ( (1.-promptrate)*(1.-promptrate) ) / ( (promptrate - fakerate)*(promptrate - fakerate) ) );
	float d2F = -(2*fakerate*(promptrate-1)*(promptrate-1)*promptrate) / ( pow((fakerate-promptrate),3) );
	float d2P =  (2*(fakerate-1)*pow(fakerate,2)*(promptrate-1)) / ( pow((fakerate-promptrate),3) ); 
	thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
      }

      // for F-P
      weightFP      = tmpWeight * thisPartWeightFP;
      weightErrorFP = tmpWeight * thisPartWeightErrFP;   
      // for F-F
      weightFF      = tmpWeight * thisPartWeightFF;
      weightErrorFF = tmpWeight * thisPartWeightErrFF;   

    } else {
      
      weightFP      = tmpWeight;
      weightErrorFP = tmpWeight;
      weightFF      = tmpWeight;
      weightErrorFF = tmpWeight;
    }


    // -------------------------------------------------------------    
    // look for PV in the event (there is always at least 1 PV)
    m_closestPV = getPV();    // fixme: si chiama closest ma e' quello a piu' alto pT. 
    
    int njets[1], nuncorrjets[1];
    float dphiLLJ[1], btag[1];
    int nsoftmu[1],nextraleptons[1];
    for(int ichan=0; ichan<1; ichan++) {

      // jet counter
      njets[ichan] = numJets(eleCands[ichan],muCands[ichan],ichan);
      nuncorrjets[ichan] = numUncorrJets(eleCands[ichan],muCands[ichan]);

      // if 1-jet bin, use deltaphi(ll-jet)
      dphiLLJ[ichan] = deltaPhiLLJet(ichan);   
      
      // b veto
      btag[ichan] = bVetoJets(eleCands[ichan],muCands[ichan]);

      // soft muon counter
      nsoftmu[ichan] = numSoftMuons(muCands[ichan]);

      // extra lepton counter
      nextraleptons[ichan] = numExtraLeptons(eleCands[ichan],muCands[ichan]);
    }




    // -------------------------------------------------------
    // filling counter for P-F

    // EE
    // float forStatErr = weightFP*weightFP;                        // chiara: this two lines to compute the statistical error
    // CutBasedHiggsSelectionEE1.SetWeight(forStatErr);             // this two lines to compute the statistical error  
    CutBasedHiggsSelectionEE1.SetWeight(weightFP);                  // this to estimate FP
    // CutBasedHiggsSelectionEE1.SetWeight(weightFF);               // this to estimate FF
    CutBasedHiggsSelectionEE1.SetMcTruth(promptEE);  
    CutBasedHiggsSelectionEE1.SetHLT(passedHLT[ee]);               
    CutBasedHiggsSelectionEE1.SetIsChannel(m_channel[ee]);     
    CutBasedHiggsSelectionEE1.SetElectronId(theReal);
    CutBasedHiggsSelectionEE1.SetPositronId(theFake);
    CutBasedHiggsSelectionEE1.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionEE1.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionEE1.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionEE1.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionEE1.SetElectronIp(theReal);
    CutBasedHiggsSelectionEE1.SetPositronIp(theFake);
    // checking if the highest pT electron at each step has pT>20
    float thisMaxPtIpEE1 = TMath::Max(GetPt(pxEle[theReal],pyEle[theReal]),GetPt(pxEle[theFake],pyEle[theFake]));
    if (thisMaxPtIpEE1<20)   { 
      CutBasedHiggsSelectionEE1.SetElectronIp(-1);
      CutBasedHiggsSelectionEE1.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionEE1.SetHighElePt(hardestLeptonPt[ee]); 
    CutBasedHiggsSelectionEE1.SetLowElePt(slowestLeptonPt[ee]);  
    CutBasedHiggsSelectionEE1.SetExtraSlowLeptonPTCut(15.0); // enforce the min pT cut only on electrons 

    CutBasedHiggsSelectionEE1.SetNJets(njets[ee]);
    CutBasedHiggsSelectionEE1.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsSelectionEE1.SetBTagJets(btag[ee]);
    CutBasedHiggsSelectionEE1.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsSelectionEE1.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsSelectionEE1.SetMet(m_theMET);
    CutBasedHiggsSelectionEE1.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsSelectionEE1.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsSelectionEE1.SetDeltaPhiLLJet(dphiLLJ[ee]);   
    CutBasedHiggsSelectionEE1.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsSelectionEE1.SetInvMass(m_mll[ee]);
    CutBasedHiggsSelectionEE1.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsSelectionEE1.SetWWInvMass(2.*m_transvMass[ee]/_massVal);

    bool isSelectedEE1           = CutBasedHiggsSelectionEE1.output();    
    bool selUpToFinalLeptonsEE1  = CutBasedHiggsSelectionEE1.outputUpToFinalLeptons();
    bool selUpToJetVetoEE1       = CutBasedHiggsSelectionEE1.outputUpToJetVeto();
    bool selUpToUncorrJetVetoEE1 = CutBasedHiggsSelectionEE1.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiEE1       = CutBasedHiggsSelectionEE1.outputPreDeltaPhi();

    bool outputStep0  = CutBasedHiggsSelectionEE1.outputStep0();
    bool outputStep1  = CutBasedHiggsSelectionEE1.outputStep1();
    bool outputStep2  = CutBasedHiggsSelectionEE1.outputStep2();
    bool outputStep3  = CutBasedHiggsSelectionEE1.outputStep3();
    bool outputStep4  = CutBasedHiggsSelectionEE1.outputStep4();
    bool outputStep5  = CutBasedHiggsSelectionEE1.outputStep5();
    bool outputStep6  = CutBasedHiggsSelectionEE1.outputStep6();
    bool outputStep7  = CutBasedHiggsSelectionEE1.outputStep7();
    bool outputStep8  = CutBasedHiggsSelectionEE1.outputStep8();
    bool outputStep9  = CutBasedHiggsSelectionEE1.outputStep9();
    bool outputStep10 = CutBasedHiggsSelectionEE1.outputStep10();
    bool outputStep11 = CutBasedHiggsSelectionEE1.outputStep11();
    bool outputStep12 = CutBasedHiggsSelectionEE1.outputStep12();
    bool outputStep13 = CutBasedHiggsSelectionEE1.outputStep13();
    bool outputStep14 = CutBasedHiggsSelectionEE1.outputStep14();
    bool outputStep15 = CutBasedHiggsSelectionEE1.outputStep15();
    bool outputStep16 = CutBasedHiggsSelectionEE1.outputStep16();
    bool outputStep17 = CutBasedHiggsSelectionEE1.outputStep17();
    bool outputStep18 = CutBasedHiggsSelectionEE1.outputStep18();
    bool outputStep19 = CutBasedHiggsSelectionEE1.outputStep19();
    bool outputStep20 = CutBasedHiggsSelectionEE1.outputStep20();
    bool outputStep21 = CutBasedHiggsSelectionEE1.outputStep21();
    bool outputStep22 = CutBasedHiggsSelectionEE1.outputStep22();
    bool outputStep23 = CutBasedHiggsSelectionEE1.outputStep23();
    bool outputStep24 = CutBasedHiggsSelectionEE1.outputStep24();

    myOutTreeEE1->fillRunInfos(runNumber, lumiBlock, eventNumber, weightFP);          // chiara
    // myOutTreeEE1->fillRunInfos(runNumber, lumiBlock, eventNumber, weightFF);
    
    myOutTreeEE1 -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[ee], m_deltaPhi[ee], m_deltaErre[ee], m_transvMass[ee], m_mll[ee], 
			   hardestLeptonPt[ee], slowestLeptonPt[ee], hardestLeptonEta[ee], slowestLeptonEta[ee], 
			   m_deltaEtaLeptons[ee], nPV,
			   selUpToFinalLeptonsEE1, selUpToJetVetoEE1, selUpToUncorrJetVetoEE1, selPreDeltaPhiEE1, isSelectedEE1);

    myOutTreeEE1 -> fillMLVars(njets[ee], nuncorrjets[ee], m_maxDxyEvt, m_maxDszEvt, m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);

    myOutTreeEE1 -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 );


      
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeEE1 -> store();

    // for errors                                    
    CutBasedHiggsErrorsSelectionEE1.SetWeight(weightErrorFP);                // chiara: this to estimate FP
    // CutBasedHiggsErrorsSelectionEE1.SetWeight(weightErrorFF);             // chiara: this to estimate FF  
    CutBasedHiggsErrorsSelectionEE1.SetMcTruth(promptEE);
    CutBasedHiggsErrorsSelectionEE1.SetHLT(passedHLT[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetIsChannel(m_channel[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetElectronId(1);
    CutBasedHiggsErrorsSelectionEE1.SetPositronId(1);
    CutBasedHiggsErrorsSelectionEE1.SetElectronIsolation(1);
    CutBasedHiggsErrorsSelectionEE1.SetPositronIsolation(1);
    CutBasedHiggsErrorsSelectionEE1.SetElectronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE1.SetPositronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE1.SetElectronIp(theReal);
    CutBasedHiggsErrorsSelectionEE1.SetPositronIp(theFake);
    if (thisMaxPtIpEE1<20) {
      CutBasedHiggsErrorsSelectionEE1.SetElectronIp(-1);
      CutBasedHiggsErrorsSelectionEE1.SetPositronIp(-1);
    }
    CutBasedHiggsErrorsSelectionEE1.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetNJets(njets[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetBTagJets(btag[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetMet(m_theMET);
    CutBasedHiggsErrorsSelectionEE1.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetInvMass(m_mll[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsErrorsSelectionEE1.SetWWInvMass(2.*m_transvMass[ee]/_massVal);
    bool isSelectedErrorEE1 = CutBasedHiggsErrorsSelectionEE1.output();    
  }
}

void LeptonPlusFakeMLSelection_fullEE::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR PF SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionEE1.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR PF SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionEE1.displayEfficiencies(datasetName);

  // simple cuts based or like based ele id                                                                                         
  if (!_selectionEE1->getSwitch("asymmetricID")) {
    std::cout << "cut based symmetric ID: " << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  } else {
    std::cout << "cut based asymmetric ID: Low pT" << std::endl;
    EgammaCutBasedIDLow.displayEfficiencies();
    std::cout << "cut based asymmetric ID: High pT" << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  }
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_acceptance() {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _acceptEleAll.clear();

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    if(_selectionEE1->getSwitch("etaElectronAcc") && !_selectionEE1->passCut("etaElectronAcc",etaEle[i]) ) continue;

    if(_selectionEE1->getSwitch("ptElectronAcc") && !_selectionEE1->passCut("ptElectronAcc",thisPt) ) continue;
    
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }

    _acceptEleAll.push_back(i);   
  }
  _acceptEleAll = sortElectronsByPt(_acceptEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_id( std::vector<int> acceptEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _idEleAll.clear();

  for (int iEle=0; iEle<acceptEle.size(); iEle++) {
    int thisEle = acceptEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    float thisPt = GetPt(pxEle[thisEle],pyEle[thisEle]);
    if (!_selectionEE1->getSwitch("asymmetricID")) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionEE1->getSwitch("asymmetricID")) {
      if (thisPt>=20) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisPt<20)  isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }

    if (!theElectronID) continue;

    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _idEleAll.push_back(thisEle);  
  }
  _idEleAll = sortElectronsByPt(_idEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_isol( std::vector<int> idEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _isolEleAll.clear();

  for (int iEle=0; iEle<idEle.size(); iEle++) {
    int thisEle = idEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);

    if (!theElectronIsol) continue;
    
    float thisPt     = GetPt(pxEle[thisEle],pyEle[thisEle]);
    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _isolEleAll.push_back(thisEle);  
  }
  _isolEleAll = sortElectronsByPt(_isolEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_conv( std::vector<int> isolEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
  
  _convEleAll.clear();

  for (int iEle=0; iEle<isolEle.size(); iEle++) {
    int thisEle = isolEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    
    if (!theElectronConvRej) continue;

    float thisPt     = GetPt(pxEle[thisEle],pyEle[thisEle]);
    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _convEleAll.push_back(thisEle);      
  }
  _convEleAll = sortElectronsByPt(_convEleAll);

  return make_pair(theLep1,theLep2);
}


std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_ip( std::vector<int> convEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _ipEleAll.clear();

  for (int iEle=0; iEle<convEle.size(); iEle++) {
    int thisEle = convEle[iEle];

    int gsfTrack = gsfTrackIndexEle[thisEle]; 
    float d3dEle = impactPar3DGsfTrack[gsfTrack];
    if (_selectionEE1->getSwitch("electronIP") && (!_selectionEE1->passCut("electronIP",d3dEle)) ) continue;   

    float thisPt     = GetPt(pxEle[thisEle],pyEle[thisEle]);
    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _ipEleAll.push_back(thisEle);  
  }
  _ipEleAll = sortElectronsByPt(_ipEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_denominator() {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _denomEleAll.clear();  

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    bool isGoodDenom = isDenomFake(i);
    if (!isGoodDenom) continue;
    
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
    
    _denomEleAll.push_back(i);   
  }
  _denomEleAll = sortElectronsByPt(_denomEleAll);

  return make_pair(theLep1,theLep2);
}

int LeptonPlusFakeMLSelection_fullEE::getBestDenominator(int realEle) {
  
  int theFake=-1;
  float maxPtFake=-1000.;
  
  for(int iele=0; iele<nEle; iele++) {
    
    if (iele==realEle) continue;
    
    if (chargeEle[iele]*chargeEle[realEle]>0) continue;

    bool isGoodDenom = isDenomFake(iele);
    if (!isGoodDenom) continue;

    float thisElePt = GetPt(pxEle[iele],pyEle[iele]);

    // removing candidates passing the tight selection
    bool isTight = true;
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    if (!_selectionEE1->getSwitch("asymmetricID")) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionEE1->getSwitch("asymmetricID")) {
      if (thisElePt>=20) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisElePt<20)  isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }
    if (!theElectronID)      isTight = false;
    if (!theElectronIsol)    isTight = false; 
    if (!theElectronConvRej) isTight = false;
    int gsfTrack = gsfTrackIndexEle[iele]; 
    float d3dEle = impactPar3DGsfTrack[gsfTrack];
    if (_selectionEE1->getSwitch("electronIP") && (!_selectionEE1->passCut("electronIP",d3dEle)) ) isTight = false;
    if (isTight) continue;

    if( thisElePt > maxPtFake ) { maxPtFake = thisElePt; theFake=iele; }
  }

  return theFake;
}

bool LeptonPlusFakeMLSelection_fullEE::isDenomFake(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance	                                                                                                                        
  if( p3Ele.Eta() > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 15. )  isGoodDenom = false;    
  
  // only ecal driven	  
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  if(!ecalDriven) isGoodDenom = false;
  
  // barrel or endcap              
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);
  
  // isolation	  
  float combinedIso;
  if (isEleEB)  combinedIso = dr03TkSumPtEle[theEle] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0) + dr03HcalTowerSumEtFullConeEle[theEle];
  if (!isEleEB) combinedIso = dr03TkSumPtEle[theEle] + dr03EcalRecHitSumEtEle[theEle] + dr03HcalTowerSumEtFullConeEle[theEle];
  float corrCombinedIso = (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / p3Ele.Pt();
  if ( corrCombinedIso>0.15 ) isGoodDenom = false;
  
    // H/E                                                                                                                                     
  bool isBarrelEle;
  if ( fabs(etaEle[theEle]) <  1.479 ) isBarrelEle = true;
  if ( fabs(etaEle[theEle]) >= 1.479 ) isBarrelEle = false;
  if ( isBarrelEle && hOverEEle[theEle]>0.15) isGoodDenom = false;
  if (!isBarrelEle && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // sigmaIetaIeta	    
  bool isBarrelSc;
  int sc = superClusterIndexEle[theEle];
  if ( sc < 0 ) isGoodDenom = false;
  if ( fabs(etaSC[sc]) <  1.479 ) isBarrelSc = true;
  if ( fabs(etaSC[sc]) >= 1.479 ) isBarrelSc = false;
  if ( isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.014 ) isGoodDenom = false;
  if (!isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.035 ) isGoodDenom = false;
  
  // spikes                            
  float theE1 = eMaxSC[sc];
  float theE4SwissCross = e4SwissCrossSC[sc];
  float theSpikeSC = 1.0 - (theE4SwissCross/theE1);
  if (theSpikeSC>0.95) isGoodDenom = false;
  
  return isGoodDenom;
}

float LeptonPlusFakeMLSelection_fullEE::getFakeRate( float fakePt, bool isFakeBarrel ) {

  for (int theBin = 0; theBin<7; theBin++) {

    if( fakePt >= m_minFakePt[theBin] && fakePt < m_maxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_fakeRateEB[theBin];
      if (!isFakeBarrel) return m_fakeRateEE[theBin];
    }
  }

  std::cout << "BIG ERROR: fakePt = " << fakePt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullEE::getFakeRateError( float fakePt, bool isFakeBarrel ) {

  for (int theBin = 0; theBin < 7; theBin++) {
    if( fakePt >= m_minFakePt[theBin] && fakePt < m_maxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_fakeRateEB_err[theBin];
      if (!isFakeBarrel) return m_fakeRateEE_err[theBin];
    }
  }

  return -1.;
}

float LeptonPlusFakeMLSelection_fullEE::getPromptRate( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin<7; theBin++) {

    if( promptPt >= m_minPromptPt[theBin] && promptPt < m_maxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_promptRateEB[theBin];
      if (!isPromptBarrel) return m_promptRateEE[theBin];
    }
  }

  std::cout << "BIG ERROR: promptPt = " << promptPt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullEE::getPromptRateError( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin < 7; theBin++) {
    if( promptPt >= m_minPromptPt[theBin] && promptPt < m_maxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_promptRateEB_err[theBin];
      if (!isPromptBarrel) return m_promptRateEE_err[theBin];
    }
  }

  return -1.;
}

void LeptonPlusFakeMLSelection_fullEE::setKinematicsEE(int myReal, int myFake) {

  if (myFake > -1 && myReal > -1) {

    eleCands[ee].push_back(myReal);
    eleCands[ee].push_back(myFake);
    hardestLeptonPt[ee] = TMath::Max(GetPt(pxEle[myReal],pyEle[myReal]),GetPt(pxEle[myFake],pyEle[myFake]));
    slowestLeptonPt[ee] = TMath::Min(GetPt(pxEle[myReal],pyEle[myReal]),GetPt(pxEle[myFake],pyEle[myFake]));
    hardestLeptonEta[ee] = etaEle[myReal];
    slowestLeptonEta[ee] = etaEle[myFake];
    m_p4LeptonMinus[ee] -> SetXYZT(pxEle[myReal], pyEle[myReal], pzEle[myReal], energyEle[myReal]);
    m_p4LeptonPlus[ee]  -> SetXYZT(pxEle[myFake],pyEle[myFake],pzEle[myFake],energyEle[myFake]);
    m_mll[ee]       = (*(m_p4LeptonMinus[ee]) + *(m_p4LeptonPlus[ee])).M();
    m_deltaPhi[ee]  = fabs(180./TMath::Pi() * m_p4LeptonMinus[ee]->Vect().DeltaPhi(m_p4LeptonPlus[ee]->Vect()));
    m_deltaErre[ee] = m_p4LeptonMinus[ee]->Vect().DeltaR(m_p4LeptonPlus[ee]->Vect());
    m_deltaEtaLeptons[ee] = etaEle[myReal]-etaEle[myFake];
    m_dilepPt[ee].SetXYZ( m_p4LeptonMinus[ee]->Vect().X()+m_p4LeptonPlus[ee]->Vect().X(),m_p4LeptonMinus[ee]->Vect().Y()+m_p4LeptonPlus[ee]->Vect().Y(),0.0 );
    m_transvMass[ee]=CalcGammaMRstar(*m_p4LeptonMinus[ee],*m_p4LeptonPlus[ee]);
    m_metOptll[ee] = m_theMET / m_dilepPt[ee].Pt();
    m_mT2[ee] = 0.;
    // m_p4MET->SetXYZT(pxPFMet[0],pyPFMet[0],pzPFMet[0],energyPFMet[0]);
    m_projectedMet[ee] = GetProjectedMet(m_p4LeptonMinus[ee]->Vect(),m_p4LeptonPlus[ee]->Vect());
  }
  
}

void LeptonPlusFakeMLSelection_fullEE::resetKinematicsStart() {

  theReal         = -1;
  theFake         = -1;
  thePreElectron  = -1;
  thePrePositron  = -1;
}

void LeptonPlusFakeMLSelection_fullEE::resetKinematics() {

  for(int theChannel=0; theChannel<1; theChannel++) {
    eleCands[theChannel].clear();
    muCands[theChannel].clear();
    m_p4LeptonMinus[theChannel] -> SetXYZT(0,0,0,0);                                                        
    m_p4LeptonPlus[theChannel]  -> SetXYZT(0,0,0,0);
    m_p3PFMET                   -> SetXYZ(0,0,0);
    m_p3TKMET                   -> SetXYZ(0,0,0);
    hardestLeptonPt[theChannel]   = 0.;
    slowestLeptonPt[theChannel]   = 0.;
    hardestLeptonEta[theChannel]  = 0.;
    slowestLeptonEta[theChannel]  = 0.;
    m_mll[theChannel]             = 0.;
    m_deltaPhi[theChannel]        = 0.;
    m_deltaErre[theChannel]       = 0.;
    m_deltaEtaLeptons[theChannel] = 0.; 
    m_dilepPt[theChannel]         = 0.;
    m_transvMass[theChannel]      = 0.;
    m_metOptll[theChannel]        = 0.;
    m_mT2[theChannel]             = 0.;
    m_projectedMet[theChannel]    = 0.;
  }
}



void LeptonPlusFakeMLSelection_fullEE::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {
  
  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  float pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxSC[sc];
      e4SwissCross = e4SwissCrossSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagSC[sc];
      seedTime = timeSC[sc];
      seedChi2 = chi2SC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }


  thisCutBasedID->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  thisCutBasedID->SetRecoFlag(recoFlagsEle[eleIndex]);
  thisCutBasedID->applyElectronIDOnPFlowElectrons(true);
  thisCutBasedID->SetHOverE( HoE );
  thisCutBasedID->SetS9S25( s9s25 );
  thisCutBasedID->SetDEta( deta );
  thisCutBasedID->SetDPhiIn( dphiin );
  thisCutBasedID->SetDPhiOut( dphiout );
  thisCutBasedID->SetBremFraction( fbrem );
  thisCutBasedID->SetSigmaEtaEta( see );
  thisCutBasedID->SetSigmaPhiPhi( spp );
  thisCutBasedID->SetEOverPout( eopout );
  thisCutBasedID->SetEOverPin( eop );
  thisCutBasedID->SetElectronClass ( classificationEle[eleIndex] );
  thisCutBasedID->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  thisCutBasedID->SetLikelihood( likelihoodRatio(eleIndex,*LH) );
  thisCutBasedID->SetNBrem( nbremsEle[eleIndex] );
  thisCutBasedID->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                
  thisCutBasedID->SetTrkIsolation ( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                        
  thisCutBasedID->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );         
  float iso = 0.0;
  if ( anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB) ) iso = dr03TkSumPtEle[eleIndex] + max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else iso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  thisCutBasedID->SetCombinedIsolation( (iso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt );
  thisCutBasedID->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  thisCutBasedID->SetConvDist( fabs(convDistEle[eleIndex]) );
  thisCutBasedID->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
  thisCutBasedID->m_cleaner->SetE1(e1);
  thisCutBasedID->m_cleaner->SetE4SwissCross(e4SwissCross);
  thisCutBasedID->m_cleaner->SetFiducialFlag(fidFlagSC);
  thisCutBasedID->m_cleaner->SetSeedFlag(seedRecHitFlag);
  thisCutBasedID->m_cleaner->SetSeedTime(seedTime);
  thisCutBasedID->m_cleaner->SetSeedChi2(seedChi2);

  //  return egammaCutBasedID.output(); // class dependent result
  *eleIdOutput = thisCutBasedID->outputNoClassEleId();
  *isolOutput = thisCutBasedID->outputNoClassIso();
  *convRejOutput = thisCutBasedID->outputNoClassConv();
}

void LeptonPlusFakeMLSelection_fullEE::isMuonID(int muonIndex, bool *muonIdOutput) {

  *muonIdOutput = true;

  Utils anaUtils; 
  bool flagGlobalMu = false;
  if(anaUtils.muonIdVal(muonIdMuon[muonIndex],AllGlobalMuons)) {
    int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
    if(trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] < 10 && 
       trackValidHitsGlobalMuonTrack[globalMuonTrack] > 0 &&
       numberOfMatchesMuon[muonIndex] > 1 ) flagGlobalMu = true; // to be used when new trees are available
  }

  bool flagTrackerMu = false;
  if( (anaUtils.muonIdVal(muonIdMuon[muonIndex],AllTrackerMuons) &&
       anaUtils.muonIdVal(muonIdMuon[muonIndex],TMLastStationTight)) ) flagTrackerMu  = true;

  if(!(flagGlobalMu || flagTrackerMu)) {
    *muonIdOutput = false;
    return;
  }
    
  int track = trackIndexMuon[muonIndex];

  if(trackValidHitsTrack[track]<=10) *muonIdOutput = false;

  if( (numberOfValidPixelBarrelHitsTrack[track]+numberOfValidPixelEndcapHitsTrack[track])<1 ) *muonIdOutput = false; 

  float ptTrack = sqrt( pxTrack[track]*pxTrack[track] + pyTrack[track]*pyTrack[track] );
  float sign = fabs(ptErrorTrack[track]/ptTrack);
  if (sign>=0.1) *muonIdOutput = false;
}


int LeptonPlusFakeMLSelection_fullEE::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  m_goodJets.clear();
  float ETMax=0.;

  theLeadingJet[theChannel]=-1;   

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionEE1->getSwitch("etaJetAcc") && !_selectionEE1->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(_selectionEE1->getSwitch("etJetAcc") && !_selectionEE1->passCut("etJetAcc", pt)) continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    
    if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
                  chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    bool foundMatch = false;

    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {
      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    m_goodJets.push_back(j);
    num++;
    
    if(pt>ETMax) {
      theLeadingJet[theChannel] = j;
      ETMax = pt;
    }

  }

  return num;
}


int LeptonPlusFakeMLSelection_fullEE::numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  int num=0;

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    float uncorrEt = uncorrEnergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrEnergyAK5PFPUcorrJet[j]);
    TVector3 p3Jet = p4Jet.Vect();

    if(_selectionEE1->getSwitch("etaJetAcc")      && !_selectionEE1->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;    
    if(_selectionEE1->getSwitch("etUncorrJetAcc") && !_selectionEE1->passCut("etUncorrJetAcc", uncorrEt))   continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    
    if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
                  chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    bool foundMatch=false;
    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {
      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;
    
    num++;
  }
  
  return num;
}

float LeptonPlusFakeMLSelection_fullEE::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  float output=-999;
  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionEE1->getSwitch("etaJetAcc") && !_selectionEE1->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/energyAK5PFPUcorrJet[j];
    
    if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
                  chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    bool foundMatch=false;
    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {
      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE1->getSwitch("jetConeWidth") && _selectionEE1->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    float tmp = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];     
    if(tmp > output) output = tmp;
    
  }

  return output;

}

float LeptonPlusFakeMLSelection_fullEE::deltaPhiLLJet(int ichan) {   
  
  int myLeadingJet = theLeadingJet[ichan];

  if(myLeadingJet > -1) {
    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);    
    return fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));                           
  } else return -999.;
}

int LeptonPlusFakeMLSelection_fullEE::numSoftMuons(std::vector<int> muonToRemove) {

  int num = 0;
  for(int i=0; i<nMuon; ++i) {

    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) { 
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;

    float pt = GetPt(pxMuon[i],pyMuon[i]);
    if(pt < 3.0) continue;

    Utils anaUtils;
    if(!anaUtils.muonIdVal(muonIdMuon[i],AllTrackerMuons) || !anaUtils.muonIdVal(muonIdMuon[i],TMLastStationTight)) continue;
       
    int track = trackIndexMuon[i];
    if(trackValidHitsTrack[track]<=10) continue;

    float dxy = transvImpactParTrack[track];
    if(dxy > 0.100) continue;   

    float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    float isoSumRel = isoSumAbs / pt;
    if(pt>20 || isoSumRel<0.1) continue;
    
    num++;
  }
  return num;
}

int LeptonPlusFakeMLSelection_fullEE::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    
    bool isSelEle=false;
    for(int eleSel=0; eleSel<(int)eleToRemove.size(); eleSel++) {
      if(i==eleToRemove[eleSel]) isSelEle=true;
    }
    if(isSelEle) continue;

    if(_selectionEE1->getSwitch("etaElectronAcc") && !_selectionEE1->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_selectionEE1->getSwitch("ptElectronAcc")  && !_selectionEE1->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])) ) continue;

    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;
    if (!_selectionEE1->getSwitch("asymmetricID")) 
      isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
    if (_selectionEE1->getSwitch("asymmetricID")) {
      float pt = GetPt(pxEle[i],pyEle[i]);	
      if(pt>=20) isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
      if(pt<20)  isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedIDLow);
    }
    if(!theId || !theIso || !theConvRej) continue;

    int track = gsfTrackIndexEle[i];
    float d3dEle = impactPar3DGsfTrack[track];
    if (_selectionEE1->getSwitch("electronIP") && (!_selectionEE1->passCut("electronIP",d3dEle)) ) continue;    

    numEle++;
  }

  int numMu = 0;
  for(int i=0; i<nMuon; ++i) {
    
    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) {
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;
    
    float ptMu = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEE1->getSwitch("etaMuonAcc") && !_selectionEE1->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_selectionEE1->getSwitch("ptMuonAcc") && !_selectionEE1->passCut("ptMuonAcc",ptMu) ) continue;

    bool theId = true;
    isMuonID(i,&theId);
    if(!theId) continue;
    float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    float isoSumRel = isoSumAbs / ptMu;
    if(_selectionEE1->getSwitch("muGlobalIso") && !_selectionEE1->passCut("muGlobalIso",isoSumRel)) continue;

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = PVzPV[m_closestPV] - trackVzTrack[track];  
    if(_selectionEE1->getSwitch("muonIP") && !_selectionEE1->passCut("muonIP",dxy)) continue;
    if(_selectionEE1->getSwitch("muonDz") && !_selectionEE1->passCut("muonDz",dz))  continue;  

    numMu++;
  }
  
  return numEle + numMu;
}

int LeptonPlusFakeMLSelection_fullEE::getPV() {
  int hardestPV = -1;
  float sumPtMax = 0.0;
  for(int v=0; v<nPV; v++) {
    if(SumPtPV[v] > sumPtMax) {
      sumPtMax = SumPtPV[v];
      hardestPV = v;
    }
  }
  return hardestPV;
}

bool LeptonPlusFakeMLSelection_fullEE::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = p3Track.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(p3Track.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

float LeptonPlusFakeMLSelection_fullEE::GetProjectedMet(TVector3 p1, TVector3 p2) {

  // calculate with PF met  
  float projMET_pf = 0.0;
  float deltaPhi1_pf = fabs(p1.DeltaPhi(*m_p3PFMET));
  float deltaPhi2_pf = fabs(p2.DeltaPhi(*m_p3PFMET));
  float deltaphi_pf = TMath::Min(deltaPhi1_pf,deltaPhi2_pf);
  if(deltaphi_pf<TMath::Pi()/2.) projMET_pf = m_p3PFMET->Mag() * sin(deltaphi_pf);
  else projMET_pf = m_p3PFMET->Mag();

  // calculate with TKMET                                                
  float projMET_tk = 0.0;
  float deltaPhi1_tk = fabs(p1.DeltaPhi(*m_p3TKMET));
  float deltaPhi2_tk = fabs(p2.DeltaPhi(*m_p3TKMET));
  float deltaphi_tk = TMath::Min(deltaPhi1_tk,deltaPhi2_tk);
  if(deltaphi_tk<TMath::Pi()/2.) projMET_tk = m_p3TKMET->Mag() * sin(deltaphi_tk);
  else projMET_tk = m_p3TKMET->Mag();

  return TMath::Min(projMET_pf,projMET_tk);
}

/// specific for HWW that has multiple channels with different HLT requirements
bool LeptonPlusFakeMLSelection_fullEE::reloadTriggerMask() {

  std::vector<int> triggerMask;
  // load the triggers required for EE
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersEE.begin();fIter!=requiredTriggersEE.end();++fIter) {
    //      std::cout << "For EE required: " << *fIter << std::endl;
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
      // nameHLT[i] has ..._vXXX
      if(nameHLT->at(i).find(*fIter) != string::npos)
	{
	  triggerMask.push_back( indexHLT[i] ) ;
	  break;
	}
    }
  }
  m_requiredTriggersEE = triggerMask;
}

bool LeptonPlusFakeMLSelection_fullEE::hasPassedHLT(){
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggersEE, firedTrg);
  return true;
}

void LeptonPlusFakeMLSelection_fullEE::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggersEE=reqTriggers;
}

void LeptonPlusFakeMLSelection_fullEE::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers){
  notRequiredTriggersEE=reqTriggers;
}

