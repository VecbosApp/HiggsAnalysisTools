#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/LeptonPlusFakeMLSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <TTree.h>

using namespace bits;

LeptonPlusFakeMLSelection::LeptonPlusFakeMLSelection(TTree *tree) 
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
  CutBasedHiggsSelectionEE.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE"); 
  CutBasedHiggsErrorsSelectionEE.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION ERRORS EE");
  _selectionEE    = CutBasedHiggsSelectionEE.GetSelection();  
  _selectionErrEE = CutBasedHiggsErrorsSelectionEE.GetSelection();  

  //  extra selection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  _selectionEE->addCut("etaElectronAcc");    
  _selectionEE->addCut("ptElectronAcc");
  _selectionEE->addCut("etaMuonAcc");
  _selectionEE->addCut("ptMuonAcc");
  _selectionEE->addCut("etUncorrJetAcc");  
  _selectionEE->addSwitch("apply_kFactor");   
  _selectionEE->addSwitch("isData");
  _selectionEE->addSwitch("goodRunLS");
  _selectionEE->addSwitch("asymmetricID");
  _selectionEE->addStringParameter("electronIDType");
  _selectionEE->addStringParameter("electronIDTypeLow");  

  // cut based electron id or likelihood 
  TString selectionString(_selectionEE->getStringParameter("electronIDType"));
  if (!_selectionEE->getSwitch("asymmetricID")) 
    cout << "=== CONFIGURING " << selectionString << " CUT BASED SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/higgs/electronId/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/higgs/electronId/");
  
  if (_selectionEE->getSwitch("asymmetricID")) {
    TString selectionStringLow (_selectionEE->getStringParameter("electronIDTypeLow"));
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
  std::cout << "[GoodRunLS]::goodRunLS is " << _selectionEE->getSwitch("goodRunLS") << " isData is " <<  _selectionEE->getSwitch("isData") << std::endl;

  // To read good run list!
  if (_selectionEE->getSwitch("goodRunLS") && _selectionEE->getSwitch("isData")) {
    std::string goodRunJsonFile       = "config/json/goodCollisions2011.json";
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

LeptonPlusFakeMLSelection::~LeptonPlusFakeMLSelection(){

  for(int theChannel=0; theChannel<1; theChannel++) {  
    delete m_p4LeptonPlus[theChannel];
    delete m_p4LeptonMinus[theChannel];
  }
  delete m_p3PFMET;
  delete m_p3TKMET;

  delete _selectionEE;
  delete _selectionErrEE;
  
  myOutTreeEE -> save();
}

void LeptonPlusFakeMLSelection::initialiseFakeRate() {

  // binning                                      
  m_minFakePt[0] = 10.;   m_maxFakePt[0] = 15.;
  m_minFakePt[1] = 15.;   m_maxFakePt[1] = 20.;
  m_minFakePt[2] = 20.;   m_maxFakePt[2] = 25.;
  m_minFakePt[3] = 25.;   m_maxFakePt[3] = 50.;
  m_minFakePt[4] = 50.;   m_maxFakePt[4] = 10000.;

  /* first attempt - Matt's talk
  // fake in the barrel                                   
  m_fakeRateEB[0] = 0.149401;
  m_fakeRateEB[1] = 0.140707;
  m_fakeRateEB[2] = 0.105181;
  m_fakeRateEB[3] = 0.0995595;
  m_fakeRateEB[4] = 0.175;

  m_fakeRateEB_err[0] = 0.00813347;
  m_fakeRateEB_err[1] = 0.00889542;
  m_fakeRateEB_err[2] = 0.0085951;
  m_fakeRateEB_err[3] = 0.00888732;
  m_fakeRateEB_err[4] = 0.0424816;

  // fake in the endcap                                                                     
  m_fakeRateEE[0] = 0.0823442;
  m_fakeRateEE[1] = 0.0733788;
  m_fakeRateEE[2] = 0.0637136;
  m_fakeRateEE[3] = 0.0508143;
  m_fakeRateEE[4] = 0.0597015;

  m_fakeRateEE_err[0] = 0.00748706;
  m_fakeRateEE_err[1] = 0.00621909;
  m_fakeRateEE_err[2] = 0.00601647;
  m_fakeRateEE_err[3] = 0.0056055;
  m_fakeRateEE_err[4] = 0.028946;
  */

  // fake in the barrel from QCD MC                                  
  m_fakeRateEB[0] = 0.318519;
  m_fakeRateEB[1] = 0.093561;
  m_fakeRateEB[2] = 0.0653243; 
  m_fakeRateEB[3] = 0.0582813;
  m_fakeRateEB[4] = 0.0457173;

  m_fakeRateEB_err[0] = 0.00598496;
  m_fakeRateEB_err[1] = 0.003919;
  m_fakeRateEB_err[2] = 0.00357566;
  m_fakeRateEB_err[3] = 0.00339644;
  m_fakeRateEB_err[4] = 0.0121774;

  // fake in the endcap from QCD MC                                                                     
  m_fakeRateEE[0] = 0.208276;
  m_fakeRateEE[1] = 0.0627951;
  m_fakeRateEE[2] = 0.0546009;
  m_fakeRateEE[3] = 0.0521685;
  m_fakeRateEE[4] = 0.0524358;

  m_fakeRateEE_err[0] = 0.00714286;
  m_fakeRateEE_err[1] = 0.00346753;
  m_fakeRateEE_err[2] = 0.00320969;
  m_fakeRateEE_err[3] = 0.00307286;
  m_fakeRateEE_err[4] = 0.0147189;

  /*
  // fake in the barrel from data (jet:30))
  m_fakeRateEB[0] = 0.154532;
  m_fakeRateEB[1] = 0.136631;
  m_fakeRateEB[2] = 0.0817844;
  m_fakeRateEB[3] = 0.0799476;
  m_fakeRateEB[4] = 0.107692;

  m_fakeRateEB_err[0] = 0.0139332;
  m_fakeRateEB_err[1] = 0.0117874;
  m_fakeRateEB_err[2] = 0.00964652;
  m_fakeRateEB_err[3] = 0.00981853;
  m_fakeRateEB_err[4] = 0.0384497;

  // fake in the endcap from data (jet:30)
  m_fakeRateEE[0] = 0.119089;
  m_fakeRateEE[1] = 0.0482625; 
  m_fakeRateEE[2] = 0.074184;
  m_fakeRateEE[3] = 0.0465328;
  m_fakeRateEE[4] = 0.0566038;

  m_fakeRateEE_err[0] = 0.0135545;
  m_fakeRateEE_err[1] = 0.00665861;
  m_fakeRateEE_err[2] = 0.00824217;
  m_fakeRateEE_err[3] = 0.00636249;
  m_fakeRateEE_err[4] = 0.0317418;
  */

  /*
  // fake in the barrel from data (jet:15))
  m_fakeRateEB[0] = 0.166922;
  m_fakeRateEB[1] = 0.128733;
  m_fakeRateEB[2] = 0.102072;
  m_fakeRateEB[3] = 0.0806324;
  m_fakeRateEB[4] = 0.108108;

  m_fakeRateEB_err[0] = 0.00551924;
  m_fakeRateEB_err[1] = 0.00610026;
  m_fakeRateEB_err[2] = 0.00680535;
  m_fakeRateEB_err[3] = 0.00765516;
  m_fakeRateEB_err[4] = 0.0360969;

  // fake in the endcap from data (jet:15)
  m_fakeRateEE[0] = 0.125552;
  m_fakeRateEE[1] = 0.0627326;
  m_fakeRateEE[2] = 0.0664093;
  m_fakeRateEE[3] = 0.0513254;
  m_fakeRateEE[4] = 0.0508475;

  m_fakeRateEE_err[0] = 0.00550408;
  m_fakeRateEE_err[1] = 0.00395338;
  m_fakeRateEE_err[2] = 0.00489263;
  m_fakeRateEE_err[3] = 0.00524047;
  m_fakeRateEE_err[4] = 0.0286007;
  */

  /*
  // fake in the barrel from data (jet:50))
  m_fakeRateEB[0] = 0.133333;
  m_fakeRateEB[1] = 0.168;
  m_fakeRateEB[2] = 0.107914;
  m_fakeRateEB[3] = 0.056338;
  m_fakeRateEB[4] = 0.0740741;

  m_fakeRateEB_err[0] = 0.0392523;
  m_fakeRateEB_err[1] = 0.0334396;
  m_fakeRateEB_err[2] = 0.0263169;
  m_fakeRateEB_err[3] = 0.0157986;
  m_fakeRateEB_err[4] = 0.0356389;

  // fake in the endcap from data (jet:50)
  m_fakeRateEE[0] = 0.0491803;
  m_fakeRateEE[1] = 0.030303;
  m_fakeRateEE[2] = 0.0902778;
  m_fakeRateEE[3] = 0.029703;
  m_fakeRateEE[4] = 0.0588235;

  m_fakeRateEE_err[0] = 0.0276873;
  m_fakeRateEE_err[1] = 0.0149202;
  m_fakeRateEE_err[2] = 0.0238816;
  m_fakeRateEE_err[3] = 0.00975284;
  m_fakeRateEE_err[4] = 0.0403526;
  */
}

void LeptonPlusFakeMLSelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;

  initialiseFakeRate();
  
  // kinematics reduced tree
  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";
  myOutTreeEE = new RedHiggsTree(reducedTreeNameEE.c_str());

  if ( _selectionEE->getSwitch("isData")) myOutTreeEE->addRunInfos();
  myOutTreeEE->addMLVars();

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
    if ( !_selectionEE->getSwitch("isData") ) tmpWeight *= fPUWeight->GetWeight(nPU);  // chiara: on by default

    // dummy MC truth
    bool promptEE = true;

    // Good Run selection
    if (_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    
    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask();
    bool passedHLT[1];
    passedHLT[ee] = hasPassedHLT(ee);

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
      if (thisMaxPt>20) m_channel[ee] = true;    // fixme: hardcoded
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

    // the two highest pT electrons at this point are those I use for my analysis since the passed the full lepton selection
    int thePositron = theBestIpEle.first;    
    int theElectron = theBestIpEle.second;
    float ptPositron = GetPt(pxEle[thePositron],pyEle[thePositron]);
    float ptElectron = GetPt(pxEle[theElectron],pyEle[theElectron]);
    if (ptPositron>ptElectron) theReal = thePositron;
    else theReal = theElectron;

    // consider all possible denominators different from the selected tight candidate and use the highest pT one as a fake          
    theFake = getBestDenominator(theReal);    

    // set of kinematics: : now I've all the final leptons 
    resetKinematics();

    // MET is an event variable. Independent o the channel        
    m_p3PFMET->SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);
    m_p3TKMET->SetXYZ(pxChMetPV[0],pyChMetPV[0],pzChMetPV[0]); // the one associated to the 0th vertex     
    m_theMET = m_p3PFMET->Pt();

    setKinematicsEE(theReal, theFake);
    
    
    // weight with the Fake -> L2 probability	
    float weight       = 1.;
    float weightError  = 1.;
    float theFakePt    = GetPt(pxEle[theFake],pyEle[theFake]);
    bool  isFakeBarrel = false;
    if ( fabs(etaEle[theFake])<1.476 ) isFakeBarrel = true;
    if ( theFake>-1 ) {
      float fakerate    = getFakeRate( theFakePt, isFakeBarrel );
      float fakerateErr = getFakeRateError( theFakePt, isFakeBarrel );
      weight      = tmpWeight * fakerate;
      weightError = tmpWeight * fakerateErr;
    } else {
      weight      = tmpWeight;
      weightError = tmpWeight;
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


    // ---------------------------------------
    // filling counters for the different final states

    // EE
    CutBasedHiggsSelectionEE.SetWeight(weight);               
    CutBasedHiggsSelectionEE.SetMcTruth(promptEE);  
    CutBasedHiggsSelectionEE.SetHLT(passedHLT[ee]);               
    CutBasedHiggsSelectionEE.SetIsChannel(m_channel[ee]);     
    CutBasedHiggsSelectionEE.SetElectronId(theReal);
    CutBasedHiggsSelectionEE.SetPositronId(theFake);
    CutBasedHiggsSelectionEE.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionEE.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionEE.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionEE.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionEE.SetElectronIp(theReal);
    CutBasedHiggsSelectionEE.SetPositronIp(theFake);
    // checking if the highest pT electron at each step has pT>20
    float thisMaxPtIpEE = TMath::Max(GetPt(pxEle[theReal],pyEle[theReal]),GetPt(pxEle[theFake],pyEle[theFake]));
    if (thisMaxPtIpEE<20)   { 
      CutBasedHiggsSelectionEE.SetElectronIp(-1);
      CutBasedHiggsSelectionEE.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionEE.SetHighElePt(hardestLeptonPt[ee]); 
    CutBasedHiggsSelectionEE.SetLowElePt(slowestLeptonPt[ee]);  
    CutBasedHiggsSelectionEE.SetExtraSlowLeptonPTCut(15.0); // enforce the min pT cut only on electrons 

    CutBasedHiggsSelectionEE.SetNJets(njets[ee]);
    CutBasedHiggsSelectionEE.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsSelectionEE.SetBTagJets(btag[ee]);
    CutBasedHiggsSelectionEE.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsSelectionEE.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsSelectionEE.SetMet(m_theMET);
    CutBasedHiggsSelectionEE.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsSelectionEE.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsSelectionEE.SetDeltaPhiLLJet(dphiLLJ[ee]);   
    CutBasedHiggsSelectionEE.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsSelectionEE.SetInvMass(m_mll[ee]);
    CutBasedHiggsSelectionEE.SetDetaLeptons(m_deltaEtaLeptons[ee]);

    bool isSelectedEE           = CutBasedHiggsSelectionEE.output();    
    bool selUpToFinalLeptonsEE  = CutBasedHiggsSelectionEE.outputUpToFinalLeptons();
    bool selUpToJetVetoEE       = CutBasedHiggsSelectionEE.outputUpToJetVeto();
    bool selUpToUncorrJetVetoEE = CutBasedHiggsSelectionEE.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiEE       = CutBasedHiggsSelectionEE.outputPreDeltaPhi();
    bool outputStep1            = CutBasedHiggsSelectionEE.outputStep1();

    myOutTreeEE->fillRunInfos(runNumber, lumiBlock, eventNumber, weight);
    
    myOutTreeEE -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[ee], m_deltaPhi[ee], m_deltaErre[ee], m_transvMass[ee], m_mll[ee], 
			   hardestLeptonPt[ee], slowestLeptonPt[ee], m_deltaEtaLeptons[ee], nPV,
			   selUpToFinalLeptonsEE, selUpToJetVetoEE, selUpToUncorrJetVetoEE, selPreDeltaPhiEE, isSelectedEE);

    myOutTreeEE -> fillMLVars(njets[ee], nuncorrjets[ee], m_maxDxyEvt, m_maxDszEvt, m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);

      
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeEE -> store();

    // for errors                                    
    CutBasedHiggsErrorsSelectionEE.SetWeight(weightError);
    CutBasedHiggsErrorsSelectionEE.SetMcTruth(promptEE);
    CutBasedHiggsErrorsSelectionEE.SetHLT(passedHLT[ee]);
    CutBasedHiggsErrorsSelectionEE.SetIsChannel(m_channel[ee]);
    CutBasedHiggsErrorsSelectionEE.SetElectronId(1);
    CutBasedHiggsErrorsSelectionEE.SetPositronId(1);
    CutBasedHiggsErrorsSelectionEE.SetElectronIsolation(1);
    CutBasedHiggsErrorsSelectionEE.SetPositronIsolation(1);
    CutBasedHiggsErrorsSelectionEE.SetElectronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE.SetPositronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE.SetElectronIp(theReal);
    CutBasedHiggsErrorsSelectionEE.SetPositronIp(theFake);
    if (thisMaxPtIpEE<20) {
      CutBasedHiggsErrorsSelectionEE.SetElectronIp(-1);
      CutBasedHiggsErrorsSelectionEE.SetPositronIp(-1);
    }
    CutBasedHiggsErrorsSelectionEE.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE.SetNJets(njets[ee]);
    CutBasedHiggsErrorsSelectionEE.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsErrorsSelectionEE.SetBTagJets(btag[ee]);
    CutBasedHiggsErrorsSelectionEE.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsErrorsSelectionEE.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsErrorsSelectionEE.SetMet(m_theMET);
    CutBasedHiggsErrorsSelectionEE.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsErrorsSelectionEE.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsErrorsSelectionEE.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsErrorsSelectionEE.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsErrorsSelectionEE.SetInvMass(m_mll[ee]);
    CutBasedHiggsErrorsSelectionEE.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    bool isSelectedErrorEE = CutBasedHiggsErrorsSelectionEE.output();    
  }
}

void LeptonPlusFakeMLSelection::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR EE SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionEE.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR EE SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionEE.displayEfficiencies(datasetName);

  // simple cuts based or like based ele id                                                                                         
  if (!_selectionEE->getSwitch("asymmetricID")) {
    std::cout << "cut based symmetric ID: " << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  } else {
    std::cout << "cut based asymmetric ID: Low pT" << std::endl;
    EgammaCutBasedIDLow.displayEfficiencies();
    std::cout << "cut based asymmetric ID: High pT" << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  }
}

std::pair<int,int> LeptonPlusFakeMLSelection::getBestElectronPair_acceptance() {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _acceptEleAll.clear();

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    if(_selectionEE->getSwitch("etaElectronAcc") && !_selectionEE->passCut("etaElectronAcc",etaEle[i]) ) continue;

    if(_selectionEE->getSwitch("ptElectronAcc") && !_selectionEE->passCut("ptElectronAcc",thisPt) ) continue;
    
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }

    _acceptEleAll.push_back(i);   
  }
  _acceptEleAll = sortElectronsByPt(_acceptEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection::getBestElectronPair_id( std::vector<int> acceptEle ) {

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
    if (!_selectionEE->getSwitch("asymmetricID")) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionEE->getSwitch("asymmetricID")) {
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

std::pair<int,int> LeptonPlusFakeMLSelection::getBestElectronPair_isol( std::vector<int> idEle ) {

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

std::pair<int,int> LeptonPlusFakeMLSelection::getBestElectronPair_conv( std::vector<int> isolEle ) {

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


std::pair<int,int> LeptonPlusFakeMLSelection::getBestElectronPair_ip( std::vector<int> convEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _ipEleAll.clear();

  for (int iEle=0; iEle<convEle.size(); iEle++) {
    int thisEle = convEle[iEle];

    int gsfTrack = gsfTrackIndexEle[thisEle]; 
    float d3dEle = impactPar3DGsfTrack[gsfTrack];
    if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",d3dEle)) ) continue;   

    float thisPt     = GetPt(pxEle[thisEle],pyEle[thisEle]);
    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _ipEleAll.push_back(thisEle);  
  }
  _ipEleAll = sortElectronsByPt(_ipEleAll);

  return make_pair(theLep1,theLep2);
}

int LeptonPlusFakeMLSelection::getBestDenominator(int realEle) {
  
  int theFake=-1;
  float maxPtFake=-1000.;
  
  for(int iele=0; iele<nEle; iele++) {
    
    if (iele==realEle) continue;
    
    if (chargeEle[iele]*chargeEle[realEle]>0) continue;

    bool isGoodDenom = isDenomFake(iele);
    if (!isGoodDenom) continue;

    float thisElePt = GetPt(pxEle[iele],pyEle[iele]);
    if( thisElePt > maxPtFake ) { maxPtFake = thisElePt; theFake=iele; }
  }

  return theFake;
}

bool LeptonPlusFakeMLSelection::isDenomFake(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance	                                                                                                                                    
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;
  
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

float LeptonPlusFakeMLSelection::getFakeRate( float fakePt, bool isFakeBarrel ) {

  for (int theBin = 0; theBin<7; theBin++) {

    if( fakePt >= m_minFakePt[theBin] && fakePt < m_maxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_fakeRateEB[theBin];
      if (!isFakeBarrel) return m_fakeRateEE[theBin];
    }
  }

  std::cout << "BIG ERROR: fakePt = " << fakePt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection::getFakeRateError( float fakePt, bool isFakeBarrel ) {

  for (int theBin = 0; theBin < 7; theBin++) {
    if( fakePt >= m_minFakePt[theBin] && fakePt < m_maxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_fakeRateEB_err[theBin];
      if (!isFakeBarrel) return m_fakeRateEE_err[theBin];
    }
  }

  return -1.;
}

void LeptonPlusFakeMLSelection::setKinematicsEE(int myReal, int myFake) {

  if (myFake > -1 && myReal > -1) {

    eleCands[ee].push_back(myReal);
    eleCands[ee].push_back(myFake);
    hardestLeptonPt[ee] = TMath::Max(GetPt(pxEle[myReal],pyEle[myReal]),GetPt(pxEle[myFake],pyEle[myFake]));
    slowestLeptonPt[ee] = TMath::Min(GetPt(pxEle[myReal],pyEle[myReal]),GetPt(pxEle[myFake],pyEle[myFake]));
    m_p4LeptonMinus[ee] -> SetXYZT(pxEle[myReal], pyEle[myReal], pzEle[myReal], energyEle[myReal]);
    m_p4LeptonPlus[ee]  -> SetXYZT(pxEle[myFake],pyEle[myFake],pzEle[myFake],energyEle[myFake]);
    m_mll[ee]       = (*(m_p4LeptonMinus[ee]) + *(m_p4LeptonPlus[ee])).M();
    m_deltaPhi[ee]  = fabs(180./TMath::Pi() * m_p4LeptonMinus[ee]->Vect().DeltaPhi(m_p4LeptonPlus[ee]->Vect()));
    m_deltaErre[ee] = m_p4LeptonMinus[ee]->Vect().DeltaR(m_p4LeptonPlus[ee]->Vect());
    m_deltaEtaLeptons[ee] = etaEle[myReal]-etaEle[myFake];
    m_dilepPt[ee].SetXYZ( m_p4LeptonMinus[ee]->Vect().X()+m_p4LeptonPlus[ee]->Vect().X(),m_p4LeptonMinus[ee]->Vect().Y()+m_p4LeptonPlus[ee]->Vect().Y(),0.0 );
    // m_transvMass[ee]=mT3(*m_p4LeptonMinus[ee],*m_p4LeptonPlus[ee],m_p4MET->Vect());
    m_transvMass[ee] = 0.;
    m_metOptll[ee] = m_theMET / m_dilepPt[ee].Pt();
    m_mT2[ee] = 0.;
    // m_p4MET->SetXYZT(pxPFMet[0],pyPFMet[0],pzPFMet[0],energyPFMet[0]);
    m_projectedMet[ee] = GetProjectedMet(m_p4LeptonMinus[ee]->Vect(),m_p4LeptonPlus[ee]->Vect());
  }
  
}

void LeptonPlusFakeMLSelection::resetKinematicsStart() {

  theReal         = -1;
  theFake         = -1;
  thePreElectron  = -1;
  thePrePositron  = -1;
}

void LeptonPlusFakeMLSelection::resetKinematics() {

  for(int theChannel=0; theChannel<1; theChannel++) {
    eleCands[theChannel].clear();
    muCands[theChannel].clear();
    m_p4LeptonMinus[theChannel] -> SetXYZT(0,0,0,0);                                                        
    m_p4LeptonPlus[theChannel]  -> SetXYZT(0,0,0,0);
    m_p3PFMET                   -> SetXYZ(0,0,0);
    m_p3TKMET                   -> SetXYZ(0,0,0);
    hardestLeptonPt[theChannel]   = 0.;
    slowestLeptonPt[theChannel]   = 0.;
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



void LeptonPlusFakeMLSelection::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {
  
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

void LeptonPlusFakeMLSelection::isMuonID(int muonIndex, bool *muonIdOutput) {

  *muonIdOutput = true;

  Utils anaUtils; 
  bool flagGlobalMu = false;
  if(anaUtils.muonIdVal(muonIdMuon[muonIndex],AllGlobalMuons)) {
    int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
    if(trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] < 10 && 
       trackValidHitsGlobalMuonTrack[globalMuonTrack] > 0 ) flagGlobalMu = true;
       //       numberOfMatchesMuon[muonIndex] > 1 ) flagGlobalMu = true; // to be used when new trees are available
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


int LeptonPlusFakeMLSelection::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  m_goodJets.clear();
  float ETMax=0.;

  theLeadingJet[theChannel]=-1;   

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(_selectionEE->getSwitch("etJetAcc") && !_selectionEE->passCut("etJetAcc", pt)) continue;

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
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
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
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
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


int LeptonPlusFakeMLSelection::numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  int num=0;

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    float uncorrEt = uncorrEnergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrEnergyAK5PFPUcorrJet[j]);
    TVector3 p3Jet = p4Jet.Vect();

    if(_selectionEE->getSwitch("etaJetAcc")      && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;    
    if(_selectionEE->getSwitch("etUncorrJetAcc") && !_selectionEE->passCut("etUncorrJetAcc", uncorrEt))   continue;

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
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;
    
    num++;
  }
  
  return num;
}

float LeptonPlusFakeMLSelection::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  float output=-999;
  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

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
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    float tmp = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];     
    if(tmp > output) output = tmp;
    
  }

  return output;

}

float LeptonPlusFakeMLSelection::deltaPhiLLJet(int ichan) {   
  
  int myLeadingJet = theLeadingJet[ichan];

  if(myLeadingJet > -1) {
    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);    
    return fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));                           
  } else return -999.;
}

int LeptonPlusFakeMLSelection::numSoftMuons(std::vector<int> muonToRemove) {

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

int LeptonPlusFakeMLSelection::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    
    bool isSelEle=false;
    for(int eleSel=0; eleSel<(int)eleToRemove.size(); eleSel++) {
      if(i==eleToRemove[eleSel]) isSelEle=true;
    }
    if(isSelEle) continue;

    if(_selectionEE->getSwitch("etaElectronAcc") && !_selectionEE->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_selectionEE->getSwitch("ptElectronAcc")  && !_selectionEE->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])) ) continue;

    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;
    if (!_selectionEE->getSwitch("asymmetricID")) 
      isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
    if (_selectionEE->getSwitch("asymmetricID")) {
      float pt = GetPt(pxEle[i],pyEle[i]);	
      if(pt>=20) isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
      if(pt<20)  isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedIDLow);
    }
    if(!theId || !theIso || !theConvRej) continue;

    int track = gsfTrackIndexEle[i];
    float d3dEle = impactPar3DGsfTrack[track];
    if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",d3dEle)) ) continue;    

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
    if(_selectionEE->getSwitch("etaMuonAcc") && !_selectionEE->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_selectionEE->getSwitch("ptMuonAcc") && !_selectionEE->passCut("ptMuonAcc",ptMu) ) continue;

    bool theId = true;
    isMuonID(i,&theId);
    if(!theId) continue;
    float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    float isoSumRel = isoSumAbs / ptMu;
    if(_selectionEE->getSwitch("muGlobalIso") && !_selectionEE->passCut("muGlobalIso",isoSumRel)) continue;

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = PVzPV[m_closestPV] - trackVzTrack[track];  
    if(_selectionEE->getSwitch("muonIP") && !_selectionEE->passCut("muonIP",dxy)) continue;
    if(_selectionEE->getSwitch("muonDz") && !_selectionEE->passCut("muonDz",dz))  continue;  

    numMu++;
  }
  
  return numEle + numMu;
}

int LeptonPlusFakeMLSelection::getPV() {
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

bool LeptonPlusFakeMLSelection::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = p3Track.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(p3Track.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

float LeptonPlusFakeMLSelection::GetProjectedMet(TVector3 p1, TVector3 p2) {

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
bool LeptonPlusFakeMLSelection::reloadTriggerMask() {

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

bool LeptonPlusFakeMLSelection::hasPassedHLT(int channel) {
  Utils anaUtils;
  if(channel==ee) return anaUtils.getTriggersOR(m_requiredTriggersEE, firedTrg);
  return true;
}

void LeptonPlusFakeMLSelection::setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) requiredTriggersEE=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

void LeptonPlusFakeMLSelection::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) notRequiredTriggersEE=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

