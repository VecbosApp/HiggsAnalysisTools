#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"
#include "HiggsAnalysisTools/include/HiggsMLSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
//#include "Mt2/SUSYPhys_Mt2_222_Calculator.h"

#include <iostream>
#include <string>

#include <TTree.h>

using namespace bits;

// WARNING: TREES PRODUCED WITHOUT KFACTOR
float evtKfactor = 1.0;

HiggsMLSelection::HiggsMLSelection(TTree *tree) 
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
  
  std::string fileCutsPreselection     = higgsConfigDir + "2l2nuCutsPreselection.txt";
  std::string fileSwitchesPreselection = higgsConfigDir + "2l2nuSwitchesPreselection.txt";

  // preselection efficiencies
  CommonHiggsPreselection.Configure(fileCutsPreselection.c_str(), fileSwitchesPreselection.c_str()); 
  _preselection = CommonHiggsPreselection.GetSelection();

  //  extra preselection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  _preselection->addSwitch("apply_kFactor");   
  _preselection->addSwitch("isData");
  _preselection->addSwitch("goodRunLS");
  _preselection->addSwitch("asymmetricID");
  _preselection->addCut("etaElectronAcc");
  _preselection->addCut("ptElectronAcc");
  _preselection->addCut("etaMuonAcc");
  _preselection->addCut("ptMuonAcc");
  _preselection->addCut("etUncorrJetAcc");
  _preselection->addStringParameter("electronIDType");
  _preselection->addStringParameter("electronIDTypeLow");
  _preselection->summary();

  // selection efficiencies
  std::string fileCutsEE     = higgsConfigDirMass + "2e2nuCuts.txt";
  std::string fileSwitchesEE = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionEE.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE"); 
  CutBasedHiggsSelectionEE.AppendPreselection(_preselection);
  _selectionEE = CutBasedHiggsSelectionEE.GetSelection();  

  std::string fileCutsMM     = higgsConfigDirMass + "2mu2nuCuts.txt";
  std::string fileSwitchesMM = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionMM.Configure(fileCutsMM.c_str(),fileSwitchesMM.c_str(),"FULL SELECTION EVENT COUNTER MM"); 
  CutBasedHiggsSelectionMM.AppendPreselection(_preselection);
  _selectionMM = CutBasedHiggsSelectionMM.GetSelection();

  std::string fileCutsEM     = higgsConfigDirMass + "emu2nuCuts.txt";
  std::string fileSwitchesEM = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionEM.Configure(fileCutsEM.c_str(),fileSwitchesEM.c_str(),"FULL SELECTION EVENT COUNTER EM"); 
  CutBasedHiggsSelectionEM.AppendPreselection(_preselection);
  _selectionEM = CutBasedHiggsSelectionEM.GetSelection();
  
  // single electron efficiency
  //  EgammaCutBasedID.Configure("config/higgs"); // this is the class dependent e-ID

  TString selectionString(_preselection->getStringParameter("electronIDType"));
  if (!_preselection->getSwitch("asymmetricID")) 
    cout << "=== CONFIGURING " << selectionString << " SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/higgs/electronId/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/higgs/electronId/");

  if (_preselection->getSwitch("asymmetricID")) {
    TString selectionStringLow (_preselection->getStringParameter("electronIDTypeLow"));
    cout << "=== CONFIGURING "  << selectionStringLow << " and " 
	 << selectionString << " for ASYMMETRIC ELECTRON ID ===" << endl;
    EgammaCutBasedIDLow.ConfigureNoClass("config/higgs/electronId/"+selectionStringLow);
    EgammaCutBasedIDLow.ConfigureEcalCleaner("config/higgs/electronId/");
  }

  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EBlt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EBgt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  LH = new ElectronLikelihood(&(*EBlt15dir), &(*EElt15dir), &(*EBgt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);
  
  //Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _preselection->getSwitch("goodRunLS") << " isData is " <<  _preselection->getSwitch("isData") << std::endl;

  //To read good run list!
  if (_preselection->getSwitch("goodRunLS") && _preselection->getSwitch("isData")) {
    std::string goodRunJsonFile       = "config/json/goodRunLS.json";
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics
  m_p4ElectronPlus  = new TLorentzVector(0.,0.,0.,0.);
  m_p4ElectronMinus = new TLorentzVector(0.,0.,0.,0.);
  m_p4MuonPlus      = new TLorentzVector(0.,0.,0.,0.);
  m_p4MuonMinus     = new TLorentzVector(0.,0.,0.,0.);
  m_p4MET           = new TLorentzVector(0.,0.,0.,0.);

  // b-veto event variables
  m_maxDxyEvt = 0.0;
  m_maxDszEvt = 0.0;

  _bestElectrons = new std::vector<int>;
  _bestMuons     = new std::vector<int>;

  // histo to study jet/electron match
  H_deltaRuncorr = new TH1F("H_deltaRuncorr","uncorrected jets",100, 0.,2*TMath::Pi());
  H_deltaRcorr   = new TH1F("H_deltaRcorr",  "corrected jets",  100, 0.,2*TMath::Pi());
}

HiggsMLSelection::~HiggsMLSelection(){

  delete m_p4ElectronPlus;
  delete m_p4ElectronMinus;
  delete m_p4MuonPlus;
  delete m_p4MuonMinus;
  delete m_p4MET;
  delete _bestElectrons;
  delete _bestMuons;
  delete _preselection;
  delete _selectionEE;
  delete _selectionMM;
  delete _selectionEM;
  myOutTreeEE   -> save();
  myOutTreeMM   -> save();
  myOutTreeEM   -> save();
  //  myTriggerTree -> save();
  myEleIdTree   -> save();
}

bool HiggsMLSelection::findMcTree(const char* processType) {

  _process = "UNDEFINED";
  _theGenEle = -1;
  _theGenPos = -1;
  
  // now we look for ee || mumu || emu
  // in the acceptance and with a loose pT threshold
  float etaEleAcc_  = 2.5;
  float ptEleAcc_   = 5.0; // GeV
  float etaMuonAcc_ = 2.4;
  float ptMuonAcc_  = 0.0; // GeV
  
  // signal: 2e2nu
  if(strcmp(processType,"HtoWWto2e2nu")==0) {
    int indeminus=999, indeplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc] == -11 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc] ==  11 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    if( indeminus<25 && indeplus<25 ) {
      _theGenPos = indeplus;
      _theGenEle = indeminus;
    }
    return ( indeplus < 25 && indeminus < 25 );
  }

  // signal: 2m2nu
  if(strcmp(processType,"HtoWWto2m2nu")==0) {
    int indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc] == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc] ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
    }
    if( indmuminus<25 && indmuplus<25 ) {
      _theGenMuPlus  = indmuplus;
      _theGenMuMinus = indmuminus;
    }
    return ( indmuplus < 25 && indmuminus < 25 );
  }

  // signal: em2nu
  if(strcmp(processType,"HtoWWtoem2nu")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeplus = imc;
      if( idMc[imc]  ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  ==  11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeminus = imc;
    }
    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos = indeplus;
      _theGenMuMinus = indmuminus;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
    }
    return ( (indeplus<25 && indmuminus<25) || (indeminus<25 && indmuplus<25) );
  }

  // signal ee excluding taus
  if(strcmp(processType,"HtoWWto2e2nu_prompt")==0) {
    int indeminus=999, indeplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc]  == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    if( indeminus<25 && indeplus<25 ) {
      _theGenPos = indeplus;
      _theGenEle = indeminus;
    }
    return ( indeplus < 25 && indeminus < 25 );
  }

  // signal mm excluding taus
  if(strcmp(processType,"HtoWWto2m2nu_prompt")==0) {
    int indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = 11;
    }
    if( indmuminus<25 && indmuplus<25 ) {
      _theGenMuPlus = indmuplus;
      _theGenMuMinus = indmuminus;
    }
    return ( indmuplus < 25 && indmuminus < 25 );
  }

  // signal em excluding taus
  if(strcmp(processType,"HtoWWtoem2nu_prompt")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc]  == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos = indeplus;
      _theGenMuMinus = indmuminus;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
    }
    return ( (indeplus<25 && indmuminus<25) || (indeminus<25 && indmuplus<25) );
  }

  // signal: 2l2nu
  if(strcmp(processType,"HtoWWto2l2nu")==0) {
    int indlminus=999, indlplus=999;
    for(int imc=6;imc<25;imc++) {
      if(idMc[imc]>10 && idMc[imc]<19 && idMc[mothMc[imc]]==-24)  indlminus=imc;
      if(idMc[imc]<-10 && idMc[imc]>-19 && idMc[mothMc[imc]]==24) indlplus=imc;
    }
    if(indlminus<25 && indlplus<25) {
      if( (idMc[indlminus]==-11) || (idMc[indlminus]==-13) || (idMc[indlminus]==-15))  
	_theGenEle = indlminus;
      if( (idMc[indlplus]==11) || (idMc[indlplus]==13) || (idMc[indlplus]==15) )
	_theGenPos = indlplus;
    }
    return (indlminus<25 && indlplus<25);
  }
  

  // WW: e / mu / tau
  else if(strcmp(processType,"WW")==0) {
    _process = "WW";
    TVector3 WminusP, WplusP;
    WminusP.SetMagThetaPhi(pMc[6],thetaMc[6],phiMc[6]);
    WplusP.SetMagThetaPhi(pMc[7],thetaMc[7],phiMc[7]);
    float pT = (WminusP+WplusP).Pt();
    _theGenEle = 6;
    _theGenPos = 7;
    _genHiggsPt[0] = pT;
    return (
	    (abs(idMc[6])==24) && (abs(idMc[7])==24) &&
	    (abs(idMc[8])>10 && abs(idMc[8])<19 && abs(idMc[mothMc[8]])==24) &&
	    (abs(idMc[10])>10 && abs(idMc[10])<19 && abs(idMc[mothMc[10]])==24)
	    );
  }
  // w+jets: e / mu / tau
  else if(strcmp(processType,"Wjets")==0) {
    _process = "Wjets";
    return ( ((abs(idMc[8])==11) && abs(idMc[9])==12) || ((abs(idMc[8])==13) && abs(idMc[9])==14) || ((abs(idMc[8])==15) && abs(idMc[9])==16));
  }
  // ttbar: e / mu / tau
  else if(strcmp(processType,"ttbar")==0) {
    _process = "ttbar";
    _theGenEle = 13;
    _theGenPos = 15;
    return ( 
	    abs(idMc[9])==24 && abs(idMc[15])>10 && abs(idMc[15])<19 &&
	    abs(idMc[11])==24 && abs(idMc[13])>10 && abs(idMc[13])<19 &&
	    (idMc[13]*idMc[15]<0)
	    );
  }
  else if(strcmp(processType,"ZZleptonic")==0) {
    _process = "ZZleptonic";
    // 8,9; 10,11 are the daughters of the Z;Z
    return (fabs(idMc[8])>10 && fabs(idMc[8])<19 &&
	    fabs(idMc[9])>10 && fabs(idMc[9])<19 &&
	    fabs(idMc[10])>10 && fabs(idMc[10])<19 &&
	    fabs(idMc[11])>10 && fabs(idMc[11])<19);
  }
  else {
    std::cout << "This processType: " << processType << "is not expected, you should put MTtruth switch off" << std::endl;
    return false;
  }
}

float HiggsMLSelection::getkFactor(std::string process) {

  float weight = 1.;
  if((process.compare("Higgs")==0)) {
    weight = evtKfactor;
  }
  else if(process.compare("WW")==0) {
    weight = 1.0; // we used MC @ NLO weight in 16X   
  }
  return weight;
}

void HiggsMLSelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;
  
  // kinematics reduced tree
  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";
  std::string reducedTreeNameMM = _datasetName+"-datasetMM.root";
  std::string reducedTreeNameEM = _datasetName+"-datasetEM.root";
  myOutTreeEE = new RedHiggsTree(reducedTreeNameEE.c_str());
  myOutTreeMM = new RedHiggsTree(reducedTreeNameMM.c_str());
  myOutTreeEM = new RedHiggsTree(reducedTreeNameEM.c_str());

  //  myOutTreeEE->addHLTElectronsInfos();
  //  myOutTreeMM->addHLTMuonsInfos();
  //  myOutTreeEM->addHLTElectronsInfos(); myOutTreeEM->addHLTMuonsInfos();

  if ( !_preselection->getSwitch("isData") && _preselection->getSwitch("apply_kFactor") ) {
    myOutTreeEE->addKFactor();
    myOutTreeMM->addKFactor();
    myOutTreeEM->addKFactor();
  }

  if(!_preselection->getSwitch("isData")) {
    myOutTreeEE->addMcTruthInfos();
    myOutTreeMM->addMcTruthInfos();
    myOutTreeEM->addMcTruthInfos();
  } else {
    myOutTreeEE->addRunInfos();
    myOutTreeMM->addRunInfos();
    myOutTreeEM->addRunInfos();
  }

  myOutTreeEE->addMLVars();
  myOutTreeMM->addMLVars();
  myOutTreeEM->addMLVars();

  myOutTreeEE->addElectronInfos();
  myOutTreeEM->addElectronInfos();

  // trigger reduced tree
  //  std::string reducedTriggerTreeName = _datasetName+"-trigger.root";
  //  myTriggerTree = new RedTriggerTree(reducedTriggerTreeName.c_str());

  // eleId reduced tree
  std::string reducedEleIdTreeName = _datasetName+"-eleId.root";
  myEleIdTree = new RedEleIDTree(reducedEleIdTreeName.c_str());

  float met, deltaPhi, transvMass, deltaErre; 
  float dileptonInvMass, maxPtEle, minPtEle, detaLeptons;

  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    resetKinematics();

    // get the kFactor of the event (for signal)
    float weight = 1;
    if (!_preselection->getSwitch("isData") && _preselection->getSwitch("apply_kFactor")) weight = getkFactor("Higgs");
    
    // look to the MC truth decay tree
    bool foundMcTree = false;
    if( !_preselection->getSwitch("isData") ) foundMcTree = findMcTree("HtoWWto2e2nu");

    //    bool decayEE = findMcTree("HtoWWto2e2nu");
    //    bool decayMM = findMcTree("HtoWWto2m2nu");
    //    bool decayEM = findMcTree("HtoWWtoem2nu");

    bool promptEE, promptMM, promptEM;
    promptEE = promptMM = promptEM = false;
    if( !_preselection->getSwitch("isData") ) {
      promptEE = findMcTree("HtoWWto2e2nu_prompt");
      promptMM = findMcTree("HtoWWto2m2nu_prompt");
      promptEM = findMcTree("HtoWWtoem2nu_prompt");
    }

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    bool newTriggerMask = false;
    if (_preselection->getSwitch("isData")) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    //Good Run selection
    if (_preselection->getSwitch("isData") && _preselection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_preselection->getSwitch("isData") && _preselection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // trigger
    bool passedHLT = hasPassedHLT();

    //    myTriggerTree->fillMcTruth(decayEE,decayMM,decayEM,promptEE,promptMM,promptEM);
    //    myTriggerTree->fillHLTElectrons( firedTrg[m_requiredTriggers[0]] );
    //    myTriggerTree->fillHLTMuons( firedTrg[m_requiredTriggers[1]] );
    //    myTriggerTree->store();

    // get the best electrons, best muons  
    std::pair<int,int> theElectrons = getBestElectronPair();
    std::pair<int,int> theMuons = getBestMuonPair();
    int tbElectron(theElectrons.second), tbPositron(theElectrons.first);    
    int tbMuonPlus(theMuons.first),      tbMuonMinus(theMuons.second);
    theElectron  = tbElectron;
    thePositron  = tbPositron;
    theMuonPlus  = tbMuonPlus;
    theMuonMinus = tbMuonMinus;
    
    // how many electrons after preselection eleID and isolation   
    // (in case we want to reproduce CSA07 conditions)
    int nPreselEle = nEle;
    
    // reconstructed channel
    m_channel[ee] = false;     
    m_channel[mm] = false;
    m_channel[em] = false;
    if (theElectron > -1 && thePositron > -1)  m_channel[ee] = true;
    if (theMuonPlus > -1 && theMuonMinus > -1) m_channel[mm] = true;
    if ((theElectron > -1 && theMuonPlus > -1) || (thePositron > -1 && theMuonMinus > -1)) m_channel[em] = true;
    if (_verbose) {
      std::cout << "nEle = " << nEle << "\tnMuon = " << nMuon << std::endl;
      std::cout << "indices: " << theElectron << " " << thePositron << " " << theMuonMinus << " " << theMuonPlus << std::endl;
      std::cout << "chargeEle = " << chargeEle[theElectron] << "\tchargePos = " << chargeEle[thePositron] 
		<< "\tchargeMuonMinus = " << chargeMuon[theMuonMinus] << "\tchargeMuonPlus = " << chargeMuon[theMuonPlus] << std::endl;
      std::cout << "ee = " << m_channel[ee] << "\tmm = " << m_channel[mm] << "\temu = " << m_channel[em] << std::endl; 
    }

    // kinematics at preselection level 
    // only few variables since here we still do not know if two good leptons are reconstructed
    setPreselKinematics();
    
    // setting preselections
    CommonHiggsPreselection.SetWeight(weight);
    CommonHiggsPreselection.SetMcTruth(foundMcTree);
    CommonHiggsPreselection.SetHLT(passedHLT);
    CommonHiggsPreselection.SetNele(nPreselEle);   
    CommonHiggsPreselection.SetNmuon(nMuon);
    CommonHiggsPreselection.SetIsEE(m_channel[ee]);
    CommonHiggsPreselection.SetIsEM(m_channel[em]);
    CommonHiggsPreselection.SetIsMM(m_channel[mm]);
    CommonHiggsPreselection.SetHighElePt(hardestElectronPt);
    CommonHiggsPreselection.SetLowElePt(slowestElectronPt);
    CommonHiggsPreselection.SetHighMuonPt(hardestMuonPt);
    CommonHiggsPreselection.SetLowMuonPt(slowestMuonPt);
    CommonHiggsPreselection.SetMet(GetPt(pxTCMet[0],pyTCMet[0]));
    CommonHiggsPreselection.SetMllEE(m_mll[ee]);
    CommonHiggsPreselection.SetMllEM(m_mll[em]);
    CommonHiggsPreselection.SetMllMM(m_mll[mm]);

    // did we pass preselections?
    bool isPreselections = CommonHiggsPreselection.output();    
    if( !isPreselections ) continue;

    // kinematics after preselections
    // now we know we have two good reconstructed leptons
    setKinematics();

    // look for PV in the event (there is always at least 1 PV)
    m_closestPV = getPV();

    // ----------------------- selection ----------------------------
    // electron ID (true by default - studied only if ee or emu channel)
    bool theElectronID, thePositronID, theElectronIsol, thePositronIsol, theElectronConvRej, thePositronConvRej;
    theElectronID = thePositronID = theElectronIsol = thePositronIsol = theElectronConvRej = thePositronConvRej = true;
    // custom electron ID cuts (tight symmetric for ee, tight for emu)

    if (!_preselection->getSwitch("asymmetricID")) {
      if (theElectron > -1) isEleID(theElectron,&theElectronID,&theElectronIsol,&theElectronConvRej,EgammaCutBasedID);
      if (thePositron > -1) isEleID(thePositron,&thePositronID,&thePositronIsol,&thePositronConvRej,EgammaCutBasedID);
    }

    if (_preselection->getSwitch("asymmetricID")) {
      if (theElectron > -1) { 
	float pt = GetPt(pxEle[theElectron],pyEle[theElectron]);
	if (pt>=20) isEleID(theElectron,&theElectronID,&theElectronIsol,&theElectronConvRej,EgammaCutBasedID);
	if (pt<20)  isEleID(theElectron,&theElectronID,&theElectronIsol,&theElectronConvRej,EgammaCutBasedIDLow);
      }
      if (thePositron > -1) { 
	float pt = GetPt(pxEle[thePositron],pyEle[thePositron]);
	if (pt>=20) isEleID(thePositron,&thePositronID,&thePositronIsol,&thePositronConvRej,EgammaCutBasedID);
	if (pt<20) isEleID(thePositron,&thePositronID,&thePositronIsol,&thePositronConvRej,EgammaCutBasedIDLow);
      }
    }

    // if (theElectron > -1) theElectronID = anaUtils.electronIdVal(eleIdCutsEle[theElectron],eleIdTight);
    // if (thePositron > -1) thePositronID = anaUtils.electronIdVal(eleIdCutsEle[thePositron],eleIdTight);

    int scElectron = superClusterIndexEle[theElectron];
    int scPositron = superClusterIndexEle[thePositron];
    // filling the three to compare the distribution before the cut
    if (theElectron > -1) myEleIdTree -> fillAll(classificationEle[theElectron], hOverEEle[theElectron], eSuperClusterOverPEle[theElectron], eSeedOverPoutEle[theElectron], deltaEtaAtVtxEle[theElectron], deltaPhiAtVtxEle[theElectron], sqrt(covIEtaIEtaSC[scElectron]));
    if (thePositron > -1) myEleIdTree -> fillAll(classificationEle[thePositron], hOverEEle[thePositron], eSuperClusterOverPEle[thePositron], eSeedOverPoutEle[thePositron], deltaEtaAtVtxEle[thePositron], deltaPhiAtVtxEle[thePositron], sqrt(covIEtaIEtaSC[scPositron]));
    myEleIdTree->store();

    // --- ele isolation ---
    // separate isolations
    float theEleTrackerPtSum = ( theElectron > -1 ) ? dr03TkSumPtEle[theElectron] : 0.0;
    float thePosTrackerPtSum = ( thePositron > -1 ) ? dr03TkSumPtEle[thePositron] : 0.0;
    float theEleEcalPtSum = ( theElectron > -1 ) ? dr03EcalRecHitSumEtEle[theElectron] : 0.0;
    float thePosEcalPtSum = ( thePositron > -1 ) ? dr03EcalRecHitSumEtEle[thePositron] : 0.0;
    float theEleHcalPtSum = ( theElectron > -1 ) ? dr03HcalTowerSumEtEle[theElectron] : 0.0;
    float thePosHcalPtSum = ( thePositron > -1 ) ? dr03HcalTowerSumEtEle[thePositron] : 0.0;

    // --- muon ID / isolation ---
    bool theMuonPlusID = true;
    bool theMuonMinusID = true;
    if ( theMuonMinus > -1 ) isMuonID(theMuonMinus, &theMuonMinusID);
    if ( theMuonPlus > -1 ) isMuonID(theMuonPlus, &theMuonPlusID);

    // jet counter
    int njets = numJets();
    int nuncorrjets = numUncorrJets();

    // soft muon counter
    int nsoftmu = numSoftMuons();
    
    // extra lepton counter
    int nextraleptons = numExtraLeptons();

    // and calculate the track impact parameters and the event b-veto variables
    calcEventBVetoVariables(m_goodJets);

    // kine variables
    float theDeltaPhiEE, theDeltaErreEE, theInvMassEE, theSTransvMassEE, theDetaLeptonsEE = 0.;
    float theDeltaPhiMM, theDeltaErreMM, theInvMassMM, theSTransvMassMM, theDetaLeptonsMM = 0.;
    float theDeltaPhiEM, theDeltaErreEM, theInvMassEM, theSTransvMassEM, theDetaLeptonsEM = 0.;


    // ---------------------------------------
    // ee final state
    if (m_channel[ee]) {

      float theEleHardTrackerPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? theEleTrackerPtSum : thePosTrackerPtSum;
      float theEleSlowTrackerPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? thePosTrackerPtSum : theEleTrackerPtSum;
      float theEleHardHcalPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? theEleHcalPtSum : thePosHcalPtSum;
      float theEleSlowHcalPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? thePosHcalPtSum : theEleHcalPtSum;
      float theEleHardEcalPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? theEleEcalPtSum : thePosEcalPtSum;
      float theEleSlowEcalPtSum = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? thePosEcalPtSum : theEleEcalPtSum;

      int gsfEle = gsfTrackIndexEle[theElectron]; 
      int gsfPos = gsfTrackIndexEle[thePositron];
      float dxyEle = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                trackVxGsfTrack[gsfEle], trackVyGsfTrack[gsfEle], trackVzTrack[gsfEle], 
                                pxGsfTrack[gsfEle], pyGsfTrack[gsfEle], pzGsfTrack[gsfEle]);
      float dxyPos = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                trackVxGsfTrack[gsfPos], trackVyGsfTrack[gsfPos], trackVzTrack[gsfPos], 
                                pxGsfTrack[gsfPos], pyGsfTrack[gsfPos], pzGsfTrack[gsfPos]);

      float theEleHardDxy = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? dxyEle : dxyPos;
      float theEleSlowDxy = (m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) ? dxyPos : dxyEle;
      theDeltaPhiEE    = m_deltaPhi[ee];
      theDeltaErreEE   = m_deltaErre[ee];
      theInvMassEE     = m_mll[ee];
      theDetaLeptonsEE = etaEle[theElectron]-etaEle[thePositron];
      theSTransvMassEE  = m_transvMass[ee];

      int hardEle, slowEle;
      if(m_p4ElectronMinus->Pt() > m_p4ElectronPlus->Pt()) {
        hardEle = theElectron;
        slowEle = thePositron;
      } else {
        hardEle = thePositron;
        slowEle = theElectron;
      }

      setEleIdVariables(hardEle, slowEle);
      
      // selections
      CutBasedHiggsSelectionEE.SetWeight(weight);
      CutBasedHiggsSelectionEE.SetHighElePt(hardestElectronPt);
      CutBasedHiggsSelectionEE.SetLowElePt(slowestElectronPt);
      CutBasedHiggsSelectionEE.SetElectronId(theElectronID);
      CutBasedHiggsSelectionEE.SetPositronId(thePositronID);
      CutBasedHiggsSelectionEE.SetElectronIsolation(theElectronIsol);
      CutBasedHiggsSelectionEE.SetPositronIsolation(thePositronIsol);
      CutBasedHiggsSelectionEE.SetElectronConvRejection(theElectronConvRej);
      CutBasedHiggsSelectionEE.SetPositronConvRejection(thePositronConvRej);

      // pass through the cuts and use the bits above from EleID cuts selector
      CutBasedHiggsSelectionEE.SetEleHardTrackerPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleSlowTrackerPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleHardHcalPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleSlowHcalPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleHardEcalPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleSlowEcalPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleHardGlobalPtSum(-999);
      CutBasedHiggsSelectionEE.SetEleSlowGlobalPtSum(-999);
      
      // CutBasedHiggsSelectionEE.SetEleHardTrackerPtSum(theEleHardTrackerPtSum);
      // CutBasedHiggsSelectionEE.SetEleSlowTrackerPtSum(theEleSlowTrackerPtSum);
      // CutBasedHiggsSelectionEE.SetEleHardHcalPtSum(theEleHardHcalPtSum);
      // CutBasedHiggsSelectionEE.SetEleSlowHcalPtSum(theEleSlowHcalPtSum);
      // CutBasedHiggsSelectionEE.SetEleHardEcalPtSum(theEleHardEcalPtSum);
      // CutBasedHiggsSelectionEE.SetEleSlowEcalPtSum(theEleSlowEcalPtSum);
      // CutBasedHiggsSelectionEE.SetEleHardGlobalPtSum(theEleHardGlobalSum);
      // CutBasedHiggsSelectionEE.SetEleSlowGlobalPtSum(theEleSlowGlobalSum);
      CutBasedHiggsSelectionEE.SetEleHardD0(theEleHardDxy);
      CutBasedHiggsSelectionEE.SetEleSlowD0(theEleSlowDxy);
      CutBasedHiggsSelectionEE.SetNJets(njets);
      CutBasedHiggsSelectionEE.SetNUncorrJets(nuncorrjets);
      CutBasedHiggsSelectionEE.SetNSoftMuons(nsoftmu);
      CutBasedHiggsSelectionEE.SetNExtraLeptons(nextraleptons);
      CutBasedHiggsSelectionEE.SetMet(GetPt(pxTCMet[0],pyTCMet[0]));					
      CutBasedHiggsSelectionEE.SetProjectedMet(m_projectedMet[ee]);
      CutBasedHiggsSelectionEE.SetDeltaPhi(theDeltaPhiEE);
      CutBasedHiggsSelectionEE.SetInvMass(theInvMassEE);
      CutBasedHiggsSelectionEE.SetDetaLeptons(theDetaLeptonsEE);
      bool isSelectedEE = CutBasedHiggsSelectionEE.output();    
      bool selUpToFinalLeptonsEE = CutBasedHiggsSelectionEE.outputUpToFinalLeptons();
      bool selUpToJetVetoEE = CutBasedHiggsSelectionEE.outputUpToJetVeto();
      bool selUpToUncorrJetVetoEE = CutBasedHiggsSelectionEE.outputUpToUncorrJetVeto();
      bool selPreDeltaPhiEE = CutBasedHiggsSelectionEE.outputPreDeltaPhi();

      if(!_preselection->getSwitch("isData")) myOutTreeEE -> fillMcTruth(promptEE);
      else myOutTreeEE->fillRunInfos(runNumber, lumiBlock, eventNumber);

//       myOutTreeEE -> fillHLTElectrons( firedTrg[m_requiredTriggers[0]], 
// 				       firedTrg[m_requiredTriggers[1]],
// 				       (firedTrg[m_requiredTriggers[0]] || firedTrg[m_requiredTriggers[1]]) );

      myOutTreeEE -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), 
			     GetPt(pxPFMet[0],pyPFMet[0]), 
			     GetPt(pxMet[0],pyMet[0]), 
                             m_projectedMet[ee],
			     theDeltaPhiEE, 
			     theDeltaErreEE, 
			     theSTransvMassEE, 
			     theInvMassEE, 
			     hardestElectronPt, 
			     slowestElectronPt, 
			     theDetaLeptonsEE,
			     selUpToFinalLeptonsEE,
			     selUpToJetVetoEE,
			     selUpToUncorrJetVetoEE,
			     selPreDeltaPhiEE,
			     isSelectedEE);

      myOutTreeEE -> fillMLVars(njets,
                                nuncorrjets,
                                m_maxDxyEvt,
                                m_maxDszEvt,
                                m_maxTrackCountingHighEffBJetTags,
                                m_maxImpactParameterMVABJetTags,
                                m_maxCombinedSecondaryVertexMVABJetTags);

      myOutTreeEE -> fillElectrons( myRecoflag, myPt, myEta, myPhi,
                                    myClassification, myNBremClusters, myDeta, myDphi, myHoe, mySee, mySpp, myEop, myFbrem,
                                    myTrackerIso, myHcalIso, myEcalJIso, myEcalGTIso, myCombinedIso, myCharge, myMissHits, myDist, myDcot, myLh, myMatched );

      if ( _preselection->getSwitch("apply_kFactor") ) {
	myOutTreeEE->fillKFactor(evtKfactor);
      }

      // dumping final tree
      myOutTreeEE -> store();
    }

    // ---------------------------------------
    // mm final state
    if (m_channel[mm]){
      theDeltaPhiMM    = m_deltaPhi[mm];
      theDeltaErreMM   = m_deltaErre[mm];
      theInvMassMM     = m_mll[mm];
      theDetaLeptonsMM = etaEle[theMuonMinus]-etaEle[theMuonPlus];
      theSTransvMassMM  = m_transvMass[mm];

      float theMuonMinusGlobalSum = muonIsoGlobalSum(theMuonMinus, theMuonPlus) / m_p4MuonMinus->Pt();
      float theMuonPlusGlobalSum  = muonIsoGlobalSum(theMuonPlus, theMuonMinus) / m_p4MuonPlus->Pt();
      float theMuonHardGlobalSum  = (m_p4MuonMinus->Pt() > m_p4MuonPlus->Pt()) ? theMuonMinusGlobalSum : theMuonPlusGlobalSum;
      float theMuonSlowGlobalSum  = (m_p4MuonMinus->Pt() > m_p4MuonPlus->Pt()) ? theMuonPlusGlobalSum : theMuonMinusGlobalSum;

      int ctfMuonMinus = trackIndexMuon[theMuonMinus]; 
      int ctfMuonPlus = trackIndexMuon[theMuonPlus]; 
      float dxyMuonMinus = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                      trackVxTrack[ctfMuonMinus], trackVyTrack[ctfMuonMinus], trackVzTrack[ctfMuonMinus], 
                                      pxTrack[ctfMuonMinus], pyTrack[ctfMuonMinus], pzTrack[ctfMuonMinus]);
      float dxyMuonPlus = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                     trackVxTrack[ctfMuonPlus], trackVyTrack[ctfMuonPlus], trackVzTrack[ctfMuonPlus], 
                                     pxTrack[ctfMuonPlus], pyTrack[ctfMuonPlus], pzTrack[ctfMuonPlus]);

      float theMuonHardDxy = (m_p4MuonMinus->Pt() > m_p4MuonPlus->Pt()) ? dxyMuonMinus : dxyMuonPlus;
      float theMuonSlowDxy = (m_p4MuonMinus->Pt() > m_p4MuonPlus->Pt()) ? dxyMuonPlus : dxyMuonMinus;

      // selections
      CutBasedHiggsSelectionMM.SetWeight(weight);
      CutBasedHiggsSelectionMM.SetHighElePt(hardestMuonPt);
      CutBasedHiggsSelectionMM.SetLowElePt(slowestMuonPt);
      CutBasedHiggsSelectionMM.SetElectronId(theMuonMinusID);
      CutBasedHiggsSelectionMM.SetPositronId(theMuonPlusID);
      // bits not used for muons: pass through and use cuts (global iso)
      CutBasedHiggsSelectionMM.SetElectronIsolation(true);
      CutBasedHiggsSelectionMM.SetPositronIsolation(true);
      CutBasedHiggsSelectionMM.SetElectronConvRejection(true);
      CutBasedHiggsSelectionMM.SetPositronConvRejection(true);
      CutBasedHiggsSelectionMM.SetEleHardTrackerPtSum(0);
      CutBasedHiggsSelectionMM.SetEleSlowTrackerPtSum(0);
      CutBasedHiggsSelectionMM.SetEleHardHcalPtSum(0);
      CutBasedHiggsSelectionMM.SetEleSlowHcalPtSum(0);
      CutBasedHiggsSelectionMM.SetEleHardEcalPtSum(0);
      CutBasedHiggsSelectionMM.SetEleSlowEcalPtSum(0);
      CutBasedHiggsSelectionMM.SetEleHardGlobalPtSum(theMuonHardGlobalSum);
      CutBasedHiggsSelectionMM.SetEleSlowGlobalPtSum(theMuonSlowGlobalSum);
      CutBasedHiggsSelectionMM.SetEleHardD0(theMuonHardDxy);
      CutBasedHiggsSelectionMM.SetEleSlowD0(theMuonSlowDxy);
      CutBasedHiggsSelectionMM.SetNJets(njets);
      CutBasedHiggsSelectionMM.SetNUncorrJets(nuncorrjets);
      CutBasedHiggsSelectionMM.SetNSoftMuons(nsoftmu);
      CutBasedHiggsSelectionMM.SetNExtraLeptons(nextraleptons);
      CutBasedHiggsSelectionMM.SetMet(GetPt(pxTCMet[0],pyTCMet[0]));					
      CutBasedHiggsSelectionMM.SetProjectedMet(m_projectedMet[mm]);
      CutBasedHiggsSelectionMM.SetDeltaPhi(theDeltaPhiMM);
      CutBasedHiggsSelectionMM.SetInvMass(theInvMassMM);
      CutBasedHiggsSelectionMM.SetDetaLeptons(theDetaLeptonsMM);
      bool isSelectedMM = CutBasedHiggsSelectionMM.output();    
      bool selUpToFinalLeptonsMM = CutBasedHiggsSelectionMM.outputUpToFinalLeptons();
      bool selUpToJetVetoMM = CutBasedHiggsSelectionMM.outputUpToJetVeto();
      bool selUpToUncorrJetVetoMM = CutBasedHiggsSelectionMM.outputUpToUncorrJetVeto();
      bool selPreDeltaPhiMM = CutBasedHiggsSelectionMM.outputPreDeltaPhi();

      if(!_preselection->getSwitch("isData")) myOutTreeMM -> fillMcTruth(promptMM);
      else myOutTreeMM->fillRunInfos(runNumber, lumiBlock, eventNumber);

//       myOutTreeMM -> fillHLTMuons( firedTrg[m_requiredTriggers[2]], 
// 				   firedTrg[m_requiredTriggers[3]],
// 				   (firedTrg[m_requiredTriggers[2]] || firedTrg[m_requiredTriggers[3]]) );
      
      myOutTreeMM -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), 
			     GetPt(pxPFMet[0],pyPFMet[0]), 
			     GetPt(pxMet[0],pyMet[0]), 
			     m_projectedMet[mm],
                             theDeltaPhiMM, 
			     theDeltaErreMM, 
			     theSTransvMassMM, 
			     theInvMassMM, 
			     hardestMuonPt, 
			     slowestMuonPt, 
			     theDetaLeptonsMM,
			     selUpToFinalLeptonsMM,
			     selUpToJetVetoMM,
			     selUpToUncorrJetVetoMM,
			     selPreDeltaPhiMM,
			     isSelectedMM);

      myOutTreeMM -> fillMLVars(njets,
                                nuncorrjets,
                                m_maxDxyEvt,
                                m_maxDszEvt,
                                m_maxTrackCountingHighEffBJetTags,
                                m_maxImpactParameterMVABJetTags,
                                m_maxCombinedSecondaryVertexMVABJetTags);
      
      if ( _preselection->getSwitch("apply_kFactor") ) {
	myOutTreeMM->fillKFactor(evtKfactor);
      }

      // dumping final tree
      myOutTreeMM -> store();
      
    }
    
    if (m_channel[em]){
      
      // electron ID, isolations for the only electron / positron
      bool theEleIDEM = true;
      bool theMuonIDEM = true;
      bool theEleIsolEM = true;
      bool theEleConvRejEM = true;
      float theEleTrackerPtSumEM = 0.0;
      float theEleHcalPtSumEM = 0.0;
      float theEleEcalPtSumEM = 0.0;
      float theMuonGlobalSumEM = 0.0;
      float theEleDxy = 0.0;
      float theMuonDxy = 0.0;

      theDeltaPhiEM    = m_deltaPhi[em];
      theDeltaErreEM   = m_deltaErre[em];
      theInvMassEM     = m_mll[em];
      theSTransvMassEM  = m_transvMass[em];
      if(theElectron>-1 && theMuonPlus>-1) {
	theDetaLeptonsEM = etaEle[theElectron]-etaEle[theMuonPlus];
        theEleIDEM = theElectronID;
        theMuonIDEM = theMuonPlusID;
        theEleIsolEM = theElectronIsol;
        theEleConvRejEM = theElectronConvRej;
	theEleTrackerPtSumEM = theEleTrackerPtSum;
	theEleHcalPtSumEM = theEleHcalPtSum;
	theEleEcalPtSumEM = theEleEcalPtSum;
        theMuonGlobalSumEM = mueleIsoGlobalSum(theMuonPlus, theElectron) / m_p4MuonPlus->Pt();

        int gsfEle = gsfTrackIndexEle[theElectron]; 
        int ctfMuonPlus = trackIndexMuon[theMuonPlus]; 
        float dxyMuonPlus = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                       trackVxTrack[ctfMuonPlus], trackVyTrack[ctfMuonPlus], trackVzTrack[ctfMuonPlus], 
                                       pxTrack[ctfMuonPlus], pyTrack[ctfMuonPlus], pzTrack[ctfMuonPlus]);
        float dxyEle = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                  trackVxGsfTrack[gsfEle], trackVyGsfTrack[gsfEle], trackVzTrack[gsfEle], 
                                  pxGsfTrack[gsfEle], pyGsfTrack[gsfEle], pzGsfTrack[gsfEle]);

        theEleDxy = dxyEle;
        theMuonDxy = dxyMuonPlus;

        setEleIdVariables(theElectron, -1);
      }
      if(thePositron>-1 && theMuonMinus>-1 ) {
	theDetaLeptonsEM = etaEle[thePositron]-etaEle[theMuonMinus];
        theEleIDEM = thePositronID;
        theMuonIDEM = theMuonMinusID;
        theEleIsolEM = thePositronIsol;
        theEleConvRejEM = thePositronConvRej;
	theEleTrackerPtSumEM = thePosTrackerPtSum;
	theEleHcalPtSumEM = thePosHcalPtSum;
	theEleEcalPtSumEM = thePosEcalPtSum;
        theMuonGlobalSumEM = mueleIsoGlobalSum(theMuonMinus, thePositron) / m_p4MuonMinus->Pt();

        int gsfPos = gsfTrackIndexEle[thePositron]; 
        int ctfMuonMinus = trackIndexMuon[theMuonMinus]; 
        float dxyMuonMinus = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                        trackVxTrack[ctfMuonMinus], trackVyTrack[ctfMuonMinus], trackVzTrack[ctfMuonMinus], 
                                        pxTrack[ctfMuonMinus], pyTrack[ctfMuonMinus], pzTrack[ctfMuonMinus]);
        float dxyPos = trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                  trackVxGsfTrack[gsfPos], trackVyGsfTrack[gsfPos], trackVzTrack[gsfPos], 
                                  pxGsfTrack[gsfPos], pyGsfTrack[gsfPos], pzGsfTrack[gsfPos]);
        theEleDxy = dxyPos;
        theMuonDxy = dxyMuonMinus;
        setEleIdVariables(thePositron, -1);
      }


      // selections
      CutBasedHiggsSelectionEM.SetWeight(weight);
      CutBasedHiggsSelectionEM.SetHighElePt(hardestElectronPt);
      CutBasedHiggsSelectionEM.SetLowElePt(slowestMuonPt);
      CutBasedHiggsSelectionEM.SetElectronId(true);
      CutBasedHiggsSelectionEM.SetElectronId(theMuonIDEM);
      CutBasedHiggsSelectionEM.SetPositronId(theEleIDEM);
      CutBasedHiggsSelectionEM.SetElectronIsolation(true); // bit: pass through and use cut
      CutBasedHiggsSelectionEM.SetPositronIsolation(theEleIsolEM);
      CutBasedHiggsSelectionEM.SetPositronConvRejection(true); // bit: pass through
      CutBasedHiggsSelectionEM.SetElectronConvRejection(theEleConvRejEM);
      CutBasedHiggsSelectionEM.SetEleHardTrackerPtSum(0);
      CutBasedHiggsSelectionEM.SetEleSlowTrackerPtSum(0);
      CutBasedHiggsSelectionEM.SetEleHardHcalPtSum(0);
      CutBasedHiggsSelectionEM.SetEleSlowHcalPtSum(0);
      CutBasedHiggsSelectionEM.SetEleHardEcalPtSum(0);
      CutBasedHiggsSelectionEM.SetEleSlowEcalPtSum(0);
      CutBasedHiggsSelectionEM.SetEleHardGlobalPtSum(theMuonGlobalSumEM); //order in pt unimportant
      CutBasedHiggsSelectionEM.SetEleSlowGlobalPtSum(0); // for electron, already applied the bit
      CutBasedHiggsSelectionEM.SetEleHardD0(theEleDxy);
      CutBasedHiggsSelectionEM.SetEleSlowD0(theMuonDxy);
      CutBasedHiggsSelectionEM.SetNJets(njets);
      CutBasedHiggsSelectionEM.SetNSoftMuons(nsoftmu);
      CutBasedHiggsSelectionEM.SetNExtraLeptons(nextraleptons);
      CutBasedHiggsSelectionEM.SetNUncorrJets(nuncorrjets);
      CutBasedHiggsSelectionEM.SetMet(GetPt(pxTCMet[0],pyTCMet[0]));					
      CutBasedHiggsSelectionEM.SetProjectedMet(m_projectedMet[em]);
      CutBasedHiggsSelectionEM.SetDeltaPhi(theDeltaPhiEM);
      CutBasedHiggsSelectionEM.SetInvMass(theInvMassEM);
      CutBasedHiggsSelectionEM.SetDetaLeptons(theDetaLeptonsEM);
      bool isSelectedEM = CutBasedHiggsSelectionEM.output();    
      bool selUpToFinalLeptonsEM = CutBasedHiggsSelectionEM.outputUpToFinalLeptons();
      bool selUpToJetVetoEM = CutBasedHiggsSelectionEM.outputUpToJetVeto();
      bool selUpToUncorrJetVetoEM = CutBasedHiggsSelectionEM.outputUpToUncorrJetVeto();
      bool selPreDeltaPhiEM = CutBasedHiggsSelectionEM.outputPreDeltaPhi();

      if(!_preselection->getSwitch("isData")) myOutTreeEM -> fillMcTruth(promptEM);
      else myOutTreeEM->fillRunInfos(runNumber, lumiBlock, eventNumber);

//       myOutTreeEM -> fillHLTElectrons( firedTrg[m_requiredTriggers[0]], 
// 				       firedTrg[m_requiredTriggers[1]],
// 				       (firedTrg[m_requiredTriggers[0]] || firedTrg[m_requiredTriggers[1]]) );

//       myOutTreeEM -> fillHLTMuons( firedTrg[m_requiredTriggers[2]], 
// 				   firedTrg[m_requiredTriggers[3]],
// 				   (firedTrg[m_requiredTriggers[2]] || firedTrg[m_requiredTriggers[3]]) );

      myOutTreeEM -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), 
			     GetPt(pxPFMet[0],pyPFMet[0]), 
			     GetPt(pxMet[0],pyMet[0]), 
                             m_projectedMet[em],
			     theDeltaPhiEM, 
			     theDeltaErreEM, 
			     theSTransvMassEM, 
			     theInvMassEM, 
			     hardestMuonPt, 
			     slowestMuonPt, 
			     theDetaLeptonsEM,
			     selUpToFinalLeptonsEM,
			     selUpToJetVetoEM,
			     selUpToUncorrJetVetoEM,
			     selPreDeltaPhiEM,
			     isSelectedEM);

      myOutTreeEM -> fillMLVars(njets,
                                nuncorrjets,
                                m_maxDszEvt,
                                m_maxDszEvt,
                                m_maxTrackCountingHighEffBJetTags,
                                m_maxImpactParameterMVABJetTags,
                                m_maxCombinedSecondaryVertexMVABJetTags);
      
      myOutTreeEM -> fillElectrons( myRecoflag, myPt, myEta, myPhi,
                                    myClassification, myNBremClusters, myDeta, myDphi, myHoe, mySee, mySpp, myEop, myFbrem,
                                    myTrackerIso, myHcalIso, myEcalJIso, myEcalGTIso, myCombinedIso, myCharge, myMissHits, myDist, myDcot, myLh, myMatched );

      if ( _preselection->getSwitch("apply_kFactor") ) {
	myOutTreeEM->fillKFactor(evtKfactor);
      }

      // dumping final tree
      myOutTreeEM -> store();

    }

  }

  fMatch = new TFile("matching.root","RECREATE");
  fMatch->cd();
  H_deltaRuncorr->Write();
  H_deltaRcorr->Write();
  fMatch->Close();
}

void HiggsMLSelection::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Common preselections: " << std::endl;
  CommonHiggsPreselection.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EE selections: " << std::endl;
  CutBasedHiggsSelectionEE.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full MM selections: " << std::endl;
  CutBasedHiggsSelectionMM.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EM selections: " << std::endl;
  CutBasedHiggsSelectionEM.displayEfficiencies(datasetName);

  if (!_preselection->getSwitch("asymmetricID")) {
    std::cout << "symmetric ID: " << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  } else {
    std::cout << "asymmetric ID: Low pT" << std::endl;
    EgammaCutBasedIDLow.displayEfficiencies();
    std::cout << "asymmetric ID: High pT" << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  }
}

std::pair<int,int> HiggsMLSelection::getBestElectronPair() {

  int theLep1=-1;
  int theLep2=-1;
  float maxPtLep1=-1000.;
  float maxPtLep2=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nEle;i++) {
    
    if(_preselection->getSwitch("etaElectronAcc") && !_preselection->passCut("etaElectronAcc",etaEle[i]) ) continue;
    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();
    if(_preselection->getSwitch("ptElectronAcc") && !_preselection->passCut("ptElectronAcc",thisPt) ) continue;
    // fixme: (in the future) put here full electron ID not to loose efficiency I think
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestElectrons->clear();
  _bestElectrons->push_back(theLep1);  _bestElectrons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

std::pair<int,int> HiggsMLSelection::getBestMuonPair() {
  int theLep1=-1;
  int theLep2=-1;
  float maxPtLep1=-1000.;
  float maxPtLep2=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nMuon;i++) {
    if(_preselection->getSwitch("etaMuonAcc") && !_preselection->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    TVector3 pLepton(pxMuon[i],pyMuon[i],pzMuon[i]);
    float thisPt=pLepton.Pt();
    if(_preselection->getSwitch("ptMuonAcc") && !_preselection->passCut("ptMuonAcc",thisPt) ) continue;
    // fixme: (in the future) put here full electron ID not to loose efficiency I think
    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestMuons->clear();
  _bestMuons->push_back(theLep1);  _bestMuons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

void HiggsMLSelection::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector thisCutBasedID) {
  
  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);
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


  thisCutBasedID.SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  thisCutBasedID.SetRecoFlag(recoFlagsEle[eleIndex]);
  thisCutBasedID.applyElectronIDOnPFlowElectrons(true);
  thisCutBasedID.SetHOverE( HoE );
  thisCutBasedID.SetS9S25( s9s25 );
  thisCutBasedID.SetDEta( deta );
  thisCutBasedID.SetDPhiIn( dphiin );
  thisCutBasedID.SetDPhiOut( dphiout );
  thisCutBasedID.SetBremFraction( fbrem );
  thisCutBasedID.SetSigmaEtaEta( see );
  thisCutBasedID.SetSigmaPhiPhi( spp );
  thisCutBasedID.SetEOverPout( eopout );
  thisCutBasedID.SetEOverPin( eop );
  thisCutBasedID.SetElectronClass ( classificationEle[eleIndex] );
  thisCutBasedID.SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  thisCutBasedID.SetLikelihood( likelihoodRatio(eleIndex,*LH) );
  thisCutBasedID.SetEcalIsolation( dr03EcalRecHitSumEtEle[eleIndex] );
  thisCutBasedID.SetTrkIsolation( dr03TkSumPtEle[eleIndex] );
  thisCutBasedID.SetHcalIsolation( dr03HcalTowerSumEtEle[eleIndex] );
  thisCutBasedID.SetCombinedIsolation( (dr03TkSumPtEle[eleIndex] + 
                                          TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + 
                                          dr03HcalTowerSumEtEle[eleIndex]) / pt );
  thisCutBasedID.SetMissingHits( expInnerLayersGsfTrack[gsf] );
  thisCutBasedID.SetConvDist( fabs(convDistEle[eleIndex]) );
  thisCutBasedID.SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
  thisCutBasedID.m_cleaner->SetE1(e1);
  thisCutBasedID.m_cleaner->SetE4SwissCross(e4SwissCross);
  thisCutBasedID.m_cleaner->SetFiducialFlag(fidFlagSC);
  thisCutBasedID.m_cleaner->SetSeedFlag(seedRecHitFlag);
  thisCutBasedID.m_cleaner->SetSeedTime(seedTime);
  thisCutBasedID.m_cleaner->SetSeedChi2(seedChi2);

  //  return egammaCutBasedID.output(); // class dependent result
  *eleIdOutput = thisCutBasedID.outputNoClassEleId();
  *isolOutput = thisCutBasedID.outputNoClassIso();
  *convRejOutput = thisCutBasedID.outputNoClassConv();
}

void HiggsMLSelection::isMuonID(int muonIndex, bool *muonIdOutput) {

  *muonIdOutput = true;
  
  Utils anaUtils; 
  bool flag = anaUtils.muonIdVal(muonIdMuon[muonIndex],AllGlobalMuons) && 
    anaUtils.muonIdVal(muonIdMuon[muonIndex],AllTrackerMuons);
  // the following cuts are based on KF and global muon track. So if the cut above has failed, return here
  if(!flag) {
    *muonIdOutput = false;
    return;
  }
  int track = trackIndexMuon[muonIndex];
  if(trackValidHitsTrack[track]<=10) *muonIdOutput = false;
  int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
  if(trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] >= 10) *muonIdOutput = false;
  if(trackValidHitsGlobalMuonTrack[globalMuonTrack] == 0) *muonIdOutput = false;
  float dxy = fabs(trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                              trackVxTrack[track], trackVyTrack[track], trackVzTrack[track], 
                              pxTrack[track], pyTrack[track], pzTrack[track]));
  if(dxy > 0.020) *muonIdOutput = false;

}

void HiggsMLSelection::setPreselKinematics() {

  // highest and lowest pt for electrons
  // + mll if 2 good electrons are reconstructed
  if (thePositron > -1 && theElectron > -1) {
    hardestElectronPt = TMath::Max(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxEle[thePositron],pyEle[thePositron]));
    slowestElectronPt = TMath::Min(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxEle[thePositron],pyEle[thePositron]));
    m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
    m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
    m_mll[ee] = (*m_p4ElectronMinus + *m_p4ElectronPlus).M();
  }
  if (thePositron <= -1 && theElectron > -1) {
    hardestElectronPt = GetPt(pxEle[theElectron],pyEle[theElectron]);
    slowestElectronPt = GetPt(pxEle[theElectron],pyEle[theElectron]);
    m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
    m_mll[ee] = -999.;
  }
  if (thePositron > -1 && theElectron <= -1) {
    hardestElectronPt = GetPt(pxEle[thePositron],pyEle[thePositron]);
    slowestElectronPt = GetPt(pxEle[thePositron],pyEle[thePositron]);
    m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
    m_mll[ee] = -999.;
  }
  if (theElectron <= -1 && thePositron <= -1) {
    hardestElectronPt = -999.;
    slowestElectronPt = -999.;
    m_mll[ee] = -999.;
  }

  // highest and lowest pt for muons
  // + mll if 2 good muons are reconstructed  
  if (theMuonPlus > -1 && theMuonMinus > -1){ 
    hardestMuonPt = TMath::Max(GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
    slowestMuonPt = TMath::Min(GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
    m_p4MuonMinus -> SetXYZT(pxMuon[theMuonMinus],pyMuon[theMuonMinus],pzMuon[theMuonMinus],energyMuon[theMuonMinus]);
    m_p4MuonPlus  -> SetXYZT(pxMuon[theMuonPlus],pyMuon[theMuonPlus],pzMuon[theMuonPlus],energyMuon[theMuonPlus]);      
    m_mll[mm] = (*m_p4MuonMinus + *m_p4MuonPlus).M();
  }
  if (theMuonMinus > -1 && theMuonPlus <= -1) {
    hardestMuonPt = GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]);
    slowestMuonPt = GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]);
    m_p4MuonMinus -> SetXYZT(pxMuon[theMuonMinus],pyMuon[theMuonMinus],pzMuon[theMuonMinus],energyMuon[theMuonMinus]);
    m_mll[mm]     = -999.;
  }
  if (theMuonPlus > -1 && theMuonMinus <= -1) {
    hardestMuonPt = GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]);
    slowestMuonPt = GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]);
    m_p4MuonPlus  -> SetXYZT(pxMuon[theMuonPlus],pyMuon[theMuonPlus],pzMuon[theMuonPlus],energyMuon[theMuonPlus]);      
    m_mll[mm]     = -999.;
  }
  if (theMuonMinus <= -1 && theMuonPlus <= -1) {
    hardestMuonPt = -999.;
    slowestMuonPt = -999.;
    m_mll[mm]     = -999.;
  }
  
  // mixed channel mll: only a lower cut used to reject fake dilepton pairs originating from a single object
  // no need to consider cases with extra pairs: events with extra leptons will be rejected by selection
  if ( theElectron > -1 && theMuonPlus > -1 ) {
    m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
    m_p4MuonPlus  -> SetXYZT(pxMuon[theMuonPlus],pyMuon[theMuonPlus],pzMuon[theMuonPlus],energyMuon[theMuonPlus]);
    m_mll[em] = (*m_p4ElectronMinus + *m_p4MuonPlus).M();
  }
  if ( thePositron > -1 && theMuonMinus > -1 ) {
    m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);    
    m_p4MuonMinus -> SetXYZT(pxMuon[theMuonMinus],pyMuon[theMuonMinus],pzMuon[theMuonMinus],energyMuon[theMuonMinus]);
    m_mll[em] = (*m_p4ElectronPlus + *m_p4MuonMinus).M();
  }
  
  // MET
  m_p4MET->SetXYZT(pxTCMet[0],pyTCMet[0],pzTCMet[0],energyTCMet[0]); 
}



void HiggsMLSelection::setKinematics( ) {

  // electron variables used for ele quality in jet veto 
  m_HoEElectronMinus     = hOverEEle[theElectron];
  m_HoEElectronPlus      = hOverEEle[thePositron];
  int scEle = superClusterIndexEle[theElectron];
  int scPos = superClusterIndexEle[thePositron];
  m_CaloEneElectronMinus = energySC[scEle];
  m_CaloEneElectronPlus  = energySC[scPos];

  // compute delta Phi in degrees, di-lepton invariant mass, transverse mass
  TVector3 dilepPt;
  if ( m_channel[ee] ) {
    m_deltaPhi[ee]  = fabs(180./TMath::Pi() * m_p4ElectronMinus->Vect().DeltaPhi(m_p4ElectronPlus->Vect()));
    m_deltaErre[ee] = m_p4ElectronMinus->Vect().DeltaR(m_p4ElectronPlus->Vect());
    dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4ElectronPlus->Vect().X(),
		    m_p4ElectronMinus->Vect().Y()+m_p4ElectronPlus->Vect().Y(),
		    0.0 );
    // m_transvMass[ee]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
    // def. 3 of http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=104213
    m_transvMass[ee]=mT3(*m_p4ElectronMinus,*m_p4ElectronPlus,m_p4MET->Vect());
    // m_mT2[ee] = mT2(m_p4ElectronMinus->Vect(),m_p4ElectronPlus->Vect(),m_p4MET->Vect());
    m_mT2[ee] = 0.;
    m_projectedMet[ee] = GetProjectedMet(m_p4ElectronMinus->Vect(),m_p4ElectronPlus->Vect()); 
  }
  else {    
    m_deltaPhi[ee]   = -1.;
    m_deltaErre[ee]  = -1.;
    m_transvMass[ee] = -1.;
    m_mT2[ee] = -1.;
    m_projectedMet[ee] = -1.;
  }

  if ( m_channel[mm] ) {    
    m_deltaPhi[mm]  = fabs(180./TMath::Pi() * m_p4MuonMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()));
    m_deltaErre[mm] = m_p4MuonMinus->Vect().DeltaR(m_p4MuonPlus->Vect());
    dilepPt.SetXYZ( m_p4MuonMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		    m_p4MuonMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		    0.0 );
    // m_transvMass[mm]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect()))));
    m_transvMass[mm]=mT3(*m_p4MuonMinus,*m_p4MuonPlus,m_p4MET->Vect());
    m_projectedMet[mm] = GetProjectedMet(m_p4MuonMinus->Vect(),m_p4MuonPlus->Vect()); 
    // m_mT2[mm] = mT2(m_p4MuonMinus->Vect(),m_p4MuonPlus->Vect(),m_p4MET->Vect());
    m_mT2[mm] = 0.;
  }
  else { 
    m_deltaPhi[mm]   = -1.;
    m_deltaErre[mm]  = -1.;
    m_transvMass[mm] = -1.;
    m_mT2[mm] = -1.;
    m_projectedMet[mm] = -1.;
  }
  
  if ( m_channel[em] ) {
    float deltaPhiEPlusMuMinus  = -1.0;
    float deltaPhiEMinusMuPlus  = -1.0;
    float deltaErreEPlusMuMinus = -1.0;
    float deltaErreEMinusMuPlus = -1.0;
    float dilepPtEPlusMuMinus   = -1.0;
    float dilepPtEMinusMuPlus   = -1.0;

    if ( thePositron > -1 && theMuonMinus > -1 ) {
      deltaPhiEPlusMuMinus  = fabs(180./TMath::Pi() * m_p4ElectronPlus->Vect().DeltaPhi(m_p4MuonMinus->Vect()));
      deltaErreEPlusMuMinus = m_p4ElectronPlus->Vect().DeltaR(m_p4MuonMinus->Vect());
      m_deltaPhi[em]  = deltaPhiEPlusMuMinus;
      m_deltaErre[em] = deltaErreEPlusMuMinus;
      dilepPt.SetXYZ( m_p4ElectronPlus->Vect().X()+m_p4MuonMinus->Vect().X(),
		      m_p4ElectronPlus->Vect().Y()+m_p4MuonMinus->Vect().Y(),
		      0.0 );
      dilepPtEPlusMuMinus = dilepPt.Mag();
      //m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect()))));
      m_transvMass[em]=mT3(*m_p4ElectronPlus,*m_p4MuonMinus,m_p4MET->Vect());
      hardestLeptonPt = TMath::Max(GetPt(pxEle[thePositron],pyEle[thePositron]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
      slowestLeptonPt = TMath::Min(GetPt(pxEle[thePositron],pyEle[thePositron]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
      m_projectedMet[em] = GetProjectedMet(m_p4ElectronPlus->Vect(),m_p4MuonMinus->Vect()); 
      // m_mT2[em] = mT2(m_p4ElectronPlus->Vect(),m_p4MuonMinus->Vect(),m_p4MET->Vect());
      m_mT2[em] = 0.;
    }
    if ( theElectron > -1 && theMuonPlus > -1 ) {
      deltaPhiEMinusMuPlus  = fabs( 180./TMath::Pi() * m_p4ElectronMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()));
      deltaErreEMinusMuPlus = m_p4ElectronMinus->Vect().DeltaR(m_p4MuonPlus->Vect());
      m_deltaPhi[em]  = deltaPhiEMinusMuPlus;
      m_deltaErre[em] = deltaErreEMinusMuPlus;
      dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		      m_p4ElectronMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		      0.0 );
      dilepPtEMinusMuPlus = dilepPt.Mag();
      // m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
      m_transvMass[em]=mT3(*m_p4ElectronMinus,*m_p4MuonPlus,m_p4MET->Vect());
      hardestLeptonPt = TMath::Max(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]));
      slowestLeptonPt = TMath::Min(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]));
      m_projectedMet[em] = GetProjectedMet(m_p4ElectronMinus->Vect(),m_p4MuonPlus->Vect()); 
      // m_mT2[em] = mT2(m_p4ElectronMinus->Vect(),m_p4MuonPlus->Vect(),m_p4MET->Vect());
      m_mT2[em] = 0.;
    }
    if ( thePositron > -1 && theMuonMinus > -1 &&
	 theElectron > -1 && theMuonPlus > -1) {
      
      // if two pairs are built we choose the one with highest di-lepton pt
      if ( dilepPtEPlusMuMinus > dilepPtEMinusMuPlus ) {
	m_deltaPhi[em]  = deltaPhiEPlusMuMinus;
	m_deltaErre[em] = deltaErreEPlusMuMinus;
	dilepPt.SetXYZ( m_p4ElectronPlus->Vect().X()+m_p4MuonMinus->Vect().X(),
			m_p4ElectronPlus->Vect().Y()+m_p4MuonMinus->Vect().Y(),
			0.0 );
	// m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
        m_transvMass[em]=mT3(*m_p4ElectronPlus,*m_p4MuonMinus,m_p4MET->Vect()); 
	hardestLeptonPt = TMath::Max(GetPt(pxEle[thePositron],pyEle[thePositron]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
	slowestLeptonPt = TMath::Min(GetPt(pxEle[thePositron],pyEle[thePositron]),GetPt(pxMuon[theMuonMinus],pyMuon[theMuonMinus]));
        m_projectedMet[em] = GetProjectedMet(m_p4ElectronPlus->Vect(),m_p4MuonMinus->Vect()); 
        // m_mT2[em] = mT2(m_p4ElectronPlus->Vect(),m_p4MuonMinus->Vect(),m_p4MET->Vect());
        m_mT2[em] = 0.;
      }
      else {
	m_deltaPhi[em]  = deltaPhiEMinusMuPlus;
	m_deltaErre[em] = deltaErreEMinusMuPlus;
	dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
			m_p4ElectronMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
			0.0 );
	//m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
        m_transvMass[em]=mT3(*m_p4ElectronMinus,*m_p4MuonPlus,m_p4MET->Vect());
	hardestLeptonPt = TMath::Max(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]));
	slowestLeptonPt = TMath::Min(GetPt(pxEle[theElectron],pyEle[theElectron]),GetPt(pxMuon[theMuonPlus],pyMuon[theMuonPlus]));
        m_projectedMet[em] = GetProjectedMet(m_p4ElectronMinus->Vect(),m_p4MuonPlus->Vect()); 
        // m_mT2[em] = mT2(m_p4ElectronMinus->Vect(),m_p4MuonPlus->Vect(),m_p4MET->Vect());
        m_mT2[em] = 0.;
      }
    }
  }
  else {
    m_deltaPhi[em]  = -1.;
    m_deltaErre[em] = -1.;
    m_transvMass[em] = -1.;
    m_mT2[em] = -1.;
    m_projectedMet[em] = -1;
  }
  
  if( !_preselection->getSwitch("isData") ) {
    // --- Higgs and electron generator level ---
    if(_theGenEle>0 && _theGenPos>0) {
      if(pMc[_theGenEle]>=pMc[_theGenPos]) {
        _highestPtGen[0]=pMc[_theGenEle]*fabs(sin(thetaMc[_theGenEle]));
        _lowestPtGen[0]=pMc[_theGenPos]*fabs(sin(thetaMc[_theGenPos]));
      }
      else {
        _highestPtGen[0]=pMc[_theGenPos]*fabs(sin(thetaMc[_theGenPos]));
        _lowestPtGen[0]=pMc[_theGenEle]*fabs(sin(thetaMc[_theGenEle]));
      }
    }
  }
}


int HiggsMLSelection::numJets() {

  int num=0;
  m_goodJets.clear();
  for(int j=0;j<nAK5PFJet;j++) {

    // check if the electron/muon falls into the jet
    TVector3 p3Jet(pxAK5PFJet[j],pyAK5PFJet[j],pzAK5PFJet[j]);
    if ( theElectron > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4ElectronMinus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      // taking from ee config file, but jets veto is the same for all the channels
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) continue;
    }

    if ( thePositron > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4ElectronPlus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR) ) continue;
    }

    if ( theMuonMinus > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4MuonMinus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR) ) continue;
    }

    if ( theMuonPlus > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4MuonPlus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR) ) continue;
    }

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFJet[j]))) continue;

    if(_selectionEE->getSwitch("etJetAcc") && !_selectionEE->passCut("etJetAcc", GetPt(pxAK5PFJet[j],pyAK5PFJet[j]))) continue;

    m_goodJets.push_back(j);
    num++;
    
  }

  return num;
}


int HiggsMLSelection::numUncorrJets() {

  int num=0;

  for(int j=0;j<nAK5PFJet;j++) {

    float uncorrEt = uncorrEnergyAK5PFJet[j]*fabs(sin(thetaAK5PFJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFJet[j],phiAK5PFJet[j],uncorrEnergyAK5PFJet[j]);
    TVector3 p3Jet = p4Jet.Vect();
    
    if ( theElectron > -1 ) {
      float deltaR = p3Jet.DeltaR( m_p4ElectronMinus->Vect() );
      H_deltaRuncorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) continue;
    }
    if ( thePositron > -1 ) {
      float deltaR = p3Jet.DeltaR( m_p4ElectronPlus->Vect() );
      H_deltaRuncorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) continue;
    }

    if ( theMuonMinus > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4MuonMinus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR) ) continue;
    }

    if ( theMuonPlus > -1 ) {
      float deltaR =  fabs( p3Jet.DeltaR( m_p4MuonPlus->Vect() ) );
      H_deltaRcorr -> Fill(deltaR);
      if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR) ) continue;
    }

    if(_selectionEE->getSwitch("etaJetAcc")      && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFJet[j]))) continue;    
    if(_selectionEE->getSwitch("etUncorrJetAcc") && !_selectionEE->passCut("etUncorrJetAcc", uncorrEt))   continue;
    
    num++;
  }
  
  return num;
}

int HiggsMLSelection::numSoftMuons() {
  int num = 0;
  for(int i=0; i<nMuon; ++i) {
    if(i==theMuonMinus || i==theMuonPlus) continue;
    if(GetPt(pxMuon[i],pyMuon[i]) < 3.0) continue;
    Utils anaUtils;
    if(!anaUtils.muonIdVal(muonIdMuon[i],AllTrackerMuons) ||
       !anaUtils.muonIdVal(muonIdMuon[i],TMLastStationTight)) continue;
    int track = trackIndexMuon[i];
    if(trackValidHitsTrack[track]<=10) continue;
    float dxy = fabs(trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                trackVxTrack[track], trackVyTrack[track], trackVzTrack[track], 
                                pxTrack[track], pyTrack[track], pzTrack[track]));
    if(dxy > 0.200) continue;
    num++;
  }
  return num;
}

int HiggsMLSelection::numExtraLeptons() {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    if(i==theElectron || i==thePositron) continue;
    if(_preselection->getSwitch("etaElectronAcc") && !_preselection->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_preselection->getSwitch("ptElectronAcc") && !_preselection->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])) ) continue;
    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;

    if (!_preselection->getSwitch("asymmetricID")) 
      isEleID(i,&theId,&theIso,&theConvRej,EgammaCutBasedID);

    if (_preselection->getSwitch("asymmetricID")) {
      float pt = GetPt(pxEle[i],pyEle[i]);	
      if(pt>=20) isEleID(i,&theId,&theIso,&theConvRej,EgammaCutBasedID);
      if(pt<20)  isEleID(i,&theId,&theIso,&theConvRej,EgammaCutBasedIDLow);
    }

    if(!theId || !theIso || !theConvRej) continue;
    int track = gsfTrackIndexEle[i];
    float dxy = fabs(trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                trackVxTrack[track], trackVyTrack[track], trackVzTrack[track], 
                                pxTrack[track], pyTrack[track], pzTrack[track]));
    if(_selectionEE->getSwitch("leptonD0") && !_selectionEE->passCut("leptonD0",dxy)) continue;
    numEle++;
  }

  int numMu = 0;
  for(int i=0; i<nMuon; ++i) {
    if(i==theMuonMinus || i==theMuonPlus) continue;
    if(_preselection->getSwitch("etaMuonAcc") && !_preselection->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_preselection->getSwitch("ptMuonAcc") && !_preselection->passCut("ptMuonAcc",GetPt(pxMuon[i],pyMuon[i])) ) continue;
    bool theId = true;
    isMuonID(i,&theId);
    if(!theId) continue;
    float isoSum = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i];
    if(_selectionMM->getSwitch("globalSum") && !_selectionMM->passCut("globalSum",isoSum)) continue;
    int track = trackIndexMuon[i];
    float dxy = fabs(trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                trackVxTrack[track], trackVyTrack[track], trackVzTrack[track], 
                                pxTrack[track], pyTrack[track], pzTrack[track]));
    if(_selectionEE->getSwitch("leptonD0") && !_selectionEE->passCut("leptonD0",dxy)) continue;
    numMu++;
  }

  return numEle + numMu;

}

void HiggsMLSelection::resetKinematics() {

  theElectron = -1;
  thePositron = -1;
  theMuonMinus = -1;
  theMuonPlus = -1;
  m_p4ElectronMinus->SetXYZT(0,0,0,0);
  m_p4ElectronPlus->SetXYZT(0,0,0,0);
  m_p4MuonMinus->SetXYZT(0,0,0,0);
  m_p4MuonPlus->SetXYZT(0,0,0,0);
  m_p4MET->SetXYZT(0,0,0,0);
  m_HoEElectronMinus = 0;
  m_HoEElectronPlus = 0;
  m_CaloEneElectronMinus = 0;
  m_CaloEneElectronPlus = 0;

  for(int ichan=0;ichan<3;ichan++) {
    m_deltaPhi[ichan] = 0;
    m_deltaErre[ichan] = 0;
    m_mll[ichan] = 0;
    m_transvMass[ichan] = 0;
    m_projectedMet[ichan] = 0;
  }
  
  hardestElectronPt = 0;
  hardestMuonPt = 0;
  slowestElectronPt = 0;
  slowestMuonPt = 0;

}

void HiggsMLSelection::setEleIdVariables(int hard, int slow) {

  int selectedElectrons[2];
  selectedElectrons[0] = hard;
  selectedElectrons[1] = slow; 

  for(int i = 0; i < 2 && selectedElectrons[i] > -1; i++) {
    int eleIndex = selectedElectrons[i];
    myRecoflag[i] = recoFlagsEle[eleIndex];
    myPt[i] = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
    myEta[i] = etaEle[eleIndex];
    myPhi[i] = phiEle[eleIndex];
    myClassification[i] = classificationEle[eleIndex];
    myNBremClusters[i] = nbremsEle[eleIndex];
    myDeta[i] = deltaEtaAtVtxEle[eleIndex];
    myDphi[i] = deltaPhiAtVtxEle[eleIndex];
    myHoe[i] = hOverEEle[eleIndex];
    int sc = superClusterIndexEle[eleIndex];
    mySee[i] = SigmaiEiE(eleIndex);
    mySpp[i] = SigmaiPiP(eleIndex);
    myEop[i] = eSuperClusterOverPEle[eleIndex];
    myFbrem[i] = fbremEle[eleIndex];
    myTrackerIso[i] = dr03TkSumPtEle[eleIndex];
    myHcalIso[i] = dr03HcalTowerSumEtEle[eleIndex];
    myEcalJIso[i] = dr03EcalRecHitSumEtEle[eleIndex];
    myEcalGTIso[i] = scBasedEcalSum04Ele[eleIndex];
    myCombinedIso[i] = (dr03TkSumPtEle[eleIndex] + 
                        TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + 
                        dr03HcalTowerSumEtEle[eleIndex]) / myPt[i];
    myCharge[i] = chargeEle[eleIndex];
    int gsf = gsfTrackIndexEle[eleIndex];
    myMissHits[i] = expInnerLayersGsfTrack[gsf];
    myDist[i] = convDistEle[eleIndex];
    myDcot[i] = convDcotEle[eleIndex];
    myLh[i] = likelihoodRatio(eleIndex,*LH);
 
    // match with MC truth
    myMatched[i] = 999;
    if ( !_preselection->getSwitch("isData") ) { 
      int matchedReco = 0;
      TVector3 pReco(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
      for (int ii=0; ii<nMc; ii++) {
	TVector3 Welegen;
 	if ( (fabs(idMc[mothMc[ii]])==24) && (fabs(idMc[ii])==11) ) {
 	  Welegen.SetMagThetaPhi(pMc[ii],thetaMc[ii],phiMc[ii]);
 	  float dRmatch = pReco.DeltaR(Welegen);  	
 	  if (fabs(dRmatch)<0.3) matchedReco = 1;
 	}}
      myMatched[i] = matchedReco;
    }
  }
}

float HiggsMLSelection::getSecondEleTkPt(TVector3 firstLepton, int second, float deltaR) {

  TVector3 secondEle(pxEle[second],pyEle[second],pzEle[second]);

  float secondEleTrackPt = 0.0;
  float dr = firstLepton.DeltaR(secondEle);

  if( dr < deltaR ) { 
    secondEleTrackPt = secondEle.Pt();
  }

  return secondEleTrackPt;

}

float HiggsMLSelection::getSecondMuonTkPt(TVector3 firstLepton, int second, float deltaR) {

  TVector3 secondMuon(pxMuon[second],pyMuon[second],pzMuon[second]);

  float secondMuonTrackPt = 0.0;
  float dr = firstLepton.DeltaR(secondMuon);

  if( dr < deltaR ) { 
    secondMuonTrackPt = secondMuon.Pt();
  }

  return secondMuonTrackPt;

}

float HiggsMLSelection::getEcalPtSum(int index) {

  float ptsum = 0.0;

  // in barrel, use isolation by SC footprint removal
  if( classificationEle[index] < 100 ) ptsum = scBasedEcalSum04Ele[index];
  // in endcaps, use jurassic isolation
  else ptsum = dr04EcalRecHitSumEtEle[index];

  return ptsum;

}

float HiggsMLSelection::muonIsoGlobalSum(int theMuon, int theOther) {

  TVector3 muonP3(pxMuon[theMuon],pyMuon[theMuon],pzMuon[theMuon]);
  TVector3 otherP3(pxMuon[theOther],pyMuon[theOther],pzMuon[theOther]);

  // global isolation
  float theMuonGlobalSum = 0.;
  float muonTrackerForGlobal = sumPt03Muon[theMuon];
  float muonEcalForGlobal    = emEt03Muon[theMuon]; 
  float muonHcalForGlobal    = hadEt03Muon[theMuon];

  // if the second muon falls into the cone, remove it from tracker only
  float theDr = muonP3.DeltaR(otherP3);
  if (theDr<0.3) muonTrackerForGlobal = muonTrackerForGlobal-otherP3.Pt(); 

  theMuonGlobalSum = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal;

  return theMuonGlobalSum;

}

float HiggsMLSelection::mueleIsoGlobalSum(int theMuon, int theOtherEle) {

  TVector3 muonP3(pxMuon[theMuon],pyMuon[theMuon],pzMuon[theMuon]);
  TVector3 otherEleP3(pxEle[theOtherEle],pyEle[theOtherEle],pzEle[theOtherEle]);

  // global isolation
  float theMuonGlobalSum = 0.;
  float muonTrackerForGlobal = sumPt03Muon[theMuon];
  float muonEcalForGlobal    = emEt03Muon[theMuon]; 
  float muonHcalForGlobal    = hadEt03Muon[theMuon];

  // if the electron falls into the cone, remove it
  float theDr = muonP3.DeltaR(otherEleP3);
  if (theDr<0.3) { 
    
    // from tracker: use the KF track closer to the GSF track
    int kftrack = trackIndexEle[theOtherEle];
    // check if there is the linked KF track (not true in 100% of cases)
    if(kftrack>=0 && kftrack<nTrack) muonTrackerForGlobal -= GetPt(pxTrack[kftrack],pyTrack[kftrack]);

    // from ECAL: use the supercluster energy
    // int sc = superClusterIndexEle[theOtherEle];
    // muonEcalForGlobal = muonEcalForGlobal-rawEnergySC[sc]; 
  }
  
  theMuonGlobalSum = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal;
  
  return theMuonGlobalSum;

}


int HiggsMLSelection::getPV() {
  // search for hardest PV of the event
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

double HiggsMLSelection::trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}

double HiggsMLSelection::trackDszPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  float eleP  = sqrt(elePx*elePx + elePy*elePy + elePz*elePz);
  return (eleVz-PVz)*elePt/eleP - ((eleVx-PVx)*elePx+(eleVy-PVy)*elePy)/elePt *elePz/eleP;
}

bool HiggsMLSelection::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = p3Track.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(p3Track.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

std::vector<float> HiggsMLSelection::jetBTagVariables(int jetIndex) {

  float nTracks   = 0;
  float sumNumDxy = 0.;
  float sumNumDsz = 0;
  float sumDen    = 0.;

  TVector3 p3Jet(pxAK5PFJet[jetIndex],pyAK5PFJet[jetIndex],pzAK5PFJet[jetIndex]);
  TLorentzVector p4TracksInJet(0.,0.,0.,0.);

  for(int iTrack=0; iTrack<nTrack; iTrack++) {
    if (!isGoodTrack(iTrack,0.5,500,20,2.4,5)) continue;
    TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
    TLorentzVector p4Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack],p3Track.Mag());      // assume mass=0...

    float deltaR = p3Jet.DeltaR(p3Track);
    if(fabs(deltaR)<0.5){

      nTracks++;

      float weight = p3Track.Pt() * p3Track.Pt() * p3Track.Pt() * p3Track.Pt();

      float dxy = fabs(trackDxyPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                  trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack], 
                                  pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]));
      float w_dxy = dxy*weight;
      float dsz = fabs(trackDszPV(PVxPV[m_closestPV], PVyPV[m_closestPV], PVzPV[m_closestPV], 
                                  trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack], 
                                  pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]));
      float w_dsz = dsz*weight;

      sumNumDxy = sumNumDxy + w_dxy;
      sumNumDsz = sumNumDsz + w_dsz;
      sumDen    = sumDen + weight;
      
      p4TracksInJet += p4Track;
    }
  }
  
  float DxyAverage = 0.0;
  float DszAverage = 0.0;

  if (sumDen != 0) {
    DxyAverage = sumNumDxy / sumDen;
    DszAverage = sumNumDsz / sumDen;
  }

  float jetMass = p4TracksInJet.M();

  std::vector<float> variables;
  variables.clear();
  variables.push_back(DxyAverage);
  variables.push_back(DszAverage);
  variables.push_back(trackCountingHighEffBJetTagsAK5PFJet[jetIndex]);
  variables.push_back(jetProbabilityBJetTagsAK5PFJet[jetIndex]);
  variables.push_back(combinedSecondaryVertexMVABJetTagsAK5PFJet[jetIndex]);
  variables.push_back(jetMass);
  variables.push_back(nTracks);

  return variables;

}

void HiggsMLSelection::calcEventBVetoVariables(std::vector<int> jets) {
  
  m_maxDxyEvt=-1000; 
  m_maxDszEvt=-1000;
  m_maxTrackCountingHighEffBJetTags=-1000;
  m_maxImpactParameterMVABJetTags=-1000;
  m_maxCombinedSecondaryVertexMVABJetTags=-1000;

  for (unsigned int iJet=0; iJet<jets.size(); iJet++){     
    int thisJet = jets[iJet];
    std::vector<float> variables = jetBTagVariables(thisJet);    
    float dxy = variables[0];
    float dsz = variables[1];
    float trackCounting = variables[2];
    float impactParameter = variables[3];
    float combinedSV = variables[4];
    if (dxy > m_maxDxyEvt) m_maxDxyEvt=dxy;
    if (dsz > m_maxDszEvt) m_maxDszEvt=dsz;
    if (trackCounting > m_maxTrackCountingHighEffBJetTags ) m_maxTrackCountingHighEffBJetTags = trackCounting;
    if (impactParameter > m_maxImpactParameterMVABJetTags ) m_maxImpactParameterMVABJetTags = impactParameter;
    if (combinedSV > m_maxCombinedSecondaryVertexMVABJetTags ) m_maxCombinedSecondaryVertexMVABJetTags = combinedSV;
  }

}

double HiggsMLSelection::mT(TVector3 plep, TVector3 pneu) {

  TVector3 pTlep(plep.X(),plep.Y(),0.0);
  TVector3 pTneu(pneu.X(),pneu.Y(),0.0);

  return sqrt(2 * (pTlep.Mag()*pTneu.Mag() - pTlep*pTneu));

}

double HiggsMLSelection::mT2(TVector3 plep1, TVector3 plep2, TVector3 ptmiss) {

  // need an external dependency: removed right now
//   Mt2::TwoVector pTlep1(plep1.X(),plep1.Y());
//   Mt2::TwoVector pTlep2(plep2.X(),plep2.Y());
//   Mt2::TwoVector pTmiss(ptmiss.X(),ptmiss.Y());
//   double invis_mass = 0.0;

//   Mt2::SUSYPhys_Mt2_222_Calculator mt2Calculator;

//   // Could tell the MT2 calculating object to be verbose, and print out
//   // debug messages while it is thinking ... but we won't:

//   mt2Calculator.setDebug(false);

//   // Now we can actually calculate MT2:
//   double mt2 = mt2Calculator.mt2_222( pTlep1, pTlep2, pTmiss, invis_mass);

  double mt2 = 0.0;
  return mt2;

}

float HiggsMLSelection::GetProjectedMet(TVector3 p1, TVector3 p2) {

  TVector3 met = m_p4MET->Vect();
  float deltaPhi1 = fabs(p1.DeltaPhi(met));
  float deltaPhi2 = fabs(p2.DeltaPhi(met));
  float deltaphi = TMath::Min(deltaPhi1,deltaPhi2);
  if(deltaphi<TMath::Pi()/2.) return met.Mag() * sin(deltaphi);
  else return met.Mag();

}
