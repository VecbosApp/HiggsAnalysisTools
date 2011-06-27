#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/LeptonPlusFakeMLSelection_fullME.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <TTree.h>

using namespace bits;

LeptonPlusFakeMLSelection_fullME::LeptonPlusFakeMLSelection_fullME(TTree *tree) 
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
  
  // selection efficiencies: counters
  std::string fileCutsME     = higgsConfigDirMass + "emu2nuCuts.txt";
  std::string fileSwitchesME = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionME.   Configure(fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION EVENT COUNTER ME FP"); 
  CutBasedHiggsSelectionME_FF.Configure(fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION EVENT COUNTER ME FF");

  // selection efficiencies: stat errors due to number of events
  CutBasedHiggsSelectionStatME.Configure   (fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION STAT ERRORS ME FP");
  CutBasedHiggsSelectionStatME_FF.Configure(fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION STAT ERRORS ME FF");

  // selection efficiencies: stat errors due to fake rate 
  CutBasedHiggsErrorsSelectionME.Configure   (fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION FAKE ERRORS ME");
  CutBasedHiggsErrorsSelectionME_FF.Configure(fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION FAKE ERRORS ME FF");

  // taking the selection                                                                                                                                                    
  _selectionME        = CutBasedHiggsSelectionME.GetSelection();
  _selectionME_FF     = CutBasedHiggsSelectionME_FF.GetSelection();
  _selectionStatME    = CutBasedHiggsSelectionStatME.GetSelection();
  _selectionStatME_FF = CutBasedHiggsSelectionStatME_FF.GetSelection();
  _selectionErrME     = CutBasedHiggsErrorsSelectionME.GetSelection();
  _selectionErrME_FF  = CutBasedHiggsErrorsSelectionME_FF.GetSelection();

  //  extra selection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  _selectionME->addCut("etaElectronAcc");    
  _selectionME->addCut("ptElectronAcc");
  _selectionME->addCut("etaMuonAcc");
  _selectionME->addCut("ptMuonAcc");
  _selectionME->addCut("etUncorrJetAcc");  
  _selectionME->addSwitch("apply_kFactor");   
  _selectionME->addSwitch("isData");
  _selectionME->addSwitch("goodRunLS");
  _selectionME->addSwitch("asymmetricID");
  _selectionME->addStringParameter("electronIDType");
  _selectionME->addStringParameter("electronIDTypeLow");  

  // cut based electron id or likelihood 
  TString selectionString(_selectionME->getStringParameter("electronIDType"));
  if (!_selectionME->getSwitch("asymmetricID")) 
    cout << "=== CONFIGURING " << selectionString << " CUT BASED SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/higgs/electronId/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/higgs/electronId/");
  
  if (_selectionME->getSwitch("asymmetricID")) {
    TString selectionStringLow (_selectionME->getStringParameter("electronIDTypeLow"));
    cout << "=== CONFIGURING "  << selectionStringLow << " and " 
	 << selectionString << " for CUT BASED ASYMMETRIC ELECTRON ID ===" << endl;
    EgammaCutBasedIDLow.ConfigureNoClass("config/higgs/electronId/"+selectionStringLow);
    EgammaCutBasedIDLow.ConfigureEcalCleaner("config/higgs/electronId/");
  }
  
  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir  = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir  = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useOneOverEMinusOneOverP = true;
  LH = new ElectronLikelihood(&(*EB0lt15dir), &(*EB1lt15dir), &(*EElt15dir), &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir), defaultSwitches, std::string("class"),std::string("class"),true,true);
  
  // Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _selectionME->getSwitch("goodRunLS") << " isData is " <<  _selectionME->getSwitch("isData") << std::endl;

  // To read good run list!
  if (_selectionME->getSwitch("goodRunLS") && _selectionME->getSwitch("isData")) {
    std::string goodRunJsonFile = "config/json/goodCollisions2011.json";   // chiara
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics
  for(int theChannel=0; theChannel<1; theChannel++) {
    m_p4LeptonPlus[theChannel]  = new TLorentzVector(0.,0.,0.,0.);
    m_p4LeptonMinus[theChannel] = new TLorentzVector(0.,0.,0.,0.);
  }

  // met
  m_p3PFMET = new TVector3(0.,0.,0.);
  m_p3TKMET = new TVector3(0.,0.,0.);

  // b-veto event variables
  m_maxDxyEvt = 0.0;
  m_maxDszEvt = 0.0;
}

LeptonPlusFakeMLSelection_fullME::~LeptonPlusFakeMLSelection_fullME(){

  for(int theChannel=0; theChannel<1; theChannel++) {  
    delete m_p4LeptonPlus[theChannel];
    delete m_p4LeptonMinus[theChannel];
  }
  delete m_p3PFMET;  
  delete m_p3TKMET;  
  
  delete _selectionME;
  delete _selectionME_FF;
  delete _selectionErrME;
  delete _selectionErrME_FF;
  delete _selectionStatME;
  delete _selectionStatME_FF;

  myOutTreeME -> save();
}

// chiara
void LeptonPlusFakeMLSelection_fullME::initialiseElectronFakeRate() {

  // binning                                      
  m_eleMinFakePt[0] = 10.;   m_eleMaxFakePt[0] = 15.;              
  m_eleMinFakePt[1] = 15.;   m_eleMaxFakePt[1] = 20.;
  m_eleMinFakePt[2] = 20.;   m_eleMaxFakePt[2] = 25.;
  m_eleMinFakePt[3] = 25.;   m_eleMaxFakePt[3] = 50.;
  m_eleMinFakePt[4] = 50.;   m_eleMaxFakePt[4] = 10000.;

  /*
  // fake in the barrel from Smurf selection, QCD
  m_eleFakeRateEB[0] = 0.00154351;
  m_eleFakeRateEB[1] = 0.00145805;   
  m_eleFakeRateEB[2] = 0.0111048;  
  m_eleFakeRateEB[3] = 0.0120757; 
  m_eleFakeRateEB[4] = 0.0163262;
  
  m_eleFakeRateEB_err[0] = 0.000558265;
  m_eleFakeRateEB_err[1] = 0.000357493;
  m_eleFakeRateEB_err[2] = 0.00101031; 
  m_eleFakeRateEB_err[3] = 0.00101347;
  m_eleFakeRateEB_err[4] = 0.00456665;
  
  // fake in the endcap from Smurf selection, QCD
  m_eleFakeRateEE[0] = 0.00149395;
  m_eleFakeRateEE[1] = 0.00230214;
  m_eleFakeRateEE[2] = 0.0127027;
  m_eleFakeRateEE[3] = 0.0136544;
  m_eleFakeRateEE[4] = 0.0163247;
  
  m_eleFakeRateEE_err[0] = 0.000509235; 
  m_eleFakeRateEE_err[1] = 0.000404248;
  m_eleFakeRateEE_err[2] = 0.000981551;
  m_eleFakeRateEE_err[3] = 0.00100148;
  m_eleFakeRateEE_err[4] = 0.00537643;
  */

  // fake in the barrel from Smurf selection, ET>30
  m_eleFakeRateEB[0] = 0.00832695;  // here we subtract W and Z with lumi corresponding to HTLT8
  m_eleFakeRateEB[1] = 0.00907589;  // here we subtract W and Z with lumi corresponding to HTLT17
  m_eleFakeRateEB[2] = 0.0140018;   // " " 
  m_eleFakeRateEB[3] = 0.0136831;   // " " 
  m_eleFakeRateEB[4] = 0.0140156;   // " " 

  m_eleFakeRateEB_err[0] = 0.00140527;
  m_eleFakeRateEB_err[1] = 0.00107052;
  m_eleFakeRateEB_err[2] = 0.00129764;
  m_eleFakeRateEB_err[3] = 0.00116528;
  m_eleFakeRateEB_err[4] = 0.00385747;

  // fake in the endcap from Smurf selection, ET>30
  m_eleFakeRateEE[0] = 0.00322436;  // here we subtract W and Z with lumi corresponding to HTLT8
  m_eleFakeRateEE[1] = 0.0050205;   // here we subtract W and Z with lumi corresponding to HTLT17
  m_eleFakeRateEE[2] = 0.0108465;   // " " 
  m_eleFakeRateEE[3] = 0.00940649;  // " " 
  m_eleFakeRateEE[4] = 0.0163363;   // " " 
  
  m_eleFakeRateEE_err[0] = 0.000972436;
  m_eleFakeRateEE_err[1] = 0.000892635;
  m_eleFakeRateEE_err[2] = 0.0011173;
  m_eleFakeRateEE_err[3] = 0.000971426;
  m_eleFakeRateEE_err[4] = 0.00557363;

  /*
  // fake in the barrel from Smurf selection, ET>15
  m_eleFakeRateEB[0] = 0.0139881;    // here we subtract W and Z with lumi corresponding to HTLT8
  m_eleFakeRateEB[1] = 0.0142939;    // here we subtract W and Z with lumi corresponding to HTLT17
  m_eleFakeRateEB[2] = 0.0229292;    // " " 
  m_eleFakeRateEB[3] = 0.0167388;
  m_eleFakeRateEB[4] = 0.014033;

  m_eleFakeRateEB_err[0] = 0.000822471;
  m_eleFakeRateEB_err[1] = 0.000839656;
  m_eleFakeRateEB_err[2] = 0.00120164;
  m_eleFakeRateEB_err[3] = 0.00111492;
  m_eleFakeRateEB_err[4] = 0.00376642;

  // fake in the endcap from Smurf selection, ET>15
  m_eleFakeRateEE[0] = 0.00536005;   // here we subtract W and Z with lumi corresponding to HTLT8
  m_eleFakeRateEE[1] = 0.00617585;   // here we subtract W and Z with lumi corresponding to HTLT17
  m_eleFakeRateEE[2] = 0.0153299;    // " " 
  m_eleFakeRateEE[3] = 0.011276;
  m_eleFakeRateEE[4] = 0.015083;
  
  m_eleFakeRateEE_err[0] = 0.000538545;
  m_eleFakeRateEE_err[1] = 0.000566034;
  m_eleFakeRateEE_err[2] = 0.00088218;
  m_eleFakeRateEE_err[3] = 0.000863319;
  m_eleFakeRateEE_err[4] = 0.00521618;
  */

  /*
  // fake in the barrel from Smurf selection, ET>50
  m_eleFakeRateEB[0] = 0.00196092;
  m_eleFakeRateEB[1] = 0.00217957;
  m_eleFakeRateEB[2] = 0.0112285;
  m_eleFakeRateEB[3] = 0.00758858;
  m_eleFakeRateEB[4] = 0.00936297;

  m_eleFakeRateEB_err[0] = 0.00198052;
  m_eleFakeRateEB_err[1] = 0.00130578;
  m_eleFakeRateEB_err[2] = 0.00239751;
  m_eleFakeRateEB_err[3] = 0.00130281;
  m_eleFakeRateEB_err[4] = 0.00335611;

  // fake in the endcap from Smurf selection, ET>50
  m_eleFakeRateEE[0] = 0.00258585;
  m_eleFakeRateEE[1] = 0.00439881;
  m_eleFakeRateEE[2] = 0.00564385;
  m_eleFakeRateEE[3] = 0.00496651;
  m_eleFakeRateEE[4] = 0.0193293;
  
  m_eleFakeRateEE_err[0] = 0.00258838;
  m_eleFakeRateEE_err[1] = 0.00221489;
  m_eleFakeRateEE_err[2] = 0.00188719;
  m_eleFakeRateEE_err[3] = 0.00123241;
  m_eleFakeRateEE_err[4] = 0.00689591;
  */
}

// numbers from Emanuele - ~150 / pb
void LeptonPlusFakeMLSelection_fullME::initialiseElectronPromptRate() {

  // binning                                      
  m_eleMinPromptPt[0] = 10.;   m_eleMaxPromptPt[0] = 15.;                    
  m_eleMinPromptPt[1] = 15.;   m_eleMaxPromptPt[1] = 20.;
  m_eleMinPromptPt[2] = 20.;   m_eleMaxPromptPt[2] = 25.;
  m_eleMinPromptPt[3] = 25.;   m_eleMaxPromptPt[3] = 50.;
  m_eleMinPromptPt[4] = 50.;   m_eleMaxPromptPt[4] = 10000.;

  // prompt in the barrel
  m_elePromptRateEB[0] = 0.917;     
  m_elePromptRateEB[1] = 0.887;
  m_elePromptRateEB[2] = 0.911;
  m_elePromptRateEB[3] = 0.94380;
  m_elePromptRateEB[4] = 0.962;
  
  m_elePromptRateEB_err[0] = 0.018;
  m_elePromptRateEB_err[1] = 0.009;
  m_elePromptRateEB_err[2] = 0.005;
  m_elePromptRateEB_err[3] = 0.00003;
  m_elePromptRateEB_err[4] = 0.002;
  
  // prompt in the endcap
  m_elePromptRateEE[0] = 0.851;
  m_elePromptRateEE[1] = 0.817;
  m_elePromptRateEE[2] = 0.8546;
  m_elePromptRateEE[3] = 0.9114;
  m_elePromptRateEE[4] = 0.9482;

  m_elePromptRateEE_err[0] = 0.021;
  m_elePromptRateEE_err[1] = 0.011;
  m_elePromptRateEE_err[2] = 0.0077;
  m_elePromptRateEE_err[3] = 0.0009;
  m_elePromptRateEE_err[4] = 0.0041;
}

void LeptonPlusFakeMLSelection_fullME::initialiseMuonFakeRate() {

  // binning                                      
  m_muonMinFakePt[0] = 10.;   m_muonMaxFakePt[0] = 15.;                
  m_muonMinFakePt[1] = 15.;   m_muonMaxFakePt[1] = 20.;
  m_muonMinFakePt[2] = 20.;   m_muonMaxFakePt[2] = 25.;
  m_muonMinFakePt[3] = 25.;   m_muonMaxFakePt[3] = 50.;
  m_muonMinFakePt[4] = 50.;   m_muonMaxFakePt[4] = 10000.;

  // fake rate for muons - barrel
  m_muonFakeRateEB[0] = 0.124862;
  m_muonFakeRateEB[1] = 0.0851988;
  m_muonFakeRateEB[2] = 0.124773; 
  m_muonFakeRateEB[3] = 0.0882528; 
  m_muonFakeRateEB[4] = 0.0256791;   
  
  m_muonFakeRateEB_err[0] = 0.00264137;
  m_muonFakeRateEB_err[1] = 0.00549935;
  m_muonFakeRateEB_err[2] = 0.0146199;
  m_muonFakeRateEB_err[3] = 0.0166693;
  m_muonFakeRateEB_err[4] = 0.0756251;
  
  // fake rate for muons - endcap
  m_muonFakeRateEE[0] = 0.166006;
  m_muonFakeRateEE[1] = 0.12626;
  m_muonFakeRateEE[2] = 0.174822;
  m_muonFakeRateEE[3] = 0.134711;    
  m_muonFakeRateEE[4] = 0.318523;
  
  m_muonFakeRateEE_err[0] = 0.00561667;
  m_muonFakeRateEE_err[1] = 0.0140856;
  m_muonFakeRateEE_err[2] = 0.0334133;
  m_muonFakeRateEE_err[3] = 0.0535601;
  m_muonFakeRateEE_err[4] = 0.;
}

void LeptonPlusFakeMLSelection_fullME::initialiseMuonPromptRate() {

  // binning                                      
  m_muonMinPromptPt[0] = 10.;   m_muonMaxPromptPt[0] = 15.;               
  m_muonMinPromptPt[1] = 15.;   m_muonMaxPromptPt[1] = 20.;
  m_muonMinPromptPt[2] = 20.;   m_muonMaxPromptPt[2] = 25.;
  m_muonMinPromptPt[3] = 25.;   m_muonMaxPromptPt[3] = 50.;
  m_muonMinPromptPt[4] = 50.;   m_muonMaxPromptPt[4] = 10000.;


  // prompt in the barrel
  m_muonPromptRateEB[0] = 0.750083;  
  m_muonPromptRateEB[1] = 0.796289; 
  m_muonPromptRateEB[2] = 0.891641; 
  m_muonPromptRateEB[3] = 0.931356; 
  m_muonPromptRateEB[4] = 0.98458;    
  
  m_muonPromptRateEB_err[0] = 0.00762918;
  m_muonPromptRateEB_err[1] = 0.00711569;
  m_muonPromptRateEB_err[2] = 0.00669801;
  m_muonPromptRateEB_err[3] = 0.00433128;
  m_muonPromptRateEB_err[4] = 0.00536175;
  
  // prompt in the endcap
  m_muonPromptRateEE[0] = 0.756677; 
  m_muonPromptRateEE[1] = 0.807518; 
  m_muonPromptRateEE[2] = 0.878232;
  m_muonPromptRateEE[3] = 0.91396;   
  m_muonPromptRateEE[4] = 0.976088;

  m_muonPromptRateEE_err[0] = 0.0117427;
  m_muonPromptRateEE_err[1] = 0.0106653;
  m_muonPromptRateEE_err[2] = 0.0117499;
  m_muonPromptRateEE_err[3] = 0.00741766;
  m_muonPromptRateEE_err[4] = 0.0110794;
}

void LeptonPlusFakeMLSelection_fullME::Loop() {

  _verbose=false;
  if(fChain == 0) return;
  
  initialiseElectronFakeRate(); 
  initialiseElectronPromptRate();

  initialiseMuonFakeRate();    
  initialiseMuonPromptRate();  

  // kinematics reduced tree
  std::string reducedTreeNameME = _datasetName+"-datasetME.root";
  myOutTreeME = new RedHiggsTree(reducedTreeNameME.c_str());
  myOutTreeME->addMLVars();
  myOutTreeME->addLatinos();

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
    if ( !_selectionME->getSwitch("isData") ) tmpWeight *= fPUWeight->GetWeight(nPU);  

    // Good Run selection
    if (_selectionME->getSwitch("isData") && _selectionME->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_selectionME->getSwitch("isData") && _selectionME->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    
    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask();
    bool passedHLT[1];
    passedHLT[me] = hasPassedHLT(); 


    // -------------------------------------------------------------
    // vertex selection - we only consider the first vertex of the list ( = highest sumPT^2)
    bool isGoodVertex = true;
    if (nPV<1) isGoodVertex = false;
    float rhoVtx = sqrt(PVxPV[0]*PVxPV[0] + PVyPV[0]*PVyPV[0]);
    if ( isFakePV[0] )       isGoodVertex = false;
    if ( ndofPV[0]<4 )       isGoodVertex = false;
    if ( fabs(PVzPV[0])>24.) isGoodVertex = false;
    if ( rhoVtx>2 )          isGoodVertex = false; 


    // -------------------------------------------------------------
    // get the best electrons and best muons ==> tu be used to select ALL the possible channels at the beginning only
    std::pair<int,int> thePreElectrons      = getBestElectronPair_acceptance();
    std::pair<int,int> thePreMuons          = getBestMuonPair_acceptance();
    std::pair<int,int> theBestAcceptMuonEle = getBestMuonElePair(_acceptEleAll,_acceptMuonsAll);
    thePreElectronME = theBestAcceptMuonEle.second;
    thePreMuonME     = theBestAcceptMuonEle.first;

    // reconstructed channel
    m_channel[me] = false;
    
    // at this level the SELECTED channel should have pT > 10 and > 20. So far, at least 2 leptons with pT >20 and 10 in the event
    if ( thePreElectronME > -1 && thePreMuonME > -1 ) {
      float thisMaxPt  = GetPt(pxMuon[thePreMuonME],pyMuon[thePreMuonME]);
      float thisMinPt  = GetPt(pxEle[thePreElectronME],pyEle[thePreElectronME]);
      if (isGoodVertex && thisMaxPt>20 && thisMinPt>10) m_channel[me] = true;    // fixme: hardcoded
    }
    
    if (_verbose) {
      std::cout << "nEle = "   << nEle << "\tnMuon = " << nMuon << std::endl;
      std::cout << "indices: " << thePreElectronME << " " << thePreMuonME << std::endl;
      std::cout << "chargeEle = " << chargeEle[thePreElectronME] << "\tchargeMuon = " << chargeMuon[thePreMuonME] << std::endl;
      std::cout << "mue = "    << m_channel[me] << std::endl; 
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

    // -------------------------------------------------------------
    // MM candidates: preparing vectors of candidates and selecting the two highest pT mu- and mu+ after each step - to check the 20 GeV cut after 

    // muID, for muons in acceptance
    std::pair<int,int> theBestIdMuon = getBestMuonPair_id(_acceptMuonsAll); 

    // isolation, for identified muons
    std::pair<int,int> theBestIsolMuon = getBestMuonPair_isol(_idMuonsAll); 

    // transverse impact parameter, for isolated muons
    std::pair<int,int> theBestIpMuon = getBestMuonPair_ip(_isolMuonsAll);     

    // muons passing the denominator definition                          
    std::pair<int,int> theBestDenomMuon = getBestMuonPair_denominator();    


    // -------------------------------------------------------------
    // ME candidates: preparing vectors of candidates and selecting the two highest pT ele+- and muon-+ after each step - to check the 20 GeV cut after

    // leptonID, for leptons in acceptance
    std::pair<int,int> theBestIdMuonEle = getBestMuonElePair(_idEleAll,_idMuonsAll);
    
    // isolation, for identified leptons
    std::pair<int,int> theBestIsolMuonEle = getBestMuonElePair(_isolEleAll,_isolMuonsAll);

    // conversion rejection, for isolated leptons
    std::pair<int,int> theBestConvMuonEle = getBestMuonElePair(_convEleAll,_isolMuonsAll);

    // transverse impact parameter, for isolated leptons
    std::pair<int,int> theBestIpMuonEle = getBestMuonElePair(_ipEleAll,_ipMuonsAll);
    
    // passing denominator, for all leptons  
    std::pair<int,int> theBestDenomMuonEle = getBestMuonElePair(_denomEleAll,_denomMuonAll);

    
    // -------------------------------------------------------------
    // the two highest pT electron and muon at this point are those I use for my analysis since they passed the full lepton selection  
    int theMuon     = theBestIpMuonEle.first;
    int theElectron = theBestIpMuonEle.second;

    // 0 - 1 - 2 tight candidates                                 
    bool is0tight = false;
    bool is1tight = false;
    bool is2tight = false;

    bool isRealMu = false;

    // here I have two candidates with opposite sign passing the tight selection: is N2T 
    if (theMuon>-1 && theElectron>-1) {
      theReal  = theMuon;
      theFake  = theElectron;
      is2tight = true;
      isRealMu = true;
    }

    // here I have only 1 candidate passing the tight selection: is N1T              
    if (theMuon>-1 && theElectron<0) {
      theReal  = theMuon;
      theFake  = getBestEleDenominator(theReal);    
      is1tight = true;
      isRealMu = true;
    }
    if (theElectron>-1 && theMuon<0) {
      theReal  = theElectron;
      theFake  = getBestMuDenominator(theReal);   
      is1tight = true;
      isRealMu = false;
    }

    // here I have zero candidates passing the tight selection: is N0T              
    if (theMuon<0 && theElectron<0) {
      theReal  = theBestDenomMuonEle.first;        
      theFake  = theBestDenomMuonEle.second;     
      is0tight = true;
      isRealMu = true;
    }

    // to be sure: I assumed above that tight => loose. If not true, I discard the event....
    // chiara: controlla che con questo dia risultati come prima
    if ( (theElectron>-1 && !isEleDenomFake(theElectron)) || (theMuon>-1 && !isMuonDenomFake(theMuon)) ) {
      theReal = -1;
      theFake = -1;
    } 
    
    // sanity check
    if ( (is0tight && is1tight) || (is0tight && is2tight) || (is1tight && is2tight) ) cout << "questo non puo' succedere mai" << endl;
    
    // set of kinematics: : now I've all the final leptons 
    resetKinematics();


    // MET is an event variable. Independent o the channel
    m_p3PFMET->SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);
    m_p3TKMET->SetXYZ(pxChMetPV[0],pyChMetPV[0],pzChMetPV[0]); // the one associated to the 0th vertex
    m_theMET = m_p3PFMET->Pt();

    setKinematicsME(theReal, theFake);

    // weight with the Fake / Prompt -> L2 probability                
    float theFakePt, theRealPt;
    bool  isFakeBarrel = false;
    bool  isRealBarrel = false;

    if (isRealMu) { 
      theFakePt = GetPt(pxEle[theFake],pyEle[theFake]);
      theRealPt = GetPt(pxMuon[theReal],pyMuon[theReal]);
      if ( fabs(etaEle[theFake])<1.476 ) isFakeBarrel = true;  
      if ( fabs(etaMuon[theReal])<1.5 )  isRealBarrel = true;
    }

    if (!isRealMu) {
      theFakePt = GetPt(pxMuon[theFake],pyMuon[theFake]);
      theRealPt = GetPt(pxEle[theReal],pyEle[theReal]);
      if ( fabs(etaMuon[theFake])<1.5 )  isFakeBarrel = true;  
      if ( fabs(etaEle[theReal])<1.476 ) isRealBarrel = true;
    }

    // do both F-P and F-F analysis
    float weightFP       = 1.;      
    float weightErrorFP  = 1.;
    float weightFF       = 1.;      
    float weightErrorFF  = 1.;

    if ( theFake>-1 && theReal>-1) { 
      
      float fakerate1, fakerateErr1, promptrate1, promptrateErr1;
      float fakerate2, fakerateErr2, promptrate2, promptrateErr2;

      if (isRealMu) {
	fakerate1      = getElectronFakeRate( theFakePt, isFakeBarrel );
	fakerateErr1   = getElectronFakeRateError( theFakePt, isFakeBarrel );
	promptrate1    = getElectronPromptRate( theFakePt, isFakeBarrel );           
	promptrateErr1 = getElectronPromptRateError( theFakePt, isFakeBarrel );      
	
	fakerate2      = getMuonFakeRate( theRealPt, isRealBarrel );
	fakerateErr2   = getMuonFakeRateError( theRealPt, isRealBarrel );
	promptrate2    = getMuonPromptRate( theRealPt, isRealBarrel );
	promptrateErr2 = getMuonPromptRateError( theRealPt, isRealBarrel );
      }

      if (!isRealMu) {
	fakerate1      = getMuonFakeRate( theFakePt, isFakeBarrel );
	fakerateErr1   = getMuonFakeRateError( theFakePt, isFakeBarrel );
	promptrate1    = getMuonPromptRate( theFakePt, isFakeBarrel );           
	promptrateErr1 = getMuonPromptRateError( theFakePt, isFakeBarrel );      
	
	fakerate2      = getElectronFakeRate( theRealPt, isRealBarrel );
	fakerateErr2   = getElectronFakeRateError( theRealPt, isRealBarrel );
	promptrate2    = getElectronPromptRate( theRealPt, isRealBarrel );
	promptrateErr2 = getElectronPromptRateError( theRealPt, isRealBarrel );
      }

      /*
      cout << "-----------------------------" << endl;
      cout << endl;
      cout << is0tight << " " << is1tight << " " << is2tight << " " << isRealMu << endl;
      cout << endl;
      cout << fakerate1   << " " << fakerateErr1   << endl;
      cout << fakerate2   << " " << fakerateErr2   << endl;
      cout << promptrate1 << " " << promptrateErr1 << endl;
      cout << promptrate2 << " " << promptrateErr2 << endl;
      cout << endl;
      */

      float thisPartWeightFP    = 1.;
      float thisPartWeightErrFP = 1.;
      float thisPartWeightFF    = 1.;
      float thisPartWeightErrFF = 1.;
      
      if (is0tight) {
	float myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));
        // for F-P 
        thisPartWeightFP = -( fakerate1*promptrate2*promptrate1*fakerate2 + fakerate2*promptrate1*promptrate2*fakerate1 )*myFactor;
        // for F-F
        thisPartWeightFF = ( promptrate1*promptrate2*fakerate1*fakerate2 )*myFactor;
      }
      
      if (is1tight) {
        double myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));
        // for F-P
        thisPartWeightFP = ( fakerate1*(1-promptrate2)*promptrate1*fakerate2 + promptrate1*(1-fakerate2)*fakerate1*promptrate2 )*myFactor;
        // for F-F
        thisPartWeightFF = -( promptrate1*(1-promptrate2)*fakerate1*fakerate2 )*myFactor;
      }
      
      if (is2tight) {
	double myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));
        // for F-P 
        thisPartWeightFP = -( (promptrate2*fakerate1*(1-promptrate1)*(1-fakerate2)) + (promptrate1*fakerate2*(1-promptrate2)*(1-fakerate1)) )*myFactor;
        // for F-F
        thisPartWeightFF = ( fakerate1*fakerate2*(1-promptrate1)*(1-promptrate2) )*myFactor;
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
    // filling counter for FP ( = W+jets)      
    CutBasedHiggsSelectionME.SetWeight(weightFP);             
    CutBasedHiggsSelectionME.SetMcTruth(1);
    CutBasedHiggsSelectionME.SetHLT(passedHLT[me]);               
    CutBasedHiggsSelectionME.SetIsChannel(m_channel[me]);     
    
    CutBasedHiggsSelectionME.SetElectronId(theReal);
    CutBasedHiggsSelectionME.SetPositronId(theFake);
    CutBasedHiggsSelectionME.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionME.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionME.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionME.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionME.SetElectronIp(theReal);
    CutBasedHiggsSelectionME.SetPositronIp(theFake);

    // checking if the highest pT lepton at the end has pT>20. MU is the hardest in ME
    float thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsSelectionME.SetElectronIp(-1);
      CutBasedHiggsSelectionME.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionME.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsSelectionME.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsSelectionME.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons

    CutBasedHiggsSelectionME.SetNJets(njets[me]);
    CutBasedHiggsSelectionME.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsSelectionME.SetBTagJets(btag[me]);
    CutBasedHiggsSelectionME.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsSelectionME.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsSelectionME.SetMet(m_theMET);					
    CutBasedHiggsSelectionME.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsSelectionME.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsSelectionME.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsSelectionME.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsSelectionME.SetInvMass(m_mll[me]);
    CutBasedHiggsSelectionME.SetDetaLeptons(m_deltaEtaLeptons[me]);
    // CutBasedHiggsSelectionME.SetWWInvMass(2.*m_transvMass[me]/_massVal);
    CutBasedHiggsSelectionME.SetWWInvMass(m_transvMass[me]);

    bool isSelectedME           = CutBasedHiggsSelectionME.output();    
    bool selUpToFinalLeptonsME  = CutBasedHiggsSelectionME.outputUpToFinalLeptons();
    bool selUpToJetVetoME       = CutBasedHiggsSelectionME.outputUpToJetVeto();
    bool selUpToUncorrJetVetoME = CutBasedHiggsSelectionME.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiME       = CutBasedHiggsSelectionME.outputPreDeltaPhi();

    // latinos
    bool outputStep0  = CutBasedHiggsSelectionME.outputStep0();
    bool outputStep1  = CutBasedHiggsSelectionME.outputStep1();
    bool outputStep2  = CutBasedHiggsSelectionME.outputStep2();
    bool outputStep3  = CutBasedHiggsSelectionME.outputStep3();
    bool outputStep4  = CutBasedHiggsSelectionME.outputStep4();
    bool outputStep5  = CutBasedHiggsSelectionME.outputStep5();
    bool outputStep6  = CutBasedHiggsSelectionME.outputStep6();
    bool outputStep7  = CutBasedHiggsSelectionME.outputStep7();
    bool outputStep8  = CutBasedHiggsSelectionME.outputStep8();
    bool outputStep9  = CutBasedHiggsSelectionME.outputStep9();
    bool outputStep10 = CutBasedHiggsSelectionME.outputStep10();
    bool outputStep11 = CutBasedHiggsSelectionME.outputStep11();
    bool outputStep12 = CutBasedHiggsSelectionME.outputStep12();
    bool outputStep13 = CutBasedHiggsSelectionME.outputStep13();
    bool outputStep14 = CutBasedHiggsSelectionME.outputStep14();
    bool outputStep15 = CutBasedHiggsSelectionME.outputStep15();
    bool outputStep16 = CutBasedHiggsSelectionME.outputStep16();
    bool outputStep17 = CutBasedHiggsSelectionME.outputStep17();
    bool outputStep18 = CutBasedHiggsSelectionME.outputStep18();
    bool outputStep19 = CutBasedHiggsSelectionME.outputStep19();
    bool outputStep20 = CutBasedHiggsSelectionME.outputStep20();
    bool outputStep21 = CutBasedHiggsSelectionME.outputStep21();
    bool outputStep22 = CutBasedHiggsSelectionME.outputStep22();
    bool outputStep23 = CutBasedHiggsSelectionME.outputStep23();
    bool outputStep24 = CutBasedHiggsSelectionME.outputStep24();


    // filling the tree
    myOutTreeME->fillRunInfos(runNumber, lumiBlock, eventNumber, weightFP);  

    myOutTreeME -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[me], m_deltaPhi[me], m_deltaErre[me], m_transvMass[me], m_mll[me], 
			   hardestLeptonPt[me], slowestLeptonPt[me], hardestLeptonEta[me], slowestLeptonEta[me],
			   m_deltaEtaLeptons[me], nPV,
			   selUpToFinalLeptonsME, selUpToJetVetoME, selUpToUncorrJetVetoME, selPreDeltaPhiME, isSelectedME);

    myOutTreeME -> fillMLVars(njets[me], nuncorrjets[me], m_maxDxyEvt, m_maxDszEvt, m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);
    
    myOutTreeME -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 );
    
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeME -> store();
    




    // ----------------------------------------------------------------
    // filling counter for FF ( = QCD)      
    CutBasedHiggsSelectionME_FF.SetWeight(weightFF);   
    CutBasedHiggsSelectionME_FF.SetMcTruth(1);
    CutBasedHiggsSelectionME_FF.SetHLT(passedHLT[me]);               
    CutBasedHiggsSelectionME_FF.SetIsChannel(m_channel[me]);   
    CutBasedHiggsSelectionME_FF.SetElectronId(theReal);
    CutBasedHiggsSelectionME_FF.SetPositronId(theFake);
    CutBasedHiggsSelectionME_FF.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionME_FF.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionME_FF.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionME_FF.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionME_FF.SetElectronIp(theReal);
    CutBasedHiggsSelectionME_FF.SetPositronIp(theFake);
    thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsSelectionME_FF.SetElectronIp(-1);
      CutBasedHiggsSelectionME_FF.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionME_FF.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsSelectionME_FF.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsSelectionME_FF.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons
    CutBasedHiggsSelectionME_FF.SetNJets(njets[me]);
    CutBasedHiggsSelectionME_FF.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsSelectionME_FF.SetBTagJets(btag[me]);
    CutBasedHiggsSelectionME_FF.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsSelectionME_FF.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsSelectionME_FF.SetMet(m_theMET);					
    CutBasedHiggsSelectionME_FF.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsSelectionME_FF.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsSelectionME_FF.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsSelectionME_FF.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsSelectionME_FF.SetInvMass(m_mll[me]);
    CutBasedHiggsSelectionME_FF.SetDetaLeptons(m_deltaEtaLeptons[me]);
    CutBasedHiggsSelectionME.SetWWInvMass(m_transvMass[me]);
    bool isSelectedME_FF = CutBasedHiggsSelectionME_FF.output();



    // ----------------------------------------------------------------
    // filling counters for statistical errors due to data / mc - FP      
    float forStatErrFP = weightFP*weightFP;
    CutBasedHiggsSelectionStatME.SetWeight(forStatErrFP);
    CutBasedHiggsSelectionStatME.SetMcTruth(1);
    CutBasedHiggsSelectionStatME.SetHLT(passedHLT[me]);               
    CutBasedHiggsSelectionStatME.SetIsChannel(m_channel[me]);   
    CutBasedHiggsSelectionStatME.SetElectronId(theReal);
    CutBasedHiggsSelectionStatME.SetPositronId(theFake);
    CutBasedHiggsSelectionStatME.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionStatME.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionStatME.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionStatME.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionStatME.SetElectronIp(theReal);
    CutBasedHiggsSelectionStatME.SetPositronIp(theFake);
    thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsSelectionStatME.SetElectronIp(-1);
      CutBasedHiggsSelectionStatME.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionStatME.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsSelectionStatME.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsSelectionStatME.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons
    CutBasedHiggsSelectionStatME.SetNJets(njets[me]);
    CutBasedHiggsSelectionStatME.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsSelectionStatME.SetBTagJets(btag[me]);
    CutBasedHiggsSelectionStatME.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsSelectionStatME.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsSelectionStatME.SetMet(m_theMET);					
    CutBasedHiggsSelectionStatME.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsSelectionStatME.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsSelectionStatME.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsSelectionStatME.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsSelectionStatME.SetInvMass(m_mll[me]);
    CutBasedHiggsSelectionStatME.SetDetaLeptons(m_deltaEtaLeptons[me]);
    CutBasedHiggsSelectionStatME.SetWWInvMass(m_transvMass[me]);
    bool isSelectedStatME = CutBasedHiggsSelectionStatME.output();





    // ----------------------------------------------------------------
    // filling counters for statistical errors due to data / mc - FF      
    float forStatErrFF = weightFF*weightFF;
    CutBasedHiggsSelectionStatME_FF.SetWeight(forStatErrFF);
    CutBasedHiggsSelectionStatME_FF.SetMcTruth(1);
    CutBasedHiggsSelectionStatME_FF.SetHLT(passedHLT[me]);               
    CutBasedHiggsSelectionStatME_FF.SetIsChannel(m_channel[me]);   
    CutBasedHiggsSelectionStatME_FF.SetElectronId(theReal);
    CutBasedHiggsSelectionStatME_FF.SetPositronId(theFake);
    CutBasedHiggsSelectionStatME_FF.SetElectronIsolation(theReal);
    CutBasedHiggsSelectionStatME_FF.SetPositronIsolation(theFake);
    CutBasedHiggsSelectionStatME_FF.SetElectronConvRejection(theReal);
    CutBasedHiggsSelectionStatME_FF.SetPositronConvRejection(theFake);
    CutBasedHiggsSelectionStatME_FF.SetElectronIp(theReal);
    CutBasedHiggsSelectionStatME_FF.SetPositronIp(theFake);
    thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsSelectionStatME_FF.SetElectronIp(-1);
      CutBasedHiggsSelectionStatME_FF.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionStatME_FF.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsSelectionStatME_FF.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsSelectionStatME_FF.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons
    CutBasedHiggsSelectionStatME_FF.SetNJets(njets[me]);
    CutBasedHiggsSelectionStatME_FF.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsSelectionStatME_FF.SetBTagJets(btag[me]);
    CutBasedHiggsSelectionStatME_FF.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsSelectionStatME_FF.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsSelectionStatME_FF.SetMet(m_theMET);					
    CutBasedHiggsSelectionStatME_FF.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsSelectionStatME_FF.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsSelectionStatME_FF.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsSelectionStatME_FF.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsSelectionStatME_FF.SetInvMass(m_mll[me]);
    CutBasedHiggsSelectionStatME_FF.SetDetaLeptons(m_deltaEtaLeptons[me]);
    CutBasedHiggsSelectionStatME_FF.SetWWInvMass(m_transvMass[me]);
    bool isSelectedStatME_FF = CutBasedHiggsSelectionStatME_FF.output();



    // ----------------------------------------------------------------------------------------                                                                              
    // to compute statistical errors due to fake rate - FP                                          
    CutBasedHiggsErrorsSelectionME.SetWeight(weightErrorFP);
    CutBasedHiggsErrorsSelectionME.SetMcTruth(1);
    CutBasedHiggsErrorsSelectionME.SetHLT(passedHLT[me]);               
    CutBasedHiggsErrorsSelectionME.SetIsChannel(m_channel[me]);   
    CutBasedHiggsErrorsSelectionME.SetElectronId(theReal);
    CutBasedHiggsErrorsSelectionME.SetPositronId(theFake);
    CutBasedHiggsErrorsSelectionME.SetElectronIsolation(theReal);
    CutBasedHiggsErrorsSelectionME.SetPositronIsolation(theFake);
    CutBasedHiggsErrorsSelectionME.SetElectronConvRejection(theReal);
    CutBasedHiggsErrorsSelectionME.SetPositronConvRejection(theFake);
    CutBasedHiggsErrorsSelectionME.SetElectronIp(theReal);
    CutBasedHiggsErrorsSelectionME.SetPositronIp(theFake);
    thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsErrorsSelectionME.SetElectronIp(-1);
      CutBasedHiggsErrorsSelectionME.SetPositronIp(-1);
    }
    CutBasedHiggsErrorsSelectionME.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsErrorsSelectionME.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsErrorsSelectionME.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons
    CutBasedHiggsErrorsSelectionME.SetNJets(njets[me]);
    CutBasedHiggsErrorsSelectionME.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsErrorsSelectionME.SetBTagJets(btag[me]);
    CutBasedHiggsErrorsSelectionME.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsErrorsSelectionME.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsErrorsSelectionME.SetMet(m_theMET);					
    CutBasedHiggsErrorsSelectionME.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsErrorsSelectionME.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsErrorsSelectionME.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsErrorsSelectionME.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsErrorsSelectionME.SetInvMass(m_mll[me]);
    CutBasedHiggsErrorsSelectionME.SetDetaLeptons(m_deltaEtaLeptons[me]);
    CutBasedHiggsSelectionME.SetWWInvMass(m_transvMass[me]);
    bool isSelectedErrorME = CutBasedHiggsErrorsSelectionME.output();



    // ----------------------------------------------------------------------------------------                                                                              
    // to compute statistical errors due to fake rate - FF                                        
    CutBasedHiggsErrorsSelectionME_FF.SetWeight(weightErrorFF);
    CutBasedHiggsErrorsSelectionME_FF.SetMcTruth(1);
    CutBasedHiggsErrorsSelectionME_FF.SetHLT(passedHLT[me]);               
    CutBasedHiggsErrorsSelectionME_FF.SetIsChannel(m_channel[me]);   
    CutBasedHiggsErrorsSelectionME_FF.SetElectronId(theReal);
    CutBasedHiggsErrorsSelectionME_FF.SetPositronId(theFake);
    CutBasedHiggsErrorsSelectionME_FF.SetElectronIsolation(theReal);
    CutBasedHiggsErrorsSelectionME_FF.SetPositronIsolation(theFake);
    CutBasedHiggsErrorsSelectionME_FF.SetElectronConvRejection(theReal);
    CutBasedHiggsErrorsSelectionME_FF.SetPositronConvRejection(theFake);
    CutBasedHiggsErrorsSelectionME_FF.SetElectronIp(theReal);
    CutBasedHiggsErrorsSelectionME_FF.SetPositronIp(theFake);
    thisMaxPtIpME = GetPt(pxMuon[theReal],pyMuon[theReal]);
    if (thisMaxPtIpME<20)   { 
      CutBasedHiggsErrorsSelectionME_FF.SetElectronIp(-1);
      CutBasedHiggsErrorsSelectionME_FF.SetPositronIp(-1);
    }
    CutBasedHiggsErrorsSelectionME_FF.SetHighElePt(hardestLeptonPt[me]); 
    CutBasedHiggsErrorsSelectionME_FF.SetLowElePt(slowestLeptonPt[me]);  
    CutBasedHiggsErrorsSelectionME_FF.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT only on electrons
    CutBasedHiggsErrorsSelectionME_FF.SetNJets(njets[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetNUncorrJets(nuncorrjets[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetBTagJets(btag[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetNSoftMuons(nsoftmu[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetNExtraLeptons(nextraleptons[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetMet(m_theMET);					
    CutBasedHiggsErrorsSelectionME_FF.SetProjectedMet(m_projectedMet[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetMetOverPtLL(m_metOptll[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetDeltaPhiLLJet(dphiLLJ[me]);   
    CutBasedHiggsErrorsSelectionME_FF.SetDeltaPhi(m_deltaPhi[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetInvMass(m_mll[me]);
    CutBasedHiggsErrorsSelectionME_FF.SetDetaLeptons(m_deltaEtaLeptons[me]);
    CutBasedHiggsSelectionME_FF.SetWWInvMass(m_transvMass[me]);
    bool isSelectedErrorME_FF = CutBasedHiggsErrorsSelectionME_FF.output();
  }
}

void LeptonPlusFakeMLSelection_fullME::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR PF SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionME.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR PF SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionME.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM STATISTICS FOR PF SELECTION ===" << std::endl;
  CutBasedHiggsSelectionStatME.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR FF SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionME_FF.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR FF SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionME_FF.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM STATISTICS FOR FF SELECTION ===" << std::endl;
  CutBasedHiggsSelectionStatME_FF.displayEfficiencies(datasetName);
  std::cout << "--------------------------------" << std::endl;


  // simple cuts based or like based ele id
  if (!_selectionME->getSwitch("asymmetricID")) {
    std::cout << "cut based symmetric ID: " << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  } else {
    std::cout << "cut based asymmetric ID: Low pT" << std::endl;
    EgammaCutBasedIDLow.displayEfficiencies();
    std::cout << "cut based asymmetric ID: High pT" << std::endl;
    EgammaCutBasedID.displayEfficiencies();
  }
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_acceptance() {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _acceptEleAll.clear();

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    if(_selectionME->getSwitch("etaElectronAcc") && !_selectionME->passCut("etaElectronAcc",etaEle[i]) ) continue;

    if(_selectionME->getSwitch("ptElectronAcc") && !_selectionME->passCut("ptElectronAcc",thisPt) ) continue;
    
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }

    _acceptEleAll.push_back(i);   
  }
  _acceptEleAll = sortElectronsByPt(_acceptEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_id( std::vector<int> acceptEle ) {

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
    if (!_selectionME->getSwitch("asymmetricID")) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionME->getSwitch("asymmetricID")) {
      if (thisPt>=20) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisPt<20)  isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }

    if (!theElectronID) continue;

    // further requests if we apply the smurf ID and pT<20
    TString stringIdLow (_selectionME->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      if ( thisPt<20  ) {
	if ( fbremEle[iEle]<0.15 && !(fabs(etaEle[iEle])<1.0 && eSuperClusterOverPEle[iEle]>0.95) ) continue;  // hardcoded
      }
    }

    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }

    _idEleAll.push_back(thisEle);  
  }
  _idEleAll = sortElectronsByPt(_idEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_isol( std::vector<int> idEle ) {

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

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_conv( std::vector<int> isolEle ) {

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


std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_ip( std::vector<int> convEle ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
  
  _ipEleAll.clear();

  for (int iEle=0; iEle<convEle.size(); iEle++) {
    int thisEle = convEle[iEle];

    int gsfTrack = gsfTrackIndexEle[thisEle]; 
    // float d3dEle = impactPar3DGsfTrack[gsfTrack];
    // if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",d3dEle)) ) continue;   
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionME->getSwitch("electronDz") && (!_selectionME->passCut("electronDz",dzEle)) )  continue;

    float thisPt     = GetPt(pxEle[thisEle],pyEle[thisEle]);
    float thisCharge = chargeEle[thisEle];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisEle; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisEle; }
    
    _ipEleAll.push_back(thisEle);  
  }
  _ipEleAll = sortElectronsByPt(_ipEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestElectronPair_denominator() {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _denomEleAll.clear();   

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    bool isGoodDenom = isEleDenomFake(i);
    if (!isGoodDenom) continue;

    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }

    _denomEleAll.push_back(i);
  }
  _denomEleAll = sortElectronsByPt(_denomEleAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonPair_acceptance() {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
  
  _acceptMuonsAll.clear();

  for(int i=0;i<nMuon;i++) {

    // only not standalone muons latinos
    Utils anaUtils;
    bool thisMuFlag1 = anaUtils.muonIdVal(muonIdMuon[i],AllGlobalMuons);
    bool thisMuFlag2 = anaUtils.muonIdVal(muonIdMuon[i],AllTrackerMuons);
    if (!thisMuFlag1 && !thisMuFlag2) continue;
    if(typeMuon[i]==8) continue;  // to be used when new trees are available
    // only not standalone muons latinos

    if(_selectionME->getSwitch("etaMuonAcc") && !_selectionME->passCut("etaMuonAcc",etaMuon[i]) ) continue;

    float thisPt = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionME->getSwitch("ptMuonAcc") && !_selectionME->passCut("ptMuonAcc",thisPt) ) continue;

    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
    
    _acceptMuonsAll.push_back(i);  
  }
  _acceptMuonsAll = sortMuonsByPt(_acceptMuonsAll);

  return make_pair(theLep1,theLep2);
}


std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonPair_id( std::vector<int> acceptMu ) {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
  
  _idMuonsAll.clear();

  for(int iMu=0; iMu<acceptMu.size(); iMu++) {

    int thisMu = acceptMu[iMu];

    bool theMuonID = true;
    isMuonID(thisMu, &theMuonID);
    if (!theMuonID) continue;

    float thisPt     = GetPt(pxMuon[thisMu],pyMuon[thisMu]);
    float thisCharge = chargeMuon[thisMu];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisMu; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisMu; }
    
    _idMuonsAll.push_back(thisMu);   
  }
  _idMuonsAll = sortMuonsByPt(_idMuonsAll);
  
  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonPair_isol( std::vector<int> idMu ) {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
   
  _isolMuonsAll.clear();
  
  for(int iMu=0; iMu<idMu.size(); iMu++) {
    
    int thisMu   = idMu[iMu];
    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);

    // fixme: diverso da prima: rimuovevo il secondo leptone....
    // float muonTrackerForGlobal = sumPt03Muon[thisMu];
    // float muonEcalForGlobal    = emEt03Muon[thisMu];
    // float muonHcalForGlobal    = hadEt03Muon[thisMu]; 
    // float theMuonGlobalSum     = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal - rhoFastjet*TMath::Pi()*0.3*0.3;
    // float theRelMuonIso        = theMuonGlobalSum/thisPt; 
    // if(_selectionME->getSwitch("muGlobalIso") && !_selectionME->passCut("muGlobalIso",theRelMuonIso)) continue;  
    if( ! isPFIsolatedMuon(thisMu) ) continue;

    float thisCharge = chargeMuon[thisMu];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisMu; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisMu; }

    _isolMuonsAll.push_back(thisMu);   
  }
  _isolMuonsAll = sortMuonsByPt(_isolMuonsAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonPair_ip( std::vector<int> isoMu ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _ipMuonsAll.clear();  

  for(int iMu=0; iMu<isoMu.size(); iMu++) {

    int thisMu = isoMu[iMu];

    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);
    
    int ctfMuon   = trackIndexMuon[thisMu]; 
    float dxyMuon = transvImpactParTrack[ctfMuon];
    float dzMuon  = PVzPV[m_closestPV] - trackVzTrack[ctfMuon];   
    // if (_selectionME->getSwitch("muonIP") && (!_selectionME->passCut("muonIP",dxyMuon)) ) continue;   
    // if (_selectionME->getSwitch("muonDz") && (!_selectionME->passCut("muonDz",dzMuon)) )  continue;   
    if (thisPt>20)    // hardcoded
      if (_selectionME->getSwitch("muonIPhighPT") && (!_selectionME->passCut("muonIPhighPT",dxyMuon)) ) continue;   
    
    if (thisPt<20)    // hardcoded
      if (_selectionME->getSwitch("muonIPlowPT")  && (!_selectionME->passCut("muonIPlowPT",dxyMuon)) ) continue;   
    
    if (_selectionME->getSwitch("muonDz") && (!_selectionME->passCut("muonDz",dzMuon)) )  continue;   

    float thisCharge = chargeMuon[thisMu];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisMu; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisMu; }
    
    _ipMuonsAll.push_back(thisMu);   
  }
  _ipMuonsAll = sortMuonsByPt(_ipMuonsAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonPair_denominator() {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
  
  _denomMuonAll.clear();   

  for(int i=0;i<nMuon;i++) {

    TVector3 pLepton(pxMuon[i],pyMuon[i],pzMuon[i]);
    float thisPt=pLepton.Pt();

    bool isGoodDenom = isMuonDenomFake(i);
    if (!isGoodDenom) continue;

    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }

    _denomMuonAll.push_back(i);
  }

  _denomMuonAll = sortMuonsByPt(_denomMuonAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> LeptonPlusFakeMLSelection_fullME::getBestMuonElePair(std::vector<int> electrons, std::vector<int> muons) {

  int theEle=-1;
  int theMuonA=-1;

  std::vector<int>::const_iterator muiter;
  for(muiter=muons.begin(); muiter<muons.end();++muiter) {
    int muCharge = chargeMuon[*muiter];
    float muPt = GetPt(pxMuon[*muiter],pyMuon[*muiter]);
    theMuonA = *muiter;
    std::vector<int>::const_iterator eleiter;
    for(eleiter=electrons.begin(); eleiter<electrons.end();++eleiter) {
      int eleCharge = chargeEle[*eleiter];
      float elePt = GetPt(pxEle[*eleiter],pyEle[*eleiter]);
      if(eleCharge*muCharge<0 && muPt > elePt) return std::make_pair(*muiter,*eleiter); 
    }
  }
  return std::make_pair(theMuonA,theEle); 
}

int LeptonPlusFakeMLSelection_fullME::getBestEleDenominator(int realMuon) {  // this requires not only denom ok, but also NOT tight!
  
  int theFake=-1;
  
  float maxPtFake=-1000.;
  
  float theRealPt = GetPt(pxMuon[realMuon],pyMuon[realMuon]);    

  for(int iele=0; iele<nEle; iele++) {
    
    float thisElePt = GetPt(pxEle[iele],pyEle[iele]);    

    if (chargeEle[iele]*chargeMuon[realMuon]>0) continue;

    if (thisElePt>theRealPt) continue;

    bool isGoodDenom = isEleDenomFake(iele);   
    if (!isGoodDenom) continue;
    
    // removing candidates passing the tight selection
    bool isTight = true;
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    if (!_selectionME->getSwitch("asymmetricID")) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionME->getSwitch("asymmetricID")) {
      if (thisElePt>=20) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisElePt<20)  isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }
    if (!theElectronID)      isTight = false;
    TString stringIdLow (_selectionME->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      if ( thisElePt<20  ) {
	if ( fbremEle[iele]<0.15 && !(fabs(etaEle[iele])<1.0 && eSuperClusterOverPEle[iele]>0.95) ) isTight = false;  // hardcoded
      }
    }
    if (!theElectronIsol)    isTight = false; 
    if (!theElectronConvRej) isTight = false;
    int gsfTrack = gsfTrackIndexEle[iele]; 
    // float d3dEle = impactPar3DGsfTrack[gsfTrack];
    // if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",d3dEle)) ) isTight = false;
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",dxyEle)) ) isTight = false;
    if (_selectionME->getSwitch("electronDz") && (!_selectionME->passCut("electronDz",dzEle)) )  isTight = false;
    if (isTight) continue;

    if( thisElePt > maxPtFake ) { maxPtFake = thisElePt; theFake=iele; }
  }
  
  return theFake;
}

bool LeptonPlusFakeMLSelection_fullME::isEleDenomFake(int theEle) {

  Utils anaUtils;
  bool isGoodDenom = true;
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance	         
  if( p3Ele.Eta() > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )  isGoodDenom = false;    
  
  // taking the supercluster     
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  if( ecalDriven) sc = superClusterIndexEle[theEle];
  if(!ecalDriven) sc = PFsuperClusterIndexEle[theEle];
  if ( sc < 0 ) isGoodDenom = false;
  
  // barrel or endcap 
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);
  
  // isolation 
  float ecalIsol = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
  float hcalIsol = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  if(ecalIsol>0.2) isGoodDenom = false;
  if(hcalIsol>0.2) isGoodDenom = false;

  // H/E 
  if ( isEleEB && hOverEEle[theEle]>0.15) isGoodDenom = false;
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;

  // sigmaIetaIeta 
  bool isBarrelSc;
  if ( fabs(etaSC[sc]) <  1.479 ) isBarrelSc = true;
  if ( fabs(etaSC[sc]) >= 1.479 ) isBarrelSc = false;
  if ( isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.013 ) isGoodDenom = false;
  if (!isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.034 ) isGoodDenom = false;

  // spikes
  float theE1 = eMaxSC[sc];
  float theE4SwissCross = e4SwissCrossSC[sc];
  float theSpikeSC = 1.0 - (theE4SwissCross/theE1);
  if (theSpikeSC>0.95) isGoodDenom = false;

  return isGoodDenom;
}

int LeptonPlusFakeMLSelection_fullME::getBestMuDenominator(int realEle) {  // this requires not only denom ok, but also NOT tight!
  
  int theFake=-1;
  
  float maxPtFake=-1000.;
  
  float theRealPt = GetPt(pxEle[realEle],pyEle[realEle]);    

  for(int imu=0; imu<nMuon; imu++) {
    
    float thisMuPt = GetPt(pxMuon[imu],pyMuon[imu]);    
    
    if (chargeEle[realEle]*chargeMuon[imu]>0) continue;

    if (theRealPt>thisMuPt) continue;

    bool isGoodDenom = isMuonDenomFake(imu);   
    if (!isGoodDenom) continue;

    // removing candidates passing the tight selection
    bool isTight = true;
    isMuonID(imu, &isTight);
    if (isTight) continue;

    if( thisMuPt > maxPtFake ) { maxPtFake = thisMuPt; theFake=imu; }
  }
  
  return theFake;
}

bool LeptonPlusFakeMLSelection_fullME::isMuonDenomFake(int theMuon) {

  bool isGoodDenom = true;
  TVector3 p3Muon(pxMuon[theMuon], pyMuon[theMuon], pzMuon[theMuon]);
  
  // acceptance	   
  if( fabs(p3Muon.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Muon.Pt() < 10. )        isGoodDenom = false;

  // muonID
  bool isTight = true;
  isMuonID(theMuon, &isTight);
  if (!isTight) isGoodDenom = false;
  
  // isolation
  float thePFMuonIso = pfCombinedIsoMuon[theMuon]/p3Muon.Pt();
  if ( thePFMuonIso > 0.4 ) isGoodDenom = false;
  
  // IP
  int ctfMuon   = trackIndexMuon[theMuon]; 
  float dxyMuon = transvImpactParTrack[ctfMuon];
  float dzMuon  = PVzPV[m_closestPV] - trackVzTrack[ctfMuon];  
  if (fabs(dxyMuon)>0.1 ) isGoodDenom = false;
  if (fabs(dzMuon)>0.1  ) isGoodDenom = false;
  
  return isGoodDenom;
}

float LeptonPlusFakeMLSelection_fullME::getElectronFakeRate( float fakePt, bool isFakeBarrel ) {
  
  for (int theBin = 0; theBin<7; theBin++) {
    
    if( fakePt >= m_eleMinFakePt[theBin] && fakePt < m_eleMaxFakePt[theBin] ) {
      if ( isFakeBarrel) return m_eleFakeRateEB[theBin];
      if (!isFakeBarrel) return m_eleFakeRateEE[theBin];
    }
  }
  
  std::cout << "BIG ERROR: fakePt = " << fakePt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getElectronFakeRateError( float fakePt, bool isFakeBarrel ) {
  
  for (int theBin = 0; theBin < 7; theBin++) {
    if( fakePt >= m_eleMinFakePt[theBin] && fakePt < m_eleMaxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_eleFakeRateEB_err[theBin];
      if (!isFakeBarrel) return m_eleFakeRateEE_err[theBin];
    }
  }
  
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getElectronPromptRate( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin<7; theBin++) {
    
    if( promptPt >= m_eleMinPromptPt[theBin] && promptPt < m_eleMaxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_elePromptRateEB[theBin];
      if (!isPromptBarrel) return m_elePromptRateEE[theBin];
    }
  }
  
  std::cout << "BIG ERROR: promptPt = " << promptPt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getElectronPromptRateError( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin < 7; theBin++) {
    if( promptPt >= m_eleMinPromptPt[theBin] && promptPt < m_eleMaxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_elePromptRateEB_err[theBin];
      if (!isPromptBarrel) return m_elePromptRateEE_err[theBin];
    }
  }

  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getMuonFakeRate( float fakePt, bool isFakeBarrel ) {
  
  for (int theBin = 0; theBin<7; theBin++) {
    
    if( fakePt >= m_muonMinFakePt[theBin] && fakePt < m_muonMaxFakePt[theBin] ) {
      if ( isFakeBarrel) return m_muonFakeRateEB[theBin];
      if (!isFakeBarrel) return m_muonFakeRateEE[theBin];
    }
  }
  
  std::cout << "BIG ERROR: fakePt = " << fakePt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getMuonFakeRateError( float fakePt, bool isFakeBarrel ) {
  
  for (int theBin = 0; theBin < 7; theBin++) {
    if( fakePt >= m_muonMinFakePt[theBin] && fakePt < m_muonMaxFakePt[theBin] ) {
      if (isFakeBarrel)  return m_muonFakeRateEB_err[theBin];
      if (!isFakeBarrel) return m_muonFakeRateEE_err[theBin];
    }
  }
  
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getMuonPromptRate( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin<7; theBin++) {
    
    if( promptPt >= m_muonMinPromptPt[theBin] && promptPt < m_muonMaxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_muonPromptRateEB[theBin];
      if (!isPromptBarrel) return m_muonPromptRateEE[theBin];
    }
  }
  
  std::cout << "BIG ERROR: promptPt = " << promptPt << endl;
  return -1.;
}

float LeptonPlusFakeMLSelection_fullME::getMuonPromptRateError( float promptPt, bool isPromptBarrel ) {

  for (int theBin = 0; theBin < 7; theBin++) {
    if( promptPt >= m_muonMinPromptPt[theBin] && promptPt < m_muonMaxPromptPt[theBin] ) {
      if (isPromptBarrel)  return m_muonPromptRateEB_err[theBin];
      if (!isPromptBarrel) return m_muonPromptRateEE_err[theBin];
    }
  }

  return -1.;
}

void LeptonPlusFakeMLSelection_fullME::setKinematicsME(int myReal, int myFake) {

  if (myFake > -1 && myReal > -1) {

    eleCands[me].push_back(myFake);
    muCands[me] .push_back(myReal);
    m_p4LeptonMinus[me] -> SetXYZT(pxEle[myFake],pyEle[myFake],pzEle[myFake],energyEle[myFake]);
    m_p4LeptonPlus[me]  -> SetXYZT(pxMuon[myReal],pyMuon[myReal],pzMuon[myReal],energyMuon[myReal]);
    hardestLeptonPt[me] = GetPt(pxMuon[myReal],pyMuon[myReal]);
    slowestLeptonPt[me] = GetPt(pxEle[myFake], pyEle[myFake]);

    if ( hardestLeptonPt[me] != GetPt(pxMuon[myReal],pyMuon[myReal]) ) {
      cout << "questo non puo succedere mai" << endl;
      cout << "myReal = " << myReal << ", myFake = " << myFake << endl;
      cout << "hardest PT = " << hardestLeptonPt[me] << endl;
      cout << "slowest PT = " << slowestLeptonPt[me] << endl;
      cout << "real PT = "    << GetPt(pxMuon[myReal],pyMuon[myReal]) << endl;
      cout << "fake PT = "    << GetPt(pxEle[myFake],pyEle[myFake]) << endl;
    }

    hardestLeptonEta[me] = etaMuon[myReal];
    slowestLeptonEta[me] = etaEle[myFake];
    m_mll[me]           = (*(m_p4LeptonMinus[me]) + *(m_p4LeptonPlus[me])).M();
    m_deltaPhi[me]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[me]->Vect().DeltaPhi(m_p4LeptonPlus[me]->Vect()));
    m_deltaErre[me]     = m_p4LeptonMinus[me]->Vect().DeltaR(m_p4LeptonPlus[me]->Vect());
    m_deltaEtaLeptons[me] = etaEle[myFake]-etaEle[myReal];
    m_dilepPt[me].SetXYZ( m_p4LeptonMinus[me]->Vect().X()+m_p4LeptonPlus[me]->Vect().X(),m_p4LeptonMinus[me]->Vect().Y()+m_p4LeptonPlus[me]->Vect().Y(),0.0 );
    // usual definition
    m_transvMass[me] = sqrt( 2.*(m_dilepPt[me].Pt())*(m_p3PFMET->Pt())*(1- cos(m_p3PFMET->DeltaPhi(m_dilepPt[me]))) );
    // chris' variable
    // m_transvMass[me]    = CalcGammaMRstar(*m_p4LeptonMinus[me],*m_p4LeptonPlus[me]);
    m_metOptll[me]      = m_theMET / m_dilepPt[me].Pt();
    m_mT2[me]           = 0.;
    m_projectedMet[me]  = GetProjectedMet(m_p4LeptonMinus[me]->Vect(),m_p4LeptonPlus[me]->Vect());
  }  
}

void LeptonPlusFakeMLSelection_fullME::resetKinematicsStart() {

  thePreElectronME = -1;
  thePreMuonME     = -1;
  theReal          = -1;
  theFake          = -1;
}

void LeptonPlusFakeMLSelection_fullME::resetKinematics() {

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

void LeptonPlusFakeMLSelection_fullME::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {
  
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
  thisCutBasedID->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  thisCutBasedID->SetNBrem( nbremsEle[eleIndex] );
  thisCutBasedID->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                
  thisCutBasedID->SetTrkIsolation ( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                        
  thisCutBasedID->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );         
  float iso = 0.0;
  if ( anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB) ) iso = dr03TkSumPtEle[eleIndex] + max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else iso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  thisCutBasedID->SetCombinedIsolation( (iso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt );
  thisCutBasedID->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pt );
  thisCutBasedID->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  thisCutBasedID->SetConvDist( fabs(convDistEle[eleIndex]) );
  thisCutBasedID->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  thisCutBasedID->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );

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

void LeptonPlusFakeMLSelection_fullME::isMuonID(int muonIndex, bool *muonIdOutput) {

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


int LeptonPlusFakeMLSelection_fullME::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  m_goodJets.clear();
  float ETMax=0.;

  theLeadingJet[theChannel]=-1;   

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionME->getSwitch("etaJetAcc") && !_selectionME->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(_selectionME->getSwitch("etJetAcc") && !_selectionME->passCut("etJetAcc", pt)) continue;

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
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
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
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
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


int LeptonPlusFakeMLSelection_fullME::numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  int num=0;

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    float uncorrEt = uncorrEnergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrEnergyAK5PFPUcorrJet[j]);
    TVector3 p3Jet = p4Jet.Vect();

    if(_selectionME->getSwitch("etaJetAcc")      && !_selectionME->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;    
    if(_selectionME->getSwitch("etUncorrJetAcc") && !_selectionME->passCut("etUncorrJetAcc", uncorrEt))   continue;

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
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;
    
    num++;
  }
  
  return num;
}

float LeptonPlusFakeMLSelection_fullME::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  float output=-999;
  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);

    if(_selectionME->getSwitch("etaJetAcc") && !_selectionME->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    // hardcoded
    float rawpt = uncorrEnergyAK5PFPUcorrJet[j] * fabs(sin(thetaAK5PFPUcorrJet[j]));
    if(rawpt < 7.0) continue;

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
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionME->getSwitch("jetConeWidth") && _selectionME->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    float tmp = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];     
    if(tmp > output) output = tmp;
    
  }

  return output;

}

float LeptonPlusFakeMLSelection_fullME::deltaPhiLLJet(int ichan) {   
  
  int myLeadingJet = theLeadingJet[ichan];

  if(myLeadingJet > -1) {
    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);    
    return fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));                           
  } else return -999.;
}

int LeptonPlusFakeMLSelection_fullME::numSoftMuons(std::vector<int> muonToRemove) {

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

    // float dxy = transvImpactParTrack[track];
    // if(dxy > 0.100) continue;   
    float dxyMuon= transvImpactParTrack[track];
    float dzMuon = fabs(PVzPV[0] - trackVzTrack[track]);   
    if(dxyMuon > 0.200) continue;     // hardcoded  
    if(dzMuon  > 0.100) continue;     // hardcoded  

    // float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    // float isoSumRel = isoSumAbs / pt;
    float isoSumRel = pfCombinedIsoMuon[i] / pt;
    if(pt>20 || isoSumRel<0.1) continue;
    
    num++;
  }
  return num;
}

int LeptonPlusFakeMLSelection_fullME::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    
    bool isSelEle=false;
    for(int eleSel=0; eleSel<(int)eleToRemove.size(); eleSel++) {
      if(i==eleToRemove[eleSel]) isSelEle=true;
    }
    if(isSelEle) continue;

    if(_selectionME->getSwitch("etaElectronAcc") && !_selectionME->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_selectionME->getSwitch("ptElectronAcc")  && !_selectionME->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])) ) continue;

    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;
    if (!_selectionME->getSwitch("asymmetricID")) 
      isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
    if (_selectionME->getSwitch("asymmetricID")) {
      float pt = GetPt(pxEle[i],pyEle[i]);	
      if(pt>=20) isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
      if(pt<20)  isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedIDLow);
    }
    if(!theId || !theIso || !theConvRej) continue;

    // further requests if we apply the smurf ID and pT<15
    TString stringIdLow (_selectionME->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      float pt = GetPt(pxEle[i],pyEle[i]);
      if ( pt<20  ) {
	if ( fbremEle[i]>0.15 || ((fabs(etaEle[i])<1.0 && eSuperClusterOverPEle[i]>0.95)) ) continue;
      }
    }

    int track = gsfTrackIndexEle[i];
    // float d3dEle = impactPar3DGsfTrack[track];
    // if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",d3dEle)) ) continue;    
    float dxyEle = transvImpactParGsfTrack[track];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[track];   
    if (_selectionME->getSwitch("electronIP") && (!_selectionME->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionME->getSwitch("electronDz") && (!_selectionME->passCut("electronDz",dzEle)) ) continue;

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
    if(_selectionME->getSwitch("etaMuonAcc") && !_selectionME->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_selectionME->getSwitch("ptMuonAcc") && !_selectionME->passCut("ptMuonAcc",ptMu) ) continue;

    bool theId = true;
    isMuonID(i,&theId);
    if(!theId) continue;
    //float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    //float isoSumRel = isoSumAbs / ptMu;
    //if(_selectionME->getSwitch("muGlobalIso") && !_selectionME->passCut("muGlobalIso",isoSumRel)) continue;
    if( ! isPFIsolatedMuon(i) ) continue; 

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = PVzPV[m_closestPV] - trackVzTrack[track];  
    //if(_selectionME->getSwitch("muonIP") && !_selectionME->passCut("muonIP",dxy)) continue;
    //if(_selectionME->getSwitch("muonDz") && !_selectionME->passCut("muonDz",dz))  continue;  
    if (ptMu>20)   // hardcoded
      if (_selectionME->getSwitch("muonIPhighPT") && (!_selectionME->passCut("muonIPhighPT",dxy)) ) continue;   
    if (ptMu<20)    // hardcoded
      if (_selectionME->getSwitch("muonIPlowPT")  && (!_selectionME->passCut("muonIPlowPT",dxy)) ) continue;   
    if (_selectionME->getSwitch("muonDz") && (!_selectionME->passCut("muonDz",dz)) )  continue;   

    numMu++;
  }
  
  return numEle + numMu;
}

float LeptonPlusFakeMLSelection_fullME::GetProjectedMet(TVector3 p1, TVector3 p2) {

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
bool LeptonPlusFakeMLSelection_fullME::reloadTriggerMask()
{
  std::vector<int> triggerMask;

  // load the triggers required for ME
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersME.begin();fIter!=requiredTriggersME.end();++fIter) {
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(*fIter) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_requiredTriggersME = triggerMask;

  // load the triggers NOT required for ME
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersME.begin();fIter!=notRequiredTriggersME.end();++fIter) {
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(*fIter) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_notRequiredTriggersME = triggerMask;
}

bool LeptonPlusFakeMLSelection_fullME::hasPassedHLT() {

  Utils anaUtils;
  bool required    = anaUtils.getTriggersOR(m_requiredTriggersME, firedTrg);
  bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersME, firedTrg);
  return (required && !notRequired);
}

void LeptonPlusFakeMLSelection_fullME::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggersME=reqTriggers;
}

void LeptonPlusFakeMLSelection_fullME::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  notRequiredTriggersME=reqTriggers;
}

bool LeptonPlusFakeMLSelection_fullME::isPFIsolatedMuon(int muonIndex) {
  float eta = etaMuon[muonIndex];
  float pt = GetPt(pxMuon[muonIndex],pyMuon[muonIndex]);
  float iso = pfCombinedIsoMuon[muonIndex]/pt;
  if( pt>=20. && fabs(eta)<1.479 ) return (iso < 0.13);
  if( pt>=20. && fabs(eta)>=1.479 ) return (iso < 0.09);
  if( pt<20. && fabs(eta)<1.479 ) return (iso < 0.06);
  if( pt<20. && fabs(eta)>=1.479 ) return (iso < 0.05);
  return true;
}
