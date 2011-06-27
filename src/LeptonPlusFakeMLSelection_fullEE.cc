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
  
  // selection efficiencies: counters
  std::string fileCutsEE     = higgsConfigDirMass + "2e2nuCuts.txt";
  std::string fileSwitchesEE = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionEE.Configure   (fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE FP");
  CutBasedHiggsSelectionEE_FF.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE FF");

  // selection efficiencies: stat errors due to number of events 
  CutBasedHiggsSelectionStatEE.Configure   (fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION STAT ERRORS EE FP");
  CutBasedHiggsSelectionStatEE_FF.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION STAT ERRORS EE FF");

  // selection efficiencies: stat errors due to fake rate
  CutBasedHiggsErrorsSelectionEE.Configure   (fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION FAKE ERRORS EE FP"); 
  CutBasedHiggsErrorsSelectionEE_FF.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION FAKE ERRORS EE FF"); 

  // taking the selection
  _selectionEE        = CutBasedHiggsSelectionEE.GetSelection();  
  _selectionEE_FF     = CutBasedHiggsSelectionEE_FF.GetSelection();  
  _selectionStatEE    = CutBasedHiggsSelectionStatEE.GetSelection();  
  _selectionStatEE_FF = CutBasedHiggsSelectionStatEE_FF.GetSelection();  
  _selectionErrEE     = CutBasedHiggsErrorsSelectionEE.GetSelection();  
  _selectionErrEE_FF  = CutBasedHiggsErrorsSelectionEE_FF.GetSelection();  

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
    std::string goodRunJsonFile = "config/json/goodCollisions2011.json";         // chiara
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics
  for(int theChannel=0; theChannel<1; theChannel++) {
    m_p4LeptonPlus[theChannel]  = new TLorentzVector(0.,0.,0.,0.);
    m_p4LeptonMinus[theChannel] = new TLorentzVector(0.,0.,0.,0.);
    m_p3PFMET = new TVector3(0.,0.,0.);
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
  
  delete _selectionEE;
  delete _selectionEE_FF;
  delete _selectionErrEE;
  delete _selectionErrEE_FF;
  delete _selectionStatEE;
  delete _selectionStatEE_FF;

  myOutTreeEE -> save();
}

// chiara
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
  m_fakeRateEB[4] = 0.0458273;ID->SetCombinedPFIsolation( (pfCombinedIsoEle[el

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

  /*
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
  */

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

  /*
  // fake in the barrel from Smurf selection, QCD
  m_fakeRateEB[0] = 0.00154351;
  m_fakeRateEB[1] = 0.00145805;   
  m_fakeRateEB[2] = 0.0111048;  
  m_fakeRateEB[3] = 0.0120757; 
  m_fakeRateEB[4] = 0.0163262;
  
  m_fakeRateEB_err[0] = 0.000558265;
  m_fakeRateEB_err[1] = 0.000357493;
  m_fakeRateEB_err[2] = 0.00101031; 
  m_fakeRateEB_err[3] = 0.00101347;
  m_fakeRateEB_err[4] = 0.00456665;
  
  // fake in the endcap from Smurf selection, QCD
  m_fakeRateEE[0] = 0.00149395;
  m_fakeRateEE[1] = 0.00230214;
  m_fakeRateEE[2] = 0.0127027;
  m_fakeRateEE[3] = 0.0136544;
  m_fakeRateEE[4] = 0.0163247;
  
  m_fakeRateEE_err[0] = 0.000509235; 
  m_fakeRateEE_err[1] = 0.000404248;
  m_fakeRateEE_err[2] = 0.000981551;
  m_fakeRateEE_err[3] = 0.00100148;
  m_fakeRateEE_err[4] = 0.00537643;
  */

  // fake in the barrel from Smurf selection, ET>30
  m_fakeRateEB[0] = 0.00832695;  // here we subtract W and Z with lumi corresponding to HTLT8
  m_fakeRateEB[1] = 0.00907589;  // here we subtract W and Z with lumi corresponding to HTLT17
  m_fakeRateEB[2] = 0.0140018;   // " " 
  m_fakeRateEB[3] = 0.0136831;   // " " 
  m_fakeRateEB[4] = 0.0140156;   // " " 

  m_fakeRateEB_err[0] = 0.00140527;
  m_fakeRateEB_err[1] = 0.00107052;
  m_fakeRateEB_err[2] = 0.00129764;
  m_fakeRateEB_err[3] = 0.00116528;
  m_fakeRateEB_err[4] = 0.00385747;

  // fake in the endcap from Smurf selection, ET>30
  m_fakeRateEE[0] = 0.00322436;  // here we subtract W and Z with lumi corresponding to HTLT8
  m_fakeRateEE[1] = 0.0050205;   // here we subtract W and Z with lumi corresponding to HTLT17
  m_fakeRateEE[2] = 0.0108465;   // " " 
  m_fakeRateEE[3] = 0.00940649;  // " " 
  m_fakeRateEE[4] = 0.0163363;   // " " 
  
  m_fakeRateEE_err[0] = 0.000972436;
  m_fakeRateEE_err[1] = 0.000892635;
  m_fakeRateEE_err[2] = 0.0011173;
  m_fakeRateEE_err[3] = 0.000971426;
  m_fakeRateEE_err[4] = 0.00557363;

  /*
  // fake in the barrel from Smurf selection, ET>15
  m_fakeRateEB[0] = 0.0139881;    // here we subtract W and Z with lumi corresponding to HTLT8
  m_fakeRateEB[1] = 0.0142939;    // here we subtract W and Z with lumi corresponding to HTLT17
  m_fakeRateEB[2] = 0.0229292;    // " " 
  m_fakeRateEB[3] = 0.0167388;
  m_fakeRateEB[4] = 0.014033;

  m_fakeRateEB_err[0] = 0.000822471;
  m_fakeRateEB_err[1] = 0.000839656;
  m_fakeRateEB_err[2] = 0.00120164;
  m_fakeRateEB_err[3] = 0.00111492;
  m_fakeRateEB_err[4] = 0.00376642;

  // fake in the endcap from Smurf selection, ET>15
  m_fakeRateEE[0] = 0.00536005;   // here we subtract W and Z with lumi corresponding to HTLT8
  m_fakeRateEE[1] = 0.00617585;   // here we subtract W and Z with lumi corresponding to HTLT17
  m_fakeRateEE[2] = 0.0153299;    // " " 
  m_fakeRateEE[3] = 0.011276;
  m_fakeRateEE[4] = 0.015083;
  
  m_fakeRateEE_err[0] = 0.000538545;
  m_fakeRateEE_err[1] = 0.000566034;
  m_fakeRateEE_err[2] = 0.00088218;
  m_fakeRateEE_err[3] = 0.000863319;
  m_fakeRateEE_err[4] = 0.00521618;
  */

  /*
  // fake in the barrel from Smurf selection, ET>50
  m_fakeRateEB[0] = 0.00196092;
  m_fakeRateEB[1] = 0.00217957;
  m_fakeRateEB[2] = 0.0112285;
  m_fakeRateEB[3] = 0.00758858;
  m_fakeRateEB[4] = 0.00936297; 

  m_fakeRateEB_err[0] = 0.00198052;
  m_fakeRateEB_err[1] = 0.00130578;
  m_fakeRateEB_err[2] = 0.00239751;
  m_fakeRateEB_err[3] = 0.00130281;
  m_fakeRateEB_err[4] = 0.00335611;

  // fake in the endcap from Smurf selection, ET>50
  m_fakeRateEE[0] = 0.00258585;
  m_fakeRateEE[1] = 0.00439881;
  m_fakeRateEE[2] = 0.00564385;
  m_fakeRateEE[3] = 0.00496651;
  m_fakeRateEE[4] = 0.0193293;
  
  m_fakeRateEE_err[0] = 0.00258838;
  m_fakeRateEE_err[1] = 0.00221489;
  m_fakeRateEE_err[2] = 0.00188719;
  m_fakeRateEE_err[3] = 0.00123241;
  m_fakeRateEE_err[4] = 0.00689591;
  */
}

// numbers provided by Emanuele, ~150/pb
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
  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";    
  myOutTreeEE = new RedHiggsTree(reducedTreeNameEE.c_str());

  if ( _selectionEE->getSwitch("isData")) myOutTreeEE->addRunInfos();
  myOutTreeEE->addMLVars();
  myOutTreeEE->addLatinos();

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
    if ( !_selectionEE->getSwitch("isData") ) tmpWeight *= fPUWeight->GetWeight(nPU);  

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
    reloadTriggerMask(runNumber);
    bool passedHLT[1];
    passedHLT[ee] = hasPassedHLT();


    // -------------------------------------------------------------
    // vertex selection - we only consider the first vertex of the list ( = highest sumPT^2)
    bool isGoodVertex = true;
    if (nPV<1) isGoodVertex = false;
    float rhoVtx = sqrt(PVxPV[0]*PVxPV[0] + PVyPV[0]*PVyPV[0]);
    if ( isFakePV[0] )       isGoodVertex = false;
    if ( ndofPV[0]<4 )       isGoodVertex = false;
    if ( fabs(PVzPV[0])>24.) isGoodVertex = false;
    if ( rhoVtx>2 )          isGoodVertex = false; 
    

    
    // ----------------------------------------------------------------------------
    // get the best electrons and best muons at acceptnace level 
    // ==> tu be used to select ALL the possible channels at the beginning only
    std::pair<int,int> thePreElectrons = getBestElectronPair_acceptance();
    thePreElectron  = thePreElectrons.second;
    thePrePositron  = thePreElectrons.first;

    // reconstructed channel
    m_channel[ee] = false;     

    // at this level the SELECTED channel should have pT > 10 and > 20. So far, at least 2 leptons with pT >20 and 10 in the event
    if (thePreElectron > -1 && thePrePositron > -1) {
      float thisMaxPt = TMath::Max(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      float thisMinPt = TMath::Min(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      if (isGoodVertex && thisMaxPt>20 && thisMinPt>10) m_channel[ee] = true;    // fixme: hardcoded
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
    
    // to be sure: I assumed above that tight => loose. If not true, I discard the event....
    if ( (thePositron>-1 && !isDenomFake(thePositron)) || (theElectron>-1 && !isDenomFake(theElectron)) ) {
      theReal = -1;
      theFake = -1;
    } 

    // sanity check
    if ( (is0tight && is1tight) || (is0tight && is2tight) || (is1tight && is2tight) ) cout << "questo non puo' succedere mai" << endl;

    // set of kinematics: : now I've all the final leptons 
    resetKinematics();

    // MET is an event variable. Independent o the channel        
    m_p3PFMET->SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);             // the one associated to the 0th vertex      
    m_theMET = m_p3PFMET->Pt();

    setKinematicsEE(theReal, theFake);
        

    // weight with the Fake / Prompt -> L2 probability	
    float theFakePt = GetPt(pxEle[theFake],pyEle[theFake]);
    float theRealPt = GetPt(pxEle[theReal],pyEle[theReal]);
    bool  isFakeBarrel = false;
    bool  isRealBarrel = false;
    if ( fabs(etaEle[theFake])<1.476 ) isFakeBarrel = true;
    if ( fabs(etaEle[theReal])<1.476 ) isRealBarrel = true;

    // do both F-P and F-F analysis
    float weightFP      = 1.;
    float weightErrorFP = 1.;
    float weightFF      = 1.;
    float weightErrorFF = 1.;

    if ( theFake>-1 && theReal>-1) {
      float fakerate1      = getFakeRate( theFakePt, isFakeBarrel );
      float fakerateErr1   = getFakeRateError( theFakePt, isFakeBarrel );
      float promptrate1    = getPromptRate( theFakePt, isFakeBarrel );
      float promptrateErr1 = getPromptRateError( theFakePt, isFakeBarrel );
      // 
      float fakerate2      = getFakeRate( theRealPt, isRealBarrel );
      float fakerateErr2   = getFakeRateError( theRealPt, isRealBarrel );
      float promptrate2    = getPromptRate( theRealPt, isRealBarrel );
      float promptrateErr2 = getPromptRateError( theRealPt, isRealBarrel );
      
      float thisPartWeightFP    = 1.;
      float thisPartWeightErrFP = 1.;
      float thisPartWeightFF    = 1.;
      float thisPartWeightErrFF = 1.;

      if (is0tight) {
	float myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));

	// for F-P
	thisPartWeightFP = -( fakerate1*promptrate2*promptrate1*fakerate2 + fakerate2*promptrate1*promptrate2*fakerate1 )*myFactor;
	// float d1F = (4*fakerate*pow(promptrate,3)) / ( pow((fakerate-promptrate),3) );
	// float d1P = (4*promptrate*pow(fakerate,3)) / ( pow((fakerate-promptrate),3) );
	// thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);

	// for F-F
	thisPartWeightFF = ( promptrate1*promptrate2*fakerate1*fakerate2 )*myFactor;     
	// float d2F = -(2*fakerate*pow(promptrate,3)) / ( pow((fakerate-promptrate),3) );
	// float d2P =  (2*promptrate*pow(fakerate,3)) / ( pow((fakerate-promptrate),3) );
	// thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
      }

      if (is1tight) {
	double myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));

	// for F-P
	thisPartWeightFP = ( fakerate1*(1-promptrate2)*promptrate1*fakerate2 + promptrate1*(1-fakerate2)*fakerate1*promptrate2 )*myFactor;  
	// float d1F = ((promptrate*promptrate)*(-promptrate + fakerate*(-3+4*promptrate))) / ( pow((fakerate-promptrate),3) );
	// float d1P = ((fakerate*fakerate)*(fakerate + promptrate*(3-4*fakerate))) / ( pow((fakerate-promptrate),3) );
	// thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);

	// for F-F 
	thisPartWeightFF = -( promptrate1*(1-promptrate2)*fakerate1*fakerate2 )*myFactor;    
	// float d2F = (2*fakerate*(promptrate-1)*pow(promptrate,2))/( pow((fakerate-promptrate),3) );
	// float d2P = (fakerate*fakerate*(fakerate+promptrate-2*fakerate*promptrate)) / ( pow((fakerate-promptrate),3) );
	// thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
      }
      
      if (is2tight) {
	double myFactor = 1/((promptrate1-fakerate1)*(promptrate2-fakerate2));

	// for F-P
	thisPartWeightFP = -( (promptrate2*fakerate1*(1-promptrate1)*(1-fakerate2)) + (promptrate1*fakerate2*(1-promptrate2)*(1-fakerate1)) )*myFactor; 
	// float d1P = -( 2*(fakerate-1)*fakerate*(-promptrate + fakerate*(-1 + 2*promptrate))) / ( pow((fakerate-promptrate),3) );
	// float d1F =  ( 2*(promptrate-1)*promptrate*(-promptrate + fakerate*(-1 + 2*promptrate))) / ( pow((fakerate-promptrate),3) );
	// thisPartWeightErrFP = sqrt(d1F*fakerateErr*d1F*fakerateErr + d1P*promptrateErr*d1P*promptrateErr);
	
	// for F-F
	thisPartWeightFF = ( fakerate1*fakerate2*(1-promptrate1)*(1-promptrate2) )*myFactor;
	// float d2F = -(2*fakerate*(promptrate-1)*(promptrate-1)*promptrate) / ( pow((fakerate-promptrate),3) );
	// float d2P =  (2*(fakerate-1)*pow(fakerate,2)*(promptrate-1)) / ( pow((fakerate-promptrate),3) ); 
	// thisPartWeightErrFF = sqrt(d2F*fakerateErr*d2F*fakerateErr + d2P*promptrateErr*d2P*promptrateErr);
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


    // -------------------------------------------------------
    // filling counter for FP ( = W+jets)
    CutBasedHiggsSelectionEE.SetWeight(weightFP);                  
    CutBasedHiggsSelectionEE.SetMcTruth(true);  
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
    // checking if the highest pT electron has pT>20
    float thisMaxPtIpEE = TMath::Max(GetPt(pxEle[theReal],pyEle[theReal]),GetPt(pxEle[theFake],pyEle[theFake]));
    if (thisMaxPtIpEE<20)   { 
      CutBasedHiggsSelectionEE.SetElectronIp(-1);
      CutBasedHiggsSelectionEE.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionEE.SetHighElePt(hardestLeptonPt[ee]); 
    CutBasedHiggsSelectionEE.SetLowElePt(slowestLeptonPt[ee]);  
    CutBasedHiggsSelectionEE.SetExtraSlowLeptonPTCut(10.0); // enforce the min pT cut only on electrons 

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
    // CutBasedHiggsSelectionEE.SetWWInvMass(2.*m_transvMass[ee]/_massVal);
    CutBasedHiggsSelectionEE.SetWWInvMass(m_transvMass[ee]);

    bool isSelectedEE           = CutBasedHiggsSelectionEE.output();    
    bool selUpToFinalLeptonsEE  = CutBasedHiggsSelectionEE.outputUpToFinalLeptons();
    bool selUpToJetVetoEE       = CutBasedHiggsSelectionEE.outputUpToJetVeto();
    bool selUpToUncorrJetVetoEE = CutBasedHiggsSelectionEE.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiEE       = CutBasedHiggsSelectionEE.outputPreDeltaPhi();

    bool outputStep0  = CutBasedHiggsSelectionEE.outputStep0();
    bool outputStep1  = CutBasedHiggsSelectionEE.outputStep1();
    bool outputStep2  = CutBasedHiggsSelectionEE.outputStep2();
    bool outputStep3  = CutBasedHiggsSelectionEE.outputStep3();
    bool outputStep4  = CutBasedHiggsSelectionEE.outputStep4();
    bool outputStep5  = CutBasedHiggsSelectionEE.outputStep5();
    bool outputStep6  = CutBasedHiggsSelectionEE.outputStep6();
    bool outputStep7  = CutBasedHiggsSelectionEE.outputStep7();
    bool outputStep8  = CutBasedHiggsSelectionEE.outputStep8();
    bool outputStep9  = CutBasedHiggsSelectionEE.outputStep9();
    bool outputStep10 = CutBasedHiggsSelectionEE.outputStep10();
    bool outputStep11 = CutBasedHiggsSelectionEE.outputStep11();
    bool outputStep12 = CutBasedHiggsSelectionEE.outputStep12();
    bool outputStep13 = CutBasedHiggsSelectionEE.outputStep13();
    bool outputStep14 = CutBasedHiggsSelectionEE.outputStep14();
    bool outputStep15 = CutBasedHiggsSelectionEE.outputStep15();
    bool outputStep16 = CutBasedHiggsSelectionEE.outputStep16();
    bool outputStep17 = CutBasedHiggsSelectionEE.outputStep17();
    bool outputStep18 = CutBasedHiggsSelectionEE.outputStep18();
    bool outputStep19 = CutBasedHiggsSelectionEE.outputStep19();
    bool outputStep20 = CutBasedHiggsSelectionEE.outputStep20();
    bool outputStep21 = CutBasedHiggsSelectionEE.outputStep21();
    bool outputStep22 = CutBasedHiggsSelectionEE.outputStep22();
    bool outputStep23 = CutBasedHiggsSelectionEE.outputStep23();
    bool outputStep24 = CutBasedHiggsSelectionEE.outputStep24();

    myOutTreeEE->fillRunInfos(runNumber, lumiBlock, eventNumber, weightFP); 
    
    myOutTreeEE -> fillAll(m_chMet[ee], GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]),
			   m_projectedMet[ee], m_deltaPhi[ee], m_deltaErre[ee], m_transvMass[ee], m_mll[ee], 
			   hardestLeptonPt[ee], slowestLeptonPt[ee], hardestLeptonEta[ee], slowestLeptonEta[ee], 
			   m_deltaEtaLeptons[ee], nPV,
			   selUpToFinalLeptonsEE, selUpToJetVetoEE, selUpToUncorrJetVetoEE, selPreDeltaPhiEE, isSelectedEE);
    
    myOutTreeEE -> fillMLVars(njets[ee], nuncorrjets[ee], m_maxDxyEvt, m_maxDszEvt, m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);
    
    myOutTreeEE -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 );
    
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeEE -> store();





    // ----------------------------------------------------------------------------------------
    // filling counters for FF ( = QCD )
    CutBasedHiggsSelectionEE_FF.SetWeight(weightFF);            
    CutBasedHiggsSelectionEE_FF.SetMcTruth(true);
    CutBasedHiggsSelectionEE_FF.SetHLT(passedHLT[ee]);
    CutBasedHiggsSelectionEE_FF.SetIsChannel(m_channel[ee]);
    CutBasedHiggsSelectionEE_FF.SetElectronId(1);
    CutBasedHiggsSelectionEE_FF.SetPositronId(1);
    CutBasedHiggsSelectionEE_FF.SetElectronIsolation(1);
    CutBasedHiggsSelectionEE_FF.SetPositronIsolation(1);
    CutBasedHiggsSelectionEE_FF.SetElectronConvRejection(1);
    CutBasedHiggsSelectionEE_FF.SetPositronConvRejection(1);
    CutBasedHiggsSelectionEE_FF.SetElectronIp(theReal);
    CutBasedHiggsSelectionEE_FF.SetPositronIp(theFake);
    if (thisMaxPtIpEE<20) {
      CutBasedHiggsSelectionEE_FF.SetElectronIp(-1);
      CutBasedHiggsSelectionEE_FF.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionEE_FF.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsSelectionEE_FF.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsSelectionEE_FF.SetExtraSlowLeptonPTCut(10.0);  // enforce the min pT cut only on electrons 
    CutBasedHiggsSelectionEE_FF.SetNJets(njets[ee]);
    CutBasedHiggsSelectionEE_FF.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsSelectionEE_FF.SetBTagJets(btag[ee]);
    CutBasedHiggsSelectionEE_FF.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsSelectionEE_FF.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsSelectionEE_FF.SetMet(m_theMET);
    CutBasedHiggsSelectionEE_FF.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsSelectionEE_FF.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsSelectionEE_FF.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsSelectionEE_FF.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsSelectionEE_FF.SetInvMass(m_mll[ee]);
    CutBasedHiggsSelectionEE_FF.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsSelectionEE_FF.SetWWInvMass(m_transvMass[ee]);
    bool isSelectedEE_FF = CutBasedHiggsSelectionEE_FF.output();    


    // ----------------------------------------------------------------------------------------
    // filling counters for statistical errors due to data / mc - FP
    float forStatErrFP = weightFP*weightFP;                    
    CutBasedHiggsSelectionStatEE.SetWeight(forStatErrFP);      
    CutBasedHiggsSelectionStatEE.SetMcTruth(true);
    CutBasedHiggsSelectionStatEE.SetHLT(passedHLT[ee]);
    CutBasedHiggsSelectionStatEE.SetIsChannel(m_channel[ee]);
    CutBasedHiggsSelectionStatEE.SetElectronId(1);
    CutBasedHiggsSelectionStatEE.SetPositronId(1);
    CutBasedHiggsSelectionStatEE.SetElectronIsolation(1);
    CutBasedHiggsSelectionStatEE.SetPositronIsolation(1);
    CutBasedHiggsSelectionStatEE.SetElectronConvRejection(1);
    CutBasedHiggsSelectionStatEE.SetPositronConvRejection(1);
    CutBasedHiggsSelectionStatEE.SetElectronIp(theReal);
    CutBasedHiggsSelectionStatEE.SetPositronIp(theFake);
    if (thisMaxPtIpEE<20) {
      CutBasedHiggsSelectionStatEE.SetElectronIp(-1);
      CutBasedHiggsSelectionStatEE.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionStatEE.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsSelectionStatEE.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsSelectionStatEE.SetExtraSlowLeptonPTCut(10.0);  // enforce the min pT cut only on electrons 
    CutBasedHiggsSelectionStatEE.SetNJets(njets[ee]);
    CutBasedHiggsSelectionStatEE.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsSelectionStatEE.SetBTagJets(btag[ee]);
    CutBasedHiggsSelectionStatEE.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsSelectionStatEE.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsSelectionStatEE.SetMet(m_theMET);
    CutBasedHiggsSelectionStatEE.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsSelectionStatEE.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsSelectionStatEE.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsSelectionStatEE.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsSelectionStatEE.SetInvMass(m_mll[ee]);
    CutBasedHiggsSelectionStatEE.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsSelectionStatEE.SetWWInvMass(m_transvMass[ee]);
    bool isSelectedStatEE = CutBasedHiggsSelectionStatEE.output();    


    // ----------------------------------------------------------------------------------------
    // filling counters for statistical errors due to data / mc - FF
    float forStatErrFF = weightFF*weightFF;                    
    CutBasedHiggsSelectionStatEE_FF.SetWeight(forStatErrFF);      
    CutBasedHiggsSelectionStatEE_FF.SetMcTruth(true);
    CutBasedHiggsSelectionStatEE_FF.SetHLT(passedHLT[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetIsChannel(m_channel[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetElectronId(1);
    CutBasedHiggsSelectionStatEE_FF.SetPositronId(1);
    CutBasedHiggsSelectionStatEE_FF.SetElectronIsolation(1);
    CutBasedHiggsSelectionStatEE_FF.SetPositronIsolation(1);
    CutBasedHiggsSelectionStatEE_FF.SetElectronConvRejection(1);
    CutBasedHiggsSelectionStatEE_FF.SetPositronConvRejection(1);
    CutBasedHiggsSelectionStatEE_FF.SetElectronIp(theReal);
    CutBasedHiggsSelectionStatEE_FF.SetPositronIp(theFake);
    if (thisMaxPtIpEE<20) {
      CutBasedHiggsSelectionStatEE_FF.SetElectronIp(-1);
      CutBasedHiggsSelectionStatEE_FF.SetPositronIp(-1);
    }
    CutBasedHiggsSelectionStatEE_FF.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetExtraSlowLeptonPTCut(10.0);  // enforce the min pT cut only on electrons 
    CutBasedHiggsSelectionStatEE_FF.SetNJets(njets[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetBTagJets(btag[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetMet(m_theMET);
    CutBasedHiggsSelectionStatEE_FF.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetInvMass(m_mll[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsSelectionStatEE_FF.SetWWInvMass(m_transvMass[ee]);
    bool isSelectedStatEE_FF = CutBasedHiggsSelectionStatEE_FF.output();    


    // ----------------------------------------------------------------------------------------
    // to compute statistical errors due to fake rate - FP                                  
    CutBasedHiggsErrorsSelectionEE.SetWeight(weightErrorFP);       
    CutBasedHiggsErrorsSelectionEE.SetMcTruth(true);
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
    CutBasedHiggsErrorsSelectionEE.SetExtraSlowLeptonPTCut(10.0);  // enforce the min pT cut only on electrons 
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
    CutBasedHiggsErrorsSelectionEE.SetWWInvMass(m_transvMass[ee]);
    bool isSelectedErrorEE = CutBasedHiggsErrorsSelectionEE.output();    


    // ----------------------------------------------------------------------------------------
    // to compute statistical errors due to fake rate - FF                                  
    CutBasedHiggsErrorsSelectionEE_FF.SetWeight(weightErrorFF);       
    CutBasedHiggsErrorsSelectionEE_FF.SetMcTruth(true);
    CutBasedHiggsErrorsSelectionEE_FF.SetHLT(passedHLT[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetIsChannel(m_channel[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetElectronId(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetPositronId(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetElectronIsolation(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetPositronIsolation(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetElectronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetPositronConvRejection(1);
    CutBasedHiggsErrorsSelectionEE_FF.SetElectronIp(theReal);
    CutBasedHiggsErrorsSelectionEE_FF.SetPositronIp(theFake);
    if (thisMaxPtIpEE<20) {
      CutBasedHiggsErrorsSelectionEE_FF.SetElectronIp(-1);
      CutBasedHiggsErrorsSelectionEE_FF.SetPositronIp(-1);
    }
    CutBasedHiggsErrorsSelectionEE_FF.SetHighElePt(hardestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetLowElePt(slowestLeptonPt[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetExtraSlowLeptonPTCut(10.0);  // enforce the min pT cut only on electrons 
    CutBasedHiggsErrorsSelectionEE_FF.SetNJets(njets[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetNUncorrJets(nuncorrjets[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetBTagJets(btag[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetNSoftMuons(nsoftmu[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetMet(m_theMET);
    CutBasedHiggsErrorsSelectionEE_FF.SetProjectedMet(m_projectedMet[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetMetOverPtLL(m_metOptll[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetDeltaPhiLLJet(dphiLLJ[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetDeltaPhi(m_deltaPhi[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetInvMass(m_mll[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetDetaLeptons(m_deltaEtaLeptons[ee]);
    CutBasedHiggsErrorsSelectionEE_FF.SetWWInvMass(m_transvMass[ee]);
    bool isSelectedErrorEE_FF = CutBasedHiggsErrorsSelectionEE_FF.output();    
  }
}

void LeptonPlusFakeMLSelection_fullEE::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR PF SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionEE.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR PF SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionEE.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM STATISTICS FOR PF SELECTION ===" << std::endl;
  CutBasedHiggsSelectionStatEE.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "=== RATE ESTIMATED FROM FAKE RATE FOR FF SELECTION ===: " << std::endl;
  CutBasedHiggsSelectionEE_FF.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM FAKE RATE FOR FF SELECTION ===" << std::endl;
  CutBasedHiggsErrorsSelectionEE_FF.displayEfficiencies(datasetName);

  std::cout << "=== RATE UNCERTAINTY ESTIMATED FROM STATISTICS FOR FF SELECTION ===" << std::endl;
  CutBasedHiggsSelectionStatEE_FF.displayEfficiencies(datasetName);
  std::cout << "--------------------------------" << std::endl;


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

std::pair<int,int> LeptonPlusFakeMLSelection_fullEE::getBestElectronPair_acceptance() {
  
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
    if (!_selectionEE->getSwitch("asymmetricID")) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionEE->getSwitch("asymmetricID")) {
      if (thisPt>=20) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisPt<20)  isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }

    if (!theElectronID) continue;

    // further requests if we apply the smurf ID and pT<20
    TString stringIdLow (_selectionEE->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      if ( thisPt<20  ) {
	if ( fbremEle[thisEle]<0.15 && !(fabs(etaEle[thisEle])<1.0 && eSuperClusterOverPEle[thisEle]>0.95) ) continue; // hardcoded
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
    // float d3dEle = impactPar3DGsfTrack[gsfTrack];
    // if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",d3dEle)) ) continue;   
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionEE->getSwitch("electronDz") && (!_selectionEE->passCut("electronDz",dzEle)) )  continue;

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
    if (!_selectionEE->getSwitch("asymmetricID")) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if ( _selectionEE->getSwitch("asymmetricID")) {
      if (thisElePt>=20) isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisElePt<20)  isEleID(iele,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
    }
    if (!theElectronID) isTight = false;
    TString stringIdLow (_selectionEE->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      if ( thisElePt<20  ) {
	if ( fbremEle[iele]<0.15 && !(fabs(etaEle[iele])<1.0 && eSuperClusterOverPEle[iele]>0.95) ) isTight = false;  // hardcoded
      }
    }
    if (!theElectronIsol)    isTight = false; 
    if (!theElectronConvRej) isTight = false;
    int gsfTrack = gsfTrackIndexEle[iele]; 
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",dxyEle)) ) isTight = false;
    if (_selectionEE->getSwitch("electronDz") && (!_selectionEE->passCut("electronDz",dzEle)) )  isTight = false;
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
    m_dilepPt[ee].SetXYZ(m_p4LeptonMinus[ee]->Vect().X()+m_p4LeptonPlus[ee]->Vect().X(),m_p4LeptonMinus[ee]->Vect().Y()+m_p4LeptonPlus[ee]->Vect().Y(),0.0);
    // chris' definition
    // m_transvMass[ee]=CalcGammaMRstar(*m_p4LeptonMinus[ee],*m_p4LeptonPlus[ee]);
    // usual definition
    m_transvMass[ee] = sqrt( 2.*(m_dilepPt[ee].Pt())*(m_p3PFMET->Pt())*(1- cos(m_p3PFMET->DeltaPhi(m_dilepPt[ee]))) );
    m_metOptll[ee] = m_theMET / m_dilepPt[ee].Pt();
    m_mT2[ee] = 0.;
    m_projectedMet[ee] = GetProjectedMet(m_p4LeptonMinus[ee]->Vect(),m_p4LeptonPlus[ee]->Vect());
    m_chMet[ee] = pfChargedMet(m_p4LeptonMinus[ee]->Vect(),m_p4LeptonPlus[ee]->Vect()).Pt();
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
    m_chMet[theChannel]           = 0.;
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

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    // if(_selectionEE->getSwitch("etJetAcc") && !_selectionEE->passCut("etJetAcc", pt)) continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    
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

    if(pt>ETMax) {
      theLeadingJet[theChannel] = j;
      ETMax = pt;
    }
    
    if(_selectionEE->getSwitch("etJetAcc") && !_selectionEE->passCut("etJetAcc", pt)) continue;

    m_goodJets.push_back(j);
    num++;

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

    if(_selectionEE->getSwitch("etaJetAcc")      && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;    
    if(_selectionEE->getSwitch("etUncorrJetAcc") && !_selectionEE->passCut("etUncorrJetAcc", uncorrEt))   continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    
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

float LeptonPlusFakeMLSelection_fullEE::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  TString JESUncertainty(_selectionEE->getStringParameter("JESUncertainty"));

  float output=-999;
  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    // no threshold is applied here on pt. Not affected by JES uncertainties
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) pt = (GetJESCorrected(p4Jet,JESUncertainty.Data())).Pt();

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    // hardcoded
    float rawpt = uncorrEnergyAK5PFPUcorrJet[j] * fabs(sin(thetaAK5PFPUcorrJet[j]));
    if(rawpt < 7.0) continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    
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

int LeptonPlusFakeMLSelection_fullEE::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

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

    // further requests if we apply the smurf ID and pT<15
    TString stringIdLow (_selectionEE->getStringParameter("electronIDTypeLow"));
    if( stringIdLow.Contains("Smurf") ) {
      float pt = GetPt(pxEle[i],pyEle[i]);
      if ( pt<20  ) {
	if ( fbremEle[i]>0.15 || ((fabs(etaEle[i])<1.0 && eSuperClusterOverPEle[i]>0.95)) ) continue;
      }
    }

    int track = gsfTrackIndexEle[i];
    // float d3dEle = impactPar3DGsfTrack[track];
    // if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",d3dEle)) ) continue;    
    float dxyEle = transvImpactParGsfTrack[track];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[track];   
    if (_selectionEE->getSwitch("electronIP") && (!_selectionEE->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionEE->getSwitch("electronDz") && (!_selectionEE->passCut("electronDz",dzEle)) ) continue;

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
    // float isoSumAbs = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i] - rhoFastjet*TMath::Pi()*0.3*0.3;
    // float isoSumRel = isoSumAbs / ptMu;
    // if(_selectionEE->getSwitch("muGlobalIso") && !_selectionEE->passCut("muGlobalIso",isoSumRel)) continue;
    if( ! isPFIsolatedMuon(i) ) continue; 

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = PVzPV[m_closestPV] - trackVzTrack[track];  
    // if(_selectionEE->getSwitch("muonIP") && !_selectionEE->passCut("muonIP",dxy)) continue;
    // if(_selectionEE->getSwitch("muonDz") && !_selectionEE->passCut("muonDz",dz))  continue;  
    if (ptMu>20)    // hardcoded
      if (_selectionEE->getSwitch("muonIPhighPT") && (!_selectionEE->passCut("muonIPhighPT",dxy)) ) continue;   
    
    if (ptMu<20)    // hardcoded
      if (_selectionEE->getSwitch("muonIPlowPT")  && (!_selectionEE->passCut("muonIPlowPT",dxy)) ) continue;   
    
    if (_selectionEE->getSwitch("muonDz") && (!_selectionEE->passCut("muonDz",dz)) )  continue;   

    numMu++;
  }
  
  return numEle + numMu;
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
  TVector3 p3tkMet = pfChargedMet(p1,p2);

  float projMET_tk = 0.0;
  float deltaPhi1_tk = fabs(p1.DeltaPhi(p3tkMet));
  float deltaPhi2_tk = fabs(p2.DeltaPhi(p3tkMet));

  float deltaphi_tk = TMath::Min(deltaPhi1_tk,deltaPhi2_tk);

  if(deltaphi_tk<TMath::Pi()/2.) projMET_tk = p3tkMet.Mag() * sin(deltaphi_tk);
  else projMET_tk = p3tkMet.Mag();

  return TMath::Min(projMET_pf,projMET_tk);
}

/// specific for HWW that has multiple channels with different HLT requirements
bool LeptonPlusFakeMLSelection_fullEE::reloadTriggerMask(int runN) {

  std::vector<int> triggerMask;
  // load the triggers required for EE
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersEE.begin();fIter!=requiredTriggersEE.end();++fIter) {
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos)
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

bool LeptonPlusFakeMLSelection_fullEE::isPFIsolatedMuon(int muonIndex) {
  float eta = etaMuon[muonIndex];
  float pt = GetPt(pxMuon[muonIndex],pyMuon[muonIndex]);
  float iso = pfCombinedIsoMuon[muonIndex]/pt;
  if( pt>=20. && fabs(eta)<1.479 ) return (iso < 0.13);
  if( pt>=20. && fabs(eta)>=1.479 ) return (iso < 0.09);
  if( pt<20. && fabs(eta)<1.479 ) return (iso < 0.06);
  if( pt<20. && fabs(eta)>=1.479 ) return (iso < 0.05);
  return true;
}
