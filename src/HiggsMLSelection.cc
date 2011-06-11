#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/HiggsMLSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"
//#include "Mt2/SUSYPhys_Mt2_222_Calculator.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <TTree.h>

using namespace bits;

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
  
  // selection efficiencies
  std::string fileCutsEE     = higgsConfigDirMass + "2e2nuCuts.txt";
  std::string fileSwitchesEE = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionEE.Configure(fileCutsEE.c_str(),fileSwitchesEE.c_str(),"FULL SELECTION EVENT COUNTER EE"); 
  _selectionEE = CutBasedHiggsSelectionEE.GetSelection();  

  std::string fileCutsMM     = higgsConfigDirMass + "2mu2nuCuts.txt";
  std::string fileSwitchesMM = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionMM.Configure(fileCutsMM.c_str(),fileSwitchesMM.c_str(),"FULL SELECTION EVENT COUNTER MM"); 
  _selectionMM = CutBasedHiggsSelectionMM.GetSelection();

  std::string fileCutsEM     = higgsConfigDirMass + "emu2nuCuts.txt";
  std::string fileSwitchesEM = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionEM.Configure(fileCutsEM.c_str(),fileSwitchesEM.c_str(),"FULL SELECTION EVENT COUNTER EM"); 
  _selectionEM = CutBasedHiggsSelectionEM.GetSelection();

  std::string fileCutsME     = higgsConfigDirMass + "emu2nuCuts.txt";
  std::string fileSwitchesME = higgsConfigDir + "2l2nuSwitches.txt";
  CutBasedHiggsSelectionME.Configure(fileCutsME.c_str(),fileSwitchesME.c_str(),"FULL SELECTION EVENT COUNTER ME"); 
  _selectionME = CutBasedHiggsSelectionME.GetSelection();


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
  _selectionEE->addStringParameter("JESUncertainty");

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
  for(int theChannel=0; theChannel<4; theChannel++) {
    m_p4LeptonPlus[theChannel]  = new TLorentzVector(0.,0.,0.,0.);
    m_p4LeptonMinus[theChannel] = new TLorentzVector(0.,0.,0.,0.);
    m_p3PFMET = new TVector3(0.,0.,0.);
    m_p3TKMET = new TVector3(0.,0.,0.);
  }    

  // b-veto event variables
  m_maxDxyEvt = 0.0;
  m_maxDszEvt = 0.0;

  // histo to study jet/electron match
  H_deltaRuncorr = new TH1F("H_deltaRuncorr","uncorrected jets",100, 0.,2*TMath::Pi());
  H_deltaRcorr   = new TH1F("H_deltaRcorr",  "corrected jets",  100, 0.,2*TMath::Pi());
}

HiggsMLSelection::~HiggsMLSelection(){

  for(int theChannel=0; theChannel<4; theChannel++) {  
    delete m_p4LeptonPlus[theChannel];
    delete m_p4LeptonMinus[theChannel];
  }
  delete m_p3PFMET;  
  delete m_p3TKMET;  

  delete _selectionEE;
  delete _selectionMM;
  delete _selectionEM;
  delete _selectionME;
  
  myOutTreeEE -> save();
  myOutTreeMM -> save();
  myOutTreeEM -> save();
  myOutTreeME -> save();

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

    bool isEM = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeplus = imc;
      if( idMc[imc]  ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  ==  11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos>ptMcMuMinus) isEM = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle>ptMcMuPlus) isEM = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isEM) || (indeminus<25 && indmuplus<25 && isEM) );
  }

  // signal: me2nu
  if(strcmp(processType,"HtoWWtome2nu")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;
    
    bool isME = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeplus = imc;
      if( idMc[imc]  ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  ==  11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos<ptMcMuMinus) isME = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle<ptMcMuPlus) isME = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isME) || (indeminus<25 && indmuplus<25 && isME) );
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

    bool isEM = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc]  == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos>ptMcMuMinus) isEM = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle>ptMcMuPlus) isEM = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isEM) || (indeminus<25 && indmuplus<25 && isEM) );
  }

  // signal me excluding taus
  if(strcmp(processType,"HtoWWtome2nu_prompt")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;

    bool isME = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc]  == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    
    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos<ptMcMuMinus) isME = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle<ptMcMuPlus) isME = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isME) || (indeminus<25 && indmuplus<25 && isME) );
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
    // if computed in AOD
    //    weight = evtKfactor;

    // if computed offline
    for(int imc=2;imc<10;imc++) {
      if(idMc[imc]==25 && statusMc[imc]==3) {
        float ptHiggs = pMc[imc]*fabs(sin(thetaMc[imc]));
        return calculator_->evaluate(ptHiggs);
      }
    }
  }
  else if(process.compare("WW")==0) weight = 1.0; // we used MC @ NLO weight in 16X   

  return 1.0;
}

void HiggsMLSelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;
  
  // kinematics reduced tree
  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";
  std::string reducedTreeNameMM = _datasetName+"-datasetMM.root";
  std::string reducedTreeNameEM = _datasetName+"-datasetEM.root";
  std::string reducedTreeNameME = _datasetName+"-datasetME.root";
  myOutTreeEE = new RedHiggsTree(reducedTreeNameEE.c_str());
  myOutTreeMM = new RedHiggsTree(reducedTreeNameMM.c_str());
  myOutTreeEM = new RedHiggsTree(reducedTreeNameEM.c_str());
  myOutTreeME = new RedHiggsTree(reducedTreeNameME.c_str());

  //  myOutTreeEE->addHLTElectronsInfos();
  //  myOutTreeMM->addHLTMuonsInfos();
  //  myOutTreeEM->addHLTElectronsInfos(); myOutTreeEM->addHLTMuonsInfos();

  if ( !_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor") ) {
    myOutTreeEE->addKFactor();
    myOutTreeMM->addKFactor();
    myOutTreeEM->addKFactor();
    myOutTreeME->addKFactor();
  }

  if(!_selectionEE->getSwitch("isData")) {
    myOutTreeEE->addMcTruthInfos();
    myOutTreeMM->addMcTruthInfos();
    myOutTreeEM->addMcTruthInfos();
    myOutTreeME->addMcTruthInfos();
  } else {
    myOutTreeEE->addRunInfos();
    myOutTreeMM->addRunInfos();
    myOutTreeEM->addRunInfos();
    myOutTreeME->addRunInfos();
  }

  myOutTreeEE->addMLVars();
  myOutTreeMM->addMLVars();
  myOutTreeEM->addMLVars();
  myOutTreeME->addMLVars();

  myOutTreeEE->addLatinos();
  myOutTreeMM->addLatinos();
  myOutTreeEM->addLatinos();
  myOutTreeME->addLatinos();

  myOutTreeEE->addElectronInfos();
  myOutTreeEM->addElectronInfos();
  myOutTreeME->addElectronInfos();

  // trigger reduced tree
  //  std::string reducedTriggerTreeName = _datasetName+"-trigger.root";
  //  myTriggerTree = new RedTriggerTree(reducedTriggerTreeName.c_str());

  // eleId reduced tree
  std::string reducedEleIdTreeName = _datasetName+"-eleId.root";
  myEleIdTree = new RedEleIDTree(reducedEleIdTreeName.c_str());


  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  /// kfactors
  std::ifstream file("config/higgs/higgsMass.txt");
  std::string var;
  int mh;
  if(!file.good()) {
    std::cout << "Error! Unable to open the mass file. Using kFactor = 1 always!" << std::endl;   
  } else {
    while(!file.eof()) {
      file >> var >> mh;
    }
  }
  
  cout << "higgs mass = " << mh << endl;

  calculator_ = new kFactorEvaluator(mh);

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

    // get the kFactor of the event (for signal)
    float weight = 1;
    float evtKfactor = 1.0;
    
    // weight for the PU observed in 2011 data
    if ( !_selectionEE->getSwitch("isData") ) weight *= fPUWeight->GetWeight(nPU);

    if (!_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor")) {
      evtKfactor = getkFactor("Higgs");
      weight *= evtKfactor;
    }

    // look to the MC truth decay tree 
    // bool decayEE = findMcTree("HtoWWto2e2nu");
    // bool decayMM = findMcTree("HtoWWto2m2nu");
    // bool decayEM = findMcTree("HtoWWtoem2nu");

    bool promptEE, promptMM, promptEM, promptME;
    promptEE = promptMM = promptEM = promptME = false;
    if( !_selectionEE->getSwitch("isData") ) {
      promptEE = findMcTree("HtoWWto2e2nu_prompt");
      promptMM = findMcTree("HtoWWto2m2nu_prompt");
      promptEM = findMcTree("HtoWWtoem2nu_prompt");  
      promptME = findMcTree("HtoWWtome2nu_prompt");  
    }

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
    //    bool passedHLT = (_selectionEE->getSwitch("isData")) ? hasPassedHLT() : true;
    bool passedHLT[4];
    passedHLT[ee] = hasPassedHLT(ee);
    passedHLT[mm] = hasPassedHLT(mm);
    passedHLT[em] = hasPassedHLT(em);
    passedHLT[me] = hasPassedHLT(em); // same triggers for em and me

    // trigger tree
    // myTriggerTree->fillMcTruth(decayEE,decayMM,decayEM,promptEE,promptMM,promptEM);
    // myTriggerTree->fillHLTElectrons( firedTrg[m_requiredTriggers[0]] );
    // myTriggerTree->fillHLTMuons( firedTrg[m_requiredTriggers[1]] );
    // myTriggerTree->store();




    // -------------------------------------------------------------
    
    // get the best electrons and best muons ==> tu be used to select ALL the possible channels at the beginning only
    std::pair<int,int> thePreElectrons = getBestElectronPair_acceptance();
    std::pair<int,int> thePreMuons     = getBestMuonPair_acceptance();
    // filling vectors with ele-muons passing the acceptance cuts   
    std::pair<int,int> theBestAcceptEleMuon = getBestEleMuonPair(_acceptEleAll,_acceptMuonsAll);
    // filling vectors with ele-muons passing the acceptance cuts
    std::pair<int,int> theBestAcceptMuonEle = getBestMuonElePair(_acceptEleAll,_acceptMuonsAll);

    thePreElectron  = thePreElectrons.second;
    thePrePositron  = thePreElectrons.first;
    thePreMuonPlus  = thePreMuons.first;
    thePreMuonMinus = thePreMuons.second;
    thePreElectronEM = theBestAcceptEleMuon.first;
    thePreMuonEM = theBestAcceptEleMuon.second;
    thePreElectronME = theBestAcceptMuonEle.second;
    thePreMuonME = theBestAcceptMuonEle.first;

    // reconstructed channel
    m_channel[ee] = false;     
    m_channel[mm] = false;
    m_channel[em] = false;
    m_channel[me] = false;
    
    // at this level the SELECTED channel should have pT > 10 and > 20. So far, at least 2 leptons with pT >20 and 10 in the event
    if (thePreElectron > -1 && thePrePositron > -1) {
      float thisMaxPt = TMath::Max(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      float thisMinPt = TMath::Min(GetPt(pxEle[thePreElectron],pyEle[thePreElectron]),GetPt(pxEle[thePrePositron],pyEle[thePrePositron]));
      if (thisMaxPt>20 && thisMinPt>10) m_channel[ee] = true;    // fixme: hardcoded
    }

    if (thePreMuonPlus > -1 && thePreMuonMinus > -1) {
      float thisMaxPt = TMath::Max(GetPt(pxMuon[thePreMuonMinus],pyMuon[thePreMuonMinus]),GetPt(pxMuon[thePreMuonPlus],pyMuon[thePreMuonPlus]));
      if (thisMaxPt>20) m_channel[mm] = true;    // fixme: hardcoded
    }

    if ( thePreElectronEM > -1 && thePreMuonEM > -1 ) {
      float thisMaxPt = GetPt(pxEle[thePreElectronEM],pyEle[thePreElectronEM]);
      if (thisMaxPt>20) m_channel[em] = true;    // fixme: hardcoded
    }

    if ( thePreElectronME > -1 && thePreMuonME > -1 ) {
      float thisMaxPt  = GetPt(pxMuon[thePreMuonME],pyMuon[thePreMuonME]);
      float thisMinPt  = GetPt(pxEle[thePreElectronME],pyEle[thePreElectronME]);
      if (thisMaxPt>20 && thisMinPt>10) m_channel[me] = true;    // fixme: hardcoded
    }
    
    if (_verbose) {
      std::cout << "nEle = "   << nEle << "\tnMuon = " << nMuon << std::endl;
      std::cout << "indices: " << thePreElectron << " " << thePrePositron << " " << thePreMuonMinus << " " << thePreMuonPlus << std::endl;
      std::cout << "chargeEle = " << chargeEle[thePreElectron] << "\tchargePos = " << chargeEle[thePrePositron] 
		<< "\tchargeMuonMinus = " << chargeMuon[thePreMuonMinus] << "\tchargeMuonPlus = " << chargeMuon[thePreMuonPlus] << std::endl;
      std::cout << "ee = " << m_channel[ee] << "\tmm = " << m_channel[mm] << "\temu = " << m_channel[em] << "\tmue = " << m_channel[me] << std::endl; 
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
    thePositron = theBestIpEle.first;    
    theElectron = theBestIpEle.second;
    
    // to be used on the following
    int theIdElectron(theBestIdEle.second);
    int theIdPositron(theBestIdEle.first);    
    int theIsolElectron(theBestIsolEle.second);
    int theIsolPositron(theBestIsolEle.first);    
    int theConvElectron(theBestConvEle.second);
    int theConvPositron(theBestConvEle.first);    
    int theIpElectron(theBestIpEle.second);
    int theIpPositron(theBestIpEle.first);    
    
    // filling the three to compare the distribution before the cut
    int scElectron = superClusterIndexEle[theElectron];
    int scPositron = superClusterIndexEle[thePositron];
    if (theElectron > -1) myEleIdTree -> fillAll(classificationEle[theElectron], hOverEEle[theElectron], eSuperClusterOverPEle[theElectron], eSeedOverPoutEle[theElectron], deltaEtaAtVtxEle[theElectron], deltaPhiAtVtxEle[theElectron], sqrt(covIEtaIEtaSC[scElectron]));
    if (thePositron > -1) myEleIdTree -> fillAll(classificationEle[thePositron], hOverEEle[thePositron], eSuperClusterOverPEle[thePositron], eSeedOverPoutEle[thePositron], deltaEtaAtVtxEle[thePositron], deltaPhiAtVtxEle[thePositron], sqrt(covIEtaIEtaSC[scPositron]));
    myEleIdTree->store();


    // -------------------------------------------------------------
    // MM candidates: preparing vectors of candidates and selecting the two highest pT mu- and mu+ after each step - to check the 20 GeV cut after 

    // muID, for muons in acceptance
    std::pair<int,int> theBestIdMuon = getBestMuonPair_id(_acceptMuonsAll); 

    // isolation, for identified muons
    std::pair<int,int> theBestIsolMuon = getBestMuonPair_isol(_idMuonsAll); 

    // transverse impact parameter, for isolated muons
    std::pair<int,int> theBestIpMuon = getBestMuonPair_ip(_isolMuonsAll);     

    // the two highest pT muons at this point are those I use for my analysis since the passed the full lepton selection
    theMuonPlus  = theBestIpMuon.first;
    theMuonMinus = theBestIpMuon.second;    

    // to be used in the following
    int theIdMuonMinus(theBestIdMuon.second);
    int theIdMuonPlus(theBestIdMuon.first);
    int theIsolMuonMinus(theBestIsolMuon.second);
    int theIsolMuonPlus(theBestIsolMuon.first);
    int theIpMuonMinus(theBestIpMuon.second);
    int theIpMuonPlus(theBestIpMuon.first);


    // -------------------------------------------------------------
    // EM candidates: preparing vectors of candidates and selecting the two highest pT ele+- and muon-+ after each step - to check the 20 GeV cut after
    
    // leptonID, for leptons in acceptance
    std::pair<int,int> theBestIdEleMuon = getBestEleMuonPair(_idEleAll,_idMuonsAll);
    
    // isolation, for identified leptons
    std::pair<int,int> theBestIsolEleMuon = getBestEleMuonPair(_isolEleAll,_isolMuonsAll);

    // conversion rejection, for isolated leptons
    std::pair<int,int> theBestConvEleMuon = getBestEleMuonPair(_convEleAll,_isolMuonsAll);

    // transverse impact parameter, for leptons passing conversion rejection
    std::pair<int,int> theBestIpEleMuon = getBestEleMuonPair(_ipEleAll,_ipMuonsAll);


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



    // -------------------------------------------------------------
    // set of kinematics: : now I've all the final leptons 
    resetKinematics();
    
    // MET is an event variable. Independent o the channel
    m_p3PFMET->SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);
    m_p3TKMET->SetXYZ(pxChMetPV[0],pyChMetPV[0],pzChMetPV[0]); // the one associated to the 0th vertex
    m_theMET = m_p3PFMET->Pt();

    setKinematicsEE(theElectron, thePositron);
    setKinematicsMM(theMuonMinus, theMuonPlus);
    setKinematicsEMME(theElectron, thePositron, theMuonPlus, theMuonMinus);

    
    // -------------------------------------------------------------    
    // look for PV in the event (there is always at least 1 PV)
    m_closestPV = getPV();    // fixme: si chiama closest ma e' quello a piu' alto pT. 
    
    int njets[4], nuncorrjets[4];
    float dphiLLJ[4];
    float btag[4];
    int nsoftmu[4],nextraleptons[4];
    for(int ichan=0; ichan<4; ichan++) {

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

    float genPtHiggs = -1.;
    if ( !_selectionEE->getSwitch("isData") ) {
      for(int imc=2;imc<10;imc++) {
        if(idMc[imc]==25 && statusMc[imc]==3) genPtHiggs = pMc[imc]*fabs(sin(thetaMc[imc]));
      }}
    

    // ---------------------------------------
    // filling counters for the different final states

    // EE
    CutBasedHiggsSelectionEE.SetWeight(weight);               
    CutBasedHiggsSelectionEE.SetMcTruth(promptEE);  
    CutBasedHiggsSelectionEE.SetHLT(passedHLT[ee]);               
    CutBasedHiggsSelectionEE.SetIsChannel(m_channel[ee]);     
    
    CutBasedHiggsSelectionEE.SetElectronId(theIdElectron);                 
    CutBasedHiggsSelectionEE.SetPositronId(theIdPositron);                 
    CutBasedHiggsSelectionEE.SetElectronIsolation(theIsolElectron);        
    CutBasedHiggsSelectionEE.SetPositronIsolation(theIsolPositron);        
    CutBasedHiggsSelectionEE.SetElectronConvRejection(theConvElectron);    
    CutBasedHiggsSelectionEE.SetPositronConvRejection(theConvPositron);    
    CutBasedHiggsSelectionEE.SetElectronIp(theIpElectron);                 
    CutBasedHiggsSelectionEE.SetPositronIp(theIpPositron);                 
    // checking if the highest pT electron at each step has pT>20
    float thisMaxPtIdEE   = TMath::Max(GetPt(pxEle[theIdElectron],pyEle[theIdElectron]),GetPt(pxEle[theIdPositron],pyEle[theIdPositron]));
    float thisMaxPtIsolEE = TMath::Max(GetPt(pxEle[theIsolElectron],pyEle[theIsolElectron]),GetPt(pxEle[theIsolPositron],pyEle[theIsolPositron]));
    float thisMaxPtConvEE = TMath::Max(GetPt(pxEle[theConvElectron],pyEle[theConvElectron]),GetPt(pxEle[theConvPositron],pyEle[theConvPositron]));
    float thisMaxPtIpEE   = TMath::Max(GetPt(pxEle[theIpElectron],pyEle[theIpElectron]),GetPt(pxEle[theIpPositron],pyEle[theIpPositron]));
    if (thisMaxPtIdEE<20)   { 
      CutBasedHiggsSelectionEE.SetElectronId(-1);
      CutBasedHiggsSelectionEE.SetPositronId(-1);
    }
    if (thisMaxPtIsolEE<20) { 
      CutBasedHiggsSelectionEE.SetElectronIsolation(-1);
      CutBasedHiggsSelectionEE.SetPositronIsolation(-1);
    }
    if (thisMaxPtConvEE<20) { 
      CutBasedHiggsSelectionEE.SetElectronConvRejection(-1);
      CutBasedHiggsSelectionEE.SetPositronConvRejection(-1);
    }
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
    CutBasedHiggsSelectionEE.SetWWInvMass(2.*m_transvMass[ee]/_massVal);

    bool isSelectedEE           = CutBasedHiggsSelectionEE.output();    
    bool selUpToFinalLeptonsEE  = CutBasedHiggsSelectionEE.outputUpToFinalLeptons();
    bool selUpToJetVetoEE       = CutBasedHiggsSelectionEE.outputUpToJetVeto();
    bool selUpToUncorrJetVetoEE = CutBasedHiggsSelectionEE.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiEE       = CutBasedHiggsSelectionEE.outputPreDeltaPhi();

    // latinos
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


    // eleID variables to fill the tree (after each cut)
    if( GetPt(pxEle[thePreElectron],pyEle[thePreElectron]) > GetPt(pxEle[thePrePositron],pyEle[thePrePositron]) ) setEleIdVariables(thePreElectron,thePrePositron);
    else setEleIdVariables(thePrePositron, thePreElectron);

    if(!_selectionEE->getSwitch("isData")) myOutTreeEE -> fillMcTruth(promptEE);

    myOutTreeEE->fillRunInfos(runNumber, lumiBlock, eventNumber, weight);

    if ( !_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor") ) {
      int theLJ  = theLeadingJet[ee];
      float ptLJ = sqrt(pxAK5PFPUcorrJet[theLJ]*pxAK5PFPUcorrJet[theLJ] + pyAK5PFPUcorrJet[theLJ]*pyAK5PFPUcorrJet[theLJ]);
      myOutTreeEE->fillKFactor(evtKfactor, genPtHiggs, ptLJ);
    }

    //       myOutTreeEE -> fillHLTElectrons( firedTrg[m_requiredTriggers[0]], 
    // 				       firedTrg[m_requiredTriggers[1]],
    // 				       (firedTrg[m_requiredTriggers[0]] || firedTrg[m_requiredTriggers[1]]) );
    
    myOutTreeEE -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[ee], m_deltaPhi[ee], m_deltaErre[ee], m_transvMass[ee], m_mll[ee], 
			   hardestLeptonPt[ee], slowestLeptonPt[ee], m_deltaEtaLeptons[ee], nPV,
			   selUpToFinalLeptonsEE, selUpToJetVetoEE, selUpToUncorrJetVetoEE, selPreDeltaPhiEE, isSelectedEE);

    myOutTreeEE -> fillMLVars(njets[ee], nuncorrjets[ee], m_maxDxyEvt, m_maxDszEvt, btag[ee], m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);

    myOutTreeEE -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 ); 
    
    myOutTreeEE -> fillElectrons( myRecoflag, myPt, myEta, myPhi,
				  myClassification, myNBremClusters, myDeta, myDphi, myHoe, mySee, mySpp, myEop, myFbrem,
				  myTrackerIso, myHcalIso, myEcalJIso, myEcalGTIso, myCombinedIso, myCharge, myMissHits, myDist, myDcot, myLh, myMatched );
    

      
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeEE -> store();


    
    // ---------------------------------------
    // MM
    CutBasedHiggsSelectionMM.SetWeight(weight);               
    CutBasedHiggsSelectionMM.SetMcTruth(promptMM); 
    CutBasedHiggsSelectionMM.SetHLT(passedHLT[mm]);               
    CutBasedHiggsSelectionMM.SetIsChannel(m_channel[mm]);     
    
    CutBasedHiggsSelectionMM.SetElectronId(theIdMuonMinus);                 
    CutBasedHiggsSelectionMM.SetPositronId(theIdMuonPlus);                 
    CutBasedHiggsSelectionMM.SetElectronIsolation(theIsolMuonMinus);        
    CutBasedHiggsSelectionMM.SetPositronIsolation(theIsolMuonPlus);        
    CutBasedHiggsSelectionMM.SetElectronConvRejection(true);    
    CutBasedHiggsSelectionMM.SetPositronConvRejection(true);    
    CutBasedHiggsSelectionMM.SetElectronIp(theIpMuonMinus);                 
    CutBasedHiggsSelectionMM.SetPositronIp(theIpMuonPlus);                 
    // checking if the highest pT electron at each step has pT>20
    float thisMaxPtIdMM   = TMath::Max(GetPt(pxMuon[theIdMuonMinus],pyMuon[theIdMuonMinus]),GetPt(pxMuon[theIdMuonPlus],pyMuon[theIdMuonPlus]));
    float thisMaxPtIsolMM = TMath::Max(GetPt(pxMuon[theIsolMuonMinus],pyMuon[theIsolMuonMinus]),GetPt(pxMuon[theIsolMuonPlus],pyMuon[theIsolMuonPlus]));
    float thisMaxPtIpMM   = TMath::Max(GetPt(pxMuon[theIpMuonMinus],pyMuon[theIpMuonMinus]),GetPt(pxMuon[theIpMuonPlus],pyMuon[theIpMuonPlus]));
    if (thisMaxPtIdMM<20)   { 
      CutBasedHiggsSelectionMM.SetElectronId(-1);
      CutBasedHiggsSelectionMM.SetPositronId(-1);
    }
    if (thisMaxPtIsolMM<20) { 
      CutBasedHiggsSelectionMM.SetElectronIsolation(-1);
      CutBasedHiggsSelectionMM.SetPositronIsolation(-1);
    }
    if (thisMaxPtIpMM<20)   { 
      CutBasedHiggsSelectionMM.SetElectronIp(-1);
      CutBasedHiggsSelectionMM.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionMM.SetHighElePt(hardestLeptonPt[mm]); 
    CutBasedHiggsSelectionMM.SetLowElePt(slowestLeptonPt[mm]);  

    CutBasedHiggsSelectionMM.SetNJets(njets[mm]);
    CutBasedHiggsSelectionMM.SetNUncorrJets(nuncorrjets[mm]);
    CutBasedHiggsSelectionMM.SetBTagJets(btag[mm]);
    CutBasedHiggsSelectionMM.SetNSoftMuons(nsoftmu[mm]);
    CutBasedHiggsSelectionMM.SetNExtraLeptons(nextraleptons[mm]);
    CutBasedHiggsSelectionMM.SetMet(m_theMET);					
    CutBasedHiggsSelectionMM.SetProjectedMet(m_projectedMet[mm]);
    CutBasedHiggsSelectionMM.SetMetOverPtLL(m_metOptll[mm]);
    CutBasedHiggsSelectionMM.SetDeltaPhiLLJet(dphiLLJ[mm]);   
    CutBasedHiggsSelectionMM.SetDeltaPhi(m_deltaPhi[mm]);
    CutBasedHiggsSelectionMM.SetInvMass(m_mll[mm]);
    CutBasedHiggsSelectionMM.SetDetaLeptons(m_deltaEtaLeptons[mm]);
    CutBasedHiggsSelectionMM.SetWWInvMass(2.*m_transvMass[mm]/_massVal);

    bool isSelectedMM           = CutBasedHiggsSelectionMM.output();    
    bool selUpToFinalLeptonsMM  = CutBasedHiggsSelectionMM.outputUpToFinalLeptons();
    bool selUpToJetVetoMM       = CutBasedHiggsSelectionMM.outputUpToJetVeto();
    bool selUpToUncorrJetVetoMM = CutBasedHiggsSelectionMM.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiMM       = CutBasedHiggsSelectionMM.outputPreDeltaPhi();

    // latinos
    outputStep0  = CutBasedHiggsSelectionMM.outputStep0();
    outputStep1  = CutBasedHiggsSelectionMM.outputStep1();
    outputStep2  = CutBasedHiggsSelectionMM.outputStep2();
    outputStep3  = CutBasedHiggsSelectionMM.outputStep3();
    outputStep4  = CutBasedHiggsSelectionMM.outputStep4();
    outputStep5  = CutBasedHiggsSelectionMM.outputStep5();
    outputStep6  = CutBasedHiggsSelectionMM.outputStep6();
    outputStep7  = CutBasedHiggsSelectionMM.outputStep7();
    outputStep8  = CutBasedHiggsSelectionMM.outputStep8();
    outputStep9  = CutBasedHiggsSelectionMM.outputStep9();
    outputStep10 = CutBasedHiggsSelectionMM.outputStep10();
    outputStep11 = CutBasedHiggsSelectionMM.outputStep11();
    outputStep12 = CutBasedHiggsSelectionMM.outputStep12();
    outputStep13 = CutBasedHiggsSelectionMM.outputStep13();
    outputStep14 = CutBasedHiggsSelectionMM.outputStep14();
    outputStep15 = CutBasedHiggsSelectionMM.outputStep15();
    outputStep16 = CutBasedHiggsSelectionMM.outputStep16();
    outputStep17 = CutBasedHiggsSelectionMM.outputStep17();
    outputStep18 = CutBasedHiggsSelectionMM.outputStep18();
    outputStep19 = CutBasedHiggsSelectionMM.outputStep19();
    outputStep20 = CutBasedHiggsSelectionMM.outputStep20();
    outputStep21 = CutBasedHiggsSelectionMM.outputStep21();
    outputStep22 = CutBasedHiggsSelectionMM.outputStep22();
    outputStep23 = CutBasedHiggsSelectionMM.outputStep23();
    outputStep24 = CutBasedHiggsSelectionMM.outputStep24();


    // filling the tree
    if(!_selectionMM->getSwitch("isData")) myOutTreeMM -> fillMcTruth(promptMM);
    myOutTreeMM->fillRunInfos(runNumber, lumiBlock, eventNumber, weight);

    if ( !_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor") ) {
      int theLJ  = theLeadingJet[mm];
      float ptLJ = sqrt(pxAK5PFPUcorrJet[theLJ]*pxAK5PFPUcorrJet[theLJ] + pyAK5PFPUcorrJet[theLJ]*pyAK5PFPUcorrJet[theLJ]);
      myOutTreeMM->fillKFactor(evtKfactor, genPtHiggs, ptLJ);
    }

    myOutTreeMM -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[mm], m_deltaPhi[mm], m_deltaErre[mm], m_transvMass[mm], m_mll[mm], 
			   hardestLeptonPt[mm], slowestLeptonPt[mm], m_deltaEtaLeptons[mm], nPV,
			   selUpToFinalLeptonsMM, selUpToJetVetoMM, selUpToUncorrJetVetoMM, selPreDeltaPhiMM, isSelectedMM);

    myOutTreeMM -> fillMLVars(njets[mm], nuncorrjets[mm], m_maxDxyEvt, m_maxDszEvt, btag[mm], m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);
    
    myOutTreeMM -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 ); 
    
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeMM -> store();




    // ---------------------------------------
    // EM
    CutBasedHiggsSelectionEM.SetWeight(weight);               
    CutBasedHiggsSelectionEM.SetMcTruth(promptEM);  
    CutBasedHiggsSelectionEM.SetHLT(passedHLT[em]);               
    CutBasedHiggsSelectionEM.SetIsChannel(m_channel[em]);     

    CutBasedHiggsSelectionEM.SetElectronId(theBestIdEleMuon.first);                 
    CutBasedHiggsSelectionEM.SetPositronId(theBestIdEleMuon.second);                 
    
    CutBasedHiggsSelectionEM.SetElectronIsolation(theBestIsolEleMuon.first);        
    CutBasedHiggsSelectionEM.SetPositronIsolation(theBestIsolEleMuon.second);        
    CutBasedHiggsSelectionEM.SetElectronConvRejection(theBestConvEleMuon.first);    
    CutBasedHiggsSelectionEM.SetPositronConvRejection(theBestConvEleMuon.second);    
    CutBasedHiggsSelectionEM.SetElectronIp(theBestIpEleMuon.first);                 
    CutBasedHiggsSelectionEM.SetPositronIp(theBestIpEleMuon.second);                 
    // checking if the highest pT electron at each step has pT>20. E is the hardest in EM
    float thisMaxPtIdEM   = (theBestIdEleMuon.first > -1) ? GetPt(pxEle[theBestIdEleMuon.first],pyEle[theBestIdEleMuon.first]) : 0.;
    float thisMaxPtIsolEM = (theBestIsolEleMuon.first > -1) ? GetPt(pxEle[theBestIsolEleMuon.first],pyEle[theBestIsolEleMuon.first]) : 0.;
    float thisMaxPtIpEM   = (theBestIpEleMuon.first > -1) ? GetPt(pxEle[theBestIpEleMuon.first],pyEle[theBestIpEleMuon.first]) : 0.;
    if (thisMaxPtIdEM<20)   { 
      CutBasedHiggsSelectionEM.SetElectronId(-1);
      CutBasedHiggsSelectionEM.SetPositronId(-1);
    }
    if (thisMaxPtIsolEM<20) { 
      CutBasedHiggsSelectionEM.SetElectronIsolation(-1);
      CutBasedHiggsSelectionEM.SetPositronIsolation(-1);
    }
    if (thisMaxPtIpEM<20)   { 
      CutBasedHiggsSelectionEM.SetElectronIp(-1);
      CutBasedHiggsSelectionEM.SetPositronIp(-1);
    }

    CutBasedHiggsSelectionEM.SetHighElePt(hardestLeptonPt[em]); 
    CutBasedHiggsSelectionEM.SetLowElePt(slowestLeptonPt[em]);  

    CutBasedHiggsSelectionEM.SetNJets(njets[em]);
    CutBasedHiggsSelectionEM.SetNUncorrJets(nuncorrjets[em]);
    CutBasedHiggsSelectionEM.SetBTagJets(btag[em]);
    CutBasedHiggsSelectionEM.SetNSoftMuons(nsoftmu[em]);
    CutBasedHiggsSelectionEM.SetNExtraLeptons(nextraleptons[em]);
    CutBasedHiggsSelectionEM.SetMet(m_theMET);					
    CutBasedHiggsSelectionEM.SetProjectedMet(m_projectedMet[em]);
    CutBasedHiggsSelectionEM.SetMetOverPtLL(m_metOptll[em]);
    CutBasedHiggsSelectionEM.SetDeltaPhiLLJet(dphiLLJ[em]);  
    CutBasedHiggsSelectionEM.SetDeltaPhi(m_deltaPhi[em]);
    CutBasedHiggsSelectionEM.SetInvMass(m_mll[em]);
    CutBasedHiggsSelectionEM.SetDetaLeptons(m_deltaEtaLeptons[em]);
    CutBasedHiggsSelectionEM.SetWWInvMass(2.*m_transvMass[em]/_massVal);

    bool isSelectedEM           = CutBasedHiggsSelectionEM.output();    
    bool selUpToFinalLeptonsEM  = CutBasedHiggsSelectionEM.outputUpToFinalLeptons();
    bool selUpToJetVetoEM       = CutBasedHiggsSelectionEM.outputUpToJetVeto();
    bool selUpToUncorrJetVetoEM = CutBasedHiggsSelectionEM.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiEM       = CutBasedHiggsSelectionEM.outputPreDeltaPhi();

    // latinos
    outputStep0  = CutBasedHiggsSelectionEM.outputStep0();
    outputStep1  = CutBasedHiggsSelectionEM.outputStep1();
    outputStep2  = CutBasedHiggsSelectionEM.outputStep2();
    outputStep3  = CutBasedHiggsSelectionEM.outputStep3();
    outputStep4  = CutBasedHiggsSelectionEM.outputStep4();
    outputStep5  = CutBasedHiggsSelectionEM.outputStep5();
    outputStep6  = CutBasedHiggsSelectionEM.outputStep6();
    outputStep7  = CutBasedHiggsSelectionEM.outputStep7();
    outputStep8  = CutBasedHiggsSelectionEM.outputStep8();
    outputStep9  = CutBasedHiggsSelectionEM.outputStep9();
    outputStep10 = CutBasedHiggsSelectionEM.outputStep10();
    outputStep11 = CutBasedHiggsSelectionEM.outputStep11();
    outputStep12 = CutBasedHiggsSelectionEM.outputStep12();
    outputStep13 = CutBasedHiggsSelectionEM.outputStep13();
    outputStep14 = CutBasedHiggsSelectionEM.outputStep14();
    outputStep15 = CutBasedHiggsSelectionEM.outputStep15();
    outputStep16 = CutBasedHiggsSelectionEM.outputStep16();
    outputStep17 = CutBasedHiggsSelectionEM.outputStep17();
    outputStep18 = CutBasedHiggsSelectionEM.outputStep18();
    outputStep19 = CutBasedHiggsSelectionEM.outputStep19();
    outputStep20 = CutBasedHiggsSelectionEM.outputStep20();
    outputStep21 = CutBasedHiggsSelectionEM.outputStep21();
    outputStep22 = CutBasedHiggsSelectionEM.outputStep22();
    outputStep23 = CutBasedHiggsSelectionEM.outputStep23();
    outputStep24 = CutBasedHiggsSelectionEM.outputStep24();

    // filling the tree
    if(!_selectionEM->getSwitch("isData")) myOutTreeEM -> fillMcTruth(promptEM);

    myOutTreeEM->fillRunInfos(runNumber, lumiBlock, eventNumber, weight);

    if ( !_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor") ) {
      int theLJ  = theLeadingJet[em];
      float ptLJ = sqrt(pxAK5PFPUcorrJet[theLJ]*pxAK5PFPUcorrJet[theLJ] + pyAK5PFPUcorrJet[theLJ]*pyAK5PFPUcorrJet[theLJ]);
      myOutTreeEM->fillKFactor(evtKfactor, genPtHiggs, ptLJ);
    }

    myOutTreeEM -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[em], m_deltaPhi[em], m_deltaErre[em], m_transvMass[em], m_mll[em], 
			   hardestLeptonPt[em], slowestLeptonPt[em], m_deltaEtaLeptons[em], nPV,
			   selUpToFinalLeptonsEM, selUpToJetVetoEM, selUpToUncorrJetVetoEM, selPreDeltaPhiEM, isSelectedEM);

    setEleIdVariables(theBestIpEleMuon.first, -1);
    myOutTreeEM -> fillElectrons( myRecoflag, myPt, myEta, myPhi,
				  myClassification, myNBremClusters, myDeta, myDphi, myHoe, mySee, mySpp, myEop, myFbrem,
				  myTrackerIso, myHcalIso, myEcalJIso, myEcalGTIso, myCombinedIso, myCharge, myMissHits, myDist, myDcot, myLh, myMatched );
    
    myOutTreeEM -> fillMLVars(njets[em], nuncorrjets[em], m_maxDxyEvt, m_maxDszEvt, btag[em], m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);
    
    myOutTreeEM -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 ); 
    
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeEM -> store();


  

    // ---------------------------------------
    // ME
    CutBasedHiggsSelectionME.SetWeight(weight);               
    CutBasedHiggsSelectionME.SetMcTruth(promptME);    
    CutBasedHiggsSelectionME.SetHLT(passedHLT[me]);               
    CutBasedHiggsSelectionME.SetIsChannel(m_channel[me]);     
    
    CutBasedHiggsSelectionME.SetElectronId(theBestIdMuonEle.first);                 
    CutBasedHiggsSelectionME.SetPositronId(theBestIdMuonEle.second);                 
    CutBasedHiggsSelectionME.SetElectronIsolation(theBestIsolMuonEle.first);        
    CutBasedHiggsSelectionME.SetPositronIsolation(theBestIsolMuonEle.second);        
    CutBasedHiggsSelectionME.SetElectronConvRejection(theBestConvMuonEle.first);    
    CutBasedHiggsSelectionME.SetPositronConvRejection(theBestConvMuonEle.second);    
    CutBasedHiggsSelectionME.SetElectronIp(theBestIpMuonEle.first);                 
    CutBasedHiggsSelectionME.SetPositronIp(theBestIpMuonEle.second);                 
    // checking if the highest pT electron at each step has pT>20. MU is the hardest in ME
    float thisMaxPtIdME = (theBestIdMuonEle.first > -1) ? GetPt(pxMuon[theBestIdMuonEle.first],pyMuon[theBestIdMuonEle.first]) : 0.;
    float thisMaxPtIsolME = (theBestIsolMuonEle.first > -1) ? GetPt(pxMuon[theBestIsolMuonEle.first],pyMuon[theBestIsolMuonEle.first]) : 0.;
    float thisMaxPtIpME = (theBestIpMuonEle.first > -1) ? GetPt(pxMuon[theBestIpMuonEle.first],pyMuon[theBestIpMuonEle.first]) : 0.;
    if (thisMaxPtIdME<20)   { 
      CutBasedHiggsSelectionME.SetElectronId(-1);
      CutBasedHiggsSelectionME.SetPositronId(-1);
    }
    if (thisMaxPtIsolME<20) { 
      CutBasedHiggsSelectionME.SetElectronIsolation(-1);
      CutBasedHiggsSelectionME.SetPositronIsolation(-1);
    }
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
    CutBasedHiggsSelectionME.SetWWInvMass(2.*m_transvMass[me]/_massVal);

    bool isSelectedME           = CutBasedHiggsSelectionME.output();    
    bool selUpToFinalLeptonsME  = CutBasedHiggsSelectionME.outputUpToFinalLeptons();
    bool selUpToJetVetoME       = CutBasedHiggsSelectionME.outputUpToJetVeto();
    bool selUpToUncorrJetVetoME = CutBasedHiggsSelectionME.outputUpToUncorrJetVeto();
    bool selPreDeltaPhiME       = CutBasedHiggsSelectionME.outputPreDeltaPhi();

    // latinos
    outputStep0  = CutBasedHiggsSelectionME.outputStep0();
    outputStep1  = CutBasedHiggsSelectionME.outputStep1();
    outputStep2  = CutBasedHiggsSelectionME.outputStep2();
    outputStep3  = CutBasedHiggsSelectionME.outputStep3();
    outputStep4  = CutBasedHiggsSelectionME.outputStep4();
    outputStep5  = CutBasedHiggsSelectionME.outputStep5();
    outputStep6  = CutBasedHiggsSelectionME.outputStep6();
    outputStep7  = CutBasedHiggsSelectionME.outputStep7();
    outputStep8  = CutBasedHiggsSelectionME.outputStep8();
    outputStep9  = CutBasedHiggsSelectionME.outputStep9();
    outputStep10 = CutBasedHiggsSelectionME.outputStep10();
    outputStep11 = CutBasedHiggsSelectionME.outputStep11();
    outputStep12 = CutBasedHiggsSelectionME.outputStep12();
    outputStep13 = CutBasedHiggsSelectionME.outputStep13();
    outputStep14 = CutBasedHiggsSelectionME.outputStep14();
    outputStep15 = CutBasedHiggsSelectionME.outputStep15();
    outputStep16 = CutBasedHiggsSelectionME.outputStep16();
    outputStep17 = CutBasedHiggsSelectionME.outputStep17();
    outputStep18 = CutBasedHiggsSelectionME.outputStep18();
    outputStep19 = CutBasedHiggsSelectionME.outputStep19();
    outputStep20 = CutBasedHiggsSelectionME.outputStep20();
    outputStep21 = CutBasedHiggsSelectionME.outputStep21();
    outputStep22 = CutBasedHiggsSelectionME.outputStep22();
    outputStep23 = CutBasedHiggsSelectionME.outputStep23();
    outputStep24 = CutBasedHiggsSelectionME.outputStep24();

    // filling the tree
    if(!_selectionME->getSwitch("isData")) myOutTreeME -> fillMcTruth(promptME);

    myOutTreeME->fillRunInfos(runNumber, lumiBlock, eventNumber, weight);

    if ( !_selectionEE->getSwitch("isData") && _selectionEE->getSwitch("apply_kFactor") ) {
      int theLJ  = theLeadingJet[me];
      float ptLJ = sqrt(pxAK5PFPUcorrJet[theLJ]*pxAK5PFPUcorrJet[theLJ] + pyAK5PFPUcorrJet[theLJ]*pyAK5PFPUcorrJet[theLJ]);
      myOutTreeME->fillKFactor(evtKfactor, genPtHiggs, ptLJ);
    }

    myOutTreeME -> fillAll(GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), GetPt(pxMet[0],pyMet[0]), 
			   m_projectedMet[me], m_deltaPhi[me], m_deltaErre[me], m_transvMass[me], m_mll[me], 
			   hardestLeptonPt[me], slowestLeptonPt[me], m_deltaEtaLeptons[me], nPV,
			   selUpToFinalLeptonsME, selUpToJetVetoME, selUpToUncorrJetVetoME, selPreDeltaPhiME, isSelectedME);

    setEleIdVariables(theBestIpMuonEle.second, -1);
    myOutTreeME -> fillElectrons( myRecoflag, myPt, myEta, myPhi,
				  myClassification, myNBremClusters, myDeta, myDphi, myHoe, mySee, mySpp, myEop, myFbrem,
				  myTrackerIso, myHcalIso, myEcalJIso, myEcalGTIso, myCombinedIso, myCharge, myMissHits, myDist, myDcot, myLh, myMatched );

    myOutTreeME -> fillMLVars(njets[me], nuncorrjets[me], m_maxDxyEvt, m_maxDszEvt, btag[me], m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags);
    
    myOutTreeME -> fillLatinos( outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16, outputStep17, outputStep18, outputStep19, outputStep20, outputStep21, outputStep22, outputStep23, outputStep24 ); 
    
    // dumping final tree, only if there are 2 leptons in the acceptance
    if(outputStep1) myOutTreeME -> store();

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
  std::cout << "Full EE selections: " << std::endl;
  CutBasedHiggsSelectionEE.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full MM selections: " << std::endl;
  CutBasedHiggsSelectionMM.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EM selections: " << std::endl;
  CutBasedHiggsSelectionEM.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full ME selections: " << std::endl;
  CutBasedHiggsSelectionME.displayEfficiencies(datasetName);


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

std::pair<int,int> HiggsMLSelection::getBestElectronPair_acceptance() {
  
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

std::pair<int,int> HiggsMLSelection::getBestElectronPair_id( std::vector<int> acceptEle ) {

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
      if (thisPt>=15) isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
      if (thisPt<15)  isEleID(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedIDLow);
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

std::pair<int,int> HiggsMLSelection::getBestElectronPair_isol( std::vector<int> idEle ) {

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

std::pair<int,int> HiggsMLSelection::getBestElectronPair_conv( std::vector<int> isolEle ) {

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


std::pair<int,int> HiggsMLSelection::getBestElectronPair_ip( std::vector<int> convEle ) {

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

std::pair<int,int> HiggsMLSelection::getBestMuonPair_acceptance() {
  
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

    if(_selectionEE->getSwitch("etaMuonAcc") && !_selectionEE->passCut("etaMuonAcc",etaMuon[i]) ) continue;

    float thisPt = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEE->getSwitch("ptMuonAcc") && !_selectionEE->passCut("ptMuonAcc",thisPt) ) continue;

    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
    
    _acceptMuonsAll.push_back(i);  
  }
  _acceptMuonsAll = sortMuonsByPt(_acceptMuonsAll);

  return make_pair(theLep1,theLep2);
}


std::pair<int,int> HiggsMLSelection::getBestMuonPair_id( std::vector<int> acceptMu ) {
  
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

std::pair<int,int> HiggsMLSelection::getBestMuonPair_isol( std::vector<int> idMu ) {
  
  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;
   
  _isolMuonsAll.clear();

  for(int iMu=0; iMu<idMu.size(); iMu++) {

    int thisMu   = idMu[iMu];
    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);

    // fixme: diverso da prima: rimuovevo il secondo leptone....
    float muonTrackerForGlobal = sumPt03Muon[thisMu];
    float muonEcalForGlobal    = emEt03Muon[thisMu];
    float muonHcalForGlobal    = hadEt03Muon[thisMu]; 
    float theMuonGlobalSum     = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal - rhoFastjet*TMath::Pi()*0.3*0.3;
    float theRelMuonIso        = theMuonGlobalSum/thisPt; 
    if(_selectionEE->getSwitch("muGlobalIso") && !_selectionEE->passCut("muGlobalIso",theRelMuonIso)) continue;  

    float thisCharge = chargeMuon[thisMu];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisMu; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisMu; }

    _isolMuonsAll.push_back(thisMu);   
  }
  _isolMuonsAll = sortMuonsByPt(_isolMuonsAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> HiggsMLSelection::getBestMuonPair_ip( std::vector<int> isoMu ) {

  int theLep1 = -1;
  int theLep2 = -1;
  float maxPtLep1 = -1000.;
  float maxPtLep2 = -1000.;

  _ipMuonsAll.clear();  

  for(int iMu=0; iMu<isoMu.size(); iMu++) {

    int thisMu = isoMu[iMu];
    
    int ctfMuon   = trackIndexMuon[thisMu]; 
    float dxyMuon = transvImpactParTrack[ctfMuon];
    float dzMuon  = PVzPV[m_closestPV] - trackVzTrack[ctfMuon];   
    if (_selectionEE->getSwitch("muonIP") && (!_selectionEE->passCut("muonIP",dxyMuon)) ) continue;   
    if (_selectionEE->getSwitch("muonDz") && (!_selectionEE->passCut("muonDz",dzMuon)) )  continue;   

    float thisPt     = GetPt(pxMuon[thisMu],pyMuon[thisMu]);
    float thisCharge = chargeMuon[thisMu];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = thisMu; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = thisMu; }

    _ipMuonsAll.push_back(thisMu);   
  }
  _ipMuonsAll = sortMuonsByPt(_ipMuonsAll);

  return make_pair(theLep1,theLep2);
}

std::pair<int,int> HiggsMLSelection::getBestEleMuonPair(std::vector<int> electrons, std::vector<int> muons) {

  int theEle=-1;
  int theMuon=-1;

  std::vector<int>::const_iterator eleiter;
  for(eleiter=electrons.begin(); eleiter<electrons.end();++eleiter) {
    int eleCharge = chargeEle[*eleiter];
    float elePt = GetPt(pxEle[*eleiter],pyEle[*eleiter]);
    theEle = *eleiter;
    std::vector<int>::const_iterator muiter;
    for(muiter=muons.begin(); muiter<muons.end();++muiter) {
      int muCharge = chargeMuon[*muiter];
      float muPt = GetPt(pxMuon[*muiter],pyMuon[*muiter]);
      if(eleCharge*muCharge<0 && elePt > muPt) return std::make_pair(*eleiter,*muiter); 
    }
  }
  
  return std::make_pair(theEle,theMuon);
}

std::pair<int,int> HiggsMLSelection::getBestMuonElePair(std::vector<int> electrons, std::vector<int> muons) {

  int theEle=-1;
  int theMuon=-1;

  std::vector<int>::const_iterator muiter;
  for(muiter=muons.begin(); muiter<muons.end();++muiter) {
    int muCharge = chargeMuon[*muiter];
    float muPt = GetPt(pxMuon[*muiter],pyMuon[*muiter]);
    theMuon = *muiter;
    std::vector<int>::const_iterator eleiter;
    for(eleiter=electrons.begin(); eleiter<electrons.end();++eleiter) {
      int eleCharge = chargeEle[*eleiter];
      float elePt = GetPt(pxEle[*eleiter],pyEle[*eleiter]);
      if(eleCharge*muCharge<0 && muPt > elePt) return std::make_pair(*muiter,*eleiter); 
    }
  }
  return std::make_pair(theMuon,theEle); 
}

void HiggsMLSelection::setKinematicsEE(int myEle, int myPosi) {

  if (myPosi > -1 && myEle > -1) {

    eleCands[ee].push_back(myEle);
    eleCands[ee].push_back(myPosi);   
    hardestLeptonPt[ee] = TMath::Max(GetPt(pxEle[myEle],pyEle[myEle]),GetPt(pxEle[myPosi],pyEle[myPosi]));
    slowestLeptonPt[ee] = TMath::Min(GetPt(pxEle[myEle],pyEle[myEle]),GetPt(pxEle[myPosi],pyEle[myPosi]));
    m_p4LeptonMinus[ee] -> SetXYZT(pxEle[myEle], pyEle[myEle], pzEle[myEle], energyEle[myEle]);
    m_p4LeptonPlus[ee]  -> SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
    m_mll[ee]       = (*(m_p4LeptonMinus[ee]) + *(m_p4LeptonPlus[ee])).M();
    m_deltaPhi[ee]  = fabs(180./TMath::Pi() * m_p4LeptonMinus[ee]->Vect().DeltaPhi(m_p4LeptonPlus[ee]->Vect()));
    m_deltaErre[ee] = m_p4LeptonMinus[ee]->Vect().DeltaR(m_p4LeptonPlus[ee]->Vect());
    m_deltaEtaLeptons[ee] = etaEle[myEle]-etaEle[myPosi];
    m_dilepPt[ee].SetXYZ( m_p4LeptonMinus[ee]->Vect().X()+m_p4LeptonPlus[ee]->Vect().X(),m_p4LeptonMinus[ee]->Vect().Y()+m_p4LeptonPlus[ee]->Vect().Y(),0.0 );
    // def. 3 of http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=104213                           
    //m_transvMass[ee]=mT3(*m_p4LeptonMinus[ee],*m_p4LeptonPlus[ee],m_p3MET);
    m_transvMass[ee]=CalcGammaMRstar(*m_p4LeptonMinus[ee],*m_p4LeptonPlus[ee]);
    m_metOptll[ee] = m_theMET / m_dilepPt[ee].Pt();
    m_mT2[ee] = 0.;
    m_projectedMet[ee] = GetProjectedMet(m_p4LeptonMinus[ee]->Vect(),m_p4LeptonPlus[ee]->Vect());
  }

}

void HiggsMLSelection::setKinematicsMM(int myMuMinus, int myMuPlus) {
  
  if (myMuPlus > -1 && myMuMinus > -1){

    muCands[mm].push_back(myMuMinus);
    muCands[mm].push_back(myMuPlus);
    hardestLeptonPt[mm] = TMath::Max(GetPt(pxMuon[myMuPlus],pyMuon[myMuPlus]),GetPt(pxMuon[myMuMinus],pyMuon[myMuMinus]));
    slowestLeptonPt[mm] = TMath::Min(GetPt(pxMuon[myMuPlus],pyMuon[myMuPlus]),GetPt(pxMuon[myMuMinus],pyMuon[myMuMinus]));
    m_p4LeptonMinus[mm] -> SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
    m_p4LeptonPlus[mm]  -> SetXYZT(pxMuon[myMuPlus], pyMuon[myMuPlus], pzMuon[myMuPlus], energyMuon[myMuPlus]);
    m_mll[mm]             = (*(m_p4LeptonMinus[mm]) + *(m_p4LeptonPlus[mm])).M();
    m_deltaPhi[mm]        = fabs(180./TMath::Pi() * m_p4LeptonMinus[mm]->Vect().DeltaPhi(m_p4LeptonPlus[mm]->Vect()));
    m_deltaErre[mm]       = m_p4LeptonMinus[mm]->Vect().DeltaR(m_p4LeptonPlus[mm]->Vect());
    m_deltaEtaLeptons[mm] = etaEle[myMuMinus]-etaEle[myMuPlus];
    m_dilepPt[mm].SetXYZ( m_p4LeptonMinus[mm]->Vect().X()+m_p4LeptonPlus[mm]->Vect().X(),m_p4LeptonMinus[mm]->Vect().Y()+m_p4LeptonPlus[mm]->Vect().Y(),0.0 );
    m_transvMass[mm]      = CalcGammaMRstar(*m_p4LeptonMinus[mm],*m_p4LeptonPlus[mm]);
    m_metOptll[mm]        = m_theMET / m_dilepPt[mm].Pt();
    m_mT2[mm]             = 0.;
    m_projectedMet[mm]    = GetProjectedMet(m_p4LeptonMinus[mm]->Vect(),m_p4LeptonPlus[mm]->Vect());
  }

}

void HiggsMLSelection::setKinematicsEMME(int myEle, int myPosi, int myMuPlus, int myMuMinus) {

  if ( myEle > -1 && myMuPlus > -1 && ( myPosi<0 || myMuMinus<0 ) ) {  // only 1 pair reconstructed                                    

    float ptE = GetPt(pxEle[myEle],pyEle[myEle]);
    float ptM = GetPt(pxMuon[myMuPlus],pyMuon[myMuPlus]);
    
    if (ptE > ptM) {
      eleCands[em].push_back(myEle);
      muCands[em].push_back(myMuPlus);
      m_p4LeptonMinus[em] -> SetXYZT(pxEle[myEle],pyEle[myEle],pzEle[myEle],energyEle[myEle]);
      m_p4LeptonPlus[em]  -> SetXYZT(pxMuon[myMuPlus],pyMuon[myMuPlus],pzMuon[myMuPlus],energyMuon[myMuPlus]);
      hardestLeptonPt[em] = ptE;
      slowestLeptonPt[em] = ptM;
      m_mll[em]           = (*(m_p4LeptonMinus[em]) + *(m_p4LeptonPlus[em])).M();
      m_deltaPhi[em]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[em]->Vect().DeltaPhi(m_p4LeptonPlus[em]->Vect()));
      m_deltaErre[em]     = m_p4LeptonMinus[em]->Vect().DeltaR(m_p4LeptonPlus[em]->Vect());
      m_deltaEtaLeptons[em] = etaEle[myEle]-etaEle[myMuPlus];
      m_dilepPt[em].SetXYZ( m_p4LeptonMinus[em]->Vect().X()+m_p4LeptonPlus[em]->Vect().X(),m_p4LeptonMinus[em]->Vect().Y()+m_p4LeptonPlus[em]->Vect().Y(),0.0 );
      m_transvMass[em]    = CalcGammaMRstar(*m_p4LeptonMinus[em],*m_p4LeptonPlus[em]);
      m_metOptll[em]      = m_theMET / m_dilepPt[em].Pt();
      m_mT2[em]           = 0.;
      m_projectedMet[em]  = GetProjectedMet(m_p4LeptonMinus[em]->Vect(),m_p4LeptonPlus[em]->Vect());

    } else {
      eleCands[me].push_back(myEle);
      muCands[me].push_back(myMuPlus);
      m_p4LeptonMinus[me] -> SetXYZT(pxEle[myEle],pyEle[myEle],pzEle[myEle],energyEle[myEle]);
      m_p4LeptonPlus[me]  -> SetXYZT(pxMuon[myMuPlus],pyMuon[myMuPlus],pzMuon[myMuPlus],energyMuon[myMuPlus]);
      hardestLeptonPt[me] = ptM;
      slowestLeptonPt[me] = ptE;
      m_mll[me]           = (*(m_p4LeptonMinus[me]) + *(m_p4LeptonPlus[me])).M();
      m_deltaPhi[me]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[me]->Vect().DeltaPhi(m_p4LeptonPlus[me]->Vect()));
      m_deltaErre[me]     = m_p4LeptonMinus[me]->Vect().DeltaR(m_p4LeptonPlus[me]->Vect());
      m_deltaEtaLeptons[me] = etaEle[myEle]-etaEle[myMuPlus];
      m_dilepPt[me].SetXYZ( m_p4LeptonMinus[me]->Vect().X()+m_p4LeptonPlus[me]->Vect().X(),m_p4LeptonMinus[me]->Vect().Y()+m_p4LeptonPlus[me]->Vect().Y(),0.0 );
      m_transvMass[me]    = CalcGammaMRstar(*m_p4LeptonMinus[me],*m_p4LeptonPlus[me]);
      m_metOptll[me]      = m_theMET / m_dilepPt[me].Pt();
      m_mT2[me]           = 0.;
      m_projectedMet[me]  = GetProjectedMet(m_p4LeptonMinus[me]->Vect(),m_p4LeptonPlus[me]->Vect());
    }
  }
  
  if ( myPosi > -1 && myMuMinus > -1 && (myEle<0 || myMuPlus<0 )) {   // only 1 pair reconstructed                                     

    float ptE = GetPt(pxEle[myPosi],pyEle[myPosi]);
    float ptM = GetPt(pxMuon[myMuMinus],pyMuon[myMuMinus]);

    if (ptE > ptM) {
      eleCands[em].push_back(myPosi);
      muCands[em].push_back(myMuMinus);
      m_p4LeptonMinus[em] -> SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
      m_p4LeptonPlus[em]  -> SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
      hardestLeptonPt[em] = ptE;
      slowestLeptonPt[em] = ptM;
      m_mll[em]           = (*(m_p4LeptonMinus[em]) + *(m_p4LeptonPlus[em])).M();
      m_deltaPhi[em]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[em]->Vect().DeltaPhi(m_p4LeptonPlus[em]->Vect()));
      m_deltaErre[em]     = m_p4LeptonMinus[em]->Vect().DeltaR(m_p4LeptonPlus[em]->Vect());
      m_deltaEtaLeptons[em] = etaEle[myMuMinus]-etaEle[myPosi];
      m_dilepPt[em].SetXYZ( m_p4LeptonMinus[em]->Vect().X()+m_p4LeptonPlus[em]->Vect().X(),m_p4LeptonMinus[em]->Vect().Y()+m_p4LeptonPlus[em]->Vect().Y(),0.0 );
      m_transvMass[em]    = CalcGammaMRstar(*m_p4LeptonMinus[em],*m_p4LeptonPlus[em]);
      m_metOptll[em]      = m_theMET / m_dilepPt[em].Pt();
      m_mT2[em]           = 0.;
      m_projectedMet[em]  = GetProjectedMet(m_p4LeptonMinus[em]->Vect(),m_p4LeptonPlus[em]->Vect());

    } else {
      eleCands[me].push_back(myPosi);
      muCands[me].push_back(myMuMinus);
      m_p4LeptonMinus[me] -> SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
      m_p4LeptonPlus[me]  -> SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
      hardestLeptonPt[me] = ptM;
      slowestLeptonPt[me] = ptE;
      m_mll[me]           = (*(m_p4LeptonMinus[me]) + *(m_p4LeptonPlus[me])).M();
      m_deltaPhi[me]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[me]->Vect().DeltaPhi(m_p4LeptonPlus[me]->Vect()));
      m_deltaErre[me]     = m_p4LeptonMinus[me]->Vect().DeltaR(m_p4LeptonPlus[me]->Vect());
      m_deltaEtaLeptons[me] = etaEle[myMuMinus]-etaEle[myPosi];
      m_dilepPt[me].SetXYZ( m_p4LeptonMinus[me]->Vect().X()+m_p4LeptonPlus[me]->Vect().X(),m_p4LeptonMinus[me]->Vect().Y()+m_p4LeptonPlus[me]->Vect().Y(),0.0 );
      m_transvMass[me]    = CalcGammaMRstar(*m_p4LeptonMinus[me],*m_p4LeptonPlus[me]);
      m_metOptll[me]      = m_theMET / m_dilepPt[me].Pt();
      m_mT2[me]           = 0.;
      m_projectedMet[me]  = GetProjectedMet(m_p4LeptonMinus[me]->Vect(),m_p4LeptonPlus[me]->Vect());
    }
  }
  
  // if two pairs are built we choose the one with highest di-lepton pt                                                                
  if ( myPosi>-1 && myEle>-1 &&  myMuPlus>-1 && myMuMinus>-1) {   // 2 pairs reconstructed                                             
    
    TLorentzVector m_posi, m_ele, m_mum, m_mup;
    m_posi.SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
    m_ele.SetXYZT(pxEle[myEle], pyEle[myEle], pzEle[myEle], energyEle[myEle]);
    m_mum.SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
    m_mup.SetXYZT(pxMuon[myMuPlus], pyMuon[myMuPlus], pzMuon[myMuPlus], energyMuon[myMuPlus]);
    TVector3 m_ep_mm( m_posi.Vect().X()+m_mum.Vect().X(), m_posi.Vect().Y()+m_mum.Vect().Y(), 0.0 );
    TVector3 m_em_mp( m_ele.Vect().X()+m_mup.Vect().X(),  m_ele.Vect().Y()+m_mup.Vect().Y(),  0.0 );
    float mod_ep_mm = m_ep_mm.Mag();
    float mod_em_mp = m_em_mp.Mag();
    
    if (mod_ep_mm>mod_em_mp) {

      float ptE = GetPt(pxEle[myPosi],pyEle[myPosi]);
      float ptM = GetPt(pxMuon[myMuMinus],pyMuon[myMuMinus]);
      
      if (ptE > ptM) {
        eleCands[em].push_back(myPosi);
        muCands[em].push_back(myMuMinus);
        m_p4LeptonMinus[em] -> SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
        m_p4LeptonPlus[em]  -> SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
        hardestLeptonPt[em] = ptE;
        slowestLeptonPt[em] = ptM;
        m_mll[em]           = (*(m_p4LeptonMinus[em]) + *(m_p4LeptonPlus[em])).M();
        m_deltaPhi[em]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[em]->Vect().DeltaPhi(m_p4LeptonPlus[em]->Vect()));
        m_deltaErre[em]     = m_p4LeptonMinus[em]->Vect().DeltaR(m_p4LeptonPlus[em]->Vect());
	m_deltaEtaLeptons[em] = etaEle[myMuMinus]-etaEle[myPosi];
        m_dilepPt[em].SetXYZ( m_p4LeptonMinus[em]->Vect().X()+m_p4LeptonPlus[em]->Vect().X(),m_p4LeptonMinus[em]->Vect().Y()+m_p4LeptonPlus[em]->Vect().Y(),0.0 );
        m_transvMass[em]    = CalcGammaMRstar(*m_p4LeptonMinus[em],*m_p4LeptonPlus[em]);
        m_metOptll[em]      = m_theMET / m_dilepPt[em].Pt();
        m_mT2[em]           = 0.;
        m_projectedMet[em]  = GetProjectedMet(m_p4LeptonMinus[em]->Vect(),m_p4LeptonPlus[em]->Vect());

      } else {
        eleCands[me].push_back(myPosi);
        muCands[me].push_back(myMuMinus);
	m_p4LeptonMinus[me] -> SetXYZT(pxMuon[myMuMinus],pyMuon[myMuMinus],pzMuon[myMuMinus],energyMuon[myMuMinus]);
        m_p4LeptonPlus[me]  -> SetXYZT(pxEle[myPosi],pyEle[myPosi],pzEle[myPosi],energyEle[myPosi]);
        hardestLeptonPt[me] = ptM;
        slowestLeptonPt[me] = ptE;
        m_mll[me]           = (*(m_p4LeptonMinus[me]) + *(m_p4LeptonPlus[me])).M();
        m_deltaPhi[me]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[me]->Vect().DeltaPhi(m_p4LeptonPlus[me]->Vect()));
        m_deltaErre[me]     = m_p4LeptonMinus[me]->Vect().DeltaR(m_p4LeptonPlus[me]->Vect());
	m_deltaEtaLeptons[me] = etaEle[myMuMinus]-etaEle[myPosi];
        m_dilepPt[me].SetXYZ( m_p4LeptonMinus[me]->Vect().X()+m_p4LeptonPlus[me]->Vect().X(),m_p4LeptonMinus[me]->Vect().Y()+m_p4LeptonPlus[me]->Vect().Y(),0.0 );
        m_transvMass[me]    = CalcGammaMRstar(*m_p4LeptonMinus[me],*m_p4LeptonPlus[me]);
        m_metOptll[me]      = m_theMET / m_dilepPt[me].Pt();
        m_mT2[me]           = 0.;
        m_projectedMet[me]  = GetProjectedMet(m_p4LeptonMinus[me]->Vect(),m_p4LeptonPlus[me]->Vect());
      }

    } else {
      float ptE = GetPt(pxEle[myEle],pyEle[myEle]);
      float ptM = GetPt(pxMuon[myMuPlus],pyMuon[myMuPlus]);
      
      if (ptE > ptM) {
        eleCands[em].push_back(myEle);
        muCands[em].push_back(myMuPlus);
        m_p4LeptonMinus[em] -> SetXYZT(pxEle[myEle],pyEle[myEle],pzEle[myEle],energyEle[myEle]);
        m_p4LeptonPlus[em]  -> SetXYZT(pxMuon[myMuPlus],pyMuon[myMuPlus],pzMuon[myMuPlus],energyMuon[myMuPlus]);
        hardestLeptonPt[em] = ptE;
        slowestLeptonPt[em] = ptM;
        m_mll[em]           = (*(m_p4LeptonMinus[em]) + *(m_p4LeptonPlus[em])).M();
        m_deltaPhi[em]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[em]->Vect().DeltaPhi(m_p4LeptonPlus[em]->Vect()));
        m_deltaErre[em]     = m_p4LeptonMinus[em]->Vect().DeltaR(m_p4LeptonPlus[em]->Vect());
	m_deltaEtaLeptons[em] = etaEle[myEle]-etaEle[myMuPlus];
        m_dilepPt[em].SetXYZ( m_p4LeptonMinus[em]->Vect().X()+m_p4LeptonPlus[em]->Vect().X(),m_p4LeptonMinus[em]->Vect().Y()+m_p4LeptonPlus[em]->Vect().Y(),0.0 );
        m_transvMass[em]    = CalcGammaMRstar(*m_p4LeptonMinus[em],*m_p4LeptonPlus[em]);
        m_metOptll[em]      = m_theMET / m_dilepPt[em].Pt();
        m_mT2[em]           = 0.;
        m_projectedMet[em]  = GetProjectedMet(m_p4LeptonMinus[em]->Vect(),m_p4LeptonPlus[em]->Vect());

      } else {
        eleCands[me].push_back(myEle);
        muCands[me].push_back(myMuPlus);
        m_p4LeptonMinus[me] -> SetXYZT(pxEle[myEle],pyEle[myEle],pzEle[myEle],energyEle[myEle]);
        m_p4LeptonPlus[me]  -> SetXYZT(pxMuon[myMuPlus],pyMuon[myMuPlus],pzMuon[myMuPlus],energyMuon[myMuPlus]);
        hardestLeptonPt[me] = ptM;
        slowestLeptonPt[me] = ptE;
        m_mll[me]           = (*(m_p4LeptonMinus[me]) + *(m_p4LeptonPlus[me])).M();
        m_deltaPhi[me]      = fabs(180./TMath::Pi() * m_p4LeptonMinus[me]->Vect().DeltaPhi(m_p4LeptonPlus[me]->Vect()));
        m_deltaErre[me]     = m_p4LeptonMinus[me]->Vect().DeltaR(m_p4LeptonPlus[me]->Vect());
	m_deltaEtaLeptons[me] = etaEle[myEle]-etaEle[myMuPlus];
        m_dilepPt[me].SetXYZ( m_p4LeptonMinus[me]->Vect().X()+m_p4LeptonPlus[me]->Vect().X(),m_p4LeptonMinus[me]->Vect().Y()+m_p4LeptonPlus[me]->Vect().Y(),0.0 );
        m_transvMass[me]    = CalcGammaMRstar(*m_p4LeptonMinus[me],*m_p4LeptonPlus[me]);
        m_metOptll[me]      = m_theMET / m_dilepPt[me].Pt();
        m_mT2[me]           = 0.;
        m_projectedMet[me]  = GetProjectedMet(m_p4LeptonMinus[me]->Vect(),m_p4LeptonPlus[me]->Vect());
      }
    }
  }
  
}

void HiggsMLSelection::resetKinematicsStart() {

  theElectron  = -1;
  thePositron  = -1;
  theMuonMinus = -1;
  theMuonPlus  = -1;

  thePreElectron  = -1;
  thePrePositron  = -1;
  thePreMuonMinus = -1;
  thePreMuonPlus  = -1;
}

void HiggsMLSelection::resetKinematics() {

  for(int theChannel=0; theChannel<4; theChannel++) {
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



void HiggsMLSelection::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {
  
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
  //  thisCutBasedID->SetLikelihood( likelihoodRatio(eleIndex,*LH) );
  thisCutBasedID->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  thisCutBasedID->SetNBrem( nbremsEle[eleIndex] );
  thisCutBasedID->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                
  thisCutBasedID->SetTrkIsolation ( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                        
  thisCutBasedID->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );         
  float iso = 0.0;
  if ( anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB) ) iso = dr03TkSumPtEle[eleIndex] + max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else iso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  thisCutBasedID->SetCombinedIsolation( (iso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt );
  thisCutBasedID->SetCombinedPFIsolation( (pfGenericChargedIsoEle[eleIndex] + pfGenericNeutralIsoEle[eleIndex] + pfGenericPhotonIsoEle[eleIndex]) / pt );
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

void HiggsMLSelection::isMuonID(int muonIndex, bool *muonIdOutput) {

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


int HiggsMLSelection::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  m_goodJets.clear();
  float ETMax=0.;

  theLeadingJet[theChannel]=-1;   

  TString JESUncertainty(_selectionEE->getStringParameter("JESUncertainty"));

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    if(_selectionEE->getSwitch("etaJetAcc") && !_selectionEE->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) pt = (GetJESCorrected(p4Jet,JESUncertainty.Data())).Pt();

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
        H_deltaRcorr -> Fill(deltaR);
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
        H_deltaRcorr -> Fill(deltaR);
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


int HiggsMLSelection::numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  int num=0;

  TString JESUncertainty(_selectionEE->getStringParameter("JESUncertainty"));

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    float uncorrEt = uncorrEnergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrEnergyAK5PFPUcorrJet[j]);
    TVector3 p3Jet = p4Jet.Vect();

    TLorentzVector p4JESJet(p3Jet, uncorrEnergyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) uncorrEt = (GetJESCorrected(p4JESJet,JESUncertainty.Data())).Pt();

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
        H_deltaRcorr -> Fill(deltaR);
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
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEE->getSwitch("jetConeWidth") && _selectionEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;
    
    num++;
  }
  
  return num;
}

float HiggsMLSelection::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove ) {

  TString JESUncertainty(_selectionEE->getStringParameter("JESUncertainty"));

  float output=-999;
  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    // no threshold is applied here on pt. Not affected by JES uncertainties
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) pt = (GetJESCorrected(p4Jet,JESUncertainty.Data())).Pt();

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
        H_deltaRcorr -> Fill(deltaR);
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
        H_deltaRcorr -> Fill(deltaR);
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

float HiggsMLSelection::deltaPhiLLJet(int ichan) {   
  
  int myLeadingJet = theLeadingJet[ichan];

  if(myLeadingJet > -1 && m_dilepPt[ichan].Pt()>0) {
    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);    
    return fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));                           
  } else return -999.;
}

int HiggsMLSelection::numSoftMuons(std::vector<int> muonToRemove) {

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

int HiggsMLSelection::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

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
      if(pt>=15) isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
      if(pt<15)  isEleID(i,&theId,&theIso,&theConvRej,&EgammaCutBasedIDLow);
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
    if(_selectionMM->getSwitch("muGlobalIso") && !_selectionMM->passCut("muGlobalIso",isoSumRel)) continue;

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = PVzPV[m_closestPV] - trackVzTrack[track];  
    if(_selectionEE->getSwitch("muonIP") && !_selectionEE->passCut("muonIP",dxy)) continue;
    if(_selectionEE->getSwitch("muonDz") && !_selectionEE->passCut("muonDz",dz))  continue;  

    numMu++;
  }
  
  return numEle + numMu;
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
    myTrackerIso[i]  = dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
    myHcalIso[i]     = dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
    myEcalJIso[i]    = dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
    myEcalGTIso[i]   = 1000.;
    float combinedIso = 0.0;
    if ( fabs(myEta[i])<1.476 ) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
    else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
    myCombinedIso[i] = ( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / myPt[i] );
    myCharge[i] = chargeEle[eleIndex];
    int gsf = gsfTrackIndexEle[eleIndex];
    myMissHits[i] = expInnerLayersGsfTrack[gsf];
    myDist[i] = convDistEle[eleIndex];
    myDcot[i] = convDcotEle[eleIndex];
    //    myLh[i] = likelihoodRatio(eleIndex,*LH);
    myLh[i] = eleIdLikelihoodEle[eleIndex];

    // match with MC truth
    myMatched[i] = 999;
    if ( !_selectionEE->getSwitch("isData") ) { 
      int matchedReco = 0;
      TVector3 pReco(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
      for (int ii=0; ii<25; ii++) {  // chiara
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

int HiggsMLSelection::getPV() {
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
bool HiggsMLSelection::reloadTriggerMask()
{
  std::vector<int> triggerMask;

  // load the triggers required for EE
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersEE.begin();fIter!=requiredTriggersEE.end();++fIter)
    {   
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
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

  // load the triggers required for MM
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersMM.begin();fIter!=requiredTriggersMM.end();++fIter)
    {   
      //      std::cout << "For MM required: " << *fIter << std::endl;
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(*fIter) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersMM = triggerMask;

  // load the triggers NOT required for MM
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersMM.begin();fIter!=notRequiredTriggersMM.end();++fIter)
    {   
      //      std::cout << "For MM not required: " << *fIter << std::endl;
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(*fIter) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_notRequiredTriggersMM = triggerMask;

  // load the triggers required for EM
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersEM.begin();fIter!=requiredTriggersEM.end();++fIter)
    {   
      //      std::cout << "For EM required: " << *fIter << std::endl;
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(*fIter) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersEM = triggerMask;

  // load the triggers NOT required for EM
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersEM.begin();fIter!=notRequiredTriggersEM.end();++fIter)
    {   
      //      std::cout << "For EM not required: " << *fIter << std::endl;
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(*fIter) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_notRequiredTriggersEM = triggerMask;

}

bool HiggsMLSelection::hasPassedHLT(int channel) {
  Utils anaUtils;
  if(channel==ee) return anaUtils.getTriggersOR(m_requiredTriggersEE, firedTrg);
  else if(channel==mm) {
    bool required = anaUtils.getTriggersOR(m_requiredTriggersMM, firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersMM, firedTrg);
    return (required && !notRequired);
  } else if(channel==em) {
    bool required = anaUtils.getTriggersOR(m_requiredTriggersEM, firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersEM, firedTrg);
    return (required && !notRequired);
  }
  return true;
}

void HiggsMLSelection::setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) requiredTriggersEE=reqTriggers;
  else if(channel==mm) requiredTriggersMM=reqTriggers;
  else if(channel==em) requiredTriggersEM=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

void HiggsMLSelection::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) notRequiredTriggersEE=reqTriggers;
  else if(channel==mm) notRequiredTriggersMM=reqTriggers;
  else if(channel==em) notRequiredTriggersEM=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

