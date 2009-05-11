#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include <iostream>
#include <math.h>

CutBasedHiggsSelector::CutBasedHiggsSelector() {
  m_finalLeptons = false;
  m_jetVeto = false;
  m_preDeltaPhi = false;
  m_processID = -1;
}

CutBasedHiggsSelector::CutBasedHiggsSelector( const CutBasedHiggsSelector& selector ) {
  m_weight = selector.m_weight;
  m_highPt = selector.m_highPt;
  m_isElectronId = selector.m_isElectronId;
  m_isPositronId = selector.m_isPositronId;
  m_invMass = selector.m_invMass;
  m_eleHardTkPtSum   = selector.m_eleHardTkPtSum;
  m_eleHardHcalPtSum = selector.m_eleHardHcalPtSum;
  m_eleHardEcalPtSum = selector.m_eleHardEcalPtSum;
  m_eleHardGlobalSum = selector.m_eleHardGlobalSum;
  m_eleSlowTkPtSum   = selector.m_eleSlowTkPtSum;
  m_eleSlowHcalPtSum = selector.m_eleSlowHcalPtSum;
  m_eleSlowEcalPtSum = selector.m_eleSlowEcalPtSum;
  m_eleSlowGlobalSum = selector.m_eleSlowGlobalSum;
  m_nJets            = selector.m_nJets;
  m_met = selector.m_met;
  m_deltaPhi = selector.m_deltaPhi;
  m_detaLeptons = selector.m_detaLeptons;
  m_maxPtElectron = selector.m_maxPtElectron;
  m_minPtElectron = selector.m_minPtElectron;
  m_processID = selector.m_processID;
  *_selection = *selector._selection;
  *globalCounter = *selector.globalCounter;
  m_finalLeptons = selector.m_finalLeptons;
  m_jetVeto = selector.m_jetVeto;
  m_preDeltaPhi = selector.m_preDeltaPhi;
  multiProcessCounter = selector.multiProcessCounter;

}

CutBasedHiggsSelector::~CutBasedHiggsSelector() {}

void CutBasedHiggsSelector::Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle) {

  _selection = new Selection(std::string(fileCuts),std::string(fileSwitches));

  // tehse cuts are applied in the HiggsSelection class, but are configured here
  _selection->addSwitch("classDepEleId");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("hardLeptonThreshold"); 
  _selection->addCut("slowLeptonThreshold"); 
  _selection->addCut("dileptonInvMassMin");  
  _selection->addCut("trackerPtSumAss");
  _selection->addCut("trackerPtSumRel");
  _selection->addCut("hcalPtSumAss");
  _selection->addCut("hcalPtSumRel");
  _selection->addCut("ecalPtSumAss");
  _selection->addCut("ecalPtSumRel");
  _selection->addCut("globalSum");
  _selection->addCut("MET");
  _selection->addCut("maxPtLepton");
  _selection->addCut("minPtLepton");
  _selection->addCut("dileptonInvMassMax");
  _selection->addCut("detaLeptons");
  _selection->addCut("deltaPhi");

  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle(theTitle);
  globalCounter->AddVar("preselected");
  globalCounter->AddVar("hardLeptonThreshold");
  globalCounter->AddVar("slowLeptonThreshold");
  globalCounter->AddVar("dileptonInvMassMin");
  globalCounter->AddVar("classDepEleId");
  globalCounter->AddVar("trackerIso");
  globalCounter->AddVar("hcalIso");
  globalCounter->AddVar("ecalIso");
  globalCounter->AddVar("globalIso");
  globalCounter->AddVar("MET");
  globalCounter->AddVar("maxPtLepton");
  globalCounter->AddVar("minPtLepton");
  globalCounter->AddVar("dileptonInvMassMax");
  globalCounter->AddVar("detaLeptons");
  globalCounter->AddVar("deltaPhi");
  globalCounter->AddVar("final");
  globalCounter->AddVar("zeroJets");
  globalCounter->AddVar("oneJet");
  globalCounter->AddVar("gt1Jets");
}

bool CutBasedHiggsSelector::output() {

  Counters *theCounter=0;

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter = multiProcessCounter.find(m_processID);

    if ( iter == multiProcessCounter.end() ) {
      
      std::cout << "First time I get process " << m_processID 
		<< ": adding a counter" << std::endl;

      char buffer[200];
      sprintf(buffer,"Event counter for process %d", m_processID);
      
      Counters *processCounter = new Counters();
      processCounter->SetTitle(buffer);
      processCounter->AddVar("preselected");
      processCounter->AddVar("hardLeptonThreshold");
      processCounter->AddVar("slowLeptonThreshold");
      processCounter->AddVar("dileptonInvMassMin");
      processCounter->AddVar("classDepEleId");
      processCounter->AddVar("trackerIso");
      processCounter->AddVar("hcalIso");
      processCounter->AddVar("ecalIso");
      processCounter->AddVar("globalIso");
      processCounter->AddVar("MET");
      processCounter->AddVar("maxPtLepton");
      processCounter->AddVar("minPtLepton");
      processCounter->AddVar("dileptonInvMassMax");
      processCounter->AddVar("detaLeptons");
      processCounter->AddVar("deltaPhi");
      processCounter->AddVar("final");
      processCounter->AddVar("zeroJets");
      processCounter->AddVar("oneJet");
      processCounter->AddVar("gt1Jets");
      
      multiProcessCounter.insert( std::make_pair(m_processID,processCounter) );      
    }

    theCounter = multiProcessCounter[m_processID];

  }
  
  else theCounter = globalCounter;

  m_finalLeptons = false;
  m_jetVeto = false;
  m_preDeltaPhi = false;

  theCounter->IncrVar("preselected",m_weight);

  // here I repeat preselection cuts only for the type of leptons I want
  if (_selection->getSwitch("hardLeptonThreshold") && !_selection->passCut("hardLeptonThreshold", m_highPt)) return false;
  theCounter->IncrVar("hardLeptonThreshold",m_weight);
  
  if (_selection->getSwitch("slowLeptonThreshold") && !_selection->passCut("slowLeptonThreshold", m_lowPt)) return false;
  theCounter->IncrVar("slowLeptonThreshold",m_weight);

  if (_selection->getSwitch("dileptonInvMassMin") && !_selection->passCut("dileptonInvMassMin", m_invMass)) return false;
  theCounter->IncrVar("dileptonInvMassMin",m_weight);

  // real selections (after preselection step)
  if ((_selection->getSwitch("classDepEleId")) && (!m_isElectronId || !m_isPositronId)) return false; 
  theCounter->IncrVar("classDepEleId",m_weight);
  
  if ((_selection->getSwitch("trackerPtSumAss"))) {
    if (m_highPt>=25 && !_selection->passCut("trackerPtSumAss",m_eleHardTkPtSum)) return false; 
    if (m_lowPt>=25  && !_selection->passCut("trackerPtSumAss",m_eleSlowTkPtSum)) return false; 
  }
  if ((_selection->getSwitch("trackerPtSumRel"))) {
    if (m_highPt<25  && !_selection->passCut("trackerPtSumRel",m_eleHardTkPtSum)) return false; 
    if (m_lowPt<25   && !_selection->passCut("trackerPtSumRel",m_eleSlowTkPtSum)) return false; 
  }
  theCounter->IncrVar("trackerIso",m_weight);

  if ((_selection->getSwitch("hcalPtSumAss"))) {
    if (m_highPt>=25 && !_selection->passCut("hcalPtSumAss",m_eleHardHcalPtSum)) return false;
    if (m_lowPt>=25  && !_selection->passCut("hcalPtSumAss",m_eleSlowHcalPtSum)) return false; 
  }
  if ((_selection->getSwitch("hcalPtSumRel"))) {
    if (m_highPt<25  && !_selection->passCut("hcalPtSumRel",m_eleHardHcalPtSum)) return false;
    if (m_lowPt<25   && !_selection->passCut("hcalPtSumRel",m_eleSlowHcalPtSum)) return false; 
  }
  theCounter->IncrVar("hcalIso",m_weight);

  if ((_selection->getSwitch("ecalPtSumAss"))) {
    if (m_highPt>=25 && !_selection->passCut("ecalPtSumAss",m_eleHardEcalPtSum)) return false; 
    if (m_lowPt>=25  && !_selection->passCut("ecalPtSumAss",m_eleSlowEcalPtSum)) return false; 
  }
  if ((_selection->getSwitch("ecalPtSumRel"))) {
    if (m_highPt<25  && !_selection->passCut("ecalPtSumRel",m_eleHardEcalPtSum)) return false; 
    if (m_lowPt<25   && !_selection->passCut("ecalPtSumRel",m_eleSlowEcalPtSum)) return false; 
  }
  theCounter->IncrVar("ecalIso",m_weight);

  if ((_selection->getSwitch("globalSum"))) {
    if (m_highPt>=25 && !_selection->passCut("globalSum",m_eleHardGlobalSum)) return false; 
    if (m_lowPt>=25  && !_selection->passCut("globalSum",m_eleSlowGlobalSum)) return false; 
    if (m_highPt<25  && (m_eleHardGlobalSum > ((m_highPt-10.)/3.) ))          return false; 
    if (m_lowPt<25   && (m_eleSlowGlobalSum > ((m_lowPt-10.)/3.) ))           return false; 
  }
  theCounter->IncrVar("globalIso",m_weight);

  m_finalLeptons = true;

  if (m_nJets==0) m_jetVeto = true;

  if(_selection->getSwitch("MET") && !_selection->passCut("MET",m_met)) return false; 
  theCounter->IncrVar("MET",m_weight);

  if (_selection->getSwitch("maxPtLepton") && !_selection->passCut("maxPtLepton", m_highPt)) return false;
  theCounter->IncrVar("maxPtLepton",m_weight);

  if (_selection->getSwitch("minPtLepton") && !_selection->passCut("minPtLepton", m_lowPt)) return false;
  theCounter->IncrVar("minPtLepton",m_weight);

  if (_selection->getSwitch("dileptonInvMassMax") && !_selection->passCut("dileptonInvMassMax", m_invMass)) return false;
  theCounter->IncrVar("dileptonInvMassMax",m_weight);

  if (_selection->getSwitch("detaLeptons") && !_selection->passCut("detaLeptons", m_detaLeptons)) return false;
  theCounter->IncrVar("detaLeptons",m_weight);

  m_preDeltaPhi = true;

  if (_selection->getSwitch("deltaPhi") && !_selection->passCut("deltaPhi", m_deltaPhi)) return false;
  theCounter->IncrVar("deltaPhi",m_weight); 


  theCounter->IncrVar("final",m_weight);

  if (m_nJets==0) theCounter->IncrVar("zeroJets",m_weight);
  if (m_nJets==1) theCounter->IncrVar("oneJet",m_weight);
  if (m_nJets>1)  theCounter->IncrVar("gt1Jets",m_weight);

  return true;
}


void CutBasedHiggsSelector::diplayEfficiencies(std::string datasetName) {

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter;
    for( iter=multiProcessCounter.begin(); iter!=multiProcessCounter.end(); ++iter ) {

      Counters *theCounter = iter->second;

      theCounter->Draw();
      theCounter->Draw("hardLeptonThreshold","preselected");
      theCounter->Draw("slowLeptonThreshold","hardLeptonThreshold");
      theCounter->Draw("dileptonInvMassMin","slowLeptonThreshold");
      theCounter->Draw("classDepEleId","dileptonInvMassMin");
      theCounter->Draw("trackerIso","classDepEleId");
      theCounter->Draw("hcalIso","trackerIso");
      theCounter->Draw("ecalIso","hcalIso");
      theCounter->Draw("globalIso","ecalIso");
      theCounter->Draw("MET","globalIso");
      theCounter->Draw("maxPtLepton","MET");   
      theCounter->Draw("minPtLepton","maxPtLepton");
      theCounter->Draw("dileptonInvMassMax","minPtLepton");
      theCounter->Draw("detaLeptons","dileptonInvMassMax");
      theCounter->Draw("deltaPhi","detaLeptons");
      theCounter->Draw("final","preselected");
      theCounter->Draw("zeroJets", "final");
      theCounter->Draw("oneJet",   "final");
      theCounter->Draw("gt1Jets",  "final");
    }
  }

  else {

    char namefile[500];
    sprintf(namefile,"%s-Counters.root",datasetName.c_str());
    
    globalCounter->Draw();
    globalCounter->Draw("hardLeptonThreshold","preselected");
    globalCounter->Draw("slowLeptonThreshold","hardLeptonThreshold");
    globalCounter->Draw("dileptonInvMassMin","slowLeptonThreshold");
    globalCounter->Draw("classDepEleId","dileptonInvMassMin");
    globalCounter->Draw("trackerIso","classDepEleId");
    globalCounter->Draw("hcalIso","trackerIso");
    globalCounter->Draw("ecalIso","hcalIso");
    globalCounter->Draw("globalIso","ecalIso");
    globalCounter->Draw("MET","globalIso");
    globalCounter->Draw("maxPtLepton","MET");
    globalCounter->Draw("minPtLepton","maxPtLepton");
    globalCounter->Draw("dileptonInvMassMax","minPtLepton");
    globalCounter->Draw("detaLeptons","dileptonInvMassMax");
    globalCounter->Draw("deltaPhi","detaLeptons");
    globalCounter->Draw("final","preselected");
    globalCounter->Draw("zeroJets", "final");
    globalCounter->Draw("oneJet",   "final");
    globalCounter->Draw("gt1Jets",  "final");

    globalCounter->Save(namefile,"update");
  }

}
