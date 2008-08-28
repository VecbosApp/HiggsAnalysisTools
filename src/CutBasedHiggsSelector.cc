#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include <iostream>
#include <math.h>

CutBasedHiggsSelector::CutBasedHiggsSelector() {
  m_finalLeptons = false;
  m_jetVeto = false;
  m_preDeltaPhi = false;
  m_processID = -1;
}

CutBasedHiggsSelector::~CutBasedHiggsSelector() {}

void CutBasedHiggsSelector::Configure(const char *fileCuts, const char* fileSwitches) {

  _selection = new Selection(std::string(fileCuts),std::string(fileSwitches));

  // tehse cuts are applied in the HiggsSelection class, but are configured here
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetLowAcc");
  _selection->addCut("etJetHighAcc");
  _selection->addCut("alphaJet");

  _selection->addSwitch("classDepEleId");
  _selection->addSwitch("jetVeto");
  _selection->addCut("hardLeptonThreshold"); 
  _selection->addCut("slowLeptonThreshold"); 
  _selection->addCut("dileptonInvMassMin");  
  _selection->addCut("trackerPtSum");
  _selection->addCut("hcalPtSum");
  _selection->addCut("ecalPtSum");
  _selection->addCut("MET");
  _selection->addCut("maxPtLepton");
  _selection->addCut("minPtLepton");
  _selection->addCut("dileptonInvMassMax");
  _selection->addCut("detaLeptons");
  _selection->addCut("deltaPhi");

  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle("FULL SELECTION EVENT COUNTER");
  globalCounter->AddVar("preselected");
  globalCounter->AddVar("hardLeptonThreshold");
  globalCounter->AddVar("slowLeptonThreshold");
  globalCounter->AddVar("dileptonInvMassMin");
  globalCounter->AddVar("classDepEleId");
  globalCounter->AddVar("trackerIso");
  globalCounter->AddVar("hcalIso");
  globalCounter->AddVar("ecalIso");
  globalCounter->AddVar("jetVeto");
  globalCounter->AddVar("MET");
  globalCounter->AddVar("deltaPhi");
  globalCounter->AddVar("maxPtLepton");
  globalCounter->AddVar("minPtLepton");
  globalCounter->AddVar("dileptonInvMassMax");
  globalCounter->AddVar("detaLeptons");
  globalCounter->AddVar("final");
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
      processCounter->AddVar("jetVeto");
      processCounter->AddVar("MET");
      processCounter->AddVar("deltaPhi");
      processCounter->AddVar("maxPtLepton");
      processCounter->AddVar("minPtLepton");
      processCounter->AddVar("dileptonInvMassMax");
      processCounter->AddVar("detaLeptons");
      processCounter->AddVar("final");
      
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

  if ((_selection->getSwitch("trackerPtSum")) && 
      (!_selection->passCut("trackerPtSum",m_eleTkPtSum) || !_selection->passCut("trackerPtSum",m_posTkPtSum))) return false; 
  theCounter->IncrVar("trackerIso",m_weight);

  if ((_selection->getSwitch("hcalPtSum")) && 
      (!_selection->passCut("hcalPtSum",m_eleHcalPtSum) || !_selection->passCut("hcalPtSum",m_posHcalPtSum))) return false; 
  theCounter->IncrVar("hcalIso",m_weight);

  if ((_selection->getSwitch("ecalPtSum")) && 
      (!_selection->passCut("ecalPtSum",m_eleEcalPtSum) || !_selection->passCut("ecalPtSum",m_posEcalPtSum))) return false; 
  theCounter->IncrVar("ecalIso",m_weight);

  m_finalLeptons = true;

  if(_selection->getSwitch("jetVeto") && !m_passedJetVeto) return false; 
  theCounter->IncrVar("jetVeto",m_weight);

  m_jetVeto = true;

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

  return true;
}


void CutBasedHiggsSelector::diplayEfficiencies() {

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
      theCounter->Draw("jetVeto","ecalIso");
      theCounter->Draw("MET","jetVeto");
      theCounter->Draw("maxPtLepton","MET");   
      theCounter->Draw("minPtLepton","maxPtLepton");
      theCounter->Draw("dileptonInvMassMax","minPtLepton");
      theCounter->Draw("detaLeptons","dileptonInvMassMax");
      theCounter->Draw("deltaPhi","detaLeptons");
      theCounter->Draw("deltaPhi","preselected");

    }

  }

  else {
    
    globalCounter->Draw();
    globalCounter->Draw("hardLeptonThreshold","preselected");
    globalCounter->Draw("slowLeptonThreshold","hardLeptonThreshold");
    globalCounter->Draw("dileptonInvMassMin","slowLeptonThreshold");
    globalCounter->Draw("classDepEleId","dileptonInvMassMin");
    globalCounter->Draw("trackerIso","classDepEleId");
    globalCounter->Draw("hcalIso","trackerIso");
    globalCounter->Draw("ecalIso","hcalIso");
    globalCounter->Draw("jetVeto","ecalIso");
    globalCounter->Draw("MET","jetVeto");
    globalCounter->Draw("maxPtLepton","MET");
    globalCounter->Draw("minPtLepton","maxPtLepton");
    globalCounter->Draw("dileptonInvMassMax","minPtLepton");
    globalCounter->Draw("detaLeptons","dileptonInvMassMax");
    globalCounter->Draw("deltaPhi","detaLeptons");
    globalCounter->Draw("deltaPhi","preselected");
    
  }

}
