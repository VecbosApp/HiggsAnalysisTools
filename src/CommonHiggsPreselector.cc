#include "HiggsAnalysisTools/include/CommonHiggsPreselector.hh"
#include <iostream>
#include <math.h>

CommonHiggsPreselector::CommonHiggsPreselector() {}

CommonHiggsPreselector::~CommonHiggsPreselector() {}

void CommonHiggsPreselector::Configure(const char *configDir) {
  
  std::string fileCuts     = std::string(configDir) + "2e2nuCuts.txt";
  std::string fileSwitches = std::string(configDir) + "2e2nuSwitches.txt";

  _selection = new Selection(fileCuts,fileSwitches);
  
  _selection->addSwitch("MCtruth");
  _selection->addSwitch("trigger");
  _selection->addSwitch("preselection");
  _selection->addSwitch("leptonAcceptance");
  _selection->addCut("nRecoLeptons");
  _selection->addCut("nRecoLeptonsMix");
  _selection->addCut("hardLeptonThreshold");
  _selection->addCut("slowLeptonThreshold");
  _selection->addCut("METpreselection");
  _selection->addCut("dileptonInvMassMin");

  _selection->summary();

  presCounter.SetTitle("PRESELECTION EVENT COUNTER");
  presCounter.AddVar("event");
  presCounter.AddVar("MCtruth");
  presCounter.AddVar("trigger");
  presCounter.AddVar("nRecoLeptons");
  presCounter.AddVar("twoGoodRec");
  presCounter.AddVar("hardLeptonThreshold");
  presCounter.AddVar("slowLeptonThreshold");
  presCounter.AddVar("METpreselection");
  presCounter.AddVar("dileptonInvMassMin");
  presCounter.AddVar("finalOURPreselection");
  presCounter.AddVar("preselection");
}



bool CommonHiggsPreselector::output() {
  
  float m_weight = m_kFactor;
  presCounter.IncrVar("event",m_weight);

  // common preselection cut: this is not really applied, should be reproduced by hand
  if( _selection->getSwitch("preselection") && m_evtPresel ) { 
    presCounter.IncrVar("preselection",m_weight);
  }
  
  // MC truth
  if(_selection->getSwitch("MCtruth") && !m_foundMcTree ) return false;
  presCounter.IncrVar("MCtruth",m_weight);
  
  // HLT
  if(_selection->getSwitch("trigger") && !m_passedHLT ) return false;
  presCounter.IncrVar("trigger",m_weight); 
  
  // nEle, nMuon
  if ( _selection->getSwitch("nRecoLeptons") && 
       ( !_selection->passCut("nRecoLeptons", m_nEle) && !_selection->passCut("nRecoLeptons", m_nMuon) ) &&
       ( !_selection->passCut("nRecoLeptonsMix", m_nEle) || !_selection->passCut("nRecoLeptonsMix", m_nMuon) ) 
       ) return false;
  presCounter.IncrVar("nRecoLeptons",m_weight);
  
  // leptons in the acceptance
  if( _selection->getSwitch("leptonAcceptance") && !m_isEE && !m_isEM && !m_isMM ) return false;
  presCounter.IncrVar("twoGoodRec",m_weight);
  
  // high pt lepton
  if ( _selection->getSwitch("hardLeptonThreshold") &&
       (!_selection->passCut("hardLeptonThreshold", m_highElePt) && 
	!_selection->passCut("hardLeptonThreshold", m_highMuonPt) )
       ) return false;
  presCounter.IncrVar("hardLeptonThreshold",m_weight);
  
  // low pt lepton
  if ( _selection->getSwitch("slowLeptonThreshold") &&
       (!_selection->passCut("slowLeptonThreshold", m_lowElePt) && 
	!_selection->passCut("slowLeptonThreshold", m_lowMuonPt) )
       ) return false;
  presCounter.IncrVar("slowLeptonThreshold",m_weight);
  
  // met cut
  if ( _selection->getSwitch("METpreselection") && 
       !_selection->passCut("METpreselection", m_met) ) return false;
  presCounter.IncrVar("METpreselection",m_weight);
  
  // m(ll) cut
  if ( _selection->getSwitch("dileptonInvMassMin")) {
    if ( ! ( (m_isEE && _selection->passCut("dileptonInvMassMin", m_mllEE)) ||
	     (m_isMM && _selection->passCut("dileptonInvMassMin", m_mllMM)) ||
	     (m_isEM && _selection->passCut("dileptonInvMassMin", m_mllEM)) ) 
	 ) return false;
  }
  presCounter.IncrVar("dileptonInvMassMin",m_weight);

  presCounter.IncrVar("finalOURPreselection",m_weight);
    
  return true;
}


void CommonHiggsPreselector::diplayEfficiencies() {

  presCounter.Draw();
  presCounter.Draw("MCtruth","event");
  presCounter.Draw("trigger","MCtruth");
  presCounter.Draw("nRecoLeptons","trigger");
  presCounter.Draw("twoGoodRec","nRecoLeptons");
  presCounter.Draw("hardLeptonThreshold","twoGoodRec");
  presCounter.Draw("slowLeptonThreshold","hardLeptonThreshold");
  presCounter.Draw("METpreselection","slowLeptonThreshold");
  presCounter.Draw("dileptonInvMassMin","METpreselection");
  presCounter.Draw("finalOURPreselection","MCtruth");
  presCounter.Draw("preselection","MCtruth");
}
