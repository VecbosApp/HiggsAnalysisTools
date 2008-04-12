#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include <iostream>
#include <math.h>

CutBasedHiggsSelector::CutBasedHiggsSelector() {}

CutBasedHiggsSelector::~CutBasedHiggsSelector() {}

void CutBasedHiggsSelector::Configure(const char *configDir) {

  std::string fileCuts     = std::string(configDir) + "2e2nuCuts.txt";
  std::string fileSwitches = std::string(configDir) + "2e2nuSwitches.txt";

  _selection = new Selection(fileCuts,fileSwitches);

  _selection->addSwitch("classDepEleId");
  _selection->addSwitch("jetVeto");
  _selection->addCut("hardLeptonThreshold"); 
  _selection->addCut("slowLeptonThreshold"); 
  _selection->addCut("dileptonInvMassMin");  
  _selection->addCut("trackerPtSum");
  _selection->addCut("hcalPtSum");
  _selection->addCut("MET");
  _selection->addCut("deltaPhi");
  _selection->addCut("dileptonInvMassMax");
  _selection->addCut("detaLeptons");

  cout << "sommario cutbasedselection" << endl;
  _selection->summary();

  higgsSelCounter.SetTitle("FULL SELECTION EVENT COUNTER");
  higgsSelCounter.AddVar("preselected");
  higgsSelCounter.AddVar("hardLeptonThreshold");
  higgsSelCounter.AddVar("slowLeptonThreshold");
  higgsSelCounter.AddVar("dileptonInvMassMin");
  higgsSelCounter.AddVar("trackerIsol");
  higgsSelCounter.AddVar("hcalIsol");
  higgsSelCounter.AddVar("jetVeto");
  higgsSelCounter.AddVar("MET");
  higgsSelCounter.AddVar("deltaPhi");
  higgsSelCounter.AddVar("dileptonInvMassMax");
  higgsSelCounter.AddVar("detaLeptons");
  higgsSelCounter.AddVar("final");
}

bool CutBasedHiggsSelector::output() {

  higgsSelCounter.IncrVar("preselected",m_weight);

  // here I repeat preselection cuts only for the type of leptons I want
  if (_selection->getSwitch("hardLeptonThreshold") && !_selection->passCut("hardLeptonThreshold", m_highPt)) return false;
  higgsSelCounter.IncrVar("hardLeptonThreshold",m_weight);
  
  if (_selection->getSwitch("slowLeptonThreshold") && !_selection->passCut("slowLeptonThreshold", m_lowPt)) return false;
  higgsSelCounter.IncrVar("slowLeptonThreshold",m_weight);

  if (_selection->getSwitch("dileptonInvMassMin") && !_selection->passCut("dileptonInvMassMin", m_invMass)) return false;
  higgsSelCounter.IncrVar("dileptonInvMassMin",m_weight);

  // real selections (after preselection step)
  if ((_selection->getSwitch("classDepEleId")) && (!m_isElectronId || !m_isPositronId)) return false; 
  higgsSelCounter.IncrVar("eleID",m_weight);

  if ((_selection->getSwitch("trackerPtSum")) && 
      (!_selection->passCut("trackerPtSum",m_eleTkPtSum) || !_selection->passCut("trackerPtSum",m_posTkPtSum))) return false; 
  higgsSelCounter.IncrVar("trackerIso",m_weight);

  if ((_selection->getSwitch("hcalPtSum")) && 
      (!_selection->passCut("hcalPtSum",m_eleCaloPtSum) || !_selection->passCut("hcalPtSum",m_posCaloPtSum))) return false; 
  higgsSelCounter.IncrVar("caloIso",m_weight);

  if(_selection->getSwitch("jetVeto") && !m_passedJetVeto) return false; 
  higgsSelCounter.IncrVar("jetVeto",m_weight);

  if(_selection->getSwitch("MET") && !_selection->passCut("MET",m_met)) return false; 
  higgsSelCounter.IncrVar("MET",m_weight);

  if (_selection->getSwitch("deltaPhi") && !_selection->passCut("deltaPhi", m_deltaPhi)) return false;
  higgsSelCounter.IncrVar("deltaPhi",m_weight); 

  if (_selection->getSwitch("dileptonInvMassMax") && !_selection->passCut("dileptonInvMassMax", m_invMass)) return false;
  higgsSelCounter.IncrVar("dileptonInvMassMax",m_weight);

  if (_selection->getSwitch("detaLeptons") && !_selection->passCut("detaLeptons", m_detaLeptons)) return false;
  higgsSelCounter.IncrVar("detaLeptons",m_weight);

  higgsSelCounter.IncrVar("final",m_weight);

  return true;
}


void CutBasedHiggsSelector::diplayEfficiencies() {

  higgsSelCounter.Draw();
  higgsSelCounter.Draw("hardLeptonThreshold","preselected");
  higgsSelCounter.Draw("slowLeptonThreshold","hardLeptonThreshold");
  higgsSelCounter.Draw("dileptonInvMassMin","slowLeptonThreshold");
  higgsSelCounter.Draw("trackerIsol","dileptonInvMassMin");
  higgsSelCounter.Draw("hcalIsol","trackerIsol");
  higgsSelCounter.Draw("jetVeto","hcalIsol");
  higgsSelCounter.Draw("MET","jetVeto");
  higgsSelCounter.Draw("deltaPhi","MET");
  higgsSelCounter.Draw("dileptonInvMassMax","deltaPhi");
  higgsSelCounter.Draw("detaLeptons","dileptonInvMassMax");
  higgsSelCounter.Draw("detaLeptons","preselected");
}
