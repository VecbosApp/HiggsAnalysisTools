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
  m_isElectronIsol = selector.m_isElectronIsol;
  m_isPositronIsol = selector.m_isPositronIsol;
  m_isElectronConvRej = selector.m_isElectronConvRej;
  m_isPositronConvRej = selector.m_isPositronConvRej;
  m_invMass = selector.m_invMass;
  m_eleHardTkPtSum   = selector.m_eleHardTkPtSum;
  m_eleHardHcalPtSum = selector.m_eleHardHcalPtSum;
  m_eleHardEcalPtSum = selector.m_eleHardEcalPtSum;
  m_eleHardGlobalPtSum = selector.m_eleHardGlobalPtSum;
  m_eleSlowTkPtSum   = selector.m_eleSlowTkPtSum;
  m_eleSlowHcalPtSum = selector.m_eleSlowHcalPtSum;
  m_eleSlowEcalPtSum = selector.m_eleSlowEcalPtSum;
  m_eleSlowGlobalPtSum = selector.m_eleSlowGlobalPtSum;
  m_eleSlowD0 = selector.m_eleSlowD0;
  m_eleHardD0 = selector.m_eleHardD0;
  m_nJets            = selector.m_nJets;
  m_nUncorrJets      = selector.m_nUncorrJets;
  m_btagJets         = selector.m_btagJets;
  m_nSoftMuons       = selector.m_nSoftMuons;
  m_nExtraLeptons    = selector.m_nExtraLeptons;
  m_met = selector.m_met;
  m_projectedMet = selector.m_projectedMet;
  m_metOverPtLL = selector.m_metOverPtLL;
  m_deltaPhiLLJet = selector.m_deltaPhiLLJet;
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
  _selection->addCut("maxPtLepton");
  _selection->addCut("minPtLepton");
  _selection->addCut("leptonD0");
  _selection->addSwitch("eleIso");
  _selection->addCut("muTrackerIso");
  _selection->addCut("muHcalIso");
  _selection->addCut("muEcalIso");
  _selection->addCut("muGlobalIso");
  _selection->addSwitch("leptonId");
  _selection->addSwitch("convRej");
  _selection->addCut("looseMET");
  _selection->addCut("mll");
  _selection->addCut("mllZPeak");
  _selection->addCut("tightMET");
  _selection->addCut("projectedMET");
  _selection->addCut("metOverPtLL");
  _selection->addCut("deltaPhiLLJet");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("bTagVeto");
  _selection->addCut("nSoftMuons");
  _selection->addCut("nExtraLeptons");
  _selection->addCut("deltaPhi");

  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle(theTitle);
  globalCounter->AddVar("preselected");
  globalCounter->AddVar("maxPtLepton");
  globalCounter->AddVar("minPtLepton");
  globalCounter->AddVar("leptonD0");
  globalCounter->AddVar("eleIso");
  globalCounter->AddVar("muTrackerIso");
  globalCounter->AddVar("muHcalIso");
  globalCounter->AddVar("muEcalIso");
  globalCounter->AddVar("muGlobalIso");
  globalCounter->AddVar("leptonId");
  globalCounter->AddVar("convRej");
  globalCounter->AddVar("looseMET");
  globalCounter->AddVar("mll");
  globalCounter->AddVar("mllZPeak");
  globalCounter->AddVar("tightMETandPrMET");
  globalCounter->AddVar("metOverPtLL");
  globalCounter->AddVar("zeroJets");
  globalCounter->AddVar("bTagVeto");
  globalCounter->AddVar("nSoftMuons");
  globalCounter->AddVar("nExtraLeptons");
  globalCounter->AddVar("deltaPhi");
  globalCounter->AddVar("final");
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
      processCounter->AddVar("maxPtLepton");
      processCounter->AddVar("minPtLepton");
      processCounter->AddVar("leptonD0");
      processCounter->AddVar("eleIso");
      processCounter->AddVar("muTrackerIso");
      processCounter->AddVar("muHcalIso");
      processCounter->AddVar("muEcalIso");
      processCounter->AddVar("muGlobalIso");
      processCounter->AddVar("leptonId");
      processCounter->AddVar("convRej");
      processCounter->AddVar("looseMET");
      processCounter->AddVar("mll");
      processCounter->AddVar("mllZPeak");
      processCounter->AddVar("tightMETandPrMET");
      processCounter->AddVar("metOverPtLL");
      processCounter->AddVar("zeroJets");
      processCounter->AddVar("bTagVeto");
      processCounter->AddVar("nSoftMuons");
      processCounter->AddVar("nExtraLeptons");
      processCounter->AddVar("deltaPhi");
      processCounter->AddVar("final");
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

  if (_selection->getSwitch("maxPtLepton") && !_selection->passCut("maxPtLepton", m_highPt)) return false;
  theCounter->IncrVar("maxPtLepton",m_weight);

  if (_selection->getSwitch("minPtLepton") && !_selection->passCut("minPtLepton", m_lowPt)) return false;
  theCounter->IncrVar("minPtLepton",m_weight);

  if (_selection->getSwitch("leptonD0") && (!_selection->passCut("leptonD0",m_eleSlowD0) || !_selection->passCut("leptonD0",m_eleHardD0)) ) return false;
  theCounter->IncrVar("leptonD0",m_weight);

  if (_selection->getSwitch("eleIso") && (!m_isElectronIsol || !m_isPositronIsol)) return false; 
  theCounter->IncrVar("eleIso",m_weight);

  if (_selection->getSwitch("muTrackerIso") && 
      (!_selection->passCut("muTrackerIso",m_eleHardTkPtSum) || !_selection->passCut("muTrackerIso",m_eleSlowTkPtSum)) ) return false;
  theCounter->IncrVar("muTrackerIso",m_weight);

  if (_selection->getSwitch("muHcalIso") && 
      (!_selection->passCut("muHcalIso",m_eleHardHcalPtSum) || !_selection->passCut("muHcalIso",m_eleSlowHcalPtSum)) ) return false;
  theCounter->IncrVar("muHcalIso",m_weight);

  if (_selection->getSwitch("muEcalIso") && 
      (!_selection->passCut("muEcalIso",m_eleHardEcalPtSum) || !_selection->passCut("muEcalIso",m_eleSlowEcalPtSum)) ) return false;
  theCounter->IncrVar("muEcalIso",m_weight);

  if (_selection->getSwitch("muGlobalIso") && 
      (!_selection->passCut("muGlobalIso",m_eleHardGlobalPtSum) || !_selection->passCut("muGlobalIso",m_eleSlowGlobalPtSum)) ) return false;
  theCounter->IncrVar("muGlobalIso",m_weight);

  if (_selection->getSwitch("leptonId") && (!m_isElectronId || !m_isPositronId) ) return false; 
  theCounter->IncrVar("leptonId",m_weight);

  if (_selection->getSwitch("convRej") && (!m_isElectronConvRej || !m_isPositronConvRej) ) return false; 
  theCounter->IncrVar("convRej",m_weight);

  m_finalLeptons = true;

  if (_selection->getSwitch("looseMET") && !_selection->passCut("looseMET",m_met)) return false; 
  theCounter->IncrVar("looseMET",m_weight);

  if (_selection->getSwitch("mll") && !_selection->passCut("mll", m_invMass)) return false;
  theCounter->IncrVar("mll",m_weight);

  if (_selection->getSwitch("mllZPeak") && !_selection->passCut("mllZPeak", fabs(m_invMass-91.1876))) return false;
  theCounter->IncrVar("mllZPeak",m_weight);

  if (_selection->getSwitch("tightMET") && !_selection->passCut("tightMET",m_met)) return false; 
  if (_selection->getSwitch("projectedMET") && !_selection->passCut("projectedMET",m_projectedMet)) return false; 
  theCounter->IncrVar("tightMETandPrMET",m_weight);
  
  if (_selection->getSwitch("metOverPtLL") && !_selection->passCut("metOverPtLL",m_metOverPtLL)) return false;
  theCounter->IncrVar("metOverPtLL",m_weight);

  if (m_nJets==0) {
    theCounter->IncrVar("zeroJets",m_weight);
    m_jetVeto = true;

    // track counting = -100 means no output. Do not reject the event.
    if(!_selection->getSwitch("bTagVeto") ||
       (_selection->getSwitch("bTagVeto") && (m_btagJets==-100 || _selection->passCut("bTagVeto",m_btagJets)))) {
      theCounter->IncrVar("bTagVeto",m_weight);
      
      if(!_selection->getSwitch("nSoftMuons") ||
         (_selection->getSwitch("nSoftMuons") && _selection->passCut("nSoftMuons",m_nSoftMuons))) {
        theCounter->IncrVar("nSoftMuons",m_weight);
        
        if(!_selection->getSwitch("nExtraLeptons") ||
           _selection->getSwitch("nExtraLeptons") && _selection->passCut("nExtraLeptons",m_nExtraLeptons)) {
          theCounter->IncrVar("nExtraLeptons",m_weight);
          
          m_preDeltaPhi = true;
          
          if (!_selection->getSwitch("deltaPhi") ||
              _selection->getSwitch("deltaPhi") && _selection->passCut("deltaPhi", m_deltaPhi)) {
            theCounter->IncrVar("deltaPhi",m_weight); 
            theCounter->IncrVar("final",m_weight);
            return true;
          }
        }
      }
    }
  } else {

    // for nJets>0 we do not need cut by cut: just final counter
    if(_selection->getSwitch("bTagVeto") && !_selection->passCut("bTagVeto",m_btagJets)) return false;

    if(_selection->getSwitch("nSoftMuons") && !_selection->passCut("nSoftMuons",m_nSoftMuons)) return false; 

    if(_selection->getSwitch("nExtraLeptons") && !_selection->passCut("nExtraLeptons",m_nExtraLeptons)) return false; 

    if(_selection->getSwitch("deltaPhiLLJet") && !_selection->passCut("deltaPhiLLJet",m_deltaPhiLLJet)) return false;
    
    if (m_nJets==1) theCounter->IncrVar("oneJet",m_weight);
    if (m_nJets>1)  theCounter->IncrVar("gt1Jets",m_weight);
    
    return true;
  }

  return false;
}


void CutBasedHiggsSelector::displayEfficiencies(std::string datasetName) {

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter;
    for( iter=multiProcessCounter.begin(); iter!=multiProcessCounter.end(); ++iter ) {

      Counters *theCounter = iter->second;

      theCounter->Draw();
      theCounter->Draw("maxPtLepton","preselected");
      theCounter->Draw("minPtLepton","maxPtLepton");
      theCounter->Draw("leptonD0","minPtLepton");
      theCounter->Draw("eleIso","leptonD0");
      theCounter->Draw("muTrackerIso","eleIso");
      theCounter->Draw("muHcalIso","muTrackerIso");
      theCounter->Draw("muEcalIso","muHcalIso");
      theCounter->Draw("muGlobalIso","muEcalIso");
      theCounter->Draw("leptonId","muGlobalIso");      
      theCounter->Draw("convRej","leptonId");
      theCounter->Draw("looseMET","convRej");
      theCounter->Draw("mll","looseMET");
      theCounter->Draw("mllZPeak","mll");
      theCounter->Draw("tightMETandPrMET","mllZPeak");
      theCounter->Draw("metOverPtLL","tightMETandPrMET");
      theCounter->Draw("zeroJets","metOverPtLL");
      theCounter->Draw("bTagVeto","zeroJets");
      theCounter->Draw("nSoftMuons","bTagVeto");
      theCounter->Draw("nExtraLeptons","nSoftMuons");
      theCounter->Draw("deltaPhi","nExtraLeptons");
      theCounter->Draw("final","preselected");
      theCounter->Draw("oneJet","preselected");
      theCounter->Draw("gt1Jets","preselected");
    }
  }

  else {

    char namefile[500];
    sprintf(namefile,"%s-Counters.root",datasetName.c_str());
    
    globalCounter->Draw();
    globalCounter->Draw("maxPtLepton","preselected");
    globalCounter->Draw("minPtLepton","maxPtLepton");
    globalCounter->Draw("leptonD0","minPtLepton");
    globalCounter->Draw("eleIso","leptonD0");
    globalCounter->Draw("muTrackerIso","eleIso");
    globalCounter->Draw("muHcalIso","muTrackerIso");
    globalCounter->Draw("muEcalIso","muHcalIso");
    globalCounter->Draw("muGlobalIso","muEcalIso");
    globalCounter->Draw("leptonId","muGlobalIso");      
    globalCounter->Draw("convRej","leptonId");
    globalCounter->Draw("looseMET","convRej");
    globalCounter->Draw("mll","looseMET");
    globalCounter->Draw("mllZPeak","mll");
    globalCounter->Draw("tightMETandPrMET","mllZPeak");
    globalCounter->Draw("metOverPtLL","tightMETandPrMET");
    globalCounter->Draw("zeroJets","metOverPtLL");
    globalCounter->Draw("bTagVeto","zeroJets");
    globalCounter->Draw("nSoftMuons","bTagVeto");
    globalCounter->Draw("nExtraLeptons","nSoftMuons");
    globalCounter->Draw("deltaPhi","nExtraLeptons");
    globalCounter->Draw("final","preselected");
    globalCounter->Draw("oneJet","preselected");
    globalCounter->Draw("gt1Jets","preselected");
    
    globalCounter->Save(namefile,"update");
  }

}
