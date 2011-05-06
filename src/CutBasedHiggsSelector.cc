#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include <iostream>
#include <math.h>

CutBasedHiggsSelector::CutBasedHiggsSelector() {

  m_finalLeptons = false;
  m_jetVeto = false;
  m_preDeltaPhi = false;
  m_extraSlowLeptonPTMin = 0.;

  // latinos
  m_step1 = false;
  m_step2 = false;
  m_step3 = false;
  m_step4 = false;
  m_step5 = false;
  m_step6 = false;
  m_step6bis = false;
  m_step7 = false;
  m_step8 = false;
  m_step9 = false;
  m_step10 = false;
  m_step11 = false;
  m_step12 = false;
  m_step13 = false;
  m_step14 = false;
  m_step15 = false;
  m_step16 = false;
  m_step17 = false;

  m_processID = -1;
}

CutBasedHiggsSelector::CutBasedHiggsSelector( const CutBasedHiggsSelector& selector ) {

  m_weight = selector.m_weight;
  m_foundMcTree = selector.m_foundMcTree;
  m_passedHLT = selector.m_passedHLT;
  m_isThisChannel = selector.m_isThisChannel;
  m_highPt = selector.m_highPt;
  m_isElectronId = selector.m_isElectronId;
  m_isPositronId = selector.m_isPositronId;
  m_isElectronIsol = selector.m_isElectronIsol;
  m_isPositronIsol = selector.m_isPositronIsol;
  m_isElectronConvRej = selector.m_isElectronConvRej;
  m_isPositronConvRej = selector.m_isPositronConvRej;
  m_isElectronIp = selector.m_isElectronIp;
  m_isPositronIp = selector.m_isPositronIp;
  m_invMass = selector.m_invMass;
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
  m_extraSlowLeptonPTMin = selector.m_extraSlowLeptonPTMin;
  m_processID = selector.m_processID;
  *_selection = *selector._selection;
  *globalCounter = *selector.globalCounter;
  m_finalLeptons = selector.m_finalLeptons;
  m_jetVeto = selector.m_jetVeto;
  m_preDeltaPhi = selector.m_preDeltaPhi;
  multiProcessCounter = selector.multiProcessCounter;

  // latinos
  m_step1  = selector.m_step1;
  m_step2  = selector.m_step2;
  m_step3  = selector.m_step3;
  m_step4  = selector.m_step4;
  m_step5  = selector.m_step5;
  m_step6  = selector.m_step6;
  m_step6bis = selector.m_step6;
  m_step7  = selector.m_step7;
  m_step8  = selector.m_step8;
  m_step9  = selector.m_step9;
  m_step10 = selector.m_step10;
  m_step11 = selector.m_step11;
  m_step12 = selector.m_step12;
  m_step13 = selector.m_step13;
  m_step14 = selector.m_step14;
  m_step15 = selector.m_step15;
  m_step16 = selector.m_step16;
  m_step17 = selector.m_step17;
}

CutBasedHiggsSelector::~CutBasedHiggsSelector() {}

void CutBasedHiggsSelector::Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle) {

  _selection = new Selection(std::string(fileCuts),std::string(fileSwitches));

  // tehse cuts are applied in the HiggsSelection class, but are configured here
  _selection->addSwitch("MCtruth");
  _selection->addSwitch("trigger");
  _selection->addSwitch("leptonId");
  _selection->addSwitch("leptonIso");
  _selection->addSwitch("leptonD0"); 
  _selection->addSwitch("convRej");
  _selection->addCut("muGlobalIso");
  _selection->addCut("electronIP");
  _selection->addCut("muonIP");
  _selection->addCut("muonDz");
  _selection->addCut("nExtraLeptons");
  _selection->addCut("looseMET");
  _selection->addCut("mll");
  _selection->addCut("mllZPeak");
  _selection->addCut("tightMET");
  _selection->addCut("projectedMET");
  _selection->addCut("metOverPtLL");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("jetConeWidth");
  _selection->addCut("nSoftMuons");
  _selection->addCut("bTagVeto");
  _selection->addCut("mll2");
  _selection->addCut("maxPtLepton");
  _selection->addCut("minPtLepton");
  _selection->addCut("deltaPhiLLJet");
  _selection->addCut("deltaPhi");

  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle(theTitle);
  globalCounter->AddVar("event");
  globalCounter->AddVar("MCtruth");
  globalCounter->AddVar("preselected");
  globalCounter->AddVar("leptonId");
  globalCounter->AddVar("leptonIso");
  globalCounter->AddVar("convRej");
  globalCounter->AddVar("leptonD0");       
  globalCounter->AddVar("nExtraLeptons");
  globalCounter->AddVar("trigger");
  globalCounter->AddVar("looseMET");
  globalCounter->AddVar("mll");
  globalCounter->AddVar("mllZPeak");
  globalCounter->AddVar("tightMETandPrMET");
  globalCounter->AddVar("metOverPtLL");
  globalCounter->AddVar("zeroJets");
  globalCounter->AddVar("nSoftMuons");
  globalCounter->AddVar("bTagVeto");
  globalCounter->AddVar("mll2");
  globalCounter->AddVar("maxPtLepton");
  globalCounter->AddVar("minPtLepton");
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
      processCounter->AddVar("event");
      processCounter->AddVar("MCtruth");
      processCounter->AddVar("preselected");
      processCounter->AddVar("leptonId");
      processCounter->AddVar("leptonIso");      
      processCounter->AddVar("convRej");
      processCounter->AddVar("leptonD0");
      processCounter->AddVar("nExtraLeptons");
      processCounter->AddVar("trigger");
      processCounter->AddVar("looseMET");
      processCounter->AddVar("mll");
      processCounter->AddVar("mllZPeak");
      processCounter->AddVar("tightMETandPrMET");
      processCounter->AddVar("metOverPtLL");
      processCounter->AddVar("zeroJets");
      processCounter->AddVar("nSoftMuons");
      processCounter->AddVar("bTagVeto");
      processCounter->AddVar("mll2");
      processCounter->AddVar("maxPtLepton");
      processCounter->AddVar("minPtLepton");
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
  
  // latinos
  m_step1 = false;
  m_step2 = false;
  m_step3 = false;
  m_step4 = false;
  m_step5 = false;
  m_step6 = false;
  m_step6bis = false;
  m_step7 = false;
  m_step8 = false;
  m_step9 = false;
  m_step10 = false;
  m_step11 = false;
  m_step12 = false;
  m_step13 = false;
  m_step14 = false;
  m_step15 = false;
  m_step16 = false;
  m_step17 = false;
  m_step18 = false;
  m_step19 = false;
  m_step20 = false;
  m_step21 = false;
  m_step22 = false;
  m_step23 = false;
  m_step24 = false;

  theCounter->IncrVar("event",m_weight);
  
  if (_selection->getSwitch("MCtruth") && !m_foundMcTree) return false;
  theCounter->IncrVar("MCtruth",m_weight);

  if(!m_isThisChannel) return false;
  theCounter->IncrVar("preselected",m_weight);
  m_step1 = true;

  if (_selection->getSwitch("leptonId") && (m_isElectronId<0 || m_isPositronId<0) ) return false; 
  theCounter->IncrVar("leptonId",m_weight);
  m_step2 = true;

  if (_selection->getSwitch("leptonIso") && (m_isElectronIsol<0 || m_isPositronIsol<0)) return false; 
  theCounter->IncrVar("leptonIso",m_weight);
  m_step3 = true;

  if (_selection->getSwitch("convRej") && (m_isElectronConvRej<0 || m_isPositronConvRej<0) ) return false; 
  theCounter->IncrVar("convRej",m_weight);
  m_step4 = true;

  if (_selection->getSwitch("leptonD0") && (m_isElectronIp<0 || m_isPositronIp<0) ) return false; 
  theCounter->IncrVar("leptonD0",m_weight);
  m_step5 = true;
  
  m_finalLeptons = true;

  if (_selection->getSwitch("nExtraLeptons") && (!_selection->passCut("nExtraLeptons",m_nExtraLeptons)) ) return false;
  theCounter->IncrVar("nExtraLeptons",m_weight);
  m_step6 = true;

  if(_selection->getSwitch("trigger") && !m_passedHLT ) return false;
  theCounter->IncrVar("trigger",m_weight); 
  m_step6bis = true;
  
  if (_selection->getSwitch("looseMET") && !_selection->passCut("looseMET",m_met)) return false; 
  theCounter->IncrVar("looseMET",m_weight);
  m_step7 = true;

  if (_selection->getSwitch("mll") && !_selection->passCut("mll", m_invMass)) return false;
  theCounter->IncrVar("mll",m_weight);
  m_step8 = true;  

  if (_selection->getSwitch("mllZPeak") && !_selection->passCut("mllZPeak", fabs(m_invMass-91.1876))) return false;
  theCounter->IncrVar("mllZPeak",m_weight);
  m_step9 = true;

  if (_selection->getSwitch("tightMET") && !_selection->passCut("tightMET",m_met)) return false; 
  if (_selection->getSwitch("projectedMET") && !_selection->passCut("projectedMET",m_projectedMet)) return false; 
  theCounter->IncrVar("tightMETandPrMET",m_weight);

  if (_selection->getSwitch("metOverPtLL") && !_selection->passCut("metOverPtLL",m_metOverPtLL)) return false;
  theCounter->IncrVar("metOverPtLL",m_weight);
  m_step10 = true;

  if (m_nJets==0) {
    theCounter->IncrVar("zeroJets",m_weight);
    m_jetVeto = true;
    m_step11 = true;

    if(!_selection->getSwitch("nSoftMuons") ||
       (_selection->getSwitch("nSoftMuons") && _selection->passCut("nSoftMuons",m_nSoftMuons))) {
      theCounter->IncrVar("nSoftMuons",m_weight);
      m_step12 = true;

      if(!_selection->getSwitch("bTagVeto") ||
	 (_selection->getSwitch("bTagVeto") && _selection->passCut("bTagVeto",m_btagJets))) {
	theCounter->IncrVar("bTagVeto",m_weight);
	m_step13 = true;
      
	if (!_selection->getSwitch("mll2") ||
	    (_selection->getSwitch("mll2") && _selection->passCut("mll2", m_invMass))) {
	  theCounter->IncrVar("mll2",m_weight);
	  m_step14 = true;

	  if (!_selection->getSwitch("maxPtLepton") || 
	      (_selection->getSwitch("maxPtLepton") && _selection->passCut("maxPtLepton", m_highPt))) {
	    theCounter->IncrVar("maxPtLepton",m_weight);
	    m_step15 = true;

	    if (!_selection->getSwitch("minPtLepton") || 
		(_selection->getSwitch("minPtLepton") && _selection->passCut("minPtLepton", m_lowPt))
                && m_lowPt >= m_extraSlowLeptonPTMin ) {
	      theCounter->IncrVar("minPtLepton",m_weight);
	      m_step16 = true;

	      m_preDeltaPhi = true;
	      
	      if (!_selection->getSwitch("deltaPhi") ||
		  _selection->getSwitch("deltaPhi") && _selection->passCut("deltaPhi", m_deltaPhi)) {
		theCounter->IncrVar("deltaPhi",m_weight); 
		theCounter->IncrVar("final",m_weight);
		m_step17 = true;
		return true;
	      }
	    }
	  }
	}
      }
    }
  } else {

    // for nJets>0 we do not need cut by cut: just final counter
    if(_selection->getSwitch("bTagVeto") && !_selection->passCut("bTagVeto",m_btagJets)) return false;
    m_step18 = true;

    if(_selection->getSwitch("nSoftMuons") && !_selection->passCut("nSoftMuons",m_nSoftMuons)) return false; 
    m_step19 = true;

    if(_selection->getSwitch("deltaPhiLLJet") && !_selection->passCut("deltaPhiLLJet",m_deltaPhiLLJet)) return false;
    m_step20 = true;

    if(_selection->getSwitch("mll2") && !_selection->passCut("mll2", m_invMass)) return false;
    m_step21 = true;

    if(_selection->getSwitch("maxPtLepton") && !_selection->passCut("maxPtLepton", m_highPt)) return false;
    m_step22 = true;

    if (_selection->getSwitch("minPtLepton") && _selection->passCut("minPtLepton", m_lowPt)) return false;
    m_step23 = true;

    if (_selection->getSwitch("deltaPhi") && _selection->passCut("deltaPhi", m_deltaPhi)) return false;
    m_step24 = true;

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
      theCounter->Draw("leptonId","preselected");     
      theCounter->Draw("leptonIso","leptonId");
      theCounter->Draw("convRej","leptonIso");
      theCounter->Draw("leptonD0","convRej");
      theCounter->Draw("nExtraLeptons","leptonD0");
      theCounter->Draw("trigger","nExtraLeptons");
      theCounter->Draw("looseMET","trigger");
      theCounter->Draw("mll","looseMET");
      theCounter->Draw("mllZPeak","mll");
      theCounter->Draw("tightMETandPrMET","mllZPeak");
      theCounter->Draw("metOverPtLL","tightMETandPrMET");
      theCounter->Draw("zeroJets","metOverPtLL");
      theCounter->Draw("nSoftMuons","zeroJets");
      theCounter->Draw("bTagVeto","nSoftMuons");
      theCounter->Draw("mll2","bTagVeto");
      theCounter->Draw("maxPtLepton","mll2");
      theCounter->Draw("minPtLepton","maxPtLepton");
      theCounter->Draw("deltaPhi","deltaPhi");
      theCounter->Draw("final","preselected");
      theCounter->Draw("oneJet","preselected");
      theCounter->Draw("gt1Jets","preselected");
    }
  }

  else {

    char namefile[500];
    sprintf(namefile,"%s-Counters.root",datasetName.c_str());
    
    globalCounter->Draw();
    globalCounter->Draw("leptonId","preselected");     
    globalCounter->Draw("leptonIso","leptonId");
    globalCounter->Draw("convRej","leptonIso");
    globalCounter->Draw("leptonD0","convRej");
    globalCounter->Draw("nExtraLeptons","leptonD0");
    globalCounter->Draw("trigger","nExtraLeptons");
    globalCounter->Draw("looseMET","trigger");
    globalCounter->Draw("mll","looseMET");
    globalCounter->Draw("mllZPeak","mll");
    globalCounter->Draw("tightMETandPrMET","mllZPeak");
    globalCounter->Draw("metOverPtLL","tightMETandPrMET");
    globalCounter->Draw("zeroJets","metOverPtLL");
    globalCounter->Draw("nSoftMuons","zeroJets");
    globalCounter->Draw("bTagVeto","nSoftMuons");
    globalCounter->Draw("mll2","bTagVeto");
    globalCounter->Draw("maxPtLepton","mll2");
    globalCounter->Draw("minPtLepton","maxPtLepton");
    globalCounter->Draw("deltaPhi","deltaPhi");
    globalCounter->Draw("final","preselected");
    globalCounter->Draw("oneJet","preselected");
    globalCounter->Draw("gt1Jets","preselected");
    
    globalCounter->Save(namefile,"update");
  }

}
