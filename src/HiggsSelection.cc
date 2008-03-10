#include <iostream>
#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"
#include "HiggsAnalysisTools/include/HiggsSelection.hh"

HiggsSelection::HiggsSelection(TTree *tree) 
  : HiggsBase(tree) {

  // choose the Higgs Mass
  std::string higgsConfigDir;
  std::ifstream setfile("config/higgs/higgsMass.txt");
  if(!setfile.good()) {
    std::cout << "Cannot read the higgsMass file to choose the selection: config/higgs/higgsMass.txt" << std::endl
	      << "using default (160 GeV)" << std::endl;
    higgsConfigDir="config/higgs/h160/";
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
	higgsConfigDir="config/higgs/h" + massVal + "/";
	std::cout << "Reading configuration for Higgs mass = " << massVal << " GeV/c^2" << std::endl;
	break;
      }
    }
  }
  

  // Initialize selection
  std::string fileCuts  = higgsConfigDir + "2e2nuCuts.txt";
  std::string fileSwitches = "config/higgs/2e2nuSwitches.txt";
  _selection = new Selection(fileCuts,fileSwitches);
  _selection->addSwitch("MCtruth");
  _selection->addSwitch("trigger");
  _selection->addSwitch("preselection");
  _selection->addSwitch("jetVeto");
  _selection->addSwitch("classDepEleId");
  _selection->addSwitch("apply_kFactor");
  _selection->addCut("nRecoLeptons");
  _selection->addCut("etaElectronAcc");
  _selection->addCut("ptElectronAcc");
  _selection->addCut("etaMuonAcc");
  _selection->addCut("ptMuonAcc");
  _selection->addCut("HardLeptonThreshold");
  _selection->addCut("METpreselection");
  _selection->addCut("dileptonInvMassMin");
  _selection->addCut("trackerPtSum");
  _selection->addCut("hcalPtSum");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetLowAcc");
  _selection->addCut("etJetHighAcc");
  _selection->addCut("alphaJet");
  _selection->addCut("MET");
  _selection->addCut("deltaPhi");
  _selection->addCut("eleInvMass");
  _selection->addCut("detaLeptons");

  Selection *_goldenSelectionEB = new Selection("config/higgs/electronIDGoldenCutsEB.txt","config/higgs/electronIDSwitches.txt");
  Selection *_bigbremSelectionEB = new Selection("config/higgs/electronIDBigBremCutsEB.txt","config/higgs/electronIDSwitches.txt");
  Selection *_narrowSelectionEB = new Selection("config/higgs/electronIDNarrowCutsEB.txt","config/higgs/electronIDSwitches.txt");
  Selection *_showeringSelectionEB = new Selection("config/higgs/electronIDShoweringCutsEB.txt","config/higgs/electronIDSwitches.txt");
  Selection *_goldenSelectionEE = new Selection("config/higgs/electronIDGoldenCutsEE.txt","config/higgs/electronIDSwitches.txt");
  Selection *_bigbremSelectionEE = new Selection("config/higgs/electronIDBigBremCutsEE.txt","config/higgs/electronIDSwitches.txt");
  Selection *_narrowSelectionEE = new Selection("config/higgs/electronIDNarrowCutsEE.txt","config/higgs/electronIDSwitches.txt");
  Selection *_showeringSelectionEE = new Selection("config/higgs/electronIDShoweringCutsEE.txt","config/higgs/electronIDSwitches.txt");

  _electronSelection.push_back(_goldenSelectionEB);
  _electronSelection.push_back(_bigbremSelectionEB);
  _electronSelection.push_back(_narrowSelectionEB);
  _electronSelection.push_back(_showeringSelectionEB);
  _electronSelection.push_back(_goldenSelectionEE);
  _electronSelection.push_back(_bigbremSelectionEE);
  _electronSelection.push_back(_narrowSelectionEE);
  _electronSelection.push_back(_showeringSelectionEE);

  std::vector<Selection*>::iterator selectionItr;
  for(selectionItr=_electronSelection.begin();selectionItr!=_electronSelection.end();selectionItr++) {
    Selection *eleSelection = *selectionItr;
    eleSelection->addCut("hOverE");
    eleSelection->addCut("s9s25");
    eleSelection->addCut("deta");
    eleSelection->addCut("dphiIn");
    eleSelection->addCut("dphiOut");
    eleSelection->addCut("covEtaEta");
    eleSelection->addCut("eOverPout");
    eleSelection->addCut("Fisher");
    eleSelection->addCut("Likelihood");
    eleSelection->summary();
  }

  std::map<std::string,std::pair<float,float> > selectionMap = _selection->getSelection();
  _selection->summary();

  // Initialize additional counters...
  _counter.SetTitle("EVENT COUNTER");
  _counter.AddVar("event");
  _counter.AddVar("MCtruth");
  _counter.AddVar("trigger");
  _counter.AddVar("nRecoLeptons");
  _counter.AddVar("twoGoodRec");
  _counter.AddVar("HardLeptonThreshold");
  _counter.AddVar("METpreselection");
  _counter.AddVar("dileptonInvMassMin");
  _counter.AddVar("preselection");
  _counter.AddVar("eleID");
  _counter.AddVar("trackerIsol");
  _counter.AddVar("hcalIsol");
  _counter.AddVar("jetVeto");
  _counter.AddVar("MET");
  _counter.AddVar("deltaPhi");
  _counter.AddVar("eleInvMass");
  _counter.AddVar("detaLeptons");
  _counter.AddVar("final");
  
  _eleCounter.SetTitle("SINGLE ELECTRON COUNTER");
  _eleCounter.AddVar("electrons");
  _eleCounter.AddVar("hOverE");
  _eleCounter.AddVar("s9s25");
  _eleCounter.AddVar("deta");
  _eleCounter.AddVar("dphiIn");
  _eleCounter.AddVar("dphiOut");
  _eleCounter.AddVar("covEtaEta");
  _eleCounter.AddVar("eOverPout");
  _eleCounter.AddVar("Fisher");
  _eleCounter.AddVar("Likelihood");
  _eleCounter.AddVar("finalEleID");

  m_p4ElectronPlus = new TLorentzVector(0.,0.,0.,0.);
  m_p4ElectronMinus = new TLorentzVector(0.,0.,0.,0.);
  m_p4MuonPlus = new TLorentzVector(0.,0.,0.,0.);
  m_p4MuonMinus = new TLorentzVector(0.,0.,0.,0.);
  m_p4MET = new TLorentzVector(0.,0.,0.,0.);

  // the Monitoring Histograms
   _monitorGenerator = new Monitor(0);
   _monitorEventAfterSelection = new Monitor(0);
   _monitorEventAfterReco = new Monitor(0);
   _monitorMet = new Monitor(&nMet);
   _bestElectrons = new std::vector<int>;
   _bestMuons = new std::vector<int>;
   _monitorElectrons = new Monitor(&nEle,_bestElectrons);
   _bestJets = new std::vector<int>;
   _excludedJets = new std::vector<int>;
   _monitorJets = new Monitor(&nJet,_bestJets);
   _monitorJets->ExcludeElements(_excludedJets);
   _bestGenJets = new std::vector<int>;
   _monitorGenJets = new Monitor(&nGenJet,_bestGenJets);

}

HiggsSelection::~HiggsSelection(){

  delete m_p4ElectronPlus;
  delete m_p4ElectronMinus;
  delete m_p4MuonPlus;
  delete m_p4MuonMinus;
  delete m_p4MET;
  delete _monitorGenerator;
  delete _monitorEventAfterSelection;
  delete _monitorEventAfterReco;
  delete _monitorMet;
  delete _bestElectrons;
  delete _bestMuons;
  delete _monitorElectrons;
  delete _bestJets;
  delete _excludedJets;
  delete _bestGenJets;
  delete _monitorGenJets;
  delete _monitorJets;
  delete _selection;
  myOutTree -> save();

}

bool HiggsSelection::findMcTree(const char* processType) {
  _process = "UNDEFINED";
  _theGenEle = -1;
  _theGenPos = -1;

  if(strcmp(processType,"HtoWWto2e2nu")==0) {
    int indlminus=999, indlplus=999;
    for(int imc=6;imc<25;imc++) {
      // look for W->l nu decay
      if(idMc[imc]>10 && idMc[imc]<19 && idMc[mothMc[imc]]==-24) indlminus=imc;
      if(idMc[imc]<-10 && idMc[imc]>-19 && idMc[mothMc[imc]]==24) indlplus=imc;
      // kFactors are different for gg->H and VBF->H production
      // in our trees this can be accessed as follows:
      // VBF: higgs is in position [8], mother = quark (lund-id < 10)
      // gg: higgs is in position [6], mother = gluon (lund-id = 21)
      if(idMc[8]==25 && idMc[mothMc[8]]<10) { 
	_process = "H_VBF";
	float pT = pMc[8]*fabs(sin(thetaMc[8]));
	_genHiggsPt[0] = pT;
      }
      else if(idMc[6]==25 && idMc[mothMc[6]]==21) {
	_process = "H_gg";
	float pT = pMc[6]*fabs(sin(thetaMc[6]));
	_genHiggsPt[0] = pT;
      }
    }
    if(indlminus<25 && indlplus<25) {
      // set electron index only if the primary lepton from W is an electron
      if(idMc[indlminus]==-11)
	_theGenEle = indlminus;
      if(idMc[indlplus]==11)
	_theGenPos = indlplus;
    }
    return (indlminus<25 && indlplus<25);
  }
  else if(strcmp(processType,"WW")==0) {
    _process = "WW";
    TVector3 WminusP, WplusP;
    WminusP.SetMagThetaPhi(pMc[6],thetaMc[6],phiMc[6]);
    WplusP.SetMagThetaPhi(pMc[7],thetaMc[7],phiMc[7]);
    float pT = (WminusP+WplusP).Pt();
    _theGenEle = 6;
    _theGenPos = 7;
    _genHiggsPt[0] = pT;
    // 8,9; 10,11 are the daughters of W;W
    // W->e nu exclusive
//     return ( 
// 	    (abs(idMc[6])==24) && (abs(idMc[7])==24) &&
// 	    (abs(idMc[8])==11 && abs(idMc[mothMc[8]])==24) &&
// 	    (abs(idMc[10])==11 && abs(idMc[mothMc[10]])==24)
// 	    );
// W -> lnu all the leptons
    return (
	    (abs(idMc[6])==24) && (abs(idMc[7])==24) &&
	    (abs(idMc[8])>10 && abs(idMc[8])<19 && abs(idMc[mothMc[8]])==24) &&
	    (abs(idMc[10])>10 && abs(idMc[10])<19 && abs(idMc[mothMc[10]])==24)
	    );
  }
  else if(strcmp(processType,"Wjets")==0) {
    _process = "Wjets";
    return ( (abs(idMc[8])==11) && abs(idMc[9])==12 );
  }
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

float HiggsSelection::getkFactor(std::string process) {
  kFactorEvaluator *kFactor;
  float pT = -1.;

  if(process.compare("H_VBF")==0) {
    kFactor = new VBFHiggsEvaluator();
    pT = pMc[8]*fabs(sin(thetaMc[8]));
  }
  else if(process.compare("H_gg")==0) {
    // kFactor depends on Higgs mass
    switch(_massVal) {
    case 120 :
      kFactor = new ggHiggs120Evaluator();
      break;
    case 130 :
      kFactor = new ggHiggs130Evaluator();
      break;
    case 140 :
      kFactor = new ggHiggs140Evaluator();
      break;
    case 150 :
      kFactor = new ggHiggs150Evaluator();
      break;
    case 160 :
      kFactor = new ggHiggs160Evaluator();
      break;
    case 170 :
      kFactor = new ggHiggs170Evaluator();
      break;
    case 180 :
      kFactor = new ggHiggs180Evaluator();
      break;
    case 190 :
      kFactor = new ggHiggs190Evaluator();
      break;
    case 200 :
      kFactor = new ggHiggs200Evaluator();
      break;
    default :
      kFactor = new flatHiggsEvaluator();
    }
    pT = pMc[6]*fabs(sin(thetaMc[6]));
  }
  else if(process.compare("WW")==0) {
    kFactor = new WWEvaluator();
    TVector3 WminusP, WplusP;
    WminusP.SetMagThetaPhi(pMc[6],thetaMc[6],phiMc[6]);
    WplusP.SetMagThetaPhi(pMc[7],thetaMc[7],phiMc[7]);
    pT = (WminusP+WplusP).Pt();
  }

  float weight = kFactor->evaluate(pT);

  delete kFactor;
  return weight;

}

void HiggsSelection::Loop() {
  _verbose=true;
  if(fChain == 0) return;


  bookHistos();
  std::string::size_type loc = _datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    _datasetName.erase(loc);
  }
  std::string recoHistogramName = _datasetName+"Histograms.root";
  TFile *file = new TFile(recoHistogramName.c_str(),"RECREATE");
  _monitorGenerator->setPath("Generator");
  _monitorMet->setPath("MET");
  _monitorElectrons->setPath("Electrons");
  _monitorJets->setPath("Jets");
  _monitorGenJets->setPath("GenJets");
  _monitorEventAfterReco->setPath("EventAfterReco");
  _monitorEventAfterSelection->setPath("EventAfterSelection");

  // kinematics reduced tree
  std::string reducedTreeName = _datasetName+"-dataset.root";
  myOutTree = new RedHiggsTree(reducedTreeName.c_str());
  // myOutTree = new RedHiggsTree("/cmsrm/pc18/crovelli/tree.root");

  float met, deltaPhi, transvMass, dileptonInvMass, maxPtEle, minPtEle, detaLeptons;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  if(_verbose) std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (_verbose && jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

      // look to the MC truth decay tree
      bool foundMcTree = findMcTree("HtoWWto2e2nu");

      // get the kFactor of the event
      float weight = 1.0;
      if(_selection->getSwitch("apply_kFactor"))
	 weight = getkFactor(_process);

      _counter.IncrVar("event",weight);

      if(_selection->getSwitch("MCtruth") && 
	 !foundMcTree ) continue;              
      _counter.IncrVar("MCtruth",weight);
      _monitorGenerator->Fill(weight);

      Utils anaUtils;
      if(_selection->getSwitch("trigger") && !( anaUtils.getTriggersAND(m_requiredTriggers, firedTrg) ) ) continue;
      _counter.IncrVar("trigger",weight); 


      if ( _selection->getSwitch("nRecoLeptons") && 
	   ( !_selection->passCut("nRecoLeptons", nEle) && 
	     !_selection->passCut("nRecoLeptons", nMuon) ) 
	   ) continue;   
      _counter.IncrVar("nRecoLeptons",weight);

      // get the best electrons, best muons
      std::pair<int,int> theElectrons = getBestElectronPair();
      int theElectron(theElectrons.first), thePositron(theElectrons.second);

      std::pair<int,int> theMuons = getBestMuonPair();
      int theMuonPlus(theMuons.first), theMuonMinus(theMuons.second);

      std::cout << "processing jevt " << jentry << std::endl;
      std::cout << "DEBUG: " << "theElectron = " << theElectron << " thePositron = " << thePositron
		<< "\ttheMuonPlus = " << theMuonPlus << " theMuonMinus = " << theMuonMinus << std::endl;
      

      m_channel[ee] = false;
      m_channel[mm] = false;
      m_channel[em] = false;
      if ( theElectron > -1 && thePositron > -1 ) m_channel[ee] = true;
      if ( theMuonPlus > -1 && theMuonMinus > -1 ) m_channel[mm] = true;
      if ( ( theElectron > -1 && theMuonPlus > -1 ) ||
	   ( thePositron > -1 && theMuonMinus > -1 ) ) m_channel[em] = true;

      std::cout << "\tDEBUG: " << "m_channel[ee] = " << m_channel[ee]
		<< "m_channel[mm] = " << m_channel[mm]
		<< "m_channel[em] = " << m_channel[em] << std::endl;
	
      if( !m_channel[ee] && !m_channel[mm] && !m_channel[em]) continue; 
      _counter.IncrVar("twoGoodRec",weight);

      float hardestElectronPt = TMath::Max(etMuon[theMuonPlus],etMuon[theMuonMinus]);
      float slowestElectronPt = TMath::Min(etMuon[theMuonPlus],etMuon[theMuonMinus]);
      float hardestMuonPt = TMath::Max(etMuon[theMuonPlus],etMuon[theMuonMinus]);
      float slowestMuonPt = TMath::Min(etMuon[theMuonPlus],etMuon[theMuonMinus]);

      if ( _selection->getSwitch("HardLeptonThreshold") &&
	   (!_selection->passCut("HardLeptonThreshold", hardestElectronPt) && 
	    !_selection->passCut("HardLeptonThreshold", hardestMuonPt) )
	   ) continue;
      _counter.IncrVar("HardLeptonThreshold",weight);

      setKinematics( );
      estimateJetMatch(0.0);
      _monitorJets->Fill(weight);
      _monitorEventAfterReco->Fill(weight);


      if ( _selection->getSwitch("METpreselection") && 
	   !_selection->passCut("METpreselection", etMet[0] ) ) continue;
      _counter.IncrVar("METpreselection",weight);


      if ( _selection->getSwitch("dileptonInvMassMin") && 
	   (!_selection->passCut("dileptonInvMassMin", m_mll[ee]) &&
	    !_selection->passCut("dileptonInvMassMin", m_mll[mm]) &&
	    !_selection->passCut("dileptonInvMassMin", m_mll[em]) )
	   ) continue;
      _counter.IncrVar("dileptonInvMassMin",weight);


      // the cut is not really applied, should be reproduced by previous ones
      if( _selection->getSwitch("preselection") && evtPresel ) { 
	_counter.IncrVar("preselection",weight);
      }


      // from here: add muons, eleID with: 1) egamma cuts 2) custom cuts 3) lik... etc
      if (theElectron < 0 || thePositron < 0 ) continue; // FIXME
      if(!isEleID(theElectron) || !isEleID(thePositron)) continue; 
      _counter.IncrVar("eleID",weight);



      if(_selection->getSwitch("trackerPtSum") && 
	 ( !_selection->passCut("trackerPtSum",eleTrackerIso_sumPtEle[theElectron]) || 
	   !_selection->passCut("trackerPtSum",eleTrackerIso_sumPtEle[thePositron]) ) ) continue; 
      _counter.IncrVar("trackerIsol",weight);



      if(_selection->getSwitch("hcalPtSum") &&
	 ( !_selection->passCut("hcalPtSum",eleCaloIso_sumPtEle[theElectron]) || 
	   !_selection->passCut("hcalPtSum",eleCaloIso_sumPtEle[thePositron]) ) ) continue; 
      _counter.IncrVar("hcalIsol",weight);



      if(_selection->getSwitch("jetVeto") && jetVeto()) continue; 
      _counter.IncrVar("jetVeto",weight);


      
      if(_selection->getSwitch("MET") && !_selection->passCut("MET",met)) continue; 
      _counter.IncrVar("MET",weight);



      if(_selection->getSwitch("deltaPhi") && 
	 !_selection->passCut("deltaPhi", m_deltaPhi[ee]) ) continue; 
      _counter.IncrVar("deltaPhi",weight);



      if(_selection->getSwitch("eleInvMass") && 
	 !_selection->passCut("eleInvMass", m_mll[ee])) continue; 
      _counter.IncrVar("eleInvMass",weight);


      detaLeptons = etaEle[theElectron]-etaEle[thePositron];
      if(_selection->getSwitch("detaLeptons") && 
	 !_selection->passCut("detaLeptons",detaLeptons) ) continue; 
      _counter.IncrVar("detaLeptons",weight);
      _counter.IncrVar("final",weight);
      
      _monitorEventAfterSelection->Fill(weight);
      _monitorMet->Fill(weight);
      _monitorElectrons->Fill(weight);
      _monitorGenJets->Fill(weight);

      // dumping the reco variables in a tree - after all the cuts
      
      myOutTree -> fillAll(etMet[0], m_deltaPhi[ee], m_transvMass[ee], m_mll[ee], hardestElectronPt, slowestElectronPt, detaLeptons);
      myOutTree -> store();

  }
  _monitorGenerator->WritePs("eventGenerator.ps");
  _monitorGenerator->WriteRoot(file);
  _monitorEventAfterReco->WritePs("eventAfterReco.ps");
  _monitorEventAfterReco->WriteRoot(file);
  _monitorEventAfterSelection->WritePs("eventAfterSelection.ps");
  _monitorEventAfterSelection->WriteRoot(file);
  _monitorMet->WritePs("met.ps");
  _monitorMet->WriteRoot(file);
  _monitorElectrons->WritePs("electrons.ps");
  _monitorElectrons->WriteRoot(file);
  _monitorJets->WritePs("jets.ps");
  _monitorJets->WriteRoot(file);
  _monitorGenJets->WritePs("jets.ps");
  _monitorGenJets->WriteRoot(file);
  file->Close();

}

void HiggsSelection::displayEfficiencies() {
  _eleCounter.Draw();
  _eleCounter.Draw("hOverE","electrons");
  _eleCounter.Draw("s9s25","hOverE");
  _eleCounter.Draw("deta","s9s25");
  _eleCounter.Draw("dphiIn","deta");
  _eleCounter.Draw("dphiOut","dphiIn");
  _eleCounter.Draw("covEtaEta","dphiOut");
  _eleCounter.Draw("eOverPout","covEtaEta");
  _eleCounter.Draw("Fisher","eOverPout");
  _eleCounter.Draw("Likelihood","Fisher");
  _eleCounter.Draw("finalEleID","electrons");

  _counter.Draw();
  _counter.Draw("MCtruth","event");
  _counter.Draw("trigger","MCtruth");
  _counter.Draw("nRecoLeptons","trigger");
  _counter.Draw("twoGoodRec","nRecoLeptons");
  _counter.Draw("HardLeptonThreshold","twoGoodRec");
  _counter.Draw("METpreselection","HardLeptonThreshold");
  _counter.Draw("dileptonInvMassMin","METpreselection");
  _counter.Draw("preselection","dileptonInvMassMin");
  _counter.Draw("eleID","twoGoodRec");
  _counter.Draw("trackerIsol","eleID");
  _counter.Draw("hcalIsol","trackerIsol");
  _counter.Draw("jetVeto","presel");
  _counter.Draw("MET","jetVeto");
  _counter.Draw("deltaPhi","MET");
  _counter.Draw("eleInvMass","deltaPhi");
  _counter.Draw("detaLeptons","minPtEle");
  _counter.Draw("final","MCtruth");

  // jet match histogram
  TH1F *MatchFracJets_pt = (TH1F*)RecoJets_pt->Clone("MatchFracJets_pt");
  MatchFracJets_pt->Sumw2();
  MatchFracJets_pt->Divide(MatchedJets_pt, RecoJets_pt, 1, 1);

  TFile jetMatchFile("jetMatchFile.root","RECREATE");
  RecoJets_pt->Write();
  MatchedJets_pt->Write();
  MatchFracJets_pt->Write();
  etHighestJet->Write();
  jetMatchFile.Close();

}

std::pair<int,int> HiggsSelection::getBestElectronPair() {
  int theLep1=-1;
  int theLep2=-1;
  float maxPtLep1=-1000.;
  float maxPtLep2=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nEle;i++) {
    if(_selection->getSwitch("etaElectronAcc") && !_selection->passCut("etaElectronAcc",etaEle[i]) ) continue;
    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();
    if(_selection->getSwitch("ptElectronAcc") && !_selection->passCut("ptElectronAcc",thisPt) ) continue;
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestElectrons->clear();
  _bestElectrons->push_back(theLep1);  _bestElectrons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

std::pair<int,int> HiggsSelection::getBestMuonPair() {
  int theLep1=-1;
  int theLep2=-1;
  float maxPtLep1=-1000.;
  float maxPtLep2=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nMuon;i++) {
    if(_selection->getSwitch("etaMuonAcc") && !_selection->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    TVector3 pLepton(pxMuon[i],pyMuon[i],pzMuon[i]);
    float thisPt=pLepton.Pt();
    if(_selection->getSwitch("ptMuonAcc") && !_selection->passCut("ptMuonAcc",thisPt) ) continue;
    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestMuons->clear();
  _bestMuons->push_back(theLep1);  _bestMuons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

bool HiggsSelection::isEleID(int eleIndex) {

  int GsfClass = eleClassEle[eleIndex];
  int offset=0;
  if(GsfClass>=100) {
    GsfClass-=100;
    offset=4;
  }

  Selection *selection;
  if(GsfClass==0) 
    selection=_electronSelection[offset];
  else if(GsfClass==10)
    selection=_electronSelection[offset+1];
  else if(GsfClass==20)
    selection=_electronSelection[offset+2];
  else if(GsfClass>=30 && GsfClass<=40) {
    selection=_electronSelection[offset+3];
  }

  _eleCounter.IncrVar("electrons");
  if(selection->getSwitch("hOverE") && !selection->passCut("hOverE",eleHoEEle[eleIndex])) return false;  _eleCounter.IncrVar("hOverE");
  if(selection->getSwitch("s9s25") && !selection->passCut("s9s25",s9s25Ele[eleIndex])) return false; _eleCounter.IncrVar("s9s25");
  if(selection->getSwitch("deta") && !selection->passCut("deta",eleDeltaEtaAtVtxEle[eleIndex])) return false; _eleCounter.IncrVar("deta");
  if(selection->getSwitch("dphiIn") && !selection->passCut("dphiIn",eleDeltaPhiAtVtxEle[eleIndex])) return false; _eleCounter.IncrVar("dphiIn");
  if(selection->getSwitch("dphiOut") && !selection->passCut("dphiOut",eleDeltaPhiAtCaloEle[eleIndex])) return false; _eleCounter.IncrVar("dphiOut");
  if(selection->getSwitch("covEtaEta") && !selection->passCut("covEtaEta",covEtaEtaEle[eleIndex])) return false; _eleCounter.IncrVar("covEtaEta");
  if(selection->getSwitch("eOverPout") && !selection->passCut("eOverPout",eleCorrEoPoutEle[eleIndex])) return false; _eleCounter.IncrVar("eOverPout");
  if(selection->getSwitch("Fisher") && !selection->passCut("Fisher",Fisher(eleIndex))) return false; _eleCounter.IncrVar("Fisher");
  if(selection->getSwitch("Likelihood") && !selection->passCut("Likelihood",eleLikelihoodEle[eleIndex])) return false; _eleCounter.IncrVar("Likelihood");
  _eleCounter.IncrVar("finalEleID");

  return true;
}



void HiggsSelection::setKinematics( ) {

  int theElectronPlus = (*_bestElectrons)[0], theElectronMinus = (*_bestElectrons)[1];
  int theMuonPlus = (*_bestMuons)[0], theMuonMinus = (*_bestMuons)[1];

  // lepton four-vectors
  m_p4ElectronMinus->SetXYZT(pxEle[theElectronMinus],pyEle[theElectronMinus],pzEle[theElectronMinus],energyEle[theElectronMinus]);
  m_p4ElectronPlus->SetXYZT(pxEle[theElectronPlus],pyEle[theElectronPlus],pzEle[theElectronPlus],energyEle[theElectronPlus]);

  m_p4MuonMinus->SetXYZT(pxEle[theMuonMinus],pyEle[theMuonMinus],pzEle[theMuonMinus],energyEle[theMuonMinus]);
  m_p4MuonPlus->SetXYZT(pxEle[theMuonPlus],pyEle[theMuonPlus],pzEle[theMuonPlus],energyEle[theMuonPlus]);

  // MET
  m_p4MET->SetXYZT(pxMet[0],pyMet[0],pzMet[0],energyMet[0]);

  // compute delta Phi in degrees, di-lepton invariant mass, transverse mass
  TVector3 dilepPt;
  if ( m_channel[ee] ) {

    m_deltaPhi[ee] = fabs( 180./TMath::Pi() * 
			   m_p4ElectronMinus->Vect().DeltaPhi(m_p4ElectronPlus->Vect()) );
    m_mll[ee] = (*m_p4ElectronMinus + *m_p4ElectronPlus).M();
    dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4ElectronPlus->Vect().X(),
		    m_p4ElectronMinus->Vect().Y()+m_p4ElectronPlus->Vect().Y(),
		    0.0 );
    m_transvMass[ee]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );

  }
  else {

    m_deltaPhi[ee] = -1.;
    m_mll[ee] = -1.;
    m_transvMass[ee] = -1.;

  }

  if ( m_channel[mm] ) {

    m_deltaPhi[mm] = fabs( 180./TMath::Pi() * 
			   m_p4MuonMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()) );
    m_mll[mm] = (*m_p4MuonMinus + *m_p4MuonPlus).M();
    dilepPt.SetXYZ( m_p4MuonMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		    m_p4MuonMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		    0.0 );
    m_transvMass[mm]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );

  }
  else { 

    m_deltaPhi[mm] = -1.;
    m_mll[mm] = -1.;
    m_transvMass[mm] = -1.;

  }

  if ( m_channel[em] ) {
    float deltaPhiEPlusMuMinus = -1.0;
    float deltaPhiEMinusMuPlus = -1.0;
    if ( theElectronPlus > -1 && theMuonMinus > -1 ) {

      deltaPhiEPlusMuMinus = fabs( 180./TMath::Pi() *
				   m_p4ElectronPlus->Vect().DeltaPhi(m_p4MuonMinus->Vect()) );

      m_deltaPhi[em] = deltaPhiEPlusMuMinus;

      dilepPt.SetXYZ( m_p4ElectronPlus->Vect().X()+m_p4MuonMinus->Vect().X(),
		      m_p4ElectronPlus->Vect().Y()+m_p4MuonMinus->Vect().Y(),
		      0.0 );
      m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );

    }
    if ( theElectronMinus > -1 && theMuonPlus > -1 ) {

      deltaPhiEMinusMuPlus = fabs( 180./TMath::Pi() *
				   m_p4ElectronMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()) );

      m_deltaPhi[em] = deltaPhiEMinusMuPlus;

      dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		      m_p4ElectronMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		      0.0 );
      m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );

    }
    if ( theElectronPlus > -1 && theMuonMinus > -1 &&
	 theElectronMinus > -1 && theMuonPlus > -1) {
      
      if ( deltaPhiEPlusMuMinus < deltaPhiEMinusMuPlus ) {
	m_deltaPhi[em] = deltaPhiEPlusMuMinus;
	dilepPt.SetXYZ( m_p4ElectronPlus->Vect().X()+m_p4MuonMinus->Vect().X(),
			m_p4ElectronPlus->Vect().Y()+m_p4MuonMinus->Vect().Y(),
			0.0 );
	m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
      }
      else {
	m_deltaPhi[em] = deltaPhiEMinusMuPlus;
	dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
			m_p4ElectronMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
			0.0 );
	m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
      }

    }
  }

  else {

    m_deltaPhi[em] = -1.;
    m_transvMass[em] = -1.;

  }

  m_mll[em] = 50.; // const value inside the interval, always accepted


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

bool HiggsSelection::jetVeto() {
  // first check that kinematics has been set
  assert(m_p4ElectronPlus->Vect().Mag() && m_p4ElectronMinus->Vect().Mag());
  bool foundJet=false;
  float maxPtJet=0.;
  for(int j=0;j<nJet;j++){
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    // check if the electron or the positron falls into the jet
    // (the electron IS a jet)
    std::vector<float> deltaR;
    deltaR.push_back( p3Jet.DeltaR(m_p4ElectronMinus->Vect()) );
    deltaR.push_back( p3Jet.DeltaR(m_p4ElectronPlus->Vect()) );
    if(_selection->getSwitch("jetConeWidth") && _selection->passCut("jetConeWidth",deltaR[0])) continue;
    if(_selection->getSwitch("jetConeWidth") && _selection->passCut("jetConeWidth",deltaR[1])) continue;

    if(etJet[j]>maxPtJet) maxPtJet=etJet[j];

    // jet veto
    if(_selection->getSwitch("etaJetAcc") && !_selection->passCut("etaJetAcc",etaJet[j])) continue;
    if((_selection->getSwitch("etJetLowAcc") && !_selection->passCut("etJetLowAcc",etJet[j]) ) ||
       ((_selection->getSwitch("etJetHighAcc") && !_selection->passCut("etJetHighAcc",etJet[j]) &&
	 (_selection->getSwitch("jetAlpha") && !_selection->passCut("jetAlpha",alphaJet[j]))))
       ) continue;
    foundJet=true;
  }
  etHighestJet->Fill(maxPtJet);
  return foundJet;
}

bool HiggsSelection::preselJetVeto() {
  // first check that kinematics has been set
  assert(m_p4ElectronPlus->Vect().Mag() && m_p4ElectronMinus->Vect().Mag());
  int numJet=0;
  bool found2Jet=false;
  for(int j=0;j<nJet;j++){
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    std::vector<float> deltaR;
    deltaR.push_back(p3Jet.DeltaR( m_p4ElectronMinus->Vect()) );
    deltaR.push_back(p3Jet.DeltaR( m_p4ElectronPlus->Vect()) );
    if(_selection->getSwitch("jetConeWidth") && _selection->passCut("jetConeWidth",deltaR[0])) continue;
    if(_selection->getSwitch("jetConeWidth") && _selection->passCut("jetConeWidth",deltaR[1])) continue;
    
    // relaxed jet veto
    if(_selection->getSwitch("etaPresJetAcc") && !_selection->passCut("etaPresJetAcc",etaJet[j])) continue;
    if(_selection->getSwitch("etPresJetAcc")  && !_selection->passCut("etPresJetAcc",etJet[j])) continue;
    numJet++;
  }
  if (numJet > 2){ found2Jet = true; } 
  return found2Jet;
}

float HiggsSelection::Fisher(int eleIndex) {
  float fisher=0;
  // CMSSW_1_3_1 coefficients, obsolete!
//   if(eleClassEle[eleIndex]<100)
//     fisher = 42.0238-3.38943*s9s25Ele[eleIndex]-794.092*sqrt(covEtaEtaEle[eleIndex])-15.3449*latEle[eleIndex]-31.1032*a20Ele[eleIndex];
//   else
//     fisher = 27.2967+2.97453*s9s25Ele[eleIndex]-169.219*sqrt(covEtaEtaEle[eleIndex])-17.0445*latEle[eleIndex]-24.8542*a20Ele[eleIndex];
  return fisher;
}

void HiggsSelection::bookHistos() {
  // Generator 
  _monitorGenerator->book1D("highestPtGenEle","Highest p_{T} generated electron (GeV/c)",_highestPtGen,50,0,165);
  _monitorGenerator->book1D("lowestPtGenEle","Lowest p_{T} generated electron (GeV/c)",_lowestPtGen,50,0,60);
  _monitorGenerator->book1D("genHiggsPt","Higgs p_{T} (GeV/c)",_genHiggsPt,80,0,200);
  _monitorGenerator->book1D("nGenJets","number of generated Jets",_nGenJet,100,0,100);

  _monitorGenJets->book1D("et","generated jet tranverse energy (GeV)",etGenJet,50,0,300,"All+Fake+Best");

  // Event quantities - FIXME
//   _monitorEventAfterReco->book1D("nEle","number of reconstructed electrons",_nEle,10,0,10);
//   _monitorEventAfterReco->book1D("nJets","number of reconstructed jets",_nJet,100,0,100);
//   _monitorEventAfterReco->book1D("deltaPhi","#Delta #phi of e^{+}e^{-} (degrees)",_deltaPhi,50,0.,180);
//   _monitorEventAfterReco->book1D("mll","e^{+}e^{-} invariant mass (GeV/c^{2})",_mll,50,0.,165.);
//   _monitorEventAfterReco->book1D("WWtrMass","W^{+}W^{-} transverse mass (GeV/c^{2})",_transvMass,50,0,250);
//   _monitorEventAfterReco->book1D("highestPtEle","Highest p_{T} electron (GeV/c)",_highestPt,50,0,165);
//   _monitorEventAfterReco->book1D("lowestPtEle","Lowest p_{T} electron (GeV/c)",_lowestPt,50,0,60);

//   _monitorEventAfterSelection->book1D("nEle","number of reconstructed electrons",_nEle,10,0,10);
//   _monitorEventAfterSelection->book1D("nJets","number of reconstructed jets",_nJet,100,0,100);
//   _monitorEventAfterSelection->book1D("deltaPhi","#Delta #phi of e^{+}e^{-} (degrees)",_deltaPhi,50,0.,180);
//   _monitorEventAfterSelection->book1D("mll","e^{+}e^{-} invariant mass (GeV/c^{2})",_mll,50,0.,50.);
//   _monitorEventAfterSelection->book1D("WWtrMass","W^{+}W^{-} transverse mass (GeV/c^{2})",_transvMass,50,0,250);
//   _monitorEventAfterSelection->book1D("highestPtEle","Highest p_{T} electron (GeV/c)",_highestPt,50,0,165);
//   _monitorEventAfterSelection->book1D("lowestPtEle","Lowest p_{T} electron (GeV/c)",_lowestPt,50,0,60);

  // Met quantities
  _monitorMet->book1D("metEt","Missing trensverse energy (GeV)",etMet,50,0,150,"All");
  _monitorMet->book1D("metPhi","Missing trensverse momentum #phi",phiMet,50,-TMath::Pi(),TMath::Pi(),"All");

  // Jet quantities
  _monitorJets->book1D("et","jet tranverse energy (GeV)",etJet,50,0,300,"All+Fake+Best");
  _monitorJets->book1D("eta","jet #eta",etaJet,50,-2.5,2.5,"All+Fake+Best");
  _monitorJets->book1D("alphaJet","jet #alpha",alphaJet,50,0.,3.,"All+Fake+Best");
  _monitorJets->book1D("emFracJet","jet e.m. energy fraction",emFracJet,50,0.,1.,"All+Fake+Best");
  _monitorJets->book1D("hadFracJet","jet hadronic energy fraction",hadFracJet,50,0.,1.,"All+Fake+Best");


  // Electron quantities
  _monitorElectrons->book1D("energy","electron energy (GeV)",energyEle,50,0,150,"All+Fake+Best");
  _monitorElectrons->book1D("et","electron tranverse energy (GeV)",etEle,50,0,150,"All+Fake+Best");
  _monitorElectrons->book1D("eta","electron #eta",etaEle,50,-6.,6.,"All+Fake+Best");
  _monitorElectrons->book1D("phi","electron #phi",phiEle,50,-TMath::Pi(),TMath::Pi(),"All+Fake+Best");
  _monitorElectrons->book1D("s1s9","electron #sum 1/#sum 9",s1s9Ele,50,0.,1.,"All+Fake+Best");
  _monitorElectrons->book1D("s9s25","electron #sum 9/#sum 25",s9s25Ele,50,0.,1.,"All+Fake+Best");

  // Jet match fraction
  RecoJets_pt = new TH1F("RecoJets_pt","reconstructed jet p_{T}",25,0.0,50.0);
  MatchedJets_pt = new TH1F("MatchedJets_pt","reconstructed matched jet p_{T}",25,0.0,50.0);
  etHighestJet = new TH1F("etHighestJet","most energetic jet p_{T} (GeV/c)",100,5.0,200.0);

}

void HiggsSelection::estimateJetMatch(float ptmin) {

  int jetNotEle=0;
  _bestJets->clear();
  _excludedJets->clear();

  for(int recojet=0;recojet<nJet;recojet++) {
    // fill the denominator: all reco jets
    if(etJet[recojet]<20.0) _excludedJets->push_back(recojet);
    if(etJet[recojet]>ptmin) {
      RecoJets_pt->Fill(etJet[recojet]);
    }
    else continue;

    TVector3 pRecoJet(pxJet[recojet],pyJet[recojet],pzJet[recojet]);

    // fill the numerator: only the matched jets
    // matching is defined as DeltaR(reco-gen)<0.3
    for(int genjet=0;genjet<nGenJet;genjet++) {
      TVector3 pGenJet(pxGenJet[genjet],pyGenJet[genjet],pzGenJet[genjet]);
      if(etGenJet[genjet]>ptmin && pRecoJet.DeltaR(pGenJet)<0.3) {
	MatchedJets_pt->Fill(etJet[recojet]);
	// the matched jets are best jets, the rest are "fake" 
	// fake jets include the electron-matched jets
	_bestJets->push_back(recojet);
	break;
      }
    }
  }

}
