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
  _selection->addSwitch("presel");
  _selection->addSwitch("jetVeto");
  _selection->addSwitch("classDepEleId");
  _selection->addSwitch("apply_kFactor");
  _selection->addCut("nRecoEle");
  _selection->addCut("etaEleAcc");
  _selection->addCut("ptEleAcc");
  _selection->addCut("trackerPtSum");
  _selection->addCut("hcalPtSum");
  _selection->addCut("etaPresJetAcc");
  _selection->addCut("etPresJetAcc");
  _selection->addCut("maxPtElePres");
  _selection->addCut("minPtElePres");
  _selection->addCut("METPres");
  _selection->addCut("eleInvMassPres");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetLowAcc");
  _selection->addCut("etJetHighAcc");
  _selection->addCut("alphaJet");
  _selection->addCut("maxPtEle");
  _selection->addCut("minPtEle");
  _selection->addCut("detaLeptons");
  _selection->addCut("MET");
  _selection->addCut("deltaPhi");
  _selection->addCut("eleInvMass");

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
  _counter.AddVar("nRecoEle");
  _counter.AddVar("twoGoodRec");
  _counter.AddVar("eleID");
  _counter.AddVar("trackerIsol");
  _counter.AddVar("hcalIsol");
  _counter.AddVar("presel");
  _counter.AddVar("jetVeto");
  _counter.AddVar("MET");
  _counter.AddVar("deltaPhi");
  _counter.AddVar("eleInvMass");
  _counter.AddVar("maxPtEle");
  _counter.AddVar("minPtEle");
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


  // the kinematic vectors
  _p3Ele = new TVector3(0.,0.,0.);
  _p3Pos = new TVector3(0.,0.,0.);
  _p3Met = new TVector3(0.,0.,0.);
  
  _p4Ele = new TLorentzVector(0.,0.,0.,0.);
  _p4Pos = new TLorentzVector(0.,0.,0.,0.);

  // the Monitoring Histograms
   _monitorGenerator = new Monitor(0);
   _monitorEventAfterSelection = new Monitor(0);
   _monitorEventAfterReco = new Monitor(0);
   _monitorMet = new Monitor(&nMet);
   _bestElectrons = new std::vector<int>;
   _monitorElectrons = new Monitor(&nEle,_bestElectrons);
   _bestJets = new std::vector<int>;
   _excludedJets = new std::vector<int>;
   _monitorJets = new Monitor(&nJet,_bestJets);
   _monitorJets->ExcludeElements(_excludedJets);
   _bestGenJets = new std::vector<int>;
   _monitorGenJets = new Monitor(&nGenJet,_bestGenJets);

}

HiggsSelection::~HiggsSelection(){
  delete _p3Ele;
  delete _p3Pos;
  delete _p3Met;
  delete _p4Ele;
  delete _p4Pos;
  delete _monitorGenerator;
  delete _monitorEventAfterSelection;
  delete _monitorEventAfterReco;
  delete _monitorMet;
  delete _bestElectrons;
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

  float met, deltaPhi, transvMass, eleInvMass, maxPtEle, minPtEle, detaLeptons;

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


      if(_selection->getSwitch("nRecoEle") && 
	 !_selection->passCut("nRecoEle",nEle)) continue;   
      _counter.IncrVar("nRecoEle",weight);
      int theEle(getBestLeptonPair().first), thePos(getBestLeptonPair().second);


      if(theEle<0 || thePos<0) continue; 
      _counter.IncrVar("twoGoodRec",weight);
      setKinematics(theEle,thePos);
      addVariables();
      estimateJetMatch(0.0);
      _monitorJets->Fill(weight);
      _monitorEventAfterReco->Fill(weight);


      if(!isEleID(theEle) || !isEleID(thePos)) continue; 
      _counter.IncrVar("eleID",weight);


      if(_selection->getSwitch("trackerPtSum") && 
	 ( !_selection->passCut("trackerPtSum",eleTrackerIso_sumPtEle[theEle]) || 
	   !_selection->passCut("trackerPtSum",eleTrackerIso_sumPtEle[thePos]) ) ) continue; 
      _counter.IncrVar("trackerIsol",weight);


      if(_selection->getSwitch("hcalPtSum") &&
	 ( !_selection->passCut("hcalPtSum",eleCaloIso_sumPtEle[theEle]) || 
	   !_selection->passCut("hcalPtSum",eleCaloIso_sumPtEle[thePos]) ) ) continue; 
      _counter.IncrVar("hcalIsol",weight);


      met=etMet[0];
      eleInvMass = (*_p4Ele+*_p4Pos).M();
      minPtEle = _minPt;
      maxPtEle = _maxPt;
      if(_selection->getSwitch("presel") && 
	 ( preselJetVeto() ||
	   !_selection->passCut("maxPtElePres",maxPtEle) ||
	   !_selection->passCut("minPtElePres",minPtEle) ||
	   !_selection->passCut("METPres",met) ||
	   !_selection->passCut("eleInvMassPres",eleInvMass)
	   )
	 ) continue; 
      _counter.IncrVar("presel",weight);

      if(_selection->getSwitch("jetVeto") && jetVeto()) continue; 
      _counter.IncrVar("jetVeto",weight);
      
      if(_selection->getSwitch("MET") && !_selection->passCut("MET",met)) continue; 
      _counter.IncrVar("MET",weight);

      deltaPhi = fabs(180./TMath::Pi()*_p3Ele->DeltaPhi(*_p3Pos));
      transvMass = _transvMass[0];
      if(_selection->getSwitch("deltaPhi") && 
	 !_selection->passCut("deltaPhi",deltaPhi) ) continue; 
      _counter.IncrVar("deltaPhi",weight);

      if(_selection->getSwitch("eleInvMass") && 
	 !_selection->passCut("eleInvMass",eleInvMass)) continue; 
      _counter.IncrVar("eleInvMass",weight);

      if(_selection->getSwitch("maxPtEle") && !_selection->passCut("maxPtEle",maxPtEle)) continue; 
      _counter.IncrVar("maxPtEle",weight);

      if(_selection->getSwitch("minPtEle") && !_selection->passCut("minPtEle",minPtEle)) continue; 
      _counter.IncrVar("minPtEle",weight);

      detaLeptons = etaEle[theEle]-etaEle[thePos];
      if(_selection->getSwitch("detaLeptons") && 
	 !_selection->passCut("detaLeptons",detaLeptons) ) continue; 
      _counter.IncrVar("detaLeptons",weight);
      _counter.IncrVar("final",weight);
      
      _monitorEventAfterSelection->Fill(weight);
      _monitorMet->Fill(weight);
      _monitorElectrons->Fill(weight);
      _monitorGenJets->Fill(weight);

      // dumping the reco variables in a tree - after all the cuts
      myOutTree -> fillAll(met, deltaPhi, transvMass, eleInvMass, maxPtEle, minPtEle, detaLeptons);
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
  _counter.Draw("nRecoEle","trigger");
  _counter.Draw("twoGoodRec","nRecoEle");
  _counter.Draw("eleID","twoGoodRec");
  _counter.Draw("trackerIsol","eleID");
  _counter.Draw("hcalIsol","trackerIsol");
  _counter.Draw("presel","hcalIsol");
  _counter.Draw("jetVeto","presel");
  _counter.Draw("MET","jetVeto");
  _counter.Draw("deltaPhi","MET");
  _counter.Draw("eleInvMass","deltaPhi");
  _counter.Draw("maxPtEle","eleInvMass");
  _counter.Draw("minPtEle","maxPtEle");
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

std::pair<int,int> HiggsSelection::getBestLeptonPair() {
  int theEle=-1;
  int thePos=-1;
  float maxPtEle=-1000.;
  float maxPtPos=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nEle;i++) {
    if(_selection->getSwitch("etaEleAcc") && !_selection->passCut("etaEleAcc",etaEle[i]) ) continue;
    TVector3 pEle(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pEle.Pt();
    if(_selection->getSwitch("ptEleAcc") && !_selection->passCut("ptEleAcc",thisPt) ) continue;
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtPos){ maxPtPos = thisPt; thePos = i; }
    if (thisCharge < 0 && thisPt> maxPtEle){ maxPtEle = thisPt; theEle = i; }
  }
  _bestElectrons->clear();
  _bestElectrons->push_back(theEle);  _bestElectrons->push_back(thePos); 
  return make_pair(theEle,thePos);
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



void HiggsSelection::addVariables() {


  // --- Variables for electrons ---
  _deltaPhi[0]=fabs(180./TMath::Pi()* _p3Ele->DeltaPhi(*_p3Pos));
  _mll[0]=(*_p4Ele+*_p4Pos).M();
  for(int i=0;i<nEle;i++) {
    _eOverP[i]=ecalEle[i]/energyEle[i];
  }

  _highestPt[0]=_maxPt;
  _lowestPt[0]=_minPt;
  _nEle[0]=(float)nEle;

  // now find the electrons inside the jets to get alpha
  int ieleRecoAsJet=0;
  for(int j=0;j<nJet;j++) {
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    if(_selection->passCut("jetConeWidth",(p3Jet.DeltaR(*_p3Pos))) ||
       _selection->passCut("jetConeWidth",(p3Jet.DeltaR(*_p3Ele))) ) {
      _alphaEle[ieleRecoAsJet]=alphaJet[j];
      _emFracEle[ieleRecoAsJet]=emFracJet[j];
      _hadFracEle[ieleRecoAsJet]=hadFracJet[j];
      ieleRecoAsJet++;
    }
  }


  TVector3 dilepPt(_p3Ele->X()+_p3Pos->X(),_p3Ele->Y()+_p3Pos->Y());
  _transvMass[0]=sqrt(2*dilepPt.Mag()*_p3Met->Mag()*(1-cos(dilepPt.Angle(*_p3Met)) ) );




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



  // --- variables for jets ---
  _nJet[0]=(float)nJet;

  int theHighestPtJet=0;
  float maxJetPt=-100;
  for(int i=0;i<nJet;i++) {
    TVector3 p3GenJet(pxGenJet[i],pyGenJet[i],pzGenJet[i]);
    if(etGenJet[i]>maxJetPt) {
      maxJetPt=etGenJet[i];
      theHighestPtJet=i;
    }
  }

  _bestGenJets->clear();
  _bestGenJets->push_back(theHighestPtJet);
  _nGenJet[0]=(float)nGenJet;


}

void HiggsSelection::setKinematics(int theEle, int thePos) {
  // three-vectors
  _p3Ele->SetXYZ(pxEle[theEle],pyEle[theEle],pzEle[theEle]);
  _p3Pos->SetXYZ(pxEle[thePos],pyEle[thePos],pzEle[thePos]);
  _p3Met->SetXYZ(pxMet[0],pyMet[0],pzMet[0]);

  // four-vectors
  _p4Ele->SetXYZT(pxEle[theEle],pyEle[theEle],pzEle[theEle],energyEle[theEle]);
  _p4Pos->SetXYZT(pxEle[thePos],pyEle[thePos],pzEle[thePos],energyEle[thePos]);

  _maxPt=TMath::Max(_p3Ele->Pt(),_p3Pos->Pt());
  _minPt=TMath::Min(_p3Ele->Pt(),_p3Pos->Pt());

}

bool HiggsSelection::jetVeto() {
  // first check that kinematics has been set
  assert(_p3Ele->Mag() && _p3Pos->Mag());
  bool foundJet=false;
  float maxPtJet=0.;
  for(int j=0;j<nJet;j++){
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    // check if the electron or the positron falls into the jet
    // (the electron IS a jet)
    std::vector<float> deltaR;
    deltaR.push_back(p3Jet.DeltaR(*_p3Ele));
    deltaR.push_back(p3Jet.DeltaR(*_p3Pos));
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
  assert(_p3Ele->Mag() && _p3Pos->Mag());
  int numJet=0;
  bool found2Jet=false;
  for(int j=0;j<nJet;j++){
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);
    std::vector<float> deltaR;
    deltaR.push_back(p3Jet.DeltaR(*_p3Ele));
    deltaR.push_back(p3Jet.DeltaR(*_p3Pos));
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
  float fisher;
  if(eleClassEle[eleIndex]<100)
    fisher = 42.0238-3.38943*s9s25Ele[eleIndex]-794.092*sqrt(covEtaEtaEle[eleIndex])-15.3449*latEle[eleIndex]-31.1032*a20Ele[eleIndex];
  else
    fisher = 27.2967+2.97453*s9s25Ele[eleIndex]-169.219*sqrt(covEtaEtaEle[eleIndex])-17.0445*latEle[eleIndex]-24.8542*a20Ele[eleIndex];
  return fisher;
}

void HiggsSelection::bookHistos() {
  // Generator 
  _monitorGenerator->book1D("highestPtGenEle","Highest p_{T} generated electron (GeV/c)",_highestPtGen,50,0,165);
  _monitorGenerator->book1D("lowestPtGenEle","Lowest p_{T} generated electron (GeV/c)",_lowestPtGen,50,0,60);
  _monitorGenerator->book1D("genHiggsPt","Higgs p_{T} (GeV/c)",_genHiggsPt,80,0,200);
  _monitorGenerator->book1D("nGenJets","number of generated Jets",_nGenJet,100,0,100);

  _monitorGenJets->book1D("et","generated jet tranverse energy (GeV)",etGenJet,50,0,300,"All+Fake+Best");

  // Event quantities
  _monitorEventAfterReco->book1D("nEle","number of reconstructed electrons",_nEle,10,0,10);
  _monitorEventAfterReco->book1D("nJets","number of reconstructed jets",_nJet,100,0,100);
  _monitorEventAfterReco->book1D("deltaPhi","#Delta #phi of e^{+}e^{-} (degrees)",_deltaPhi,50,0.,180);
  _monitorEventAfterReco->book1D("mll","e^{+}e^{-} invariant mass (GeV/c^{2})",_mll,50,0.,165.);
  _monitorEventAfterReco->book1D("WWtrMass","W^{+}W^{-} transverse mass (GeV/c^{2})",_transvMass,50,0,250);
  _monitorEventAfterReco->book1D("highestPtEle","Highest p_{T} electron (GeV/c)",_highestPt,50,0,165);
  _monitorEventAfterReco->book1D("lowestPtEle","Lowest p_{T} electron (GeV/c)",_lowestPt,50,0,60);

  _monitorEventAfterSelection->book1D("nEle","number of reconstructed electrons",_nEle,10,0,10);
  _monitorEventAfterSelection->book1D("nJets","number of reconstructed jets",_nJet,100,0,100);
  _monitorEventAfterSelection->book1D("deltaPhi","#Delta #phi of e^{+}e^{-} (degrees)",_deltaPhi,50,0.,180);
  _monitorEventAfterSelection->book1D("mll","e^{+}e^{-} invariant mass (GeV/c^{2})",_mll,50,0.,50.);
  _monitorEventAfterSelection->book1D("WWtrMass","W^{+}W^{-} transverse mass (GeV/c^{2})",_transvMass,50,0,250);
  _monitorEventAfterSelection->book1D("highestPtEle","Highest p_{T} electron (GeV/c)",_highestPt,50,0,165);
  _monitorEventAfterSelection->book1D("lowestPtEle","Lowest p_{T} electron (GeV/c)",_lowestPt,50,0,60);

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
  _monitorElectrons->book1D("eOverP","electron E/P",_eOverP,50,0.75,1.1,"All+Fake+Best");
  _monitorElectrons->book1D("alphaEle","electron #alpha",_alphaEle,50,0.,3.,"All+Fake+Best");
  _monitorElectrons->book1D("emFracEle","electron e.m. energy fraction",_emFracEle,50,0.,1.,"All+Fake+Best");
  _monitorElectrons->book1D("hadFracEle","electron hadronic energy fraction",_hadFracEle,50,0.,1.,"All+Fake+Best");


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
