#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"
#include "HiggsAnalysisTools/include/HiggsSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"

#include <iostream>
#include <string>

#include <TTree.h>

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


  // preselection efficiencies
  CommonHiggsPreselection.Configure(higgsConfigDir.c_str()); 
  
  // extra preselection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  std::string fileCuts  = higgsConfigDir + "2e2nuCuts.txt";
  std::string fileSwitches = higgsConfigDir + "2e2nuSwitches.txt";
  _addedPres = new Selection(fileCuts,fileSwitches);
  _addedPres->addCut("etaElectronAcc");
  _addedPres->addCut("ptElectronAcc");
  _addedPres->addCut("etaMuonAcc");
  _addedPres->addCut("ptMuonAcc");

  // selection efficiencies
  CutBasedHiggsSelectionEE.Configure(higgsConfigDir.c_str()); 
  CutBasedHiggsSelectionMM.Configure(higgsConfigDir.c_str()); 

  // extra selection efficiencies  - to be put here not to pass the full list of jets to the selection class
  _addedSel = new Selection(fileCuts,fileSwitches);
  _addedSel->addCut("jetConeWidth");
  _addedSel->addCut("etaJetAcc");
  _addedSel->addCut("etJetLowAcc");
  _addedSel->addCut("etJetHighAcc");
  _addedSel->addCut("alphaJet");
  
  // single electron efficiency
  EgammaCutBasedID.Configure("../EgammaAnalysisTools/config/tightEleId/"); 

  
  // kinematics
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
  _monitorMuons = new Monitor(&nMuon,_bestMuons);
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
  delete _monitorMuons;
  delete _bestJets;
  delete _excludedJets;
  delete _bestGenJets;
  delete _monitorGenJets;
  delete _monitorJets;
  delete _addedPres;
  delete _addedSel;
  myOutTreeEE -> save();
  myOutTreeMM -> save();
  myOutTreeEM -> save();

}

bool HiggsSelection::findMcTree(const char* processType) {

  _process = "UNDEFINED";
  _theGenEle = -1;
  _theGenPos = -1;
  
  // now we look for ee || mumu || emu
  
  // signal: 2e2nu
  if(strcmp(processType,"HtoWWto2e2nu")==0) {
    if( idMc[9]  == -11) _theGenEle = 9;
    if( idMc[11] ==  11) _theGenPos = 11;
    return (_theGenEle > -1 && _theGenPos > -1 );
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

float HiggsSelection::getkFactor(std::string process) {

  float weight = 1.;
  if((process.compare("H_VBF")==0) || (process.compare("H_gg")==0)) {
    weight = evtKfactor;   
  }
  else if(process.compare("WW")==0) {
    weight = evtMcAtNlo;   
  }
  return weight;
}

void HiggsSelection::Loop() {

  _verbose=false;
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
  _monitorMuons->setPath("Muons");   
  _monitorJets->setPath("Jets");
  _monitorGenJets->setPath("GenJets");
  _monitorEventAfterReco->setPath("EventAfterReco");
  _monitorEventAfterSelection->setPath("EventAfterSelection");
  
  // kinematics reduced tree
  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";
  std::string reducedTreeNameMM = _datasetName+"-datasetMM.root";
  std::string reducedTreeNameEM = _datasetName+"-datasetEM.root";
  myOutTreeEE = new RedHiggsTree(reducedTreeNameEE.c_str());
  myOutTreeMM = new RedHiggsTree(reducedTreeNameMM.c_str());
  myOutTreeEM = new RedHiggsTree(reducedTreeNameEM.c_str());
  
  float met, deltaPhi, transvMass; 
  float dileptonInvMass, maxPtEle, minPtEle, detaLeptons;
  
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // get the kFactor of the event
    float weight = getkFactor(_process);

    // look to the MC truth decay tree
    bool foundMcTree = findMcTree("HtoWWto2e2nu");

    // trigger
    Utils anaUtils;
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

    // get the best electrons, best muons  
    std::pair<int,int> theElectrons = getBestElectronPair();
    std::pair<int,int> theMuons = getBestMuonPair();
    int tbElectron(theElectrons.second), tbPositron(theElectrons.first);    
    int tbMuonPlus(theMuons.first),     tbMuonMinus(theMuons.second);
    theElectron  = tbElectron;
    thePositron  = tbPositron;
    theMuonPlus  = tbMuonPlus;
    theMuonMinus = tbMuonMinus;

    // reconstructed channel
    m_channel[ee] = false;     
    m_channel[mm] = false;
    m_channel[em] = false;
    if (theElectron > -1 && thePositron > -1)  m_channel[ee] = true;
    if (theMuonPlus > -1 && theMuonMinus > -1) m_channel[mm] = true;
    if ((theElectron > -1 && theMuonPlus > -1) || (thePositron > -1 && theMuonMinus > -1)) m_channel[em] = true;
    if (_verbose) {
      std::cout << "nEle = " << nEle << "\tnMuon = " << nMuon << std::endl;
      std::cout << "indices: " << theElectron << " " << thePositron << " " << theMuonMinus << " " << theMuonPlus << std::endl;
      std::cout << "chargeEle = " << chargeEle[theElectron] << "\tchargePos = " << chargeEle[thePositron] 
	 << "\tchargeMuonMinus = " << chargeMuon[theMuonMinus] << "\tchargeMuonPlus = " << chargeMuon[theMuonPlus] << std::endl;
      std::cout << "ee = " << m_channel[ee] << "\tmm = " << m_channel[mm] << "\temu = " << m_channel[em] << std::endl; 
    }

    // kinematics at preselection level 
    // only few variables since here we still do not know if two good leptons are reconstructed
    setPreselKinematics();

    // setting preselections
    CommonHiggsPreselection.SetWeight(weight);
    CommonHiggsPreselection.SetMcTruth(foundMcTree);
    CommonHiggsPreselection.SetHLT(passedHLT);
    CommonHiggsPreselection.SetNele(nEle);
    CommonHiggsPreselection.SetNmuon(nMuon);
    CommonHiggsPreselection.SetIsEE(m_channel[ee]);
    CommonHiggsPreselection.SetIsEM(m_channel[em]);
    CommonHiggsPreselection.SetIsMM(m_channel[mm]);
    CommonHiggsPreselection.SetHighElePt(hardestElectronPt);
    CommonHiggsPreselection.SetLowElePt(slowestElectronPt);
    CommonHiggsPreselection.SetHighMuonPt(hardestMuonPt);
    CommonHiggsPreselection.SetLowMuonPt(slowestMuonPt);
    CommonHiggsPreselection.SetMet(etMet[0]);
    CommonHiggsPreselection.SetMllEE(m_mll[ee]);
    CommonHiggsPreselection.SetMllEM(m_mll[em]);
    CommonHiggsPreselection.SetMllMM(m_mll[mm]);
    CommonHiggsPreselection.SetCommonPres(evtPresel);

    // preselection histos
    _monitorGenerator->Fill(weight);
    estimateJetMatch(0.0);
    _monitorJets->Fill(weight);
    _monitorEventAfterReco->Fill(weight);
    
    // did we pass preselections?
    bool isPreselections = CommonHiggsPreselection.output();    
    if( !isPreselections ) continue;

    // kinematics after preselections
    // now we know we have two good reconstructed leptons
    setKinematics();

    // ----------------------- selection ----------------------------
    // electron ID (true by default - studied only if ee or emu channel)
    bool theElectronID = true;
    bool thePositronID = true;
    if (theElectron > -1) theElectronID = isEleID(theElectron);
    if (thePositron > -1) thePositronID = isEleID(thePositron);
    // loose egamma electron ID
//     if (theElectron > -1) theElectronID = eleIdCutBasedEle[theElectron];
//     if (thePositron > -1) thePositronID = eleIdCutBasedEle[thePositron];

    // extra tracker isolation for electrons
    float theEleTrackerPtSum = 0.;
    float thePosTrackerPtSum = 0.;
    if (theElectron > -1) theEleTrackerPtSum = eleTrackerIso_sumPtEle[theElectron];
    if (thePositron > -1) thePosTrackerPtSum = eleTrackerIso_sumPtEle[thePositron];

    // hcal isolation for electrons
    float theEleCaloPtSum = 0.;
    float thePosCaloPtSum = 0.;
    if (theElectron > -1) theEleCaloPtSum = eleCaloIso_sumPtEle[theElectron];
    if (thePositron > -1) thePosCaloPtSum = eleCaloIso_sumPtEle[thePositron];

    // jet veto
    bool passedJetVeto = jetVeto();

    // kine variables
    float theDeltaPhiEE, theInvMassEE, theTransvMassEE, theDetaLeptonsEE = 0.;
    float theDeltaPhiMM, theInvMassMM, theTransvMassMM, theDetaLeptonsMM = 0.;
    float theDeltaPhiEM, theInvMassEM, theTransvMassEM, theDetaLeptonsEM = 0.;


    // ---------------------------------------
    // ee final state
    if (m_channel[ee]){
      theDeltaPhiEE    = m_deltaPhi[ee];
      theInvMassEE     = m_mll[ee];
      theDetaLeptonsEE = etaEle[theElectron]-etaEle[thePositron];
      theTransvMassEE  = m_transvMass[ee];

      // selections
      CutBasedHiggsSelectionEE.SetWeight(weight);
      CutBasedHiggsSelectionEE.SetHighElePt(hardestElectronPt);
      CutBasedHiggsSelectionEE.SetLowElePt(slowestElectronPt);
      CutBasedHiggsSelectionEE.SetElectronId(theElectronID);
      CutBasedHiggsSelectionEE.SetPositronId(thePositronID);
      CutBasedHiggsSelectionEE.SetEleTrackerPtSum(theEleTrackerPtSum);
      CutBasedHiggsSelectionEE.SetPosTrackerPtSum(thePosTrackerPtSum);
      CutBasedHiggsSelectionEE.SetEleCaloPtSum(theEleCaloPtSum);
      CutBasedHiggsSelectionEE.SetPosCaloPtSum(thePosCaloPtSum);
      CutBasedHiggsSelectionEE.SetJetVeto(passedJetVeto);
      CutBasedHiggsSelectionEE.SetMet(etMet[0]);					
      CutBasedHiggsSelectionEE.SetDeltaPhi(theDeltaPhiEE);
      CutBasedHiggsSelectionEE.SetInvMass(theInvMassEE);
      CutBasedHiggsSelectionEE.SetDetaLeptons(theDetaLeptonsEE);
      bool isSelectedEE = CutBasedHiggsSelectionEE.output();    

      // dumping final tree
      if(isSelectedEE) myOutTreeEE -> fillAll(etMet[0], theDeltaPhiEE, theTransvMassEE, theInvMassEE, hardestElectronPt, slowestElectronPt, theDetaLeptonsEE);
    }


    // ---------------------------------------
    // mm final state
    if (m_channel[mm]){
      theDeltaPhiMM    = m_deltaPhi[mm];
      theInvMassMM     = m_mll[mm];
      theDetaLeptonsMM = etaEle[theMuonMinus]-etaEle[theMuonPlus];
      theTransvMassMM  = m_transvMass[mm];

      // selections
      CutBasedHiggsSelectionMM.SetWeight(weight);
      CutBasedHiggsSelectionMM.SetHighElePt(hardestMuonPt);
      CutBasedHiggsSelectionMM.SetLowElePt(slowestMuonPt);
      CutBasedHiggsSelectionMM.SetElectronId(true);
      CutBasedHiggsSelectionMM.SetPositronId(true);
      CutBasedHiggsSelectionMM.SetEleTrackerPtSum(0);
      CutBasedHiggsSelectionMM.SetPosTrackerPtSum(0);
      CutBasedHiggsSelectionMM.SetEleCaloPtSum(0);
      CutBasedHiggsSelectionMM.SetPosCaloPtSum(0);
      CutBasedHiggsSelectionMM.SetJetVeto(passedJetVeto);
      CutBasedHiggsSelectionMM.SetMet(etMet[0]);					
      CutBasedHiggsSelectionMM.SetDeltaPhi(theDeltaPhiMM);
      CutBasedHiggsSelectionMM.SetInvMass(theInvMassMM);
      CutBasedHiggsSelectionMM.SetDetaLeptons(theDetaLeptonsMM);
      bool isSelectedMM = CutBasedHiggsSelectionMM.output();    

      // dumping final tree
      if(isSelectedMM) myOutTreeMM -> fillAll(etMet[0], theDeltaPhiMM, theTransvMassMM, theInvMassMM, hardestMuonPt, slowestMuonPt, theDetaLeptonsMM);
    }
    
    // ancora da finire
    if (m_channel[em]){
      theDeltaPhiEM    = m_deltaPhi[em];
      theInvMassEM     = m_mll[em];
      theTransvMassEM  = m_transvMass[em];
      if(theElectron>-1 && theMuonPlus>-1){
	theDetaLeptonsEM = etaEle[theElectron]-etaEle[theMuonPlus];
      }
      if(thePositron>-1 && theMuonMinus>-1){
	theDetaLeptonsEM = etaEle[thePositron]-etaEle[theMuonMinus];
      }
    }

    /*
    // ancora da finire -- sistema --- 
    // monitor element
    _monitorEventAfterSelection->Fill(weight);
    _monitorMet->Fill(weight);
    _monitorElectrons->Fill(weight);
    _monitorMuons->Fill(weight);   
    _monitorGenJets->Fill(weight);
    */

    myOutTreeEE -> store(); 
    myOutTreeMM -> store(); 
    myOutTreeEM -> store(); 
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
  _monitorMuons->WritePs("muons.ps");
  _monitorMuons->WriteRoot(file);
  _monitorJets->WritePs("jets.ps");
  _monitorJets->WriteRoot(file);
  _monitorGenJets->WritePs("jets.ps");
  _monitorGenJets->WriteRoot(file);
  file->Close();  
}

void HiggsSelection::displayEfficiencies() {

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Common preselections: " << std::endl;
  CommonHiggsPreselection.diplayEfficiencies();

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EE selections: " << std::endl;
  CutBasedHiggsSelectionEE.diplayEfficiencies();

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full MM selections: " << std::endl;
  CutBasedHiggsSelectionMM.diplayEfficiencies();

  EgammaCutBasedID.diplayEfficiencies();

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
    if(_addedPres->getSwitch("etaElectronAcc") && !_addedPres->passCut("etaElectronAcc",etaEle[i]) ) continue;
    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();
    if(_addedPres->getSwitch("ptElectronAcc") && !_addedPres->passCut("ptElectronAcc",thisPt) ) continue;
    // fixme: (in the future) put here full electron ID not to loose efficiency I think
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
    if(_addedPres->getSwitch("etaMuonAcc") && !_addedPres->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    TVector3 pLepton(pxMuon[i],pyMuon[i],pzMuon[i]);
    float thisPt=pLepton.Pt();
    if(_addedPres->getSwitch("ptMuonAcc") && !_addedPres->passCut("ptMuonAcc",thisPt) ) continue;
    // fixme: (in the future) put here full electron ID not to loose efficiency I think
    float thisCharge = chargeMuon[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestMuons->clear();
  _bestMuons->push_back(theLep1);  _bestMuons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

bool HiggsSelection::isEleID(int eleIndex) {

  TVector3 pTrkAtOuter(pxAtOuterEle[eleIndex],pyAtOuterEle[eleIndex],pzAtOuterEle[eleIndex]);

  EgammaCutBasedID.SetHOverE( eleHoEEle[eleIndex] );
  EgammaCutBasedID.SetS9S25( s9s25Ele[eleIndex] );
  EgammaCutBasedID.SetDEta( eleDeltaEtaAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiIn( eleDeltaPhiAtVtxEle[eleIndex] );
  EgammaCutBasedID.SetDPhiOut( eleDeltaPhiAtCaloEle[eleIndex] );
  EgammaCutBasedID.SetInvEminusInvP( 1./eleCaloCorrEEle[eleIndex]-1./eleTrackerPEle[eleIndex] );
  EgammaCutBasedID.SetBremFraction( fabs(eleTrackerPEle[eleIndex]-pTrkAtOuter.Mag())/eleTrackerPEle[eleIndex] );
  EgammaCutBasedID.SetSigmaEtaEta( sqrt(covEtaEtaEle[eleIndex]) );
  EgammaCutBasedID.SetSigmaPhiPhi( sqrt(covPhiPhiEle[eleIndex]) );
  EgammaCutBasedID.SetEOverPout( eleCorrEoPoutEle[eleIndex] );
  EgammaCutBasedID.SetEOverPin( eleCorrEoPEle[eleIndex] );
  EgammaCutBasedID.SetElectronClass ( eleClassEle[eleIndex] );
  EgammaCutBasedID.SetEgammaCutBasedID ( eleIdCutBasedEle[eleIndex] );
  EgammaCutBasedID.SetLikelihood( eleLikelihoodEle[eleIndex] );

  bool isIdentified = EgammaCutBasedID.output();

  return isIdentified;
}



void HiggsSelection::setPreselKinematics() {

  // highest and lowest pt for electrons
  // + mll if 2 good electrons are reconstructed
  if (thePositron > -1 && theElectron > -1) {
    hardestElectronPt = TMath::Max(etEle[theElectron],etEle[thePositron]);
    slowestElectronPt = TMath::Min(etEle[theElectron],etEle[thePositron]);
    m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
    m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
    m_mll[ee] = (*m_p4ElectronMinus + *m_p4ElectronPlus).M();
  }
  if (thePositron <= -1 && theElectron > -1) {
    hardestElectronPt = etEle[theElectron];
    slowestElectronPt = etEle[theElectron];
    m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
    m_mll[ee] = -999.;
  }
  if (thePositron > -1 && theElectron <= -1) {
    hardestElectronPt = etEle[thePositron];
    slowestElectronPt = etEle[thePositron];
    m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
    m_mll[ee] = -999.;
  }
  if (theElectron <= -1 && thePositron <= -1) {
    hardestElectronPt = -999.;
    slowestElectronPt = -999.;
    m_mll[ee] = -999.;
  }

  // highest and lowest pt for muons
  // + mll if 2 good muons are reconstructed  
  if (theMuonPlus > -1 && theMuonMinus > -1){ 
    hardestMuonPt = TMath::Max(etMuon[theMuonPlus],etMuon[theMuonMinus]);
    slowestMuonPt = TMath::Min(etMuon[theMuonPlus],etMuon[theMuonMinus]);
    m_p4MuonMinus -> SetXYZT(pxMuon[theMuonMinus],pyMuon[theMuonMinus],pzMuon[theMuonMinus],energyMuon[theMuonMinus]);
    m_p4MuonPlus  -> SetXYZT(pxMuon[theMuonPlus],pyMuon[theMuonPlus],pzMuon[theMuonPlus],energyMuon[theMuonPlus]);      
    m_mll[mm] = (*m_p4MuonMinus + *m_p4MuonPlus).M();
  }
  if (theMuonMinus > -1 && theMuonPlus <= -1) {
    hardestMuonPt = etMuon[theMuonMinus];
    slowestMuonPt = etMuon[theMuonMinus];
    m_p4MuonMinus -> SetXYZT(pxMuon[theMuonMinus],pyMuon[theMuonMinus],pzMuon[theMuonMinus],energyMuon[theMuonMinus]);
    m_mll[mm]     = -999.;
  }
  if (theMuonPlus > -1 && theMuonMinus <= -1) {
    hardestMuonPt = etMuon[theMuonPlus];
    slowestMuonPt = etMuon[theMuonPlus];
    m_p4MuonPlus  -> SetXYZT(pxMuon[theMuonPlus],pyMuon[theMuonPlus],pzMuon[theMuonPlus],energyMuon[theMuonPlus]);      
    m_mll[mm]     = -999.;
  }
  if (theMuonMinus <= -1 && theMuonPlus <= -1) {
    hardestMuonPt = -999.;
    slowestMuonPt = -999.;
    m_mll[mm]     = -999.;
  }
  
  // mixed channel mll --> always passed
  m_mll[em]       = 50.;
  
  // MET
  m_p4MET->SetXYZT(pxMet[0],pyMet[0],pzMet[0],energyMet[0]); 
}



void HiggsSelection::setKinematics( ) {

  // electron variables used for ele quality in jet veto 
  m_HoEElectronMinus     = eleHoEEle[theElectron];
  m_HoEElectronPlus      = eleHoEEle[thePositron];
  m_CaloEneElectronMinus = eleCaloCorrEEle[theElectron];
  m_CaloEneElectronPlus  = eleCaloCorrEEle[thePositron];

  // compute delta Phi in degrees, di-lepton invariant mass, transverse mass
  TVector3 dilepPt;
  if ( m_channel[ee] ) {
    m_deltaPhi[ee] = fabs(180./TMath::Pi() * m_p4ElectronMinus->Vect().DeltaPhi(m_p4ElectronPlus->Vect()));
    dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4ElectronPlus->Vect().X(),
		    m_p4ElectronMinus->Vect().Y()+m_p4ElectronPlus->Vect().Y(),
		    0.0 );
    m_transvMass[ee]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
  }
  else {    
    m_deltaPhi[ee]   = -1.;
    m_transvMass[ee] = -1.;
  }

  if ( m_channel[mm] ) {    
    m_deltaPhi[mm] = fabs(180./TMath::Pi() * m_p4MuonMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()));
    dilepPt.SetXYZ( m_p4MuonMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		    m_p4MuonMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		    0.0 );
    m_transvMass[mm]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect()))));
  }
  else { 
    m_deltaPhi[mm]   = -1.;
    m_transvMass[mm] = -1.;
  }
  
  if ( m_channel[em] ) {
    float deltaPhiEPlusMuMinus = -1.0;
    float deltaPhiEMinusMuPlus = -1.0;
    float dilepPtEPlusMuMinus  = -1.0;
    float dilepPtEMinusMuPlus  = -1.0;

    if ( thePositron > -1 && theMuonMinus > -1 ) {
      deltaPhiEPlusMuMinus = fabs(180./TMath::Pi() * m_p4ElectronPlus->Vect().DeltaPhi(m_p4MuonMinus->Vect()));
      m_deltaPhi[em] = deltaPhiEPlusMuMinus;
      dilepPt.SetXYZ( m_p4ElectronPlus->Vect().X()+m_p4MuonMinus->Vect().X(),
		      m_p4ElectronPlus->Vect().Y()+m_p4MuonMinus->Vect().Y(),
		      0.0 );
      dilepPtEPlusMuMinus = dilepPt.Mag();
      m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect()))));
    }
    if ( theElectron > -1 && theMuonPlus > -1 ) {
      deltaPhiEMinusMuPlus = fabs( 180./TMath::Pi() * m_p4ElectronMinus->Vect().DeltaPhi(m_p4MuonPlus->Vect()));
      m_deltaPhi[em] = deltaPhiEMinusMuPlus;
      dilepPt.SetXYZ( m_p4ElectronMinus->Vect().X()+m_p4MuonPlus->Vect().X(),
		      m_p4ElectronMinus->Vect().Y()+m_p4MuonPlus->Vect().Y(),
		      0.0 );
      dilepPtEMinusMuPlus = dilepPt.Mag();
      m_transvMass[em]=sqrt(2*dilepPt.Mag() * m_p4MET->Vect().Mag() * (1-cos(dilepPt.Angle(m_p4MET->Vect())) ) );
    }
    if ( thePositron > -1 && theMuonMinus > -1 &&
	 theElectron > -1 && theMuonPlus > -1) {
      
      // if two pairs are built we choose the one with highest di-lepton pt
      if ( dilepPtEPlusMuMinus > dilepPtEMinusMuPlus ) {
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
  bool foundJet=false;
  float maxPtJet=0.;
  for(int j=0;j<nJet;j++){
    
    // check if the electron or the positron falls into the jet
    // common cleaning class
    TVector3 p3Jet(pxJet[j],pyJet[j],pzJet[j]);

    if ( m_p4ElectronMinus->Vect().Mag() != 0 ) {
      float deltaR =  p3Jet.DeltaR( m_p4ElectronMinus->Vect() );
      if(_addedSel->getSwitch("jetConeWidth") && 
	 _addedSel->passCut("jetConeWidth",deltaR) &&
	 (m_HoEElectronMinus < 0.2) &&  
	 (m_CaloEneElectronMinus/energyJet[j] > 0.9)
	 ) continue;
    }
    if ( m_p4ElectronPlus->Vect().Mag() != 0 ) {
      float deltaR =  p3Jet.DeltaR( m_p4ElectronPlus->Vect() );
      if(_addedSel->getSwitch("jetConeWidth") && 
	 _addedSel->passCut("jetConeWidth",deltaR) &&
	 (m_HoEElectronPlus < 0.2) &&  
	 (m_CaloEneElectronPlus/energyJet[j] > 0.9)
	 ) continue;
    }

    // jet veto
    if(etJet[j]>maxPtJet) maxPtJet=etJet[j];
    if(_addedSel->getSwitch("etaJetAcc") && !_addedSel->passCut("etaJetAcc",etaJet[j])) continue;
    if((_addedSel->getSwitch("etJetLowAcc") && !_addedSel->passCut("etJetLowAcc",etJet[j]) ) ||
       ((_addedSel->getSwitch("etJetHighAcc") && !_addedSel->passCut("etJetHighAcc",etJet[j]) &&
	 (_addedSel->getSwitch("jetAlpha") && !_addedSel->passCut("jetAlpha",alphaJet[j]))))
       ) continue;
    foundJet=true;
  }

  etHighestJet->Fill(maxPtJet);
  return foundJet;

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

  // Muon quantities
  _monitorMuons->book1D("energy","muon energy (GeV)",energyMuon,50,0,150,"All+Fake+Best");
  _monitorMuons->book1D("et","muon tranverse energy (GeV)",etMuon,50,0,150,"All+Fake+Best");
  _monitorMuons->book1D("eta","muon #eta",etaMuon,50,-6.,6.,"All+Fake+Best");
  _monitorMuons->book1D("phi","muon #phi",phiMuon,50,-TMath::Pi(),TMath::Pi(),"All+Fake+Best");

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
