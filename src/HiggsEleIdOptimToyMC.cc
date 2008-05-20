#include <string>
#include <iostream>

#include <TTree.h>
#include "TRandom.h"

#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/HiggsEleIdOptimToyMC.hh"

HiggsEleIdOptimToyMC::HiggsEleIdOptimToyMC(TTree *tree) 
  : HiggsBase(tree) {

  // to be changed:
  // 1) signal or background
  nVar = 6;         // 2) number of variables
                    // 3) electron class
  theClass = 0;     // 0 = golden, EB
                    // 1 = golden, EE
                    // 2 = showering, EB
                    // 3 = showering, EE

  // kinematics
  _p4ElectronPlus  = new TLorentzVector(0.,0.,0.,0.);
  _p4ElectronMinus = new TLorentzVector(0.,0.,0.,0.);
  _bestElectrons   = new std::vector<int>;
  
  // event weight
  theWeight = 1;

  // counters
  allEvents    = 0.;
  passedMc     = 0.;
  triggered    = 0.; 
  commonPresel = 0.; 
  passedReco   = 0.;
  elePresel    = 0.; 
  looseId      = 0.; 
  passedIsol   = 0.;
  twoGoodClass = 0.;
  fullKine     = 0.;

  // 1 dimension histos - high pt electron
  HH_dEta   = new TH1F("HH_dEta",   "HH_dEta",   35,  0.0, 0.009);
  HH_dPhi   = new TH1F("HH_dPhi",   "HH_dPhi",   35,  0.0, 0.09);
  HH_HoE    = new TH1F("HH_HoE",    "HH_HoE",    35, -0.1, 0.12);
  HH_S9S25  = new TH1F("HH_S9s25",  "HH_S9S25",  35,  0.5, 1.0);
  HH_See    = new TH1F("HH_See",    "HH_See",    35,  0.0, 0.02);
  HH_EoPout = new TH1F("HH_EoPout", "HH_EoPout", 35,  0.5, 4.0);

  // 1 dimension histos - low pt electron
  HL_dEta   = new TH1F("HL_dEta",   "HL_dEta",   35,  0.0, 0.009);
  HL_dPhi   = new TH1F("HL_dPhi",   "HL_dPhi",   35,  0.0, 0.09);
  HL_HoE    = new TH1F("HL_HoE",    "HL_HoE",    35, -0.1, 0.12);
  HL_S9S25  = new TH1F("HL_S9s25",  "HL_S9S25",  35,  0.5, 1.0);
  HL_See    = new TH1F("HL_See",    "HL_See",    35,  0.0, 0.02);
  HL_EoPout = new TH1F("HL_EoPout", "HL_EoPout", 35,  0.5, 4.0);

  // N dimensions histo
  theBins[0] = 35;
  theBins[1] = 35;
  theBins[2] = 35;
  theBins[3] = 35;
  theBins[4] = 35;
  theBins[5] = 35;
  // 
  theMin[0] =  0.0;
  theMin[1] =  0.0; 
  theMin[2] = -0.1; 
  theMin[3] =  0.5;
  theMin[4] =  0.0;
  theMin[5] =  0.5;
  //
  theMax[0] = 0.009;
  theMax[1] = 0.09; 
  theMax[2] = 0.12; 
  theMax[3] = 1.0;
  theMax[4] = 0.02;
  theMax[5] = 4.0;
  // 
  HH_NVarDim = new THnSparseF("HH_NVarDim", "HH_NVarDim", nVar, theBins, theMin, theMax);
  HL_NVarDim = new THnSparseF("HL_NVarDim", "HL_NVarDim", nVar, theBins, theMin, theMax);

  // output tree
  outRootTree = new RedEleIDOptimTree("myOutTree.root");

  // output root file
  outRootFile = new TFile("outHistos.root","RECREATE");
}

HiggsEleIdOptimToyMC::~HiggsEleIdOptimToyMC(){
  
  // output txt file with efficiencies:
  ofstream *outTxtFile = new ofstream("outputFile.txt",ios::app);
  *outTxtFile << "all events:    "      << allEvents    << endl;
  *outTxtFile << "passedMC:      "      << passedMc     << endl;
  *outTxtFile << "triggered:     "      << triggered    << endl;
  *outTxtFile << "common presel: "      << commonPresel << endl;
  *outTxtFile << "reco:          "      << passedReco   << endl;
  *outTxtFile << "presel x ele:  "      << elePresel    << endl;
  *outTxtFile << "loose eleID:   "      << looseId      << endl;
  *outTxtFile << "isolation:     "      << passedIsol   << endl;
  *outTxtFile << "ok class for both: "  << twoGoodClass << endl;
  *outTxtFile << "full kinematics ok: " << fullKine     << endl; 
  
  // saving histos and tree
  outRootFile->Close();
  outRootTree->save();

  // deleting
  delete _p4ElectronPlus;
  delete _p4ElectronMinus;
  delete _bestElectrons;
}

bool HiggsEleIdOptimToyMC::findMcTree(const char* processType) {

  _process = "UNDEFINED";
  _theGenEle = -1;
  _theGenPos = -1;
  
  // signal: 2e2nu
  if(strcmp(processType,"HtoWWto2e2nu")==0) {
    if( idMc[9]  == -11) _theGenEle = 9;
    if( idMc[11] ==  11) _theGenPos = 11;
    theWeight=1; 
    return (_theGenEle > -1 && _theGenPos > -1 );
  }
  // w+jets: e / mu / tau
  else if(strcmp(processType,"Wjets")==0) {
    _process = "Wjets";
    theWeight=genWeight; 
    return ( genAlpgenID>=1000 && genAlpgenID<2000);
  }
  else {
    std::cout << "This processType: " << processType << "is not expected, you should put MTtruth switch off" << std::endl;
    return false;
  }
}


void HiggsEleIdOptimToyMC::Loop() {

  _verbose=false;
  if(fChain == 0) return;
  
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    resetKinematics();

    // before any cut
    allEvents=allEvents+theWeight; 

    // look to the MC truth decay tree
    // bool foundMcTree = findMcTree("HtoWWto2e2nu");
    bool foundMcTree = findMcTree("Wjets");
    if ( !foundMcTree ) continue;
    passedMc=passedMc+theWeight;    

    // trigger
    Utils anaUtils;
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
    if (!passedHLT) continue; 
    triggered=triggered+theWeight;   

    // did we pass preselections?
    if( !evtPresel ) continue;
    commonPresel=commonPresel+theWeight;

    // did we reconstruct two electrons?
    std::pair<int,int> theElectrons = getBestElectronPair();
    int tbElectron(theElectrons.second), tbPositron(theElectrons.first);    
    theElectron = tbElectron;
    thePositron = tbPositron;
    if (theElectron<0 || thePositron<0)  continue;  
    passedReco=passedReco+theWeight;  

    // did we pass preselections specific for electrons?
    setKinematics();
    if (hardestElectronPt < 20) continue; 
    if (slowestElectronPt < 10) continue; 
    if (etMet[0] < 30)          continue; 
    if (_mll < 12)             continue; 
    elePresel=elePresel+theWeight; 

    // did we pass the loose electronId?
    if (!eleIdCutBasedEle[theElectron] || !eleIdCutBasedEle[thePositron]) continue;
    looseId=looseId+theWeight;   

    // did we pass loose electron tracker based isolation?
    if (eleTrackerIso_sumPtEle[theElectron]>0.1 || eleTrackerIso_sumPtEle[thePositron]>0.1) continue;
    passedIsol=passedIsol+theWeight;

    // is this electron belonging to the good class?
    bool isEle = false;
    bool isPos = false;
    if (theClass==0 && (eleClassEle[theElectron]==0)) isEle = true;
    if (theClass==0 && (eleClassEle[thePositron]==0)) isPos = true;
    if (theClass==1 && (eleClassEle[theElectron]==100)) isEle = true;
    if (theClass==1 && (eleClassEle[thePositron]==100)) isPos = true;
    if (theClass==2 && (eleClassEle[theElectron]>=30 && eleClassEle[theElectron]<=40)) isEle = true;
    if (theClass==2 && (eleClassEle[thePositron]>=30 && eleClassEle[thePositron]<=40)) isPos = true;
    if (theClass==3 && (eleClassEle[theElectron]>=130 && eleClassEle[theElectron]<=140)) isEle = true;
    if (theClass==3 && (eleClassEle[thePositron]>=130 && eleClassEle[thePositron]<=140)) isPos = true;
    if (isEle && isPos) { twoGoodClass=twoGoodClass+theWeight; }

    // ordering in pt
    bool isHigherEle = true; 
    if(etEle[theElectron] >= etEle[thePositron]) isHigherEle = true; 
    if(etEle[theElectron] < etEle[thePositron])  isHigherEle = false; 
    
    // filling histos: 1 dim and N dim histos
    if(isEle){
      double toFill[nVar];
      toFill[0] = fabs(eleDeltaEtaAtVtxEle[theElectron]);
      toFill[1] = fabs(eleDeltaPhiAtVtxEle[theElectron]);
      toFill[2] = eleHoEEle[theElectron];
      toFill[3] = s9s25Ele[theElectron];
      toFill[4] = sqrt(covEtaEtaEle[theElectron]);
      toFill[5] = eleCorrEoPoutEle[theElectron];

      // golden electrons: fill HH and HL together 
      if(theClass==0 || theClass==1){
	HH_dEta    -> Fill(toFill[0],theWeight);
	HH_dPhi    -> Fill(toFill[1],theWeight);
	HH_HoE     -> Fill(toFill[2],theWeight);
	HH_S9S25   -> Fill(toFill[3],theWeight);
	HH_See     -> Fill(toFill[4],theWeight);
	HH_EoPout  -> Fill(toFill[5],theWeight);
	HH_NVarDim -> Fill(toFill,theWeight);
	HL_dEta    -> Fill(toFill[0],theWeight);
	HL_dPhi    -> Fill(toFill[1],theWeight);
	HL_HoE     -> Fill(toFill[2],theWeight);
	HL_S9S25   -> Fill(toFill[3],theWeight);
	HL_See     -> Fill(toFill[4],theWeight);
	HL_EoPout  -> Fill(toFill[5],theWeight);
	HL_NVarDim -> Fill(toFill,theWeight);
      }

      // showering electrons: separating high and low pt
      if(theClass==2 || theClass==3){
	if(isHigherEle){  	// high pt lepton
	  HH_dEta    -> Fill(toFill[0],theWeight);
	  HH_dPhi    -> Fill(toFill[1],theWeight);
	  HH_HoE     -> Fill(toFill[2],theWeight);
	  HH_S9S25   -> Fill(toFill[3],theWeight);
	  HH_See     -> Fill(toFill[4],theWeight);
	  HH_EoPout  -> Fill(toFill[5],theWeight);
	  HH_NVarDim -> Fill(toFill,theWeight);
	}
	if(!isHigherEle){   	// low pt lepton
	  HL_dEta    -> Fill(toFill[0],theWeight);
	  HL_dPhi    -> Fill(toFill[1],theWeight);
	  HL_HoE     -> Fill(toFill[2],theWeight);
	  HL_S9S25   -> Fill(toFill[3],theWeight);
	  HL_See     -> Fill(toFill[4],theWeight);
	  HL_EoPout  -> Fill(toFill[5],theWeight);
	  HL_NVarDim -> Fill(toFill,theWeight);
	}
      }
    }
    if(isPos){
      double toFill[nVar];
      toFill[0] = fabs(eleDeltaEtaAtVtxEle[thePositron]);
      toFill[1] = fabs(eleDeltaPhiAtVtxEle[thePositron]);
      toFill[2] = eleHoEEle[thePositron];
      toFill[3] = s9s25Ele[thePositron];
      toFill[4] = sqrt(covEtaEtaEle[thePositron]);
      toFill[5] = eleCorrEoPoutEle[thePositron];

      // golden electrons: fill HH and HL together
      if(theClass==0 || theClass==1){
	HH_dEta    -> Fill(toFill[0],theWeight);
	HH_dPhi    -> Fill(toFill[1],theWeight);
	HH_HoE     -> Fill(toFill[2],theWeight);
	HH_S9S25   -> Fill(toFill[3],theWeight);
	HH_See     -> Fill(toFill[4],theWeight);
	HH_EoPout  -> Fill(toFill[5],theWeight);
	HH_NVarDim -> Fill(toFill,theWeight);
	HL_dEta    -> Fill(toFill[0],theWeight);
	HL_dPhi    -> Fill(toFill[1],theWeight);
	HL_HoE     -> Fill(toFill[2],theWeight);
	HL_S9S25   -> Fill(toFill[3],theWeight);
	HL_See     -> Fill(toFill[4],theWeight);
	HL_EoPout  -> Fill(toFill[5],theWeight);
	HL_NVarDim -> Fill(toFill,theWeight);
      }

      // showering electrons: separating high and low pt
      if(theClass==2 || theClass==3){ 	  
	if(!isHigherEle){       // high pt lepton
	  HH_dEta    -> Fill(toFill[0],theWeight);
	  HH_dPhi    -> Fill(toFill[1],theWeight);
	  HH_HoE     -> Fill(toFill[2],theWeight);
	  HH_S9S25   -> Fill(toFill[3],theWeight);
	  HH_See     -> Fill(toFill[4],theWeight);
	  HH_EoPout  -> Fill(toFill[5],theWeight);
	  HH_NVarDim -> Fill(toFill,theWeight);
	}
	if(isHigherEle){	// low pt lepton
	  HL_dEta    -> Fill(toFill[0],theWeight);
	  HL_dPhi    -> Fill(toFill[1],theWeight);
	  HL_HoE     -> Fill(toFill[2],theWeight);
	  HL_S9S25   -> Fill(toFill[3],theWeight);
	  HL_See     -> Fill(toFill[4],theWeight);
	  HL_EoPout  -> Fill(toFill[5],theWeight);
	  HL_NVarDim -> Fill(toFill,theWeight);
	}
      }
    }

    // full kine analysis - to be modified according the mass hypothesis
    if (!isEle || !isPos)       continue;
    if (goodJetFound())         continue;
    if (etMet[0] < 50)          continue;
    if (_deltaPhi > 45)         continue;
    if (hardestElectronPt < 25) continue;
    if (hardestElectronPt > 50) continue;
    if (slowestElectronPt < 15) continue;
    if (_mll > 50)              continue;
    fullKine=fullKine+theWeight;

  } // end loop over entries

  // generating a random distribution according to the sampled one
  TDatime *now = new TDatime();  
  for(int scan =0; scan<100000; scan++){
    if(scan%1000==0) cout << "random: " << scan << endl;

    // randomly generated electron and positron 
    int today = now->GetDate();  
    int clock = now->GetTime();  
    int seed  = today+clock+scan;
    gRandom->SetSeed(seed);
    double theRndHigh[6], theRndLow[6];
    HH_NVarDim->GetRandom(theRndHigh);
    HL_NVarDim->GetRandom(theRndLow);

    // filling the reduced tree with toy MC electrons
    outRootTree->fillAll(theRndHigh[0],theRndHigh[1],theRndHigh[2],theRndHigh[3],theRndHigh[4],theRndHigh[5],theRndLow[0],theRndLow[1],theRndLow[2],theRndLow[3],theRndLow[4],theRndLow[5]);
  outRootTree -> store();

  } // mc generation
  delete now;

  // saving the histos
  outRootFile->cd();
  HH_dEta         -> Write();
  HH_dPhi         -> Write();
  HH_HoE          -> Write();
  HH_S9S25        -> Write();
  HH_See          -> Write();
  HH_EoPout       -> Write();
  HH_NVarDim      -> Write();
  HL_dEta         -> Write();
  HL_dPhi         -> Write();
  HL_HoE          -> Write();
  HL_S9S25        -> Write();
  HL_See          -> Write();
  HL_EoPout       -> Write();
  HL_NVarDim      -> Write();
}

std::pair<int,int> HiggsEleIdOptimToyMC::getBestElectronPair() {
  int theLep1=-1;          
  int theLep2=-1;
  float maxPtLep1=-1000.;  
  float maxPtLep2=-1000.;
  std::vector<int> goodRecoLeptons;
  for(int i=0;i<nEle;i++) {
    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();
    if(fabs(etaEle[i])>2.5) continue;
    if(thisPt<10)           continue;
    float thisCharge = chargeEle[i];
    if (thisCharge > 0 && thisPt> maxPtLep1){ maxPtLep1 = thisPt; theLep1 = i; }
    if (thisCharge < 0 && thisPt> maxPtLep2){ maxPtLep2 = thisPt; theLep2 = i; }
  }
  _bestElectrons->clear();
  _bestElectrons->push_back(theLep1);  
  _bestElectrons->push_back(theLep2); 
  return make_pair(theLep1,theLep2);
}

void HiggsEleIdOptimToyMC::setKinematics() {
  hardestElectronPt = TMath::Max(etEle[theElectron],etEle[thePositron]);
  slowestElectronPt = TMath::Min(etEle[theElectron],etEle[thePositron]);
  _p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
  _p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
  _mll = (*_p4ElectronMinus + *_p4ElectronPlus).M();
  _deltaPhi = fabs(180./TMath::Pi() * _p4ElectronMinus->Vect().DeltaPhi(_p4ElectronPlus->Vect()));
  _HoEElectronMinus     = eleHoEEle[theElectron];
  _HoEElectronPlus      = eleHoEEle[thePositron];
  _CaloEneElectronMinus = eleCaloCorrEEle[theElectron];
  _CaloEneElectronPlus  = eleCaloCorrEEle[thePositron];
}

void HiggsEleIdOptimToyMC::resetKinematics() {
  theElectron = -1;
  thePositron = -1;
  _p4ElectronMinus->SetXYZT(0,0,0,0);
  _p4ElectronPlus->SetXYZT(0,0,0,0);
  _mll = 0;
  _deltaPhi = 0;
  hardestElectronPt = 0;
  slowestElectronPt = 0;
  _HoEElectronMinus = 0;
  _HoEElectronPlus  = 0;
  _CaloEneElectronMinus = 0;
  _CaloEneElectronPlus  = 0;
}

bool HiggsEleIdOptimToyMC::goodJetFound() {

  bool foundJet=false;
  for(int jj=0;jj<nJet;jj++) {
    // check if the electron falls into the jet
    TVector3 p3Jet(pxJet[jj],pyJet[jj],pzJet[jj]);
    if (_p4ElectronMinus->Vect().Mag() != 0) {
      float deltaR = p3Jet.DeltaR(_p4ElectronMinus->Vect());
      if( (deltaR<0.2) && (_HoEElectronMinus < 0.2) && (_CaloEneElectronMinus/energyJet[jj] > 0.9) ) continue;
    }
    // check if the positron falls into the jet
    if (_p4ElectronPlus->Vect().Mag() != 0) {
      float deltaR =  p3Jet.DeltaR(_p4ElectronPlus->Vect());
      if( (deltaR<0.2) && (_HoEElectronPlus < 0.2) && (_CaloEneElectronPlus/energyJet[jj] > 0.9) ) continue;
    }
    // checking the jet kinematics
    if(fabs(etaJet[jj])>2.5) continue;
    if(etJet[jj]<15)         continue;
    if(etJet[jj]>=15 && etJet[jj]<=20 && alphaJet[jj]<0.2) continue; 

    // this is a good jet
    foundJet=true;
    break;
  }
  return foundJet;
}
