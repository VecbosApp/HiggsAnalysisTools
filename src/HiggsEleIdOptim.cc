#include <string>

#include <TTree.h>

#include "CommonTools/include/Utils.hh"
#include "HiggsAnalysisTools/include/HiggsEleIdOptim.hh"

#include <iostream>
#include <string>

#include <TTree.h>

HiggsEleIdOptim::HiggsEleIdOptim(TTree *tree) 
  : HiggsBase(tree) {

  // kinematics
  m_p4ElectronPlus  = new TLorentzVector(0.,0.,0.,0.);
  m_p4ElectronMinus = new TLorentzVector(0.,0.,0.,0.);
  _bestElectrons    = new std::vector<int>;

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
  for(int iiDeta=0; iiDeta<6; iiDeta++){
    for(int iiDphi=0; iiDphi<5; iiDphi++){
      for(int iiHoE=0; iiHoE<5; iiHoE++){
	for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	  for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
	    for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
	      for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=0.;
	      }}}}}}}}
}

HiggsEleIdOptim::~HiggsEleIdOptim(){

  // output file:
  ofstream *outFile  = new ofstream("outputFile.txt",ios::app);
  *outFile << allEvents    << endl;
  *outFile << passedMc     << endl;
  *outFile << triggered    << endl;
  *outFile << commonPresel << endl;
  *outFile << passedReco   << endl;
  *outFile << elePresel    << endl;
  *outFile << looseId      << endl;
  *outFile << passedIsol   << endl;
  for(int iiDeta=0; iiDeta<6; iiDeta++){
    for(int iiDphi=0; iiDphi<5; iiDphi++){
      for(int iiHoE=0; iiHoE<5; iiHoE++){
	for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	  for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
	    for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
	      for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		  *outFile << iiDeta      << " " << iiDphi      << " " << iiHoE    << " " << iiS9S25  << " " 
			   << iiEoPoutmin << " " << iiEoPoutmax << " " << iiSeemin << " " << iiSeemax << " " 
			   << passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax] << endl;
	      }}}}}}}}

  delete m_p4ElectronPlus;
  delete m_p4ElectronMinus;
  delete _bestElectrons;
}

bool HiggsEleIdOptim::findMcTree(const char* processType) {

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


void HiggsEleIdOptim::Loop() {

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
    bool foundMcTree = findMcTree("HtoWWto2e2nu");
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
    if (m_mll < 12)             continue; 
    elePresel=elePresel+theWeight; 
  
    // did we pass the loose electronId?
    if (!eleIdCutBasedEle[theElectron] || !eleIdCutBasedEle[thePositron]) continue;
    looseId=looseId+theWeight;   

    // did we pass loose electron tracker based isolation?
    if (eleTrackerIso_sumPtEle[theElectron]>0.1 || eleTrackerIso_sumPtEle[thePositron]>0.1) continue;
    passedIsol=passedIsol+theWeight;

    // real analysis to study the best eleID cuts
    bool theElectronID = false;
    bool thePositronID = false;
    bool theEleIDScan[6][5][5][3][3][3][3][3], thePosIDScan[6][5][5][3][3][3][3][3]; 
    for(int iiDeta=0; iiDeta<6; iiDeta++){
      for(int iiDphi=0; iiDphi<5; iiDphi++){
	for(int iiHoE=0; iiHoE<5; iiHoE++){
	  for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	    for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
	      for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		  for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		    theEleIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=false;
		    thePosIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=false;
		  }}}}}}}}

    // electron is not a EB golden - apply loose egamma cuts 
    if (eleClassEle[theElectron]!=0){
      theElectronID=isEleID(theElectron);
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
		for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      theEleIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=theElectronID;
		    }}}}}}}}
    } 

    // positron is not a EB golden - apply loose egamma cuts
    if (eleClassEle[thePositron]!=0){
      thePositronID=isEleID(thePositron);
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
		for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      thePosIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=thePositronID;
		    }}}}}}}}
    } 
    
    // electron is a EB golden - perform the scan
    if (eleClassEle[theElectron]==0){
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
		for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      theEleIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=isEleIDScan(theElectron, iiDeta, iiDphi, iiHoE, iiS9S25, iiEoPoutmin, iiEoPoutmax, iiSeemin, iiSeemax);
		    }}}}}}}}
    }

    // positron is a EB golden - perform the scan
    if (eleClassEle[thePositron]==0){
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
		for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      thePosIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]=isEleIDScan(thePositron, iiDeta, iiDphi, iiHoE, iiS9S25, iiEoPoutmin, iiEoPoutmax, iiSeemin, iiSeemax);
		    }}}}}}}}
    }

    // incrementing the counter for each category
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPoutmin=0; iiEoPoutmin<3; iiEoPoutmin++){
		for(int iiEoPoutmax=0; iiEoPoutmax<3; iiEoPoutmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      if(theEleIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax] && thePosIDScan[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]){ 
			passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax] = 
			  passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPoutmin][iiEoPoutmax][iiSeemin][iiSeemax]+theWeight; }
		    }}}}}}}}
  }
}

std::pair<int,int> HiggsEleIdOptim::getBestElectronPair() {
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

bool HiggsEleIdOptim::isEleID(int eleIndex) {

  double HOverEMaxCut, S9S25MaxCut, DEtaMaxCut, DPhiMaxCut, SeeMaxCut, EoPoutMaxCut;
  double HOverEMinCut, S9S25MinCut, DEtaMinCut, DPhiMinCut, SeeMinCut, EoPoutMinCut;
  
  // default loose golden barrel
  if (eleClassEle[eleIndex]==0){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.09;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.7; 
    EoPoutMaxCut = 2.5;    EoPoutMinCut = 0.6;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose golden endcap
  if (eleClassEle[eleIndex]==100){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.09;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.8; 
    EoPoutMaxCut = 2.5;    EoPoutMinCut = 0.6;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
 
  // default loose big brem barrel
  if (eleClassEle[eleIndex]==10){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.06;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.7; 
    EoPoutMaxCut = 999.;   EoPoutMinCut = 1.7;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose big brem endcap
  if (eleClassEle[eleIndex]==110){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.06;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.8; 
    EoPoutMaxCut = 999.;   EoPoutMinCut = 1.7;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose narrow barrel
  if (eleClassEle[eleIndex]==20){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.07;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.7; 
    EoPoutMaxCut = 2.2;    EoPoutMinCut = 0.9;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose narrow endcap
  if (eleClassEle[eleIndex]==120){
    DEtaMaxCut   = 0.008;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.06;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.07;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.8; 
    EoPoutMaxCut = 2.2;    EoPoutMinCut = 0.9;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose  showering/crack barrel 
  if (eleClassEle[eleIndex]>=30 && eleClassEle[eleIndex]<=40){
    DEtaMaxCut   = 0.009;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.08;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.12;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.5; 
    EoPoutMaxCut = 999.;   EoPoutMinCut = 0.5;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  // default loose  showering/crack endcap
  if (eleClassEle[eleIndex]>=130 && eleClassEle[eleIndex]<=140){
    DEtaMaxCut   = 0.009;  DEtaMinCut   = 0.0; 
    DPhiMaxCut   = 0.09;   DPhiMinCut   = 0.0;     
    HOverEMaxCut = 0.12;   HOverEMinCut = 0.0;   
    S9S25MaxCut  = 1.0;    S9S25MinCut  = 0.5; 
    EoPoutMaxCut = 999.;   EoPoutMinCut = 0.5;  
    SeeMaxCut    = 999.;   SeeMinCut    = 0.0;     
  }
  
  bool idPassed = false;
  if( (fabs(eleDeltaEtaAtVtxEle[eleIndex])<=DEtaMaxCut) && (fabs(eleDeltaEtaAtVtxEle[eleIndex])>=DEtaMinCut)){
    if( (fabs(eleDeltaPhiAtVtxEle[eleIndex])<=DPhiMaxCut) && (fabs(eleDeltaPhiAtVtxEle[eleIndex])>=DPhiMinCut)){
      if( (eleHoEEle[eleIndex]<=HOverEMaxCut) && (eleHoEEle[eleIndex]>=HOverEMinCut) ){
	if( (s9s25Ele[eleIndex]<=S9S25MaxCut) && (s9s25Ele[eleIndex]>=S9S25MinCut) ){
	  if( (sqrt(covEtaEtaEle[eleIndex])<=SeeMaxCut) && (sqrt(covEtaEtaEle[eleIndex])>=SeeMinCut) ){
	    if( (eleCorrEoPoutEle[eleIndex]<=EoPoutMaxCut) && (eleCorrEoPoutEle[eleIndex]>=EoPoutMinCut) ){
	      idPassed = true;
	    }}}}}}

  return idPassed; 
}

bool HiggsEleIdOptim::isEleIDScan(int eleIndex, int iiDeta, int iiDphi, int iiHoE, int iiS9S25, int iiEoPoutmin, int iiEoPoutmax, int iiSeemin, int iiSeemax) {

  double HOverEMaxCut,  S9S25MaxCut,  DEtaMaxCut,  DPhiMaxCut,  SeeMaxCut,  EoPoutMaxCut;
  double HOverEMinCut,  S9S25MinCut,  DEtaMinCut,  DPhiMinCut,  SeeMinCut,  EoPoutMinCut;
  double HOverEMaxStep, S9S25MinStep, DEtaMaxStep, DPhiMaxStep, SeeMaxStep, SeeMinStep, EoPoutMaxStep, EoPoutMinStep;
  double HOverEMaxInit, S9S25MinInit, DEtaMaxInit, DPhiMaxInit, SeeMaxInit, SeeMinInit, EoPoutMaxInit, EoPoutMinInit;
  
  // fixed values
  DEtaMinCut   = 0.0; 
  DPhiMinCut   = 0.0;     
  HOverEMinCut = 0.0;   
  S9S25MaxCut  = 1.0;    

  // golden barrel
  if (eleClassEle[eleIndex]==0){
    // starting point   
    DEtaMaxInit   = 0.004;  
    DPhiMaxInit   = 0.02;   
    HOverEMaxInit = 0.05;   
    S9S25MinInit  = 0.6;    
    EoPoutMaxInit = 2.0;    
    EoPoutMinInit = 0.5;  
    SeeMaxInit    = 0.009;   
    SeeMinInit    = 0.004;   
    // step
    DEtaMaxStep   = 0.001;
    DPhiMaxStep   = 0.01;   
    HOverEMaxStep = 0.01;   
    S9S25MinStep  = 0.1; 
    EoPoutMaxStep = 1.0;    
    EoPoutMinStep = 0.1;  
    SeeMaxStep    = 0.002;   
    SeeMinStep    = 0.001;   
    // the current cut
    DEtaMaxCut   = DEtaMaxInit + iiDeta*DEtaMaxStep;
    DPhiMaxCut   = DPhiMaxInit + iiDphi*DPhiMaxStep;
    HOverEMaxCut = HOverEMaxInit + iiHoE*HOverEMaxStep;
    S9S25MinCut  = S9S25MinInit + iiS9S25*S9S25MinStep;
    EoPoutMaxCut = EoPoutMaxInit + iiEoPoutmax*EoPoutMaxStep;
    EoPoutMinCut = EoPoutMinInit + iiEoPoutmin*EoPoutMinStep;
    SeeMaxCut    = SeeMaxInit + iiSeemax*SeeMaxStep;
    SeeMinCut    = SeeMinInit + iiSeemin*SeeMinStep;
  }

  // golden endcap
  if (eleClassEle[eleIndex]==100){
    // starting point   
    DEtaMaxInit   = 0.004;  
    DPhiMaxInit   = 0.03;   
    HOverEMaxInit = 0.06;   
    S9S25MinInit  = 0.7;    
    EoPoutMaxInit = 2.0;    
    EoPoutMinInit = 0.5;  
    SeeMaxInit    = 0.025;   
    SeeMinInit    = 0.007;   
    // step
    DEtaMaxStep   = 0.001;
    DPhiMaxStep   = 0.01;   
    HOverEMaxStep = 0.01;   
    S9S25MinStep  = 0.1; 
    EoPoutMaxStep = 1.0;    
    EoPoutMinStep = 0.1;  
    SeeMaxStep    = 0.005;   
    SeeMinStep    = 0.001;   
    // the current cut
    DEtaMaxCut   = DEtaMaxInit + iiDeta*DEtaMaxStep;
    DPhiMaxCut   = DPhiMaxInit + iiDphi*DPhiMaxStep;
    HOverEMaxCut = HOverEMaxInit + iiHoE*HOverEMaxStep;
    S9S25MinCut  = S9S25MinInit + iiS9S25*S9S25MinStep;
    EoPoutMaxCut = EoPoutMaxInit + iiEoPoutmax*EoPoutMaxStep;
    EoPoutMinCut = EoPoutMinInit + iiEoPoutmin*EoPoutMinStep;
    SeeMaxCut    = SeeMaxInit + iiSeemax*SeeMaxStep;
    SeeMinCut    = SeeMinInit + iiSeemin*SeeMinStep;
  }

  // showering barrel
  if (eleClassEle[eleIndex]>=30 && eleClassEle[eleIndex]<=40){
    // starting point   
    DEtaMaxInit   = 0.005;  
    DPhiMaxInit   = 0.04;   
    HOverEMaxInit = 0.04;   
    S9S25MinInit  = 0.5;    
    EoPoutMaxInit = 999.;    
    EoPoutMinInit = 0.4;  
    SeeMaxInit    = 0.009;   
    SeeMinInit    = 0.004;   
    // step
    DEtaMaxStep   = 0.001;
    DPhiMaxStep   = 0.01;   
    HOverEMaxStep = 0.02;   
    S9S25MinStep  = 0.1; 
    EoPoutMaxStep = 1.0;    
    EoPoutMinStep = 0.2;  
    SeeMaxStep    = 0.002;   
    SeeMinStep    = 0.001;   
    // the current cut
    DEtaMaxCut   = DEtaMaxInit + iiDeta*DEtaMaxStep;
    DPhiMaxCut   = DPhiMaxInit + iiDphi*DPhiMaxStep;
    HOverEMaxCut = HOverEMaxInit + iiHoE*HOverEMaxStep;
    S9S25MinCut  = S9S25MinInit + iiS9S25*S9S25MinStep;
    EoPoutMaxCut = EoPoutMaxInit + iiEoPoutmax*EoPoutMaxStep;
    EoPoutMinCut = EoPoutMinInit + iiEoPoutmin*EoPoutMinStep;
    SeeMaxCut    = SeeMaxInit + iiSeemax*SeeMaxStep;
    SeeMinCut    = SeeMinInit + iiSeemin*SeeMinStep;
  }

  // showering endcap
  if (eleClassEle[eleIndex]>=130 && eleClassEle[eleIndex]<=140){
    // starting point   
    DEtaMaxInit   = 0.005;  
    DPhiMaxInit   = 0.05;   
    HOverEMaxInit = 0.06;   
    S9S25MinInit  = 0.5;    
    EoPoutMaxInit = 999.;    
    EoPoutMinInit = 0.4;  
    SeeMaxInit    = 0.018;   
    SeeMinInit    = 0.007;   
    // step
    DEtaMaxStep   = 0.001;
    DPhiMaxStep   = 0.01;   
    HOverEMaxStep = 0.02;   
    S9S25MinStep  = 0.1; 
    EoPoutMaxStep = 1.0;    
    EoPoutMinStep = 0.2;  
    SeeMaxStep    = 0.005;   
    SeeMinStep    = 0.001;   
    // the current cut
    DEtaMaxCut   = DEtaMaxInit + iiDeta*DEtaMaxStep;
    DPhiMaxCut   = DPhiMaxInit + iiDphi*DPhiMaxStep;
    HOverEMaxCut = HOverEMaxInit + iiHoE*HOverEMaxStep;
    S9S25MinCut  = S9S25MinInit + iiS9S25*S9S25MinStep;
    EoPoutMaxCut = EoPoutMaxInit + iiEoPoutmax*EoPoutMaxStep;
    EoPoutMinCut = EoPoutMinInit + iiEoPoutmin*EoPoutMinStep;
    SeeMaxCut    = SeeMaxInit + iiSeemax*SeeMaxStep;
    SeeMinCut    = SeeMinInit + iiSeemin*SeeMinStep;
  }

  bool idPassed = false;
  if( (fabs(eleDeltaEtaAtVtxEle[eleIndex])<=DEtaMaxCut) && (fabs(eleDeltaEtaAtVtxEle[eleIndex])>=DEtaMinCut)){
    if( (fabs(eleDeltaPhiAtVtxEle[eleIndex])<=DPhiMaxCut) && (fabs(eleDeltaPhiAtVtxEle[eleIndex])>=DPhiMinCut)){
      if( (eleHoEEle[eleIndex]<=HOverEMaxCut) && (eleHoEEle[eleIndex]>=HOverEMinCut) ){
	if( (s9s25Ele[eleIndex]<=S9S25MaxCut) && (s9s25Ele[eleIndex]>=S9S25MinCut) ){
	  if( (sqrt(covEtaEtaEle[eleIndex])<=SeeMaxCut) && (sqrt(covEtaEtaEle[eleIndex])>=SeeMinCut) ){
	    if( (eleCorrEoPoutEle[eleIndex]<=EoPoutMaxCut) && (eleCorrEoPoutEle[eleIndex]>=EoPoutMinCut) ){
	      idPassed = true;
	    }}}}}}

  return idPassed; 
}

void HiggsEleIdOptim::setKinematics() {
  hardestElectronPt = TMath::Max(etEle[theElectron],etEle[thePositron]);
  slowestElectronPt = TMath::Min(etEle[theElectron],etEle[thePositron]);
  m_p4ElectronMinus -> SetXYZT(pxEle[theElectron],pyEle[theElectron],pzEle[theElectron],energyEle[theElectron]);
  m_p4ElectronPlus  -> SetXYZT(pxEle[thePositron],pyEle[thePositron],pzEle[thePositron],energyEle[thePositron]);      
  m_mll      = (*m_p4ElectronMinus + *m_p4ElectronPlus).M();
}

void HiggsEleIdOptim::resetKinematics() {
  theElectron = -1;
  thePositron = -1;
  m_p4ElectronMinus->SetXYZT(0,0,0,0);
  m_p4ElectronPlus->SetXYZT(0,0,0,0);
  m_mll = 0;
  hardestElectronPt = 0;
  slowestElectronPt = 0;
}

