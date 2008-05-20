// ! c++ includes
#include <string>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <fstream.h>
#include <math.h>

//! ROOT includes
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TApplication.h"
#include "TBranch.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

// variables to be used later on
float HOverEMaxCut,  S9S25MaxCut,  DEtaMaxCut,  DPhiMaxCut,  SeeMaxCut,  EoPoutMaxCut;
float HOverEMinCut,  S9S25MinCut,  DEtaMinCut,  DPhiMinCut,  SeeMinCut,  EoPoutMinCut;
float HOverEMaxStep, S9S25MinStep, DEtaMaxStep, DPhiMaxStep, SeeMaxStep, SeeMinStep, EoPoutMaxStep, EoPoutMinStep;
float HOverEMaxInit, S9S25MinInit, DEtaMaxInit, DPhiMaxInit, SeeMaxInit, SeeMinInit, EoPoutMaxInit, EoPoutMinInit;

// methods
void setScanValue(int theClass);
bool isEleIDScan(float thisDeta, float thisDphi, float thisHoE, float thisS9s25, float thisEopOut, float thisSee, int deta, int dphi, int hoe, int s9s25, int eopmin, int eopmax, int seemin, int seemax);

int main ( int argc, char **argv)
{
  if (argc < 8){ 
    cout << "Argument missing! Insert: "        << std::endl; 
    cout << "1) inputFile - root tree for sgn " << std::endl;
    cout << "2) inputFile - root tree for bkg " << std::endl;
    cout << "3) class:"                         << std::endl;
    cout << "   0 --> golden EB"                << std::endl;
    cout << "   1 --> golden EE"                << std::endl;
    cout << "   2 --> showering EB"             << std::endl;
    cout << "   3 --> showering EE"             << std::endl;
    cout << "4) signal preEleID efficiency"     << std::endl;
    cout << "5) signal kine efficiency"         << std::endl;
    cout << "6) background preEleID # events"   << std::endl;
    cout << "7) background kine efficiency"     << std::endl;
    return 0;
  }


  // reading the input trees --------------------------
  TChain *T[2];
  T[0]= new TChain("T1");
  T[1]= new TChain("T1");
  T[0]->Add(argv[1]);   // signal
  T[1]->Add(argv[2]);   // background
  float H_deltaEta, L_deltaEta;
  float H_deltaPhi, L_deltaPhi;
  float H_hoe,      L_hoe;
  float H_s9s25,    L_s9s25;
  float H_see,      L_see;
  float H_eopOut,   L_eopOut;  
  for(int ii=0; ii<2; ii++){
    T[ii]->SetMakeClass(1);
    T[ii]->SetBranchStatus("*",0);
    T[ii]->SetBranchStatus("H_deltaEta",1);
    T[ii]->SetBranchStatus("H_deltaPhi",1);
    T[ii]->SetBranchStatus("H_hoe",1);
    T[ii]->SetBranchStatus("H_s9s25",1);
    T[ii]->SetBranchStatus("H_eopOut",1);
    T[ii]->SetBranchStatus("H_see",1);
    T[ii]->SetBranchStatus("L_deltaEta",1);
    T[ii]->SetBranchStatus("L_deltaPhi",1);
    T[ii]->SetBranchStatus("L_hoe",1);
    T[ii]->SetBranchStatus("L_s9s25",1);
    T[ii]->SetBranchStatus("L_see",1);
    T[ii]->SetBranchStatus("L_eopOut",1);
    T[ii]->SetBranchAddress("H_deltaEta",&H_deltaEta);
    T[ii]->SetBranchAddress("H_deltaPhi",&H_deltaPhi);
    T[ii]->SetBranchAddress("H_hoe",&H_hoe);
    T[ii]->SetBranchAddress("H_s9s25",&H_s9s25);
    T[ii]->SetBranchAddress("H_see",&H_see);
    T[ii]->SetBranchAddress("H_eopOut",&H_eopOut);
    T[ii]->SetBranchAddress("L_deltaEta",&L_deltaEta);
    T[ii]->SetBranchAddress("L_deltaPhi",&L_deltaPhi);
    T[ii]->SetBranchAddress("L_hoe",&L_hoe);
    T[ii]->SetBranchAddress("L_s9s25",&L_s9s25);
    T[ii]->SetBranchAddress("L_see",&L_see);
    T[ii]->SetBranchAddress("L_eopOut",&L_eopOut);
  }
  
  // electron class
  int theClass = atoi(argv[3]);
  setScanValue(theClass);

  // kinematical / preselection efficiencies
  float sgnPreEleIDEff = atof(argv[4]);    
  float sgnKineEff     = atof(argv[5]);   
  float bkgPreEleIDEvt = atof(argv[6]);    
  float bkgKineEff     = atof(argv[7]);     

  // counters
  float passedEleID[6][5][5][3][3][3][3][3][2];
  for(int ii=0; ii<2; ii++){
    for(int iiDeta=0; iiDeta<6; iiDeta++){
      for(int iiDphi=0; iiDphi<5; iiDphi++){
	for(int iiHoE=0; iiHoE<5; iiHoE++){
	  for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	    for(int iiEoPmin=0; iiEoPmin<3; iiEoPmin++){
	      for(int iiEoPmax=0; iiEoPmax<3; iiEoPmax++){
		for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		  for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		    passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPmin][iiEoPmax][iiSeemin][iiSeemax][ii]=0.;
		  }}}}}}}}}

  
  // loop: signal / background samples
  for(int ii=0; ii<2; ii++){
    
    // reading the tree
    float nEnt = T[ii]->GetEntries();
    cout << endl;
    cout << "Total number of events in loop for sample " << ii << " is " << nEnt << endl; 
    for (int entry=0; entry<nEnt; entry++) { 
      if (entry%1000==0) cout << "sample " << ii << ", entry " << entry << endl;
      T[ii] -> GetEntry(entry);
      
      // scan to compute the efficiencies for each bin
      bool theHighScan = false;
      bool theLowScan  = false;
      for(int iiDeta=0; iiDeta<6; iiDeta++){
	for(int iiDphi=0; iiDphi<5; iiDphi++){
	  for(int iiHoE=0; iiHoE<5; iiHoE++){
	    for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	      for(int iiEoPmin=0; iiEoPmin<3; iiEoPmin++){
		for(int iiEoPmax=0; iiEoPmax<3; iiEoPmax++){
		  for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		    for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		      theHighScan=isEleIDScan(H_deltaEta, H_deltaPhi, H_hoe, H_s9s25, H_eopOut, H_see, iiDeta, iiDphi, iiHoE, iiS9S25, iiEoPmin, iiEoPmax, iiSeemin, iiSeemax);
		      theLowScan=isEleIDScan(L_deltaEta, L_deltaPhi, L_hoe, L_s9s25, L_eopOut, L_see, iiDeta, iiDphi, iiHoE, iiS9S25, iiEoPmin, iiEoPmax, iiSeemin, iiSeemax);
      		      if(theHighScan && theLowScan){ 
			passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPmin][iiEoPmax][iiSeemin][iiSeemax][ii]=passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPmin][iiEoPmax][iiSeemin][iiSeemax][ii]+1.; 
		      }
		    }}}}}}}}
    } // loop over entries
  } // loop over signal / background 


  // maximization:
  ofstream *outTxtFile = new ofstream("outputFileScan.txt",ios::app);

  float signPunziMax = -999.;
  float effMax       = -999.;
  int signBinMax[8];
  for(int yy=0; yy<8; yy++){signBinMax[yy]=-1;}  

  for(int iiDeta=0; iiDeta<6; iiDeta++){
    for(int iiDphi=0; iiDphi<5; iiDphi++){
      for(int iiHoE=0; iiHoE<5; iiHoE++){
	for(int iiS9S25=0; iiS9S25<3; iiS9S25++){
	  for(int iiEoPmin=0; iiEoPmin<3; iiEoPmin++){
	    for(int iiEoPmax=0; iiEoPmax<3; iiEoPmax++){
	      for(int iiSeemin=0; iiSeemin<3; iiSeemin++){
		for(int iiSeemax=0; iiSeemax<3; iiSeemax++){
		  float thisBinSgnEff = passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPmin][iiEoPmax][iiSeemin][iiSeemax][0]/((float)T[0]->GetEntries());
		  float thisBinBkgEff = passedEleID[iiDeta][iiDphi][iiHoE][iiS9S25][iiEoPmin][iiEoPmax][iiSeemin][iiSeemax][1]/((float)T[1]->GetEntries());
		  float effSgn        = sgnPreEleIDEff*thisBinSgnEff*sgnKineEff;
		  float bkgEvents     = bkgPreEleIDEvt*thisBinBkgEff*bkgKineEff;
		  float sqrtB         = sqrt(bkgEvents);
		  float signPunzi     = effSgn/(0.5+sqrtB);

		  // saving the full output
		  *outTxtFile << iiDeta << " " << iiDphi << " " << iiHoE << " " << iiS9S25 << " " << iiEoPmin << " " << iiEoPmax << " " << iiSeemin << " " << iiSeemax << " " << thisBinSgnEff << " " << thisBinBkgEff << " " << signPunzi << endl;
		  
		  // looking for the maximum
		  if (signPunzi>signPunziMax) { 
		    signPunziMax  = signPunzi; 
		    effMax        = thisBinSgnEff;
		    signBinMax[0] = iiDeta;
		    signBinMax[1] = iiDphi;
		    signBinMax[2] = iiHoE;
		    signBinMax[3] = iiS9S25;
		    signBinMax[4] = iiEoPmin;
		    signBinMax[5] = iiEoPmax;
		    signBinMax[6] = iiSeemin;
		    signBinMax[7] = iiSeemax;
		  }
		}}}}}}}}
  
  // max significance bin
  float DEtaMaxCut   = DEtaMaxInit   + signBinMax[0]*DEtaMaxStep;
  float DPhiMaxCut   = DPhiMaxInit   + signBinMax[1]*DPhiMaxStep;
  float HOverEMaxCut = HOverEMaxInit + signBinMax[2]*HOverEMaxStep;
  float S9S25MinCut  = S9S25MinInit  + signBinMax[3]*S9S25MinStep;
  float EoPoutMinCut = EoPoutMinInit + signBinMax[4]*EoPoutMinStep;
  float EoPoutMaxCut = EoPoutMaxInit + signBinMax[5]*EoPoutMaxStep;
  float SeeMinCut    = SeeMinInit    + signBinMax[6]*SeeMinStep;
  float SeeMaxCut    = SeeMaxInit    + signBinMax[7]*SeeMaxStep;

  // output
  cout << endl;
  cout << "highest significance (Punzi) = " << signPunziMax << endl;
  cout << "eff eleID signal = " << passedEleID[signBinMax[0]][signBinMax[1]][signBinMax[2]][signBinMax[3]][signBinMax[4]][signBinMax[5]][signBinMax[6]][signBinMax[7]][0]/((float)T[0]->GetEntries()) << endl;
  cout << "in the following bin: " << endl;
  cout << "|deta| = " << "0 - " << DEtaMaxCut   << endl; 
  cout << "|dphi| = " << "0 - " << DPhiMaxCut   << endl;
  cout << "H/E  = "   << "0 - " << HOverEMaxCut << endl;
  cout << "S9/S25 = " << S9S25MinCut  << " - 1" << endl; 
  cout << "E/Pout = " << EoPoutMinCut << " - " << EoPoutMaxCut << endl;
  cout << "See = "    << SeeMinCut    << " - " << SeeMaxCut    << endl;
}


void setScanValue(int theClass){

  // golden barrel
  if (theClass==0){     
    // starting point           // step
    DEtaMaxInit   = 0.002;      DEtaMaxStep   = 0.001;
    DPhiMaxInit   = 0.002;      DPhiMaxStep   = 0.003;   
    HOverEMaxInit = 0.02;       HOverEMaxStep = 0.01;   
    S9S25MinInit  = 0.7;        S9S25MinStep  = 0.1; 
    SeeMaxInit    = 0.01;       SeeMaxStep    = 0.003;   
    SeeMinInit    = 0.0;        SeeMinStep    = 0.003;   
    EoPoutMaxInit = 1.9;        EoPoutMaxStep = 0.03;    
    EoPoutMinInit = 0.7;        EoPoutMinStep = 0.2;  
  }

  // golden endcap
  if (theClass==1){     
    DEtaMaxInit   = 0.002;      DEtaMaxStep   = 0.001;
    DPhiMaxInit   = 0.002;      DPhiMaxStep   = 0.003;   
    HOverEMaxInit = 0.02;       HOverEMaxStep = 0.01;   
    S9S25MinInit  = 0.7;        S9S25MinStep  = 0.1; 
    SeeMaxInit    = 0.01;       SeeMaxStep    = 0.003;   
    SeeMinInit    = 0.0;        SeeMinStep    = 0.003;   
    EoPoutMaxInit = 1.9;        EoPoutMaxStep = 0.03;    
    EoPoutMinInit = 0.7;        EoPoutMinStep = 0.2;  
  }

  // showering barrel
  if (theClass==2){     
    DEtaMaxInit   = 0.002;      DEtaMaxStep   = 0.001;
    DPhiMaxInit   = 0.025;      DPhiMaxStep   = 0.01;   
    HOverEMaxInit = 0.02;       HOverEMaxStep = 0.01;   
    S9S25MinInit  = 0.7;        S9S25MinStep  = 0.1; 
    SeeMaxInit    = 0.01;       SeeMaxStep    = 0.003;   
    SeeMinInit    = 0.0;        SeeMinStep    = 0.003;   
    EoPoutMaxInit = 999.;       EoPoutMaxStep = 1.0;    
    EoPoutMinInit = 0.7;        EoPoutMinStep = 0.2;  
  }

  // showering endcap
  if (theClass==3){     
    DEtaMaxInit   = 0.002;      DEtaMaxStep   = 0.001;
    DPhiMaxInit   = 0.025;      DPhiMaxStep   = 0.01;   
    HOverEMaxInit = 0.02;       HOverEMaxStep = 0.01;   
    S9S25MinInit  = 0.7;        S9S25MinStep  = 0.1; 
    SeeMaxInit    = 0.014;      SeeMaxStep    = 0.003;   
    SeeMinInit    = 0.002;      SeeMinStep    = 0.003;   
    EoPoutMaxInit = 999.;       EoPoutMaxStep = 1.0;    
    EoPoutMinInit = 0.7;        EoPoutMinStep = 0.2;  
  }
}

bool isEleIDScan(float thisDeta, float thisDphi, float thisHoE, float thisS9s25, float thisEopOut, float thisSee, int deta, int dphi, int hoe, int s9s25, int eopmin, int eopmax, int seemin, int seemax) {

  // the current cut
  float DEtaMaxCut   = DEtaMaxInit   + deta*DEtaMaxStep;
  float DPhiMaxCut   = DPhiMaxInit   + dphi*DPhiMaxStep;
  float HOverEMaxCut = HOverEMaxInit + hoe*HOverEMaxStep;
  float S9S25MinCut  = S9S25MinInit  + s9s25*S9S25MinStep;
  float EoPoutMaxCut = EoPoutMaxInit + eopmax*EoPoutMaxStep;
  float EoPoutMinCut = EoPoutMinInit + eopmin*EoPoutMinStep;
  float SeeMaxCut    = SeeMaxInit    + seemax*SeeMaxStep;
  float SeeMinCut    = SeeMinInit    + seemin*SeeMinStep;
  
  bool isPassed = false;
  if ( thisDeta<=DEtaMaxCut ){
    if ( thisDphi<=DPhiMaxCut ){
      if ( fabs(thisHoE)<=HOverEMaxCut ){
	if ( thisS9s25>=S9S25MinCut ) {
	  if( (thisEopOut<=EoPoutMaxCut) && (thisEopOut>=EoPoutMinCut) ){
	    if( (thisSee<=SeeMaxCut) && (thisSee>=SeeMinCut) ){
	      isPassed = true;
	    }}}}}}

    return isPassed;
}
