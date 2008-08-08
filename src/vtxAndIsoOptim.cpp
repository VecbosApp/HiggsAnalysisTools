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

// ancora da fare: nella zuppa selezionare gli eventi secondo il processo e pesarli 

// variables to be used later on
float trackerCut,   hcalCut,   ecalCut,  dzCut,  dxySignCut;
float trackerStep,  hcalStep,  ecalStep, dzStep, dxySignStep;
float trackerInit,  hcalInit,  ecalInit, dzInit, dxySignInit;

// methods
void setScanValue();
bool isScan(float thisTracker, float thisHcal, float thisEcal, float thisDz, float thisDxySign, int trackerIsol, int hcalIsol, int ecalIsol, int dzVtx, int dxyVtx);

int main ( int argc, char **argv)
{
  if (argc < 8){
    cout << "Argument missing! Insert: "                                                          << std::endl; 
    cout << "1) inputFile - root tree for sgn"                                                    << std::endl;
    cout << "2) signal pre vtx/isolation efficiency"                                              << std::endl;
    cout << "3) inputFile - text file with chowder backgrounds root file and nevents and kineEff" << std::endl;
    cout << "4) inputFile - text file with other bkgs root file, preEff, xsec and kineEff"        << std::endl;
    cout << "5) signal kine efficiency"                                                           << std::endl; 
    cout << "6) discovery (1) or exclusion (0) limits"                                            << std::endl;
    cout << "7) output scan file name"                                                            << std::endl;
    return 0;
  }

  // the luminosity for the optimization
  float lumi = 1000; // pb-1

  std::string rootFileAlpgen[2], rootFileOthers[1];
  float nEventsAlpgenPreKin[2], alpgenKineEff[2];
  float effPreKinOthers[1], xsecOthers[1], othersKineEff[1];

  // reading the text file for ALPGEN - backgrounds in chowder 
  ifstream inFileAlpgen(argv[3]);
  if( inFileAlpgen.is_open() ) {
    for(int ii=0; ii<2; ii++) {
      inFileAlpgen >> rootFileAlpgen[ii] >> nEventsAlpgenPreKin[ii] >> alpgenKineEff[ii];
      std::cout << "Expected events pre-Vtx/Iso in 1fb-1 of process " << rootFileAlpgen[ii] 
		<< " = " << nEventsAlpgenPreKin[ii] << std::endl;
      std::cout << "Kine efficiency " << alpgenKineEff[ii] << std::endl; 
    }
  }
  else {
    cout << "Unable to open file " << inFileAlpgen << std::endl;
  }

  // reading the text file for other backgrounds
  ifstream inFileOthers(argv[4]);
  if( inFileOthers.is_open() ) {
    for(int ii=0; ii<1; ii++) {
      inFileOthers >> rootFileOthers[ii] >> effPreKinOthers[ii] >> xsecOthers[ii];
      std::cout << "Sample " << rootFileOthers[ii] << " has eff pre vtx/iso = " << effPreKinOthers[ii] 
		<< ", kine efficiency = " << othersKineEff[ii]
		<< " and x-sec = " << xsecOthers[ii] << " pb" << std::endl;
    }
  }
  else {
    cout << "Unable to open file " << inFileOthers << std::endl;
  }

  // reading the input trees --------------------------
  TChain *T[4];
  T[0]= new TChain("T1"); // signal: W+jets, W -> e nu
  T[1]= new TChain("T1"); // W+jets, W -> tau nu
  T[2]= new TChain("T1"); // ttbar
  T[3]= new TChain("T1"); // QCD di-jets

  T[0]->Add(argv[1]);   // signal
  T[1]->Add(rootFileAlpgen[0].c_str());   // background
  T[2]->Add(rootFileAlpgen[1].c_str());   // background
  T[3]->Add(rootFileOthers[0].c_str());   // background

  int WToENuDecay;
  float trackerIsol, hcalIsol, ecalIsol, dzVtx, dxyVtx;
  float CSA07lumi;
  double CSA07weight, CSA07processId;

  for(int ii=0; ii<4; ii++){
    T[ii]->SetMakeClass(1);
    T[ii]->SetBranchStatus("*",0);
    T[ii]->SetBranchStatus("trackerIsol,",1);
    T[ii]->SetBranchStatus("hcalIsol,",1);
    T[ii]->SetBranchStatus("ecalIsol",1);
    T[ii]->SetBranchStatus("dzVtx",1);
    T[ii]->SetBranchStatus("dxyVtx",1);
    T[ii]->SetBranchStatus("WToENuDecay",1);
    T[ii]->SetBranchStatus("CSA07weight",1);
    T[ii]->SetBranchStatus("CSA07processId",1);
    T[ii]->SetBranchStatus("CSA07lumi",1);
    T[ii]->SetBranchAddress("trackerIsol",&trackerIsol);
    T[ii]->SetBranchAddress("hcalIsol", &hcalIsol);
    T[ii]->SetBranchAddress("ecalIsol", &ecalIsol);
    T[ii]->SetBranchAddress("dzVtx", &dzVtx);
    T[ii]->SetBranchAddress("dxyVtx", &dxyVtx);
    T[ii]->SetBranchAddress("WToENuDecay", &WToENuDecay);
    T[ii]->SetBranchAddress("CSA07weight", &CSA07weight);
    T[ii]->SetBranchAddress("CSA07processId", &CSA07processId);
    T[ii]->SetBranchAddress("CSA07lumi", &CSA07lumi);
  }

  // setting the scan
  setScanValue();
  
  // kinematical / preselection efficiencies
  float sgnPreVtxEff = atof(argv[2]);
  float sgnKineEff   = atof(argv[5]);
  
  float nBkgPreVtx[3];
  float bkgKineEff[3];
  for(int ii=0; ii<2; ii++) {      // w+jets (in taus) and ttbar in chowder
    nBkgPreVtx[ii] = nEventsAlpgenPreKin[ii] / 1000. * lumi;
    bkgKineEff[ii] = alpgenKineEff[ii];
    std::cout << "Expected nBkgPreVtx[" << ii << "]= " << nBkgPreVtx[ii] << std::endl;
    std::cout << "Kine efficiency[" << ii << "]= " << bkgKineEff[ii] << std::endl;
  }
  for(int ii=0; ii<1; ii++) {
    nBkgPreVtx[ii+2] = effPreKinOthers[ii] * xsecOthers[ii] * lumi;
    bkgKineEff[ii+2] = othersKineEff[ii];
    std::cout << "Expected nBkgPreVtx[" << ii << "]= " << nBkgPreVtx[ii+2] << std::endl;
    std::cout << "Kine efficiency[" << ii << "]= " << bkgKineEff[ii+2] << std::endl;
  }

  // discovery or exclusion
  int discovery = atoi(argv[6]); 

  // counters
  float passedVtx[10][10][10][10][10][4];
  for(int ii=0; ii<4; ii++){
    for(int iiTracker=0; iiTracker<10; iiTracker++){
      for(int iiHcal=0; iiHcal<10; iiHcal++){
	for(int iiEcal=0; iiEcal<10; iiEcal++){
	  for(int iiDz=0; iiDz<10; iiDz++){
	    for(int iiDxySign=0; iiDxySign<10; iiDxySign++){
	      passedVtx[iiTracker][iiHcal][iiEcal][iiDz][iiDxySign][ii]=0;
	    }}}}}}
  
  // loop: signal / background samples
  for(int ii=0; ii<4; ii++){
    
    // reading the tree
    float nEnt = T[ii]->GetEntries();
    cout << endl;
    cout << "Total number of events in loop for sample " << ii << " is " << nEnt << endl; 
    for (int entry=0; entry<nEnt; entry++) { 
      if (entry%1000==0) cout << "sample " << ii << ", entry " << entry << endl;
      T[ii] -> GetEntry(entry);
      
      // scan to compute the efficiencies for each bin
      bool theScan  = false;
      for(int iiTracker=0; iiTracker<10; iiTracker++){
	for(int iiHcal=0; iiHcal<10; iiHcal++){
	  for(int iiEcal=0; iiEcal<10; iiEcal++){
	    for(int iiDz=0; iiDz<10; iiDz++){
	      for(int iiDxySign=0; iiDxySign<10; iiDxySign++){
		theScan=isScan(trackerIsol, hcalIsol, ecalIsol, dzVtx, dxyVtx, iiTracker, iiHcal, iiEcal, iiDz, iiDxySign);
		if(theScan){ 
		  passedVtx[iiTracker][iiHcal][iiEcal][iiDz][iiDxySign][ii]+=1;
		}
	      }}}}}
    } // loop over entries
  } // loop over signal / backgrounds 
  

  // maximization:
  ofstream *outTxtFile = new ofstream(argv[7],ios::app);
  float signPunziMax = -999.;
  float effMax       = -999.;
  int signBinMax[5];
  for(int yy=0; yy<5; yy++){signBinMax[yy]=-1;}  
  
  for(int iiTracker=0; iiTracker<10; iiTracker++){
    for(int iiHcal=0; iiHcal<10; iiHcal++){
      for(int iiEcal=0; iiEcal<10; iiEcal++){
	for(int iiDz=0; iiDz<10; iiDz++){
	  for(int iiDxySign=0; iiDxySign<10; iiDxySign++){
	    
	    float thisBinSgnEff = passedVtx[iiTracker][iiHcal][iiEcal][iiDz][iiDxySign][0]/((float)T[0]->GetEntries());
	    float thisBinBkgEff[3];
	    for(int ibkg=0; ibkg<3; ibkg++) {
	      thisBinBkgEff[ibkg] = passedVtx[iiTracker][iiHcal][iiEcal][iiDz][iiDxySign][ibkg+1]/((float)T[ibkg+1]->GetEntries());
	    }
		    
	    float effSgn = sgnPreVtxEff*thisBinSgnEff*sgnKineEff;  
 	    float bkgEvents = 0;
	    for(int ibkg=0; ibkg<3; ibkg++) {
	      bkgEvents += nBkgPreVtx[ibkg]*thisBinBkgEff[ibkg]*bkgKineEff[ibkg];  
	    }
		    
	    float sqrtB = sqrt(bkgEvents);
	    float signPunzi;
	    if(discovery==1){  // 5 sigma 
	      signPunzi = effSgn/(2.5+sqrtB);
	    }
	    if(discovery==0){ // 2 sigma
	      signPunzi = effSgn/(1.+sqrtB);
	    }
	    
	    // saving the full output
	    *outTxtFile << iiTracker << " " << iiHcal << " " << iiEcal << " " << iiDz << " " << iiDxySign << " ";
	    
	    for(int ibkg=0; ibkg<3; ibkg++) {
	      *outTxtFile << nBkgPreVtx[ibkg]*thisBinBkgEff[ibkg] << " ";     // Aggiungi la postVtx 
	    }

	    *outTxtFile <<  signPunzi << endl;
	    
	    // looking for the maximum
	    if (signPunzi>signPunziMax) { 
	      signPunziMax  = signPunzi; 
	      effMax        = thisBinSgnEff;
	      signBinMax[0] = iiTracker;
	      signBinMax[1] = iiHcal;
	      signBinMax[2] = iiEcal;
	      signBinMax[3] = iiDz;
	      signBinMax[4] = iiDxySign;
	    }
	  }}}}}
  
  // max significance bin
  float trackerCut = trackerInit  + signBinMax[0]*trackerStep;
  float hcalCut    = hcalInit  + signBinMax[1]*hcalStep;
  float ecalCut    = ecalInit  + signBinMax[2]*ecalStep;
  float dzCut      = dzInit  + signBinMax[3]*dzStep;
  float dxySignCut = dxySignInit  + signBinMax[4]*dxySignStep;
  
  // output
  cout << endl;
  cout << "highest significance (Punzi) = " << signPunziMax << endl;
  cout << "eff vtx/isol signal = "     << passedVtx[signBinMax[0]][signBinMax[1]][signBinMax[2]][signBinMax[3]][signBinMax[4]][0]/((float)T[0]->GetEntries()) << endl;
  for(int i=0; i<3; i++) {
    cout << "eff vtx/iso background[" << i << "] = " << passedVtx[signBinMax[0]][signBinMax[1]][signBinMax[2]][signBinMax[3]][signBinMax[4]][i+1]/((float)T[i+1]->GetEntries()) << endl;
  }
  cout << "in the following bin: "   << endl;
  cout << "tracker = " << trackerCut << endl;
  cout << "hcal = "    << hcalCut    << endl;
  cout << "ecal = "    << ecalCut    << endl;
  cout << "dz = "      << dzCut      << endl;
  cout << "dxySign = " << dxySignCut << endl;

  
}


void setScanValue(){

  // starting point         
  trackerInit = 0.02;
  hcalInit    = 0.02;
  ecalInit    = 0.02;
  dzInit      = 0.;
  dxySignInit = 0.;

  // step for the scan
  trackerStep = 0.005;
  hcalStep    = 0.02;
  ecalStep    = 0.02;
  dzStep      = 0.1;
  dxySignStep = 2.;
}

bool isScan(float thisTracker, float thisHcal, float thisEcal, float thisDz, float thisDxySign, int trackerIsol, int hcalIsol, int ecalIsol, int dzVtx, int dxyVtx) {

  // the current cut
  float trackerCut = trackerInit + trackerIsol*trackerStep;
  float hcalCut    = hcalInit + hcalIsol*hcalStep;
  float ecalCut    = ecalInit + ecalIsol*ecalStep;
  float dzCut      = dzInit + dzVtx*dzStep;
  float dxySignCut = dxySignInit + dxyVtx*dxySignStep;

  bool isPassed = false;
  if ( (thisTracker<=trackerCut) && (thisTracker>=trackerCut) ) {
    if ( (thisHcal<=hcalCut) && (thisHcal>=hcalCut) ) {
      if ( (thisEcal<=ecalCut) && (thisEcal>=ecalCut) ) {
	if ( (thisDz<=dzCut) && (thisDz>=dzCut) ) {
	  if ( (thisDxySign<=dxySignCut) && (thisDxySign>=dxySignCut) ) {
      	    isPassed = true;
	  }}}}}

  return isPassed;
}
