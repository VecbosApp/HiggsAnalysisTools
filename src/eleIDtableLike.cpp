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

// methods
bool isEleIDScan(float thisLike, float cutLike);

int main ( int argc, char **argv)
{
  if (argc < 7){ 
    cout << "Argument missing! Insert: "               << std::endl; 
    cout << "1) inputFile - root tree for sgn "        << std::endl;
    cout << "2) inputFile - root tree for bkg "        << std::endl;
    cout << "3) signal preEleID efficiency"            << std::endl;
    cout << "4) signal kine efficiency"                << std::endl;
    cout << "5) background preEleID # events"          << std::endl;
    cout << "6) background kine efficiency"            << std::endl;
    return 0;
  }


  // reading the input trees --------------------------
  TChain *T[2];
  T[0]= new TChain("T1");
  T[1]= new TChain("T1");
  T[0]->Add(argv[1]);   // signal
  T[1]->Add(argv[2]);   // background
  float H_Likelihood, L_Likelihood;
  for(int ii=0; ii<2; ii++){
    T[ii]->SetMakeClass(1);
    T[ii]->SetBranchStatus("*",0);
    T[ii]->SetBranchStatus("H_Likelihood",1);
    T[ii]->SetBranchStatus("L_Likelihood",1);
    T[ii]->SetBranchAddress("H_Likelihood",&H_Likelihood);
    T[ii]->SetBranchAddress("L_Likelihood",&L_Likelihood);
  }

  // kinematical / preselection efficiencies
  float sgnPreEleIDEff = atof(argv[3]);    
  float sgnKineEff     = atof(argv[4]);   
  float bkgPreEleIDEvt = atof(argv[5]);    
  float bkgKineEff     = atof(argv[6]);     

  // counters
  float passedEleID[100][100][2];
  for(int ii=0; ii<2; ii++){    //  signal/background
    for(int iiHL=0; iiHL<100; iiHL++){
      for(int iiLL=0; iiLL<100; iiLL++){
	passedEleID[iiHL][iiLL][ii]=0.;
      }}}
  
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
      for(int iiHL=0; iiHL<100; iiHL++){
	for(int iiLL=0; iiLL<100; iiLL++){
	  theHighScan=isEleIDScan(H_Likelihood, iiHL);
	  theLowScan =isEleIDScan(L_Likelihood, iiLL);
	  if(theHighScan && theLowScan) { passedEleID[iiHL][iiLL][ii]++; }
	}}
    } // loop over entries
  } // loop over signal / background 

  // maximization:
  ofstream *outTxtFile = new ofstream("outputFileScan.txt",ios::app);
  
  float signPunziMax = -999.;
  float effMax       = -999.;
  int signBinMax[2];
  for(int yy=0; yy<2; yy++){signBinMax[yy]=-1;}  

  for(int iiHL=0; iiHL<100; iiHL++){
    for(int iiLL=0; iiLL<100; iiLL++){
      float thisBinSgnEff = passedEleID[iiHL][iiLL][0]/((float)T[0]->GetEntries());
      float thisBinBkgEff = passedEleID[iiHL][iiLL][1]/((float)T[1]->GetEntries());
      float effSgn        = sgnPreEleIDEff*thisBinSgnEff*sgnKineEff;
      float bkgEvents     = bkgPreEleIDEvt*thisBinBkgEff*bkgKineEff;
      float sqrtB         = sqrt(bkgEvents);
      float signPunzi     = effSgn/(0.5+sqrtB);

      // saving the full output
      *outTxtFile << iiHL << " " << iiLL << " " << thisBinSgnEff << " " << thisBinBkgEff << " " << signPunzi << endl;
		  
      // looking for the maximum
      if (signPunzi>signPunziMax) { 
	signPunziMax  = signPunzi; 
	effMax        = thisBinSgnEff;
	signBinMax[0] = iiHL;
	signBinMax[1] = iiLL;
      }
    }}
  
  // max significance bin
  float highCut = signBinMax[0]*0.01;
  float lowCut  = signBinMax[1]*0.01;

  // output
  cout << endl;
  cout << "highest significance (Punzi) = " << signPunziMax << endl;
  cout << "eff eleID signal = " << passedEleID[signBinMax[0]][signBinMax[1]][0]/((float)T[0]->GetEntries()) << endl;
  cout << "in the following bin: "  << endl;
  cout << "high Pt electron cut = " << highCut << endl; 
  cout << "low Pt electron cut = "  << lowCut  << endl; 
}

bool isEleIDScan(float thisLike, float cutLike) {

  bool isPassed = false;
  float theCut = cutLike*0.01;
  if ( thisLike>=theCut ){ isPassed = true; }
  return isPassed;
}
