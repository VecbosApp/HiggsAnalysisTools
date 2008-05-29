//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Original code:
//    CutAnaHiggs_2e2nu.cpp
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "HiggsAnalysisTools/include/Application.hh"
#include "CommonTools/include/TriggerMask.hh"
#if Application == 1
#include "HiggsAnalysisTools/src/HiggsSelection.cc"
#endif
#if Application == 2
#include "HiggsAnalysisTools/src/ElectronID.cc"
#endif
#if Application == 3
#include "HiggsAnalysisTools/src/plotsEleID.cc"
#endif
#if Application == 4
#include "HiggsAnalysisTools/src/ZSelection.cc"
#endif
#if Application == 5
#include "HiggsAnalysisTools/src/WSelection.cc"
#endif
#if Application == 6
#include "HiggsAnalysisTools/src/ClassEfficiencyStudy.cc"
#endif
#if Application == 7
#include "HiggsAnalysisTools/src/WplusJets.cc"
#endif
#if Application == 8
#include "HiggsAnalysisTools/src/ZTandProbe.cc"
#endif
#if Application == 9
#include "HiggsAnalysisTools/src/WewkSelection.cc"
#endif
#if Application == 10
#include "HiggsAnalysisTools/src/ZewkSelection.cc"
#endif
#if Application == 11
#include "HiggsAnalysisTools/src/HiggsEleIdOptimToyMC.cc"
#endif
#if Application == 12
#include "HiggsAnalysisTools/src/ZplusJetsSelection.cc"
#endif
#if Application == 13
#include "HiggsAnalysisTools/src/HiggsIsolationOptimToyMC.cc"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert inputFile with list of root files" << std::endl; 
    return 1;
  }
  strcpy(inputFileName,argv[1]);

  // -------------------------
  // loading file:
  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];  
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  // get the tree with the conditions from the first file
  TTree *treeCond = new TTree();

  int nfiles=1;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
	if ( nfiles==1 ) {
	  TFile *firstfile = TFile::Open(MyRootFile);
	  treeCond = (TTree*)firstfile->Get("Conditions");
	}
	std::cout << "chaining " << MyRootFile << std::endl;
	nfiles++;
      }
  }
  inputFile->close();
  delete inputFile;

#if Application == 1

  HiggsSelection htoww(theChain);
  std::string outFileName(inputFileName);
  outFileName+=".root";
  htoww.SetDatasetName(outFileName);

  TriggerMask mask(treeCond);

  // require triggers for ee channel
  mask.requireTrigger("HLT1Electron");
  mask.requireTrigger("HLT1ElectronRelaxed");
  mask.requireTrigger("HLT2Electron");
  mask.requireTrigger("HLT2ElectronRelaxed");
  mask.requireTrigger("HLT1MuonIso");
  mask.requireTrigger("HLT1MuonNonIso");
  mask.requireTrigger("HLT2MuonNonIso");
  mask.requireTrigger("HLTXElectronMuon");
  mask.requireTrigger("HLTXElectronMuonRelaxed");

  std::vector<int> requiredTriggers = mask.getBits();
  htoww.requireTrigger(requiredTriggers);

  htoww.Loop();
  htoww.displayEfficiencies();

#endif

#if Application == 2
  ElectronID eidanalyzer(theChain);
  eidanalyzer.Loop();
  eidanalyzer.displayEfficiencies();
#endif

#if Application == 3
  plotsEleID eidplots(theChain);
  eidplots.Loop();
#endif

#if Application == 4  
  ZSelection ztoee(theChain);
  ztoee.Loop();
  ztoee.displayEfficiencies();
#endif

#if Application == 5  
  WSelection wtoenu(theChain);
  wtoenu.Loop();
  wtoenu.displayEfficiencies();
#endif

#if Application == 6  
  ClassEfficiencyStudy classEff(theChain);
  classEff.Loop();
  classEff.displayEfficiencies();
#endif

#if Application == 7  
  WplusJets wjets(theChain);
  wjets.Loop();
  wjets.displayEfficiencies();
#endif

#if Application == 8
  ZTandProbe zTaP(theChain);
  zTaP.Loop();
  zTaP.displayEfficiencies();
  zTaP.writeHistos();
#endif

#if Application == 9  
  WewkSelection wtoenu(theChain);
  wtoenu.Loop();
  wtoenu.displayEfficiencies();
#endif

#if Application == 10  
  ZewkSelection ztoee(theChain);

  TriggerMask mask(treeCond);
  mask.requireTrigger("HLT1Electron");
  mask.requireTrigger("HLT1ElectronRelaxed");
  mask.requireTrigger("HLT2Electron");
  mask.requireTrigger("HLT2ElectronRelaxed");

  std::vector<int> requiredTriggers = mask.getBits();

  ztoee.requireTrigger(requiredTriggers);

  ztoee.Loop();
  ztoee.displayEfficiencies();
#endif

#if Application == 11

  HiggsEleIdOptimToyMC heleIdtoy(theChain);
  std::string outFileName(inputFileName);
  outFileName+=".root";

  TriggerMask mask(treeCond);

  // require triggers for ee channel
  mask.requireTrigger("HLT1Electron");
  mask.requireTrigger("HLT1ElectronRelaxed");
  mask.requireTrigger("HLT2Electron");
  mask.requireTrigger("HLT2ElectronRelaxed");

  std::vector<int> requiredTriggers = mask.getBits();
  heleIdtoy.requireTrigger(requiredTriggers);

  heleIdtoy.Loop();

#endif

#if Application == 12  
  ZplusJetsSelection zjetstoee(theChain);

  TriggerMask mask(treeCond);
  mask.requireTrigger("HLT1Electron");
  mask.requireTrigger("HLT1ElectronRelaxed");
  mask.requireTrigger("HLT2Electron");
  mask.requireTrigger("HLT2ElectronRelaxed");

  std::vector<int> requiredTriggers = mask.getBits();

  zjetstoee.requireTrigger(requiredTriggers);

  zjetstoee.Loop();
  zjetstoee.writeHistos();
  zjetstoee.displayEfficiencies();
#endif

#if Application == 13

  HiggsIsolationOptimToyMC hisoltoy(theChain);
  std::string outFileName(inputFileName);
  outFileName+=".root";
  TriggerMask mask(treeCond);
  mask.requireTrigger("HLT1Electron");
  mask.requireTrigger("HLT1ElectronRelaxed");
  mask.requireTrigger("HLT2Electron");
  mask.requireTrigger("HLT2ElectronRelaxed");
  std::vector<int> requiredTriggers = mask.getBits();
  hisoltoy.requireTrigger(requiredTriggers);
  hisoltoy.Loop();

#endif


  return 0;

}
