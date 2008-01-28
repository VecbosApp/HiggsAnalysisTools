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
#include "Application.hh"
#if Application == 1
#include "HiggsSelection.cc"
#endif
#if Application == 2
#include "ElectronID.cc"
#endif
#if Application == 3
#include "plotsEleID.cc"
#endif
#if Application == 4
#include "ZSelection.cc"
#endif
#if Application == 5
#include "WSelection.cc"
#endif
#if Application == 6
#include "ClassEfficiencyStudy.cc"
#endif
#if Application == 7
#include "WplusJets.cc"
#endif
#if Application == 8
#include "ZTandProbe.cc"
#endif
#if Application == 9
#include "WewkSelection.cc"
#endif
#if Application == 10
#include "ZewkSelection.cc"
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
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	theChain->Add(MyRootFile);
	std::cout << "chaining " << MyRootFile << std::endl;
      }
  }
  inputFile->close();
  delete inputFile;

#if Application == 1  
  HiggsSelection htoww(theChain);
  std::string outFileName(inputFileName);
  outFileName+=".root";
  htoww.SetDatasetName(outFileName);
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
  ztoee.Loop();
  ztoee.displayEfficiencies();
#endif

  return 0;

}
