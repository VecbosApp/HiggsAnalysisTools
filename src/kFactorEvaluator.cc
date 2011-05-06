#include <fstream>
#include <iostream>
#include "TString.h"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"

kFactorEvaluator::kFactorEvaluator(int mH) {
  HiggsMass_ = mH;
  
  ranges.clear();
  kfactors.clear();
  
  char mass[5];
  sprintf(mass,"%d",mH);

  TString path = TString("/afs/cern.ch/user/g/gpetrucc/scratch0/higgs/HqT/HqT/spectra/scalefactor.mh")+
    TString(mass)+TString(".txt");
  std::ifstream file(path.Data());

  if(!file.good()) {
    std::cout << "kFactorEvaluator::Error!   Unable to open the file " << path.Data() << std::endl;
  }
  else {
    std::cout << "KFactors read from " << path.Data() << "..." << std::endl;
    while(!file.eof()) {
      float min, max, val;
      file >> min >> max >> val;
      range bin;
      bin.min = min;
      bin.max = max;
      ranges.push_back(bin);
      kfactors.push_back(val);
    }
  }

  std::cout << kfactors.size() << " values read." << std::endl;

}

// new k-factor for WW, computed by Giovanni Petrucciani 
float kFactorEvaluator::evaluate(float ptH) {
  for(int i=0; i<(int)kfactors.size(); i++) {
    float min = (ranges[i]).min;
    float max = (ranges[i]).max;
    if(ptH>=min && ptH<max) return kfactors[i];
  }
  //  std::cout << "H pT = " << ptH << ": out of range, use 1" << std::endl;
  return 1;
}

