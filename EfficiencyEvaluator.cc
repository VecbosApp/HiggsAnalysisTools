#include <vector>
#include <iostream>
#include <TH1F.h>
#include "EfficiencyEvaluator.hh"

EfficiencyEvaluator::EfficiencyEvaluator(const char* namefile) {
  _file = new TFile(namefile,"RECREATE");
}

EfficiencyEvaluator::~EfficiencyEvaluator() {
  delete _file;
}

void EfficiencyEvaluator::ComputeEfficiencies() {
  
  // total efficiency
  std::vector<TH1F*>::const_iterator numItr;
  for(numItr=_numerators.begin(); numItr!=_numerators.end(); ++numItr) {
    TH1F *eff = (TH1F*) (*numItr)->Clone((std::string((*numItr)->GetName())+"_Eff").c_str());
    eff->Sumw2();
    eff->Divide(*numItr, _denominator, 1, 1);
    for(int i=1;i<=eff->GetNbinsX();++i) {
      float effVal = eff->GetBinContent(i);
      float efferrVal = sqrt(effVal*(1-effVal)/(*numItr)->GetBinContent(i));
      eff->SetBinError(i,efferrVal);
    }
    _efficiencies.push_back(eff);
  }
  
  // partial efficiency
  std::vector<TH1F*>::const_iterator numPreviousItr=_numerators.begin();
  for(numItr=_numerators.begin(); numItr!=_numerators.end(); ++numItr) {
    TH1F *effPartial = (TH1F*) (*numItr)->Clone((std::string((*numItr)->GetName())+"_EffWrtPrevious").c_str());
    effPartial->Sumw2();
    effPartial->Divide(*numItr, *numPreviousItr, 1, 1);
    for(int i=1;i<=effPartial->GetNbinsX();++i) {
      float effVal = effPartial->GetBinContent(i);
      float efferrVal = sqrt(effVal*(1-effVal)/(*numItr)->GetBinContent(i));
      effPartial->SetBinError(i,efferrVal);
    }
    _efficiencies.push_back(effPartial);
    if(numItr!=_numerators.begin()) ++numPreviousItr;
  }

}

void EfficiencyEvaluator::Write() {
  _file->cd();
  _denominator->Write();
  std::vector<TH1F*>::const_iterator numItr;
  for(numItr=_numerators.begin(); numItr!=_numerators.end(); ++numItr) {
    (*numItr)->Write();
  }
  std::vector<TH1F*>::const_iterator effItr;
  for(effItr=_efficiencies.begin(); effItr!=_efficiencies.end(); ++effItr) {
    (*effItr)->Write();
  }
  _file->Close();
}
