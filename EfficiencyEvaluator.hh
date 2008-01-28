#ifndef EFFICIENCYEVALUATOR_H
#define EFFICIENCYEVALUATOR_H

#include <vector>
#include <TH1F.h>
#include <TFile.h>

class EfficiencyEvaluator {
private:
  std::vector<TH1F*> _numerators;
  TH1F *_denominator;
  std::vector<TH1F*> _efficiencies;
  TFile *_file;
  
public:
  EfficiencyEvaluator(const char* namefile);
  ~EfficiencyEvaluator();
  void AddNumerator(TH1F *numerator) {_numerators.push_back(numerator);};
  void SetDenominator(TH1F *denominator) {_denominator=denominator;};
  void ComputeEfficiencies();
  void Write();
};
#endif
