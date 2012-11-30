#include "TMath.h"

class FitSelection {
public:
  FitSelection();
  ~FitSelection() {}
  float mrmin, mrmax, dphimin, dphimax, mtmin, mtmax;
};

FitSelection::FitSelection() {
  mrmin=0.0;
  mrmax=10000.;
  dphimin=0.0;
  dphimax=TMath::Pi();
  mtmin=80.;
  mtmax=10000.;
}
