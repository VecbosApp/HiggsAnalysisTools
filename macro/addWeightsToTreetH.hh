#ifndef ADDWEIGHTSTOTREETH_H
#define ADDWEIGHTSTOTREETH_H

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH2F.h>

class addWeightsToTreetH {

public:

  //! constructor
  addWeightsToTreetH(const char* filename, float baseW, int processId, int finalstate, int release);
  //! destructor
  virtual ~addWeightsToTreetH() {}
  //! the main method
  void addWeights();

private:

  float GetProjectedMet(TVector3 met, TVector3 p1, TVector3 p2);
  float calcMT(TVector3 met, TVector3 lepton);
  float getOfflineEff(float pT, float eta, TH2F *myH);

  std::string filename_;
  float baseW_;
  int processId_, finalstate_, release_;

};
#endif // ADDWEIGHTSTOTREETH_H

