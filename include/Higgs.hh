/// The Higgs class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from HiggsBase.

#ifndef Higgs_h
#define Higgs_h

#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "HiggsAnalysisTools/include/HiggsBase.h"
// ROOT includes
#include <TLorentzVector.h>
#include <TVector3.h>
// std includes
#include <string>
#include <vector>
#include <map>

class Higgs : public HiggsBase{

public:
  typedef std::pair<unsigned int,unsigned int> aLSSegment;
  typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
  typedef unsigned int aRun;
  typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
  typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

  /// Class Constructor
  Higgs(TTree *tree=0);
  /// Class Destructor
  virtual ~Higgs();

  /// Fill RunLSMap according to json file
  void fillRunLSMap();
  /// Set Good Run LS
  void setJsonGoodRunList(const std::string& jsonFilePath);
  /// check if Run/LS is a good one
  bool isGoodRunLS();
  /// reload TriggerMask if necessary (data file is changed). Should be called for each event inside the event loop
  bool reloadTriggerMask(bool newVersion=false);
  /// set the list of required trigger to produce the bitmask
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers);
  //check if the event passed HLT. To be called per event
  bool hasPassedHLT();
  //get the value of the requested bits
  std::vector<int> getHLTOutput();

  /// Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }

  // useful electron functions
  /// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiEiE(int electron);
  /// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
  float SigmaiPiP(int electron);
  // get the likelihood electron ID
  float likelihoodRatio(int eleIndex, ElectronLikelihood &lh);
  // get the PFjet ID
  bool isPFJetID(double eta, double nHFrac, double nEmFrac, int nConst, double chHFrac, double chMult, double chEmFrac, int WP);

  enum jetIdWP { none=0, loose=1, medium=2, tight=3 };

private:
  ///goodRUN/LS list
  runsLSSegmentsMap goodRunLS; 
  std::string jsonFile;

  std::string lastFile;
  std::vector<std::string> requiredTriggers;

protected:
  //! the list of required triggers
  std::vector<int> m_requiredTriggers;

  /// calculate transverse mass
  /// definitions in http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=104213
  float mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met);

};

#endif
