/// The Higgs class is an auxiliary class which contains basic
/// functionality useful for any analysis of Vecbos+jets events.
/// It derives from HiggsBase.

#ifndef Higgs_h
#define Higgs_h

#include "HiggsAnalysisTools/include/HiggsBase.h"
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

private:
  ///goodRUN/LS list
  runsLSSegmentsMap goodRunLS; 
  std::string jsonFile;

  std::string lastFile;
  std::vector<std::string> requiredTriggers;

protected:
  //! the list of required triggers
  std::vector<int> m_requiredTriggers;


};

#endif
