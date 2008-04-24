//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Original code:
//    CutAnaHiggs_2e2nu.cpp
//-------------------------------------------------------

#ifndef HiggsEleIdOptim_h
#define HiggsEleIdOptim_h

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

#include "HiggsBase.h"
//#include "RedHiggsTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class HiggsEleIdOptim : public HiggsBase{
public:
  
  //! constructor
  HiggsEleIdOptim(TTree *tree=0);
  //! destructor
  virtual ~HiggsEleIdOptim();
  //! loop over events
  void Loop();
  //! set the list of the required triggers
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  
private:

  bool findMcTree(const char* processType);

  //! get the two hardest electrons with opposite charge
  std::pair<int,int> getBestElectronPair();
  //! set the 4 vectors, invariant mass, etc. after preselections
  void setKinematics();  
  //! returns the output of the custom cut electron ID
  bool isEleID(int eleIndex);
  bool isEleIDScan(int eleIndex,int iiDeta,int iiDphi,int iiHoE,int iiS9S25min,int iiS9S25max,int iiEoPoutmin,int iiEoPoutmax,int iiSee);
  //! reset the kinematic quantities at the beginning of event
  void resetKinematics();

  //! be verbose during runtime
  bool _verbose;
  //! the list of required triggers
  vector<int> m_requiredTriggers;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  //! an integer defining the sub-channel
  enum { ee=0 };

  //! kinematics of the event
  int theElectron,  thePositron;
  TLorentzVector *m_p4ElectronPlus, *m_p4ElectronMinus;
  float m_mll, hardestElectronPt, slowestElectronPt;
  int _theGenEle, _theGenPos;

  //! vectors to store indices of best candidates
  std::vector<int> *_bestElectrons;

  //! counters
  float theWeight;
  float allEvents, passedMc, triggered, commonPresel; 
  float passedReco, elePresel, looseId, passedIsol;
  float passedEleID[6][5][5][3][3][3][3][3];


};
#endif
