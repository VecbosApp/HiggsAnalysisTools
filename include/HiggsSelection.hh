//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Original code:
//    CutAnaHiggs_2e2nu.cpp
//-------------------------------------------------------

#ifndef HiggsSelection_h
#define HiggsSelection_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "HiggsBase.h"
#include "RedHiggsTree.h"

#include <TVector3.h>
#include <TLorentzVector.h>


class HiggsSelection : public HiggsBase{
public:
  
  //! constructor
  HiggsSelection(TTree *tree=0);
  //! destructor
  virtual ~HiggsSelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies();
  //! set the list of the required triggers
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  
private:

  bool findMcTree(const char* processType);

  void bookHistos();
  std::pair<int,int> getBestLeptonPair();
  void setKinematics(int,int);  
  void addVariables();
  bool jetVeto();
  bool preselJetVeto();
  bool isEleID(int eleIndex);
  float Fisher(int eleIndex);
  void estimateJetMatch(float ptmin=10.0);
  float getkFactor(std::string process);


  //! counters event based
  Counters _counter;
  //! counters electron based
  Counters _eleCounter;
  //! initialise event selection from file
  Selection* _selection;
  //! initialise selection for electrons from file (class by class)
  std::vector<Selection*> _electronSelection;
  //! be verbose during runtime
  bool _verbose;
  //! the list of required triggers
  vector<int> m_requiredTriggers;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  TVector3 *_p3Ele, *_p3Pos, *_p3Met;
  TLorentzVector *_p4Ele, *_p4Pos;
  float _maxPt,_minPt;
  int _theGenEle, _theGenPos;

  //! vectors to store indices of best candidates
  std::vector<int> *_bestElectrons;
  std::vector<int> *_bestJets;
  std::vector<int> *_bestGenJets;

  //! vector to store indices of candidate to exclude
  std::vector<int> *_excludedJets;

  //! monitoring tools
  Monitor *_monitorEventAfterSelection, *_monitorMet, *_monitorElectrons, *_monitorJets;
  Monitor *_monitorEventAfterReco;
  Monitor *_monitorGenerator, *_monitorGenJets;

  //! simple histograms
  TH1F *RecoJets_pt, *MatchedJets_pt;
  TH1F *etHighestJet;
  
  //! reduced tree
  RedHiggsTree *myOutTree;

  //! new variables
  float _eOverP[100];
  float _deltaPhi[1];
  float _mll[1];
  float _transvMass[1];
  float _highestPt[1], _lowestPt[1];
  float _nEle[1], _nJet[1];

  float _highestPtGen[1], _lowestPtGen[1];
  float _genHiggsPt[1];
  float _nGenJet[1];
  float _alphaEle[50], _emFracEle[50], _hadFracEle[50];

  //! name of rootfile with dataset
  std::string _datasetName;
};
#endif
