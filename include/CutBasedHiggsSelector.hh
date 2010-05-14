#ifndef CutBasedHiggsSelector_h
#define CutBasedHiggsSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CutBasedHiggsSelector {

public:

  //! constructor
  CutBasedHiggsSelector();

  //! copy constructor
  CutBasedHiggsSelector( const CutBasedHiggsSelector& selector );

  //! destructor
  virtual ~CutBasedHiggsSelector();   

  //! configure from files
  void Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle);

  //! configure pre-selection (it is not necessarily applied)
  void AppendPreselection(Selection *preselection) { _selection->append(preselection); }

  //! get the applied selection
  Selection* GetSelection() { return _selection; }

  //! set event by event observables
  void SetProcessID(int processID)       { m_processID     = processID; }
  void SetWeight(float weight)           { m_weight        = weight;    }
  void SetHighElePt(float highPt)        { m_highPt        = highPt;    }
  void SetLowElePt(float lowPt)          { m_lowPt         = lowPt;     }
  void SetInvMass(float mll)             { m_invMass       = mll;       }
  void SetElectronId(bool isEleId)       { m_isElectronId  = isEleId; }
  void SetPositronId(bool isPosId)       { m_isPositronId  = isPosId; }
  void SetElectronIsolation(bool isEleIsol)   { m_isElectronIsol  = isEleIsol; }
  void SetPositronIsolation(bool isPosIsol)   { m_isPositronIsol  = isPosIsol; }
  void SetElectronConvRejection(bool isEleConvRej)   { m_isElectronConvRej  = isEleConvRej; }
  void SetPositronConvRejection(bool isPosConvRej)   { m_isPositronConvRej  = isPosConvRej; }
  void SetEleHardTrackerPtSum(float sum) { m_eleHardTkPtSum    = sum; }
  void SetEleSlowTrackerPtSum(float sum) { m_eleSlowTkPtSum    = sum; }
  void SetEleHardHcalPtSum(float sum)    { m_eleHardHcalPtSum  = sum; }
  void SetEleSlowHcalPtSum(float sum)    { m_eleSlowHcalPtSum  = sum; }
  void SetEleHardEcalPtSum(float sum)    { m_eleHardEcalPtSum  = sum; }
  void SetEleSlowEcalPtSum(float sum)    { m_eleSlowEcalPtSum  = sum; }
  void SetEleHardGlobalSum(float sum)    { m_eleHardGlobalSum  = sum; }
  void SetEleSlowGlobalSum(float sum)    { m_eleSlowGlobalSum  = sum; }
  void SetEleSlowD0(float d0)            { m_eleSlowD0 = d0;}
  void SetEleHardD0(float d0)            { m_eleHardD0 = d0;}
  void SetJetVeto(bool passedCJV)        { m_passedJetVeto = passedCJV; }
  void SetUncorrJetVeto(bool passedUncorrCJV) { m_passedUncorrJetVeto = passedUncorrCJV; }
  void SetNJets(int njets)               { m_nJets         = njets; }
  void SetNUncorrJets(int nuncorrjets)   { m_nUncorrJets   = nuncorrjets; }
  void SetMet(float met)                 { m_met           = met;}
  void SetDeltaPhi(float deltaPhi)       { m_deltaPhi      = deltaPhi;}
  void SetDetaLeptons(float deltaEta)    { m_detaLeptons   = deltaEta;}

  //! get output of the selector
  bool output();
  //! the methods outputUpToFinalLeptons and outputUpToJetVeto need 
  //! output() run before using them
  //! get output of the selector until lepton ID + isolation
  bool outputUpToFinalLeptons() { return m_finalLeptons; }
  //! get output of the selector until jet veto
  bool outputUpToJetVeto() {return m_jetVeto; }
  bool outputUpToUncorrJetVeto() {return m_uncorrJetVeto; }
  //! get output of the selector previous to deltaPhi cut
  bool outputPreDeltaPhi() { return m_preDeltaPhi; }

  //! display the electron efficiency
  void diplayEfficiencies(std::string datasetName);

private:
  
  float m_weight;
  float m_highPt, m_lowPt;
  bool m_isElectronId, m_isPositronId, m_isElectronIsol, m_isPositronIsol, m_isElectronConvRej, m_isPositronConvRej;
  float m_invMass;
  float m_eleHardTkPtSum, m_eleHardHcalPtSum, m_eleHardEcalPtSum, m_eleHardGlobalSum;
  float m_eleSlowTkPtSum, m_eleSlowHcalPtSum, m_eleSlowEcalPtSum, m_eleSlowGlobalSum;
  float m_eleSlowD0, m_eleHardD0;
  bool m_passedJetVeto;
  bool m_passedUncorrJetVeto;
  int m_nJets;
  int m_nUncorrJets;
  float m_met, m_deltaPhi, m_detaLeptons;
  float m_maxPtElectron, m_minPtElectron;
  int m_processID;

  //! contains the preselection cuts
  Selection* _selection;

  //! counters for the efficiencies display, based on electron candidates
  Counters* globalCounter;

  //! true if the selection arrived to lepton ID and isolation
  bool m_finalLeptons;
  //! true if the selection arrived to jet veto
  bool m_jetVeto;
  //! true if the selection arrived to jet veto with uncorrected jets
  bool m_uncorrJetVeto;
  //! true if passed all the selections previous to deltaPhi
  bool m_preDeltaPhi;

  //! this is to do an efficiency for each process in the sample 
  //! (if more than one is present)
  //! to turn on it, use SetProcessID(int processID) with processID=!-1
  std::map<int, Counters*> multiProcessCounter;

};

#endif
