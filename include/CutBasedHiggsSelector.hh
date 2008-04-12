#ifndef CutBasedHiggsSelector_h
#define CutBasedHiggsSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CutBasedHiggsSelector {

public:

  //! constructor
  CutBasedHiggsSelector();

  //! destructor
  virtual ~CutBasedHiggsSelector();   

  //! configure from files
  void Configure(const char *configDir);

  //! set event by event observables
  void SetWeight(float weight)                { m_weight        = weight;    }
  void SetHighElePt(float highPt)             { m_highPt        = highPt;    }
  void SetLowElePt(float lowPt)               { m_lowPt         = lowPt;     }
  void SetInvMass(float mll)                  { m_invMass       = mll;       }
  void SetElectronId(bool isEleId)            { m_isElectronId  = isEleId; }
  void SetPositronId(bool isPosId)            { m_isPositronId  = isPosId; }
  void SetEleTrackerPtSum(float eleTkPtSum)   { m_eleTkPtSum    = eleTkPtSum; }
  void SetPosTrackerPtSum(float posTkPtSum)   { m_posTkPtSum    = posTkPtSum; }
  void SetEleCaloPtSum(float eleCaloPtSum)    { m_eleCaloPtSum  = eleCaloPtSum; }
  void SetPosCaloPtSum(float posCaloPtSum)    { m_posCaloPtSum  = posCaloPtSum; }
  void SetJetVeto(float passedCJV)            { m_passedJetVeto = passedCJV; }
  void SetMet(float met)                      { m_met           = met;}
  void SetDeltaPhi(float deltaPhi)            { m_deltaPhi      = deltaPhi;}
  void SetDetaLeptons(float deltaEta)         { m_detaLeptons   = deltaEta;}

  //! get output of the selector
  bool output();

  //! display the electron efficiency
  void diplayEfficiencies();

private:
  
  float m_weight;
  float m_highPt, m_lowPt;
  bool m_isElectronId, m_isPositronId;
  float m_invMass;
  float m_eleTkPtSum, m_eleCaloPtSum;
  float m_posTkPtSum, m_posCaloPtSum;
  bool m_passedJetVeto;
  float m_met, m_deltaPhi, m_detaLeptons;
  float m_maxPtElectron, m_minPtElectron;

  //! contains the preselection cuts
  Selection* _selection;

  //! counters for the efficiencies display, based on electron candidates
  Counters higgsSelCounter;

};

#endif
