#ifndef CommonHiggsPreselector_h
#define CommonHiggsPreselector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CommonHiggsPreselector {

public:

  //! constructor
  CommonHiggsPreselector();

  //! destructor
  virtual ~CommonHiggsPreselector();   

  //! configure from files
  void Configure(const char *configDir);

  //! set event by event observables
  void SetWeight(float kFactor)       { m_kFactor     = kFactor;     }
  void SetMcTruth(bool foundMcTree)   { m_foundMcTree = foundMcTree; }
  void SetHLT(bool passedHLT)         { m_passedHLT   = passedHLT;   }
  void SetNele(int nEle)              { m_nEle        = nEle;        } 
  void SetNmuon(int nMuon)            { m_nMuon       = nMuon;       } 
  void SetIsEE(bool isEE)             { m_isEE        = isEE;        } 
  void SetIsMM(bool isMM)             { m_isMM        = isMM;        } 
  void SetIsEM(bool isEM)             { m_isEM        = isEM;        } 
  void SetHighElePt(float highElePt)  { m_highElePt   = highElePt;   }
  void SetLowElePt(float lowElePt)    { m_lowElePt    = lowElePt;    }
  void SetHighMuonPt(float highMuonPt){ m_highMuonPt  = highMuonPt;  }
  void SetLowMuonPt(float lowMuonPt)  { m_lowMuonPt   = lowMuonPt;   }
  void SetMet(float met)              { m_met         = met;         }
  void SetMllEE(float mllEE)          { m_mllEE       = mllEE;       }
  void SetMllEM(float mllEM)          { m_mllEM       = mllEM;       }
  void SetMllMM(float mllMM)          { m_mllMM       = mllMM;       }
  void SetCommonPres(bool evtPresel)  { m_evtPresel   = evtPresel;   }

  //! get output of the selector
  bool output();

  //! display the electron efficiency
  void diplayEfficiencies();

private:
  
  float m_weight;
  float m_kFactor;
  bool m_foundMcTree;
  bool m_passedHLT;
  int m_nEle, m_nMuon;
  bool m_isEE, m_isEM, m_isMM;
  float m_highElePt,  m_lowElePt;
  float m_highMuonPt, m_lowMuonPt;
  float m_met;
  float m_mllEE, m_mllEM, m_mllMM;
  bool m_evtPresel;
  
  //! contains the preselection cuts
  Selection* _selection;

  //! counters for the efficiencies display, based on electron candidates
  Counters presCounter;

};

#endif
