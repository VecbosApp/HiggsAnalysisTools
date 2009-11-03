//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef HiggsMLSelection_h
#define HiggsMLSelection_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "HiggsAnalysisTools/include/CommonHiggsPreselector.hh"
#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include "HiggsAnalysisTools/include/HiggsBase.h"
#include "HiggsAnalysisTools/include/RedHiggsTree.h"
#include "HiggsAnalysisTools/include/RedTriggerTree.hh"
#include "HiggsAnalysisTools/include/RedEleIDTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>


class HiggsMLSelection : public HiggsBase{
public:
  
  //! constructor
  HiggsMLSelection(TTree *tree=0);
  //! destructor
  virtual ~HiggsMLSelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies(std::string filename);
  //! set the list of the required triggers
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  
private:

  bool findMcTree(const char* processType);

  //! get the two hardest electrons with opposite charge
  std::pair<int,int> getBestElectronPair();
  //! get the two hardest muons with opposite charge 
  std::pair<int,int> getBestMuonPair();
  //! set the 4 vectors, invariant mass, etc. after preselections
  void setKinematics();  
  //! set the 4 vectors, invariant mass, etc. before preselections
  void setPreselKinematics();  
  //! count jet multiplicity
  int numJets();
  int numUncorrJets();
  //! returns the output of the custom cut electron ID
  bool isEleID(int eleIndex);
  //! if the 2nd ele falls in deltaR from first, get its Pt in tracker
  float getSecondEleTkPt(TVector3 firstLepton, int second, float deltaR);
  //! if the 2nd muon falls in deltaR from first, get its Pt in tracker
  float getSecondMuonTkPt(TVector3 firstLepton, int second, float deltaR);
  //! get ECAL isolation sum
  float getEcalPtSum(int index);
  //! get global isolation sum for 2 muons case
  float muonIsoGlobalSum(int theMuon, int theOther);
  //! get global isolation sum for 2 electrons case
  float electronIsoGlobalSum(int theElectron, int theOther);
  //! get global isolation sum for electrons in the electron-muon case
  float elemuIsoGlobalSum(int theElectron, int theOtherMu); 
  //! get global isolation sum for muons in the electron-muon case
  float mueleIsoGlobalSum(int theMuon, int theOtherEle); 
  //! get the kFactor of the event
  float getkFactor(std::string process);
  //! reset the kinematic quantities at the beginning of event
  void resetKinematics();
  //! search for the hardest lepton vertex
  int getPV();
  //! dxy parameter with respect to PV for electron tracks
  double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
  //! dsz parameter with respect to PV for electron tracks
  double trackDszPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
  //! methods for the jet veto: track quality
  bool isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits);
  //! methods for the jet veto: per jet variable
  std::vector<float> jetBTagVariables(int jetIndex);
  //! methods for the jet veto: event based variable
  void calcEventBVetoVariables(std::vector<int> jets);
  //! method to evaluate Mt from the lepton and neutrino pT's
  double mT(TVector3 plep, TVector3 pneu);
  //! method to evaluate Mt2 from the 2 leptons pT's and pTmiss
  double mT2(TVector3 plep1, TVector3 plep2, TVector3 ptmiss);

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  //! to evaluate preselection efficiency
  Selection *_preselection;
  CommonHiggsPreselector CommonHiggsPreselection;
  //! to evaluate full selection efficiency
  Selection *_selectionEE, *_selectionMM, *_selectionEM;
  CutBasedHiggsSelector CutBasedHiggsSelectionEE;
  CutBasedHiggsSelector CutBasedHiggsSelectionMM;
  CutBasedHiggsSelector CutBasedHiggsSelectionEM;
  //! be verbose during runtime
  bool _verbose;
  //! the list of required triggers
  vector<int> m_requiredTriggers;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  //! an integer defining the sub-channel
  enum { ee=0, mm=1, em=2 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[3];
  bool isOk[3];

  //! kinematics of the event
  int theElectron,  thePositron;
  int theMuonMinus, theMuonPlus;
  TLorentzVector *m_p4ElectronPlus, *m_p4ElectronMinus;
  TLorentzVector *m_p4MuonPlus, *m_p4MuonMinus;
  TLorentzVector *m_p4MET;
  float m_HoEElectronMinus, m_HoEElectronPlus;
  float m_CaloEneElectronMinus, m_CaloEneElectronPlus;
  float m_deltaPhi[3];
  float m_deltaErre[3];
  float m_mll[3];
  float m_transvMass[3];
  float m_mT2[3];
  //! used for ee final state
  float hardestElectronPt, hardestMuonPt;
  //! used for mm final state
  float slowestElectronPt, slowestMuonPt;    
  //! used for mixed final state
  float hardestLeptonPt, slowestLeptonPt;
  
  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  int m_closestPV;
  
  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  //! vectors to store indices of best candidates
  std::vector<int> *_bestElectrons;
  std::vector<int> *_bestMuons;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;

  //! reduced tree for event selection (on at least 2leptons events)
  RedHiggsTree *myOutTreeEE;
  RedHiggsTree *myOutTreeMM;
  RedHiggsTree *myOutTreeEM;

  //! reduced tree for trigger studies (on all events)
  RedTriggerTree *myTriggerTree;

  //! reduced tree for eleId studies
  RedEleIDTree *myEleIdTree;

  //! new variables
  float m_eOverP[100];

  float _highestPtGen[1], _lowestPtGen[1];
  float _genHiggsPt[1];
  float _nGenJet[1];
  float _emFracEle[50], _hadFracEle[50];

  //! name of rootfile with dataset
  std::string _datasetName;

  //! to check the electron/jet matching
  TFile *fMatch;
  TH1F *H_deltaRcorr;
  TH1F *H_deltaRuncorr;
};
#endif
