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
#include "HiggsAnalysisTools/include/Higgs.hh"
#include "HiggsAnalysisTools/include/CommonHiggsPreselector.hh"
#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include "HiggsAnalysisTools/include/RedHiggsTree.h"
#include "HiggsAnalysisTools/include/RedTriggerTree.hh"
#include "HiggsAnalysisTools/include/RedEleIDTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>


class HiggsMLSelection : public Higgs{
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
  //! count the soft muons
  int numSoftMuons();
  //! count the extra leptons (id, iso, d0,acceptance etc) with pt>10 GeV
  int numExtraLeptons();
  //! returns the output of the custom cut electron ID
  void isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);
  //! returns the output of the custom muon ID
  void isMuonID(int muonIndex, bool *muonIdOutput);
  //! if the 2nd ele falls in deltaR from first, get its Pt in tracker
  float getSecondEleTkPt(TVector3 firstLepton, int second, float deltaR);
  //! if the 2nd muon falls in deltaR from first, get its Pt in tracker
  float getSecondMuonTkPt(TVector3 firstLepton, int second, float deltaR);
  //! get ECAL isolation sum
  float getEcalPtSum(int index);
  //! get global isolation sum for 2 muons case
  float muonIsoGlobalSum(int theMuon, int theOther);
  //! get global isolation sum for muons in the electron-muon case
  float mueleIsoGlobalSum(int theMuon, int theOtherEle); 
  //! get the kFactor of the event
  float getkFactor(std::string process);
  //! reset the kinematic quantities at the beginning of event
  void resetKinematics();
  //! set the electron ID variables to dump
  void setEleIdVariables(int hard, int slow);
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
  //! Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }
  //! for the two leptons, find the closest one to MET in phi. if(deltaPhi(lepton-MET)<pi/2) projectedMET = MET * sin(deltaPhi(lepton-MET)); else projectedMET = MET
  float GetProjectedMet(TVector3 p1, TVector3 p2);

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
  TVector3 *m_p3ProjectedMET;
  float m_HoEElectronMinus, m_HoEElectronPlus;
  float m_CaloEneElectronMinus, m_CaloEneElectronPlus;
  float m_deltaPhi[3];
  float m_deltaErre[3];
  float m_mll[3];
  float m_transvMass[3];
  float m_mT2[3];
  float m_projectedMet[3];
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
  
  //! variables for EleID
  int myRecoflag[2];
  float myPt[2], myEta[2], myPhi[2];
  int myClassification[2], myNBremClusters[2];
  float myDeta[2], myDphi[2], myHoe[2], mySee[2], mySpp[2], myEop[2], myFbrem[2];
  float myTrackerIso[2], myHcalIso[2], myEcalJIso[2], myEcalGTIso[2], myCombinedIso[2];
  int myCharge[2];
  int myMissHits[2];
  float myDist[2], myDcot[2];
  float myLh[2];

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
