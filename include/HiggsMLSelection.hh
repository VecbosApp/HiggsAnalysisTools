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
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
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
  //! set the required triggers masks (one per channel)
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  //! set the not-required triggers masks (one per channel)
  void setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  
private:

  bool findMcTree(const char* processType);

  //! get the two hardest electrons with opposite charge after different steps
  std::pair<int,int> getBestElectronPair_acceptance();
  std::pair<int,int> getBestElectronPair_id( std::vector<int> acceptEle );
  std::pair<int,int> getBestElectronPair_isol( std::vector<int> idEle );
  std::pair<int,int> getBestElectronPair_conv( std::vector<int> isolEle );
  std::pair<int,int> getBestElectronPair_ip( std::vector<int> convEle );

  //! get the two hardest muons with opposite charge after different steps
  std::pair<int,int> getBestMuonPair_acceptance();
  std::pair<int,int> getBestMuonPair_id( std::vector<int> acceptMu ); 
  std::pair<int,int> getBestMuonPair_isol( std::vector<int> idMu );
  std::pair<int,int> getBestMuonPair_ip( std::vector<int> isoMu ); 

  //! get the two hardest ele-muon with opposite charge
  std::pair<int,int> getBestEleMuonPair(std::vector<int> electrons, std::vector<int> muons);
  std::pair<int,int> getBestMuonElePair(std::vector<int> electrons, std::vector<int> muons);

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsEE(int myEle, int myPosi);
  void setKinematicsMM(int myMuMinus, int myMuPlus);
  void setKinematicsEMME(int myEle, int myPosi, int myMuPlus, int myMuMinus);
  //! reset the kinematic quantities at the beginning of event and after the selection if needed
  void resetKinematicsStart();
  void resetKinematics();

  //! count jet multiplicity
  int numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) ;
  int numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove );
  //! give the highest b-tag of calojets in the event
  float bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove );
  //! in the 1-jet bin, deltaphi between ll system and leading jet
  float deltaPhiLLJet(int ichan);
  //! count the soft muons
  int numSoftMuons(std::vector<int> muonToRemove);
  //! count the extra leptons (id, iso, d0,acceptance etc) with pt>10 GeV
  int numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove );
  //! returns the output of the custom cut electron ID with WPXX
  void isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID);
  //! returns the output of the custom muon ID
  void isMuonID(int muonIndex, bool *muonIdOutput);
  //! if the 2nd ele falls in deltaR from first, get its Pt in tracker
  float getSecondEleTkPt(TVector3 firstLepton, int second, float deltaR);
  //! if the 2nd muon falls in deltaR from first, get its Pt in tracker
  float getSecondMuonTkPt(TVector3 firstLepton, int second, float deltaR);
  //! get global isolation sum for 2 muons case
  float muonIsoGlobalSum(int theMuon, int theOther);
  //! get global isolation sum for muons in the electron-muon case
  float mueleIsoGlobalSum(int theMuon, int theOtherEle); 
  //! get the kFactor of the event
  float getkFactor(std::string process);
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
  //! reload the trigger mask_s_ (one per channel)
  bool reloadTriggerMask();
  //! get the trigger answer depending on the channel
  bool hasPassedHLT(int channel);

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  CutBasedEleIDSelector EgammaCutBasedIDLow;
  ElectronLikelihood *LH;
  //! to evaluate preselection efficiency
  Selection *_preselection;
  CommonHiggsPreselector CommonHiggsPreselection;
  //! to evaluate full selection efficiency
  Selection *_selectionEE, *_selectionMM, *_selectionEM, *_selectionME;
  CutBasedHiggsSelector CutBasedHiggsSelectionEE;
  CutBasedHiggsSelector CutBasedHiggsSelectionMM;
  CutBasedHiggsSelector CutBasedHiggsSelectionEM;
  CutBasedHiggsSelector CutBasedHiggsSelectionME;
  //! be verbose during runtime
  bool _verbose;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  //! an integer defining the sub-channel
  enum { ee=0, mm=1, em=2, me=3 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[4];
  bool isOk[4];
  
  //! trigger masks
  std::vector<int>  m_requiredTriggersEE, m_requiredTriggersMM, m_notRequiredTriggersMM,
    m_requiredTriggersEM, m_notRequiredTriggersEM;
  std::vector<std::string> requiredTriggersEE, requiredTriggersMM, requiredTriggersEM;
  std::vector<std::string> notRequiredTriggersEE, notRequiredTriggersMM, notRequiredTriggersEM;

  //! kinematics of the event
  int theElectron,  thePositron;
  int theMuonMinus, theMuonPlus;
  int thePreElectron,  thePrePositron, thePreElectronME, thePreElectronEM;
  int thePreMuonMinus, thePreMuonPlus, thePreMuonME, thePreMuonEM;
  int theLeadingJet[4];
  std::vector<int> eleCands[4], muCands[4];
  TLorentzVector *m_p4LeptonPlus[4], *m_p4LeptonMinus[4];
  TVector3 *m_p3PFMET, *m_p3TKMET;
  float m_theMET;

  TVector3 m_dilepPt[4];
  float m_deltaPhi[4];
  float m_deltaErre[4];
  float m_deltaEtaLeptons[4];
  float m_mll[4];
  float m_transvMass[4];
  float m_mT2[4];
  float m_projectedMet[4];
  float m_metOptll[4];
  float hardestLeptonPt[4], slowestLeptonPt[4];
  
  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  int m_closestPV;
  
  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  //! vectors to store indices of best candidates
  std::vector<int> _acceptEleAll, _idEleAll, _isolEleAll, _convEleAll, _ipEleAll;
  std::vector<int> _acceptMuonsAll, _idMuonsAll, _isolMuonsAll, _ipMuonsAll;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;

  //! reduced tree for event selection (on at least 2leptons events)
  RedHiggsTree *myOutTreeEE;
  RedHiggsTree *myOutTreeMM;
  RedHiggsTree *myOutTreeEM;
  RedHiggsTree *myOutTreeME;

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
  int myMatched[2];

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
