//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef LeptonPlusFakeMLSelection_fullEE_h
#define LeptonPlusFakeMLSelection_fullEE_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "HiggsAnalysisTools/include/Higgs.hh"
#include "HiggsAnalysisTools/include/CommonHiggsPreselector.hh"
#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include "HiggsAnalysisTools/include/RedHiggsTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>


class LeptonPlusFakeMLSelection_fullEE : public Higgs{
public:
  
  //! constructor
  LeptonPlusFakeMLSelection_fullEE(TTree *tree=0);
  //! destructor
  virtual ~LeptonPlusFakeMLSelection_fullEE();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies(std::string filename);
  //! set the required triggers masks (one per channel)
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers);
  //! set the not-required triggers masks (one per channel)
  void setNotRequiredTriggers(const std::vector<std::string>& reqTriggers);
  
private:

  //! get the two hardest electrons with opposite charge after different steps
  std::pair<int,int> getBestElectronPair_acceptance();
  std::pair<int,int> getBestElectronPair_id( std::vector<int> acceptEle );
  std::pair<int,int> getBestElectronPair_isol( std::vector<int> idEle );
  std::pair<int,int> getBestElectronPair_conv( std::vector<int> isolEle );
  std::pair<int,int> getBestElectronPair_ip( std::vector<int> convEle );
  std::pair<int,int> getBestElectronPair_denominator();

  //! fake rates initialization                                         
  void initialiseFakeRate();
  float getFakeRate( float fakePt, bool isEB );
  float getFakeRateError( float fakePt, bool isEE );

  //! prompt rates initialization                                         
  void initialisePromptRate();
  float getPromptRate( float promptPt, float promptEta );
  float getPromptRateError( float promptPt, float promptEta );

  //! fake related selection                
  int getBestDenominator(int realEle);
  bool isDenomFake(int theEle);

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsEE(int myEle, int myPosi);
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
  float deltaPhiLLJet15(int ichan);
  //! count the soft muons
  int numSoftMuons(std::vector<int> muonToRemove);
  //! count the extra leptons (id, iso, d0,acceptance etc) with pt>10 GeV
  int numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove );
  //! returns the output of the custom cut electron ID with WPXX
  void isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID);
  //! PF isolation for muons (depending on pT)
  bool isPFIsolatedMuon(int muonIndex);
  //! returns the output of the custom muon ID
  void isMuonID(int muonIndex, bool *muonIdOutput);
  //! search for the hardest lepton vertex
  int getPV();
  //! methods for the jet veto: track quality
  bool isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits);
  //! methods for the jet veto: per jet variable
  std::vector<float> jetBTagVariables(int jetIndex);
  //! methods for the jet veto: event based variable
  void calcEventBVetoVariables(std::vector<int> jets);
  //! Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }
  //! for the two leptons, find the closest one to MET in phi. if(deltaPhi(lepton-MET)<pi/2) projectedMET = MET * sin(deltaPhi(lepton-MET)); else projectedMET = MET
  float GetProjectedMet(TVector3 p1, TVector3 p2);
  //! reload the trigger mask_s_ (one per channel)
  bool reloadTriggerMask(int runN);
  //! get the trigger answer depending on the channel
  bool hasPassedHLT();

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  CutBasedEleIDSelector EgammaCutBasedIDLow;
  ElectronLikelihood *LH;

  //! to evaluate full selection efficiency
  Selection *_selectionEE,     *_selectionEE_FF,     *_selectionEE_PP;     
  Selection *_selectionStatEE, *_selectionStatEE_FF, *_selectionStatEE_PP;
  Selection *_selectionErrEE,  *_selectionErrEE_FF;

  CutBasedHiggsSelector CutBasedHiggsSelectionEE;
  CutBasedHiggsSelector CutBasedHiggsSelectionStatEE;
  CutBasedHiggsSelector CutBasedHiggsErrorsSelectionEE;
  CutBasedHiggsSelector CutBasedHiggsSelectionEE_FF;
  CutBasedHiggsSelector CutBasedHiggsSelectionStatEE_FF;
  CutBasedHiggsSelector CutBasedHiggsErrorsSelectionEE_FF;
  CutBasedHiggsSelector CutBasedHiggsSelectionEE_PP;
  CutBasedHiggsSelector CutBasedHiggsSelectionStatEE_PP;

  //! be verbose during runtime
  bool _verbose;

  //! an integer defining the sub-channel
  enum { ee=0 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[1];
  bool isOk[1];
  
  //! trigger masks
  std::vector<int>  m_requiredTriggersEE;
  std::vector<std::string> requiredTriggersEE;
  std::vector<std::string> notRequiredTriggersEE;

  //! kinematics of the event
  int theReal,  theFake;
  int thePreElectron,  thePrePositron;
  int theLeadingJet[1];
  std::vector<int> eleCands[1], muCands[1];
  TLorentzVector *m_p4LeptonPlus[1], *m_p4LeptonMinus[1];
  TVector3 *m_p3PFMET;
  float m_theMET;

  TVector3 m_dilepPt[1];
  float m_deltaPhi[1];
  float m_deltaErre[1];
  float m_deltaEtaLeptons[1];
  float m_mll[1];
  float m_transvMass[1];
  float m_GammaMR[1], m_MR[1], m_MTR[1];
  float m_mT2[1];
  float m_projectedMet[1];
  float m_chMet[1];
  float m_metOptll[1];
  float hardestLeptonPt[1],  slowestLeptonPt[1];
  float hardestLeptonEta[1], slowestLeptonEta[1];
  float leadJetBtag[1];
  
  //! fake rates                                                                     
  float m_minFakePt[5],  m_maxFakePt[5];
  float m_fakeRateEB[5], m_fakeRateEB_err[5];
  float m_fakeRateEE[5], m_fakeRateEE_err[5];

  //! prompt rates                                                                     
  float m_minPromptPt[5],  m_maxPromptPt[5];
  float m_promptRateEB[5], m_promptRateEB_err[5];
  float m_promptRateEE[5], m_promptRateEE_err[5];
  float m_promptRateCr[5], m_promptRateCr_err[5];

  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  
  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  //! vectors to store indices of best candidates
  std::vector<int> _acceptEleAll, _idEleAll, _isolEleAll, _convEleAll, _ipEleAll, _denomEleAll;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;

  //! reduced tree for event selection (on at least 2leptons events)
  RedHiggsTree *myOutTreeEE;

  //! new variables
  int _massVal;
  float _highestPtGen[1], _lowestPtGen[1];
  float _genHiggsPt[1];
  float _nGenJet[1];
  float _emFracEle[50], _hadFracEle[50];

  //! name of rootfile with dataset
  std::string _datasetName;
};
#endif
