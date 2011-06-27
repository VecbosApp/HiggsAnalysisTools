//-------------------------------------------------------
// Description:
//  Class for selection of reconstructed H->WW->2e2nu
// Authors:
//   Chiara Rovelli & Emanuele Di Marco
//   Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef LeptonPlusFakeMLSelection_fullME_h
#define LeptonPlusFakeMLSelection_fullME_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "HiggsAnalysisTools/include/Higgs.hh"
#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include "HiggsAnalysisTools/include/RedHiggsTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>


class LeptonPlusFakeMLSelection_fullME : public Higgs{
public:
  
  //! constructor
  LeptonPlusFakeMLSelection_fullME(TTree *tree=0);
  //! destructor
  virtual ~LeptonPlusFakeMLSelection_fullME();
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

  //! get the two hardest muons with opposite charge after different steps
  std::pair<int,int> getBestMuonPair_acceptance();
  std::pair<int,int> getBestMuonPair_id( std::vector<int> acceptMu ); 
  std::pair<int,int> getBestMuonPair_isol( std::vector<int> idMu );
  std::pair<int,int> getBestMuonPair_ip( std::vector<int> isoMu ); 
  std::pair<int,int> getBestMuonPair_denominator(); 

  //! get the two hardest ele-muon with opposite charge
  std::pair<int,int> getBestMuonElePair(std::vector<int> electrons, std::vector<int> muons);

  //! fake rates initialization for electrons
  void initialiseElectronFakeRate();
  float getElectronFakeRate( float fakePt, bool isEB );
  float getElectronFakeRateError( float fakePt, bool isEE );

  //! prompt rates initialization for electrons
  void initialiseElectronPromptRate();
  float getElectronPromptRate( float fakePt, bool isEB );
  float getElectronPromptRateError( float fakePt, bool isEE );

  //! fake rates initialization for muons
  void initialiseMuonFakeRate();
  float getMuonFakeRate( float fakePt, bool isEB );
  float getMuonFakeRateError( float fakePt, bool isEE );

  //! prompt rates initialization for muons
  void initialiseMuonPromptRate();
  float getMuonPromptRate( float fakePt, bool isEB );
  float getMuonPromptRateError( float fakePt, bool isEE );

  //! fake related selection
  int getBestEleDenominator(int realMuon);
  int getBestMuDenominator(int realEle);
  bool isEleDenomFake(int theEle);
  bool isMuonDenomFake(int theMuon);

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsME(int myReal, int myFake);
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
  //! retrns the pT dependent isolation 
  bool isPFIsolatedMuon(int muonIndex);
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
  bool reloadTriggerMask();
  //! get the trigger answer depending on the channel
  bool hasPassedHLT();

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  CutBasedEleIDSelector EgammaCutBasedIDLow;
  ElectronLikelihood *LH;

  //! to evaluate full selection efficiency
  Selection *_selectionME,     *_selectionME_FF;
  Selection *_selectionStatME, *_selectionStatME_FF;
  Selection *_selectionErrME,  *_selectionErrME_FF;
  CutBasedHiggsSelector CutBasedHiggsSelectionME;
  CutBasedHiggsSelector CutBasedHiggsSelectionME_FF;
  CutBasedHiggsSelector CutBasedHiggsErrorsSelectionME;
  CutBasedHiggsSelector CutBasedHiggsErrorsSelectionME_FF;
  CutBasedHiggsSelector CutBasedHiggsSelectionStatME;
  CutBasedHiggsSelector CutBasedHiggsSelectionStatME_FF;

  //! be verbose during runtime
  bool _verbose;

  //! an integer defining the sub-channel
  enum { me=0 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[1];
  bool isOk[1];
  
  //! trigger masks
  std::vector<int> m_requiredTriggersME, m_notRequiredTriggersME;
  std::vector<std::string> requiredTriggersME;
  std::vector<std::string> notRequiredTriggersME;

  //! kinematics of the event
  int theReal,  theFake;
  int thePreElectronME,  thePreMuonME;
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
  float m_mT2[1];
  float m_projectedMet[1];
  float m_chMet[1];
  float m_metOptll[1];
  float hardestLeptonPt[1], slowestLeptonPt[1];
  float hardestLeptonEta[1], slowestLeptonEta[1];
  
  //! electron fake rates                             
  float m_eleMinFakePt[5],  m_eleMaxFakePt[5];
  float m_eleFakeRateEB[5], m_eleFakeRateEB_err[5];
  float m_eleFakeRateEE[5], m_eleFakeRateEE_err[5];

  //! electron prompt rates                             
  float m_eleMinPromptPt[5],  m_eleMaxPromptPt[5];
  float m_elePromptRateEB[5], m_elePromptRateEB_err[5];
  float m_elePromptRateEE[5], m_elePromptRateEE_err[5];

  //! muon fake rates                             
  float m_muonMinFakePt[5],  m_muonMaxFakePt[5];
  float m_muonFakeRateEB[5], m_muonFakeRateEB_err[5];
  float m_muonFakeRateEE[5], m_muonFakeRateEE_err[5];

  //! muon prompt rates                             
  float m_muonMinPromptPt[5],  m_muonMaxPromptPt[5];
  float m_muonPromptRateEB[5], m_muonPromptRateEB_err[5];
  float m_muonPromptRateEE[5], m_muonPromptRateEE_err[5];

  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  int m_closestPV;
  
  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  //! vectors to store indices of best candidates
  std::vector<int> _acceptEleAll, _idEleAll, _isolEleAll, _convEleAll, _ipEleAll, _denomEleAll;
  std::vector<int> _acceptMuonsAll, _idMuonsAll, _isolMuonsAll, _ipMuonsAll, _denomMuonAll;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;

  //! reduced tree for event selection (on at least 2leptons events)
  RedHiggsTree *myOutTreeME;

  int _massVal;
  float _highestPtGen[1], _lowestPtGen[1];
  float _nGenJet[1];
  float _emFracEle[50], _hadFracEle[50];

  //! name of rootfile with dataset
  std::string _datasetName;
};
#endif
