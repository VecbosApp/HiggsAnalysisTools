#ifndef RedHiggsTree_h
#define RedHiggsTree_h

#include "TLorentzVector.h"
#include "TVector3.h"

class TFile;
class TTree;

class G3EventProxy;

class RedHiggsTree {
public:
   RedHiggsTree(const char * filename = "eleID.root");
  ~RedHiggsTree();

  //! add more informations for analysis not cut based
  void addMLVars();
  //! add the electron ID+iso variables for the selected best electrons
  void addElectronInfos();
  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! add the k-Factor (used for signal only)
  void addKFactor();
  //! add the MC truth informations
  void addMcTruthInfos();
  //! add the HLT electron triggers informations
  void addHLTElectronsInfos();
  //! add the HLT muon triggers informations
  void addHLTMuonsInfos();
  //! add run,lumi, event number (for data)
  void addRunInfos();
  //! add latinos
  void addLatinos();
  //! add kinematics
  void addKinematics();
  //! add razor variables
  void addRazor();
  //! add variables for W+jets
  void addFake();
  //! add variables to study shape systematics
  void addSystematics();


  //! event by event final dataset fill
  void fillAll(float met, float pfmet, float cmet, float projmet, 
	       float dphi, float derre, float tmass, float mee, float max, float min, float deta, int nvtx,
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);

  void fillAll(float met, float pfmet, float cmet, float projmet, 
	       float dphi, float derre, float tmass, float mee, 
	       float max, float min, float maxEta, float minEta,
	       float deta, int nvtx,
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);

  void fillKinematics(float pxTkMet, float pyTkMet, float pzTkMet,
                      float pxLeadJet[3], float pyLeadJet[3], float pzLeadJet[3],
                      float pxSecJet[3], float pySecJet[3], float pzSecJet[3],
                      float pxL1, float pyL1, float pzL1,
                      float pxL2, float pyL2, float pzL2,
                      int ch[2], float lh[2], float iso[2],
                      TLorentzVector *jetSum, TLorentzVector *uncorrJetSum, TVector3 *pfmet);
  
  void fillSystematics(float scE[2], float r9[2], float ene1, float ene2, int ty1, int ty2,
                       TVector3 *metFromJets, TVector3 *pfMetUp, TVector3 *pfMetDown, float mtUp, float mtDown); 
                  
  void fillRazor(float MTR, float mR, float gammaMR);

  void fillFake(int ntigh, float puwst );

  //! fill more informations for analysis not cut based
  void fillMLVars(int njets, int nuncorrjets, float dxyEVT, float dszEVT,
                  float bTagTrackCount, float bTagImpPar, float bTagSecVertex, int nSoftMu, float leadJetBTagTrackCount, float subleadJetBTagTrackCount, int numExtraLep);
  //! fill electron ID variables
  void fillElectrons(int recoflag[2], float pt[2], float eta[2], float phi[2],
                     int classification[2], int nbrems[2], float deta[2], float dphi[2], float hoe[2], float see[2], float spp[2], float eop[2], float fbrem[2],
                     float trackerIso[2], float hcalIso[2], float ecalJIso[2], float ecalGTIso[2], float combinedIso[2], int charge[2],
                     int missHits[2], float dist[2], float dcot[2], float lh[2], int matched[2]);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! fill with the k-Factor (used for signal only)
  void fillKFactor(float kfactor, float genh, float ptlj );
  //! fill the MC truth informations
  void fillMcTruth(bool prompt);
  //! fill the HLT electron triggers informations
  void fillHLTElectrons(bool singleEle, bool singleEleRelaxed, bool singleEleOR);
  //! fill the HLT muons triggers informations
  void fillHLTMuons(bool singleMuon, bool singleMuonRelaxed, bool singleMuonOR);
  //! fill the run,lumi, event number
  void fillRunInfos(int run, int lumi, int event, float puweight, float puwst, bool HLT);   // used in fake estimate
  //! latinos 
  void fillLatinos(bool s0, bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16, bool s17,
                   bool s18, bool s19, bool s20, bool s21, bool s22, bool s23, bool s24);

  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  bool myHLT;
  bool myPromptDecay;
  bool myHLTSingleElectron;
  bool myHLTSingleElectronRelaxed;
  bool myHLTSingleElectronOR;
  bool myHLTSingleMuon;
  bool myHLTSingleMuonRelaxed;
  bool myHLTSingleMuonOR;
  float myMet;       
  float myPFMet;       
  float myCaloMet;
  float myProjectedMet;
  float myDeltaPhi;  
  float myDeltaR;  
  float myTransvMass;
  float myEleInvMass;
  float maxPtEle;  
  float minPtEle;  
  float maxEtaEle;  
  float minEtaEle;  
  float myDetaLeptons;
  int myNVtx;
  int myNjets;
  int myNuncorrjets;
  float myDxyEVT;
  float myDszEVT;
  float myBTagTrackCount;
  float myBTagImpPar;
  float myBTagSecVertex;
  int myNSoftMu;
  float myLeadingJetBTagTrackCount, mySubleadingJetBTagTrackCount;
  float myNumExtraLep;
  bool myFinalLeptons;
  bool myJetVeto;
  bool myUncorrJetVeto;
  bool myPreDeltaPhi;
  bool myFinalSelection;
  double myWeight;
  double myProcesId;
  float myLumi;
  float myKFactor, myPUWeight;
  float myGenHPt;
  float myLeadingJetPt;
  int myRun, myLS, myEvent;
  float myPxTkMet, myPyTkMet, myPzTkMet;
  float myPxLeadJet[3], myPyLeadJet[3], myPzLeadJet[3];
  float myPxSecondJet[3], myPySecondJet[3], myPzSecondJet[3];
  float myPxL1, myPyL1, myPzL1;
  float myPxL2, myPyL2, myPzL2;
  float myEneL1, myEneL2;
  int myTypeL1, myTypeL2;

  float myMTR, myMR, myGammaMR;
  
  TLorentzVector *myJetsSum, *myUncorrJetsSum;
  TVector3 *myPfMet;

  //! for W+jets
  int myTight;
  float myPUWeightSt;

  // latinos
  bool mySteps[25];

  // electron variables
  int myRecoflag[2];
  float myPt[2], myEta[2], myPhi[2];
  int myClassification[2], myNBremClusters[2];
  float myDeta[2], myDphi[2], myHoe[2], mySee[2], mySpp[2], myEop[2], myFbrem[2];
  float myTrackerIso[2], myHcalIso[2], myEcalJIso[2], myEcalGTIso[2], myCombinedIso[2];
  int myCharge[2];
  int myMissHits[2];
  float myDist[2], myDcot[2];
  float myLh[2];
  float myScEnergy[2], myR9[2];
  int myMatched[2];
  TVector3 *myMetFromJets, *myPfMetUp, *myPfMetDown;
  float myMtUp, myMtDown;

  // lepton variables
  int myLepCharge[2];
  float myIso[2], myEleLh[2];

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedHiggsTree_h
