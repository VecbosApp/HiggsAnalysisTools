#ifndef RedHiggsTree_h
#define RedHiggsTree_h
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

  //! event by event final dataset fill
  void fillAll(float met, float pfmet, float cmet, float projmet, 
	       float dphi, float derre, float tmass, float mee, float max, float min, float deta, 
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);
  //! fill more informations for analysis not cut based
  void fillMLVars(int njets, int nuncorrjets, float dxyEVT, float dszEVT,
                  float bTagTrackCount, float bTagImpPar, float bTagSecVertex);
  //! fill electron ID variables
  void fillElectrons(int recoflag[2], float pt[2], float eta[2], float phi[2],
                     int classification[2], int nbrems[2], float deta[2], float dphi[2], float hoe[2], float see[2], float spp[2], float eop[2], float fbrem[2],
                     float trackerIso[2], float hcalIso[2], float ecalJIso[2], float ecalGTIso[2], float combinedIso[2], int charge[2],
                     int missHits[2], float dist[2], float dcot[2], float lh[2], int matched[2]);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! fill with the k-Factor (used for signal only)
  void fillKFactor(double kfactor);
  //! fill the MC truth informations
  void fillMcTruth(bool prompt);
  //! fill the HLT electron triggers informations
  void fillHLTElectrons(bool singleEle, bool singleEleRelaxed, bool singleEleOR);
  //! fill the HLT muons triggers informations
  void fillHLTMuons(bool singleMuon, bool singleMuonRelaxed, bool singleMuonOR);
  //! fill the run,lumi, event number
  void fillRunInfos(int run, int lumi, int event);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
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
  float myDetaLeptons;
  int myNjets;
  int myNuncorrjets;
  float myDxyEVT;
  float myDszEVT;
  float myBTagTrackCount;
  float myBTagImpPar;
  float myBTagSecVertex;
  bool myFinalLeptons;
  bool myJetVeto;
  bool myUncorrJetVeto;
  bool myPreDeltaPhi;
  bool myFinalSelection;
  double myWeight;
  double myProcesId;
  float myLumi;
  float myKFactor;
  int myRun, myLS, myEvent;

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
  int myMatched[2];

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedHiggsTree_h
