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
  //! event by event final dataset fill
  void fillAll(float mt, float dphi, float tmass, float mee, float max, float min, float deta, 
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);
  //! fill more informations for analysis not cut based
  void fillMLVars(float maxlh, float minlh, int njets, int nuncorrjets, float dxyEVT, float dszEVT,
                  float bTagTrackCount, float bTagImpPar, float bTagSecVertex);
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
  float myDeltaPhi;  
  float myTransvMass;
  float myEleInvMass;
  float maxPtEle;  
  float minPtEle;  
  float myDetaLeptons;
  float myMaxPtLh;
  float myMinPtLh;
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

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedHiggsTree_h
