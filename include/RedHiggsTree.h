#ifndef RedHiggsTree_h
#define RedHiggsTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedHiggsTree {
public:
   RedHiggsTree(const char * filename = "eleID.root");
  ~RedHiggsTree();

  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! add the k-Factor (used for signal only)
  void addKFactor();
  //! event by event final dataset fill
  void fillAll(float mt, float dphi, float tmass, float mee, float max, float min, float deta, 
	       bool finalLeptons, bool jetVeto, bool preDeltaPhi, bool finalSelection);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! fill with the k-Factor (used for signal only)
  void fillKFactor(double kfactor);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  float myMet;       
  float myDeltaPhi;  
  float myTransvMass;
  float myEleInvMass;
  float maxPtEle;  
  float minPtEle;  
  float myDetaLeptons;
  bool myFinalLeptons;
  bool myJetVeto;
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
