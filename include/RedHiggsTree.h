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
  //! event by event final dataset fill
  void fillAll(float mt, float dphi, float tmass, float mee, float max, float min, float deta, 
	       bool finalLeptons, bool jetVeto, bool finalSelection);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
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
  bool myFinalSelection;
  float myWeight;
  int myProcesId;
  float myLumi;

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedHiggsTree_h
