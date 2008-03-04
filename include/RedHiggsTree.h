#ifndef RedHiggsTree_h
#define RedHiggsTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedHiggsTree {
public:
   RedHiggsTree(const char * filename = "eleID.root");
  ~RedHiggsTree();

  void fillAll(float mt, float dphi, float tmass, float mee, float max, float min, float deta);

  void store();
  void save();

private:
  float myMet;       
  float myDeltaPhi;  
  float myTransvMass;
  float myEleInvMass;
  float maxPtEle;  
  float minPtEle;  
  float myDetaLeptons;

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedHiggsTree_h
