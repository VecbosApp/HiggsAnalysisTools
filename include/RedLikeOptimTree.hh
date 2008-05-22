#ifndef RedLikeOptimTree_h
#define RedLikeOptimTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedLikeOptimTree {
public:
   RedLikeOptimTree(const char * filename = "eleID.root");
  ~RedLikeOptimTree();

  //! event by event final dataset fill
  void fillAll(float fhlike, float fllike);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  float myH_Likelihood;
  float myL_Likelihood;

  TFile* myFile;
  TTree* myTree;
};

#endif // RedLikeOptimTree_h
