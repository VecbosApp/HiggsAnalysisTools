class TFile;
class TTree;

class G3EventProxy;

class RedEleIDTree {
public:
   RedEleIDTree(const char * filename = "eleID.root");
  ~RedEleIDTree();

  void fillAll(int isok, int seHLT, int serHLT, int deHLT, int derHLT, int charge, float ene, float et, float mom, float theta, float eta, float phi, float lat, float a20, float s9s25, float covEE, int theclass, float hoe, float eop, float eopout, float deta, float dphi, float iso, float like, float fis);

  void store();
  void save();

private:
  int mySampleOk;
  int mySingleElePassedTrg, mySingleEleRelaxPassedTrg;
  int myDoubleElePassedTrg, myDoubleEleRelaxPassedTrg;
  int myChargeEle, myEleClassEle;
  float myEnergyEle, myEtEle, myMomentumEle;            
  float myThetaEle, myEtaEle, myPhiEle;                 
  float myLatEle, myA20Ele, myS9s25Ele, myCovEtaEtaEle;         
  float myEleHoEEle, myEleCorrEoPEle;          
  float myEleCorrEoPoutEle, myEleDeltaEtaAtVtxEle, myEleDeltaPhiAtVtxEle;    
  float myEleTrackerIso_SumPtEle; 
  float myEleLikelihoodEle;       
  float myEleFisherEle;       

  TFile* myFile;
  TTree* myTree;
 
};

