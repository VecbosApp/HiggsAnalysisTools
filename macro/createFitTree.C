#define createFitTree_cxx
#include "createFitTree.h"

#include <vector>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooRealVar.h>

using namespace std;

void createFitTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L createFitTree.C
//      Root > createFitTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


  int myNJets;
  float myMET, myDeltaPhi, myMaxPt, myMinPt, myInvMass;
  float myDxyEVT, myDszEVT;
  float myWeight;

  TTree *myTree = new TTree("data","higgs fit dataset");

  myTree->Branch("nJets",               &myNJets,               "nJets/I");
  myTree->Branch("MET",                 &myMET,                 "MET/F");
  myTree->Branch("deltaPhi",            &myDeltaPhi,            "deltaPhi/F");
  myTree->Branch("maxPt",               &myMaxPt,               "maxPt/F");
  myTree->Branch("minPt",               &myMinPt,               "minPt/F");
  myTree->Branch("invMass",             &myInvMass,             "invMass/F");
  myTree->Branch("dxyEVT",              &myDxyEVT,              "dxyEVT/F");
  myTree->Branch("dszEVT",              &myDszEVT,              "dszEVT/F");
  myTree->Branch("weight",              &myWeight,              "weight/F");
  
  if (fChain == 0) return;
   
  Long64_t nentries = fChain->GetEntriesFast();
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    myNJets = njets;
    myMET = met;
    myMaxPt = maxPtEle;
    myMinPt = minPtEle;
    myInvMass = eleInvMass;
    myDxyEVT = dxyEVT;
    myDszEVT = dszEVT;

    // weight is given as ratio N_sampleX / N_signal
    // at last selection cut
    float n_H = 52.0;
    float n_Wj = 33.0;
    float n_Zj = 284.0; // estimate
    float n_ttj = 548.0;
    float n_WW_2l = 181.0;
    float n_ZZ_4l = 5;
    float n_WZ_3l = 4.; // estimate
    float n_WZ_incl = 6.; // estimate

    if(strcmp(_datasetname,"Higgs")==0)           myWeight = 1.0;
    if(strcmp(_datasetname,"WjetsMADGRAPH")==0)   myWeight = n_Wj/n_H;
    if(strcmp(_datasetname,"ZjetsMADGRAPH")==0)   myWeight = n_Zj/n_H;
    if(strcmp(_datasetname,"ttjetsMADGRAPH")==0)  myWeight = n_ttj/n_H;
    if(strcmp(_datasetname,"WW_2l")==0)           myWeight = n_WW_2l/n_H;
    if(strcmp(_datasetname,"ZZ_4l")==0)         myWeight = n_ZZ_4l/n_H;
    if(strcmp(_datasetname,"WZ_3l")==0)         myWeight = n_WZ_3l/n_H;
    if(strcmp(_datasetname,"WZ_incl")==0)       myWeight = n_WZ_incl/n_H;

    myTree->Fill();

  }

  RooRealVar *nJetsVar = new RooRealVar("nJets","nJets",0);
  RooRealVar *METVar = new RooRealVar("MET","MET",0,200,"GeV");
  RooRealVar *deltaPhiVar = new RooRealVar("deltaPhi","deltaPhi",0,180,"#deg");
  RooRealVar *maxPtVar = new RooRealVar("MaxPt","maxPt",0,200,"GeV");
  RooRealVar *minPtVar = new RooRealVar("minPt","minPt",0,200,"GeV");
  RooRealVar *invMassVar = new RooRealVar("invMass","invMass",0,200,"GeV");
  RooRealVar *dxyEVTVar = new RooRealVar("dxyEVT","dxyEVT",0,1000,"#mum");
  RooRealVar *dszEVTVar = new RooRealVar("dszEVT","dszEVT",0,1000,"#mum");
  RooRealVar *weightVar = new RooRealVar("weight","weight",0,20);
   
  RooArgSet setWenu(*nJetsVar,*METVar,*deltaPhiVar,*maxPtVar,*minPtVar,*invMassVar);
  setWenu.add(*dxyEVTVar);
  setWenu.add(*dszEVTVar);
  setWenu.add(*weightVar);

  // to create the RooDataSets with the merged species as in the fit
  if(!_saveTree) {

    RooDataSet *njetsData = new RooDataSet(_datasetname,_datasetname,myTree,setWenu);
    njetsData->SetName(_datasetname);
    njetsData->SetTitle(_datasetname);
      
    char namefile[200];
    sprintf(namefile,"hww_2e_%s_21X.root",_datasetname);
      
    TFile file(namefile,_option);
    njetsData->Write();
    file.Close();
    
    delete njetsData;
    
  }

  delete nJetsVar;
  delete METVar;
  delete deltaPhiVar;
  delete maxPtVar;
  delete minPtVar;
  delete invMassVar;
  delete dxyEVTVar;
  delete dszEVTVar;
  delete weightVar;

  if(_saveTree) {
    char nameRootTree[200];
    sprintf(nameRootTree,"hww_2e_Tree_%s_21X.root",_datasetname);
    TFile file(nameRootTree,"recreate");
    myTree->Write();
    file.Close();
  }

  delete myTree;

}
