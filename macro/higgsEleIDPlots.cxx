// macro to get the number of exected events 
// ----
// efficiencies are hardcoded for 2 selections:
//   1) full selection
//   2) after electronID + isolation 
// ----
// lumi to normalize is hardcoded (100 pb-1)
// ----
// usage: 
// root -b
// .L macro/higgsEleIDPlots.cxx++
// expectedEvents("eleID"):          after lepton selections
// expectedEvents("finalSelection"): after the full selection 

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

std::vector<float> expectedEvents(const char *selection) {

  // normalization lumi (pb-1)
  float lumi = 100.0;

  // for the normal processes, exp ev = L_int * xsec * efficiency
  vector<float> efficiency;         // the final one - tight eleID from egamma
  vector<float> efficiencyAfterEle; // after CJV - tight eleID from egamma

  /* 
  // Mh = 160 
  // here put signal efficiency after full selection (denominator: all 2l2nu)
  // efficiency.push_back(147.000/51566.000);       // loose eGamma eleID
  // efficiency.push_back(112.000/51566.000);       // tight eGamma eleID
  // efficiency.push_back(102.000/51566.000);       // H->WW AN eleID
  // efficiency.push_back(124.000/51566.000);       // new, loose for BB & narrow
  // efficiency.push_back(121.000/51566.000);       // new, tight for BB & narrow

  // here put signal efficiency after ele selection (denominator: all 2l2nu)
  // efficiencyAfterEle.push_back(2732.000/51566.000);       // loose eGamma eleID
  // efficiencyAfterEle.push_back(2034.000/51566.000);       // tight eGamma eleID
  // efficiencyAfterEle.push_back(1704.000/51566.000);       // H->WW AN eleID
  // efficiencyAfterEle.push_back(2157.000/51566.000);       // new, loose for BB & narrow
  // efficiencyAfterEle.push_back(2133.000/51566.000);       // new, tight for BB & narrow
  */

  // Mh = 190 
  // here put signal efficiency after full selection (denominator: all 2l2nu)
  // efficiency.push_back(33.000/19238.000);       // loose eGamma eleID
  // efficiency.push_back(27.000/19238.000);       // tight eGamma eleID
  // efficiency.push_back(27.000/19238.000);       // H->WW AN eleID
  efficiency.push_back(30.000/19238.000);       // new, loose for BB & narrow
  // efficiency.push_back(30.000/19238.000);       // new, tight for BB & narrow

  // here put signal efficiency after ele selection (denominator: all 2l2nu)
  // efficiencyAfterEle.push_back(1102.000/19238.000);       // loose eGamma eleID
  // efficiencyAfterEle.push_back(827.000/19238.000);        // tight eGamma eleID
  // efficiencyAfterEle.push_back(708.000/19238.000);        // H->WW AN eleID
  efficiencyAfterEle.push_back(860.000/19238.000);        // new, loose for BB & narrow
  // efficiencyAfterEle.push_back(854.000/19238.000);        // new, tight for BB & narrow

  /*
  // Mh = 160 with eleID optimized for Mh = 190
  // here put signal efficiency after full selection (denominator: all 2l2nu)
  // efficiency.push_back(116.000/51566.000);      
  // efficiencyAfterEle.push_back(2062.000/51566.000);       
  */

  // xsecs in pb
  vector<float> xsec;
  // xsec.push_back(2.36);                      // signal H160 pb-1 -> WW -> 2l2nu
  xsec.push_back(1.51);                      // signal H190 pb-1 -> WW -> 2l2nu - taken from the note

  // expected signal events
  vector<float> expEv;
  for(int i=0; i< (int) xsec.size(); i++) {
    if(strcmp(selection,"finalSelection")==0)  expEv.push_back( efficiency[i] * xsec[i] * lumi );
    if(strcmp(selection,"eleID")==0)           expEv.push_back( efficiencyAfterEle[i] * xsec[i] * lumi );
  }
  
  // now evaluate the expected events from Chowder CSA07
  TFile *fileChowderPDElectronSkim = 0;
  // fileChowderPDElectronSkim = TFile::Open("../results_mh190/higgsPDElectronChowder_okPresel-datasetEE-eGammaLooseEleID.root");       // loose eGamma eleID
  // fileChowderPDElectronSkim = TFile::Open("../results_mh190/higgsPDElectronChowder_okPresel-datasetEE-eGammaTightEleID.root");       // tight eGamma eleID
  // fileChowderPDElectronSkim = TFile::Open("../results_mh190/higgsPDElectronChowder_okPresel-datasetEE-hwwAnEleID.root");             // H->WW AN eleID
  fileChowderPDElectronSkim = TFile::Open("../results_mh190/higgsPDElectronChowder_okPresel-datasetEE-newEleIDothersLoose.root");    // new, loose for BB & narrow
  // fileChowderPDElectronSkim = TFile::Open("../results_mh190/higgsPDElectronChowder_okPresel-datasetEE-newEleIDothersTight.root");    // new, tight for BB & narrow
  //
  // fileChowderPDElectronSkim = TFile::Open("../results_mh160/conEleIdOptimA190/higgsPDElectronChowder_okPresel-datasetEE.root");

  TTree *treeChowderPDElectronSkim = (TTree*) fileChowderPDElectronSkim->Get("T1");

  // ALPGEN procees id:
  // 1000 + jet multiplicity for W+jets
  // 2000 + jet multiplicity for Z+jets
  // 3000 + jet multiplicity for ttbar

  // evaluate W+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  TH1F *dummyVar = new TH1F("dummyVar","dummyVar",10,0,10000);
  if(strcmp(selection,"finalSelection")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=1000 && CSA07processId<2000 && finalSelection)*CSA07weight");
  if(strcmp(selection,"eleID")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=1000 && CSA07processId<2000 && finalLeptons)*CSA07weight");
  float expEvWj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  // evaluate Z+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  if(strcmp(selection,"finalSelection")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=2000 && CSA07processId<3000 && finalSelection)*CSA07weight");
  if(strcmp(selection,"eleID")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=2000 && CSA07processId<3000 && finalLeptons)*CSA07weight");
  float expEvZj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  expEv.push_back(expEvWj);
  expEv.push_back(expEvZj);

  std::cout << "Summary after the selection: " << selection << std::endl;
  for (int i=0; i< (int) expEv.size(); i++) {
    std::cout << "process n. " << i << "\texpected events in " << lumi << " pb-1:\t"
	      << expEv[i] << std::endl;
  }
  
  return expEv;
  
}

