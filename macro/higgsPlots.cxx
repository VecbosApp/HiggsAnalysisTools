#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

std::vector<float> expectedEvents() {

  // normalization lumi (pb-1)
  float lumi = 100.0;

  // for the normal processes, exp ev = L_int * xsec * efficiency
  vector<float> efficiency;
  vector<float> xsec;

  efficiency.push_back(0.0037); // signal H165 (denominator: all 2l2nu)
  efficiency.push_back(0.000039); // WW inclusive (lep + had)
  efficiency.push_back(0.000016); // WZ (incl) -- very poor stat: resubmit CRAB 
  efficiency.push_back(0.000008);  // ZZ (incl)
  efficiency.push_back(0.000010); // tW (incl)

  // xsecs in pb
  // http://ceballos.web.cern.ch/ceballos/hwwlnln/cross-sections_csa07analysis.txt
  xsec.push_back(2.36); // signal H165 pb-1 -> WW -> 2l2nu
  xsec.push_back(114.3); // WW inclusive (lep + had)
  xsec.push_back(49.9); // WZ (incl) 
  xsec.push_back(15.3); // ZZ (incl)
  xsec.push_back(62.0); // tW (incl)

  vector<float> expEv;
  for(int i=0; i< (int) xsec.size(); i++) {
    expEv.push_back( efficiency[i] * xsec[i] * lumi );
  }

  // now evaluate the expected events from Chowder CSA07
  TFile *fileChowderPDElectronSkim = TFile::Open("/cmsrm/pc18/emanuele/releases/Higgs169_new/src/OfflineAnalysis/HiggsAnalysisTools/data/ChowderPDElectronSkim-datasetEE.root");
  TTree *treeChowderPDElectronSkim = (TTree*) fileChowderPDElectronSkim->Get("T1");

  // ALPGEN procees id:
  // 1000 + jet multiplicity for W+jets
  // 2000 + jet multiplicity for Z+jets
  // 3000 + jet multiplicity for ttbar

  TH1F *dummyVar = new TH1F("dummyVar","dummyVar",10,0,10000);

  // evaluate W+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=1000 && CSA07processId<2000 && finalSelection)*CSA07weight");
  float expEvWj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  // evaluate Z+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=2000 && CSA07processId<3000 && finalSelection)*CSA07weight");
  float expEvZj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  // evaluate ttbar expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=3000 && CSA07processId<4000 && finalSelection)*CSA07weight");
  float expEvttbar = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  expEv.push_back(expEvWj);
  expEv.push_back(expEvZj);
  expEv.push_back(expEvttbar);

  for (int i=0; i< (int) expEv.size(); i++) {
    std::cout << "process n. " << i << "\texpected events in " << lumi << " pb-1:\t"
	      << expEv[i] << std::endl;
  }
  
  return expEv;

}
