// macro to get the number of exected events 
// and to draw normalized distributions
// ----
// efficiencies are hardcoded for 2 selections:
//   1) full selection
//   2) after the CJV (trigger+reco+iso+ID+CJV)
// ----
// lumi to normalize is hardcoded (100 pb-1)
// ----
// usage: 
// root -b
// .L macros/higgsPlots.cxx++
// drawKinematics("jetVeto"): draws distributions after CJV
// drawKinematics("finalSelection"): draw distributions after the full selection (but stat is poor)
// ----

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
  vector<float> efficiency; // the final one
  vector<float> efficiencyAfterCJV;
  vector<float> xsec;

  efficiency.push_back(0.0037); // signal H165 (denominator: all 2l2nu)
  efficiency.push_back(0.000039); // WW inclusive (lep + had)
  efficiency.push_back(0.000016); // WZ (incl) -- very poor stat: resubmit CRAB 
  efficiency.push_back(0.000008);  // ZZ (incl)
  efficiency.push_back(0.000010); // tW (incl)

  efficiencyAfterCJV.push_back(0.016); // signal H165 (denominator: all 2l2nu)
  efficiencyAfterCJV.push_back(0.0016); // WW inclusive (lep + had)
  efficiencyAfterCJV.push_back(0.0050); // WZ (incl) -- very poor stat: resubmit CRAB 
  efficiencyAfterCJV.push_back(0.0026);  // ZZ (incl)
  efficiencyAfterCJV.push_back(0.00045); // tW (incl)


  // xsecs in pb
  // http://ceballos.web.cern.ch/ceballos/hwwlnln/cross-sections_csa07analysis.txt
  xsec.push_back(2.36); // signal H165 pb-1 -> WW -> 2l2nu
  xsec.push_back(114.3); // WW inclusive (lep + had)
  xsec.push_back(49.9); // WZ (incl) 
  xsec.push_back(15.3); // ZZ (incl)
  xsec.push_back(62.0); // tW (incl)

  vector<float> expEv;
  for(int i=0; i< (int) xsec.size(); i++) {
    if(strcmp(selection,"finalSelection")==0)
      expEv.push_back( efficiency[i] * xsec[i] * lumi );
    if(strcmp(selection,"jetVeto")==0)
      expEv.push_back( efficiencyAfterCJV[i] * xsec[i] * lumi );
  }

  // now evaluate the expected events from Chowder CSA07
  TFile *fileChowderPDElectronSkim = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/ChowderPDElectronSkim-datasetEE.root");
  TTree *treeChowderPDElectronSkim = (TTree*) fileChowderPDElectronSkim->Get("T1");

  // ALPGEN procees id:
  // 1000 + jet multiplicity for W+jets
  // 2000 + jet multiplicity for Z+jets
  // 3000 + jet multiplicity for ttbar

  TH1F *dummyVar = new TH1F("dummyVar","dummyVar",10,0,10000);

  // evaluate W+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  if(strcmp(selection,"finalSelection")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=1000 && CSA07processId<2000 && finalSelection)*CSA07weight");
  if(strcmp(selection,"jetVeto")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=1000 && CSA07processId<2000 && jetVeto)*CSA07weight");
  float expEvWj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  // evaluate Z+jets expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  if(strcmp(selection,"finalSelection")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=2000 && CSA07processId<3000 && finalSelection)*CSA07weight");
  if(strcmp(selection,"jetVeto")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=2000 && CSA07processId<3000 && jetVeto)*CSA07weight");
  float expEvZj = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  // evaluate ttbar expected events (weights were evaluated for 1000pb-1, the equivalent lumi of CSA07 sample)
  if(strcmp(selection,"finalSelection")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=3000 && CSA07processId<4000 && finalSelection)*CSA07weight");
  if(strcmp(selection,"jetVeto")==0)
    treeChowderPDElectronSkim->Project("dummyVar","CSA07lumi","(CSA07processId>=3000 && CSA07processId<4000 && jetVeto)*CSA07weight");
  float expEvttbar = float( dummyVar->Integral() ) * lumi / 1000. ;
  dummyVar->Reset();

  expEv.push_back(expEvWj);
  expEv.push_back(expEvZj);
  expEv.push_back(expEvttbar);

  std::cout << "Summary after the selection: " << selection << std::endl;
  for (int i=0; i< (int) expEv.size(); i++) {
    std::cout << "process n. " << i << "\texpected events in " << lumi << " pb-1:\t"
	      << expEv[i] << std::endl;
  }
  
  return expEv;

}


void drawKinematics(const char* selection) {

  // get the expected events for each process considered
  std::vector<float> expEvents = expectedEvents(selection);

  std::vector<TFile*> datasets;
  TFile* H165 = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/HiggsH165_CMSSW_1_6_9-datasetEE.root");
  datasets.push_back(H165);
  TFile* WW_incl = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/WW_incl-datasetEE.root");
  datasets.push_back(WW_incl);
  TFile* WZ = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/WZ-datasetEE.root");
  datasets.push_back(WZ);
  TFile* ZZ_incl = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/ZZ_incl-datasetEE.root");
  datasets.push_back(ZZ_incl);
  TFile* tW_incl = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/tW_incl-datasetEE.root");
  datasets.push_back(tW_incl);
  TFile* Chowder = TFile::Open("rfio:/castor/cern.ch/user/e/emanuele/Higgs169_new/SelectedDatasets/ChowderPDElectronSkim-datasetEE.root");
  datasets.push_back(Chowder);

  std::vector<TH1F*> met;
  std::vector<TH1F*> mll;
  std::vector<TH1F*> deltaphi;
  std::vector<TH1F*> ptmax;
  std::vector<TH1F*> ptmin;

  TH1F* metH = new TH1F("metH","metH",50,0,200);
  TH1F* mllH = new TH1F("mllH", "mllH",50,0,200);
  TH1F* deltaphiH = new TH1F("deltaphiH", "deltaphiH",50,0,180);
  TH1F* ptmaxH = new TH1F("ptmaxH","ptmaxH",50,0,200);
  TH1F* ptminH = new TH1F("ptminH", "ptminH",50,0,200);

  for(int i=0; i<(int)expEvents.size(); i++) {

    TTree *tree;
    if(i<5)
      tree = (TTree*)datasets[i]->Get("T1");
    else // ALPGEN (Chowder)
      tree = (TTree*)datasets[5]->Get("T1");

    char buf[50];
    sprintf(buf,"met_%d",i);
    TH1F* metProcessX = (TH1F*) metH->Clone(buf);
    metProcessX->Sumw2();
    met.push_back(metProcessX);

    if(i < 5) {
      tree->Project(buf,"met",selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",selection);
      tree->Project(buf,"met",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",selection);
      tree->Project(buf,"met",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",selection);
      tree->Project(buf,"met",extracut);
    }

    met[i]->Scale( expEvents[i]/met[i]->Integral() );

    sprintf(buf,"mll_%d",i);
    TH1F* mllProcessX = (TH1F*) mllH->Clone(buf);
    mllProcessX->Sumw2();
    mll.push_back(mllProcessX);

    if(i < 5) {
      tree->Project(buf,"eleInvMass",selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    mll[i]->Scale( expEvents[i]/mll[i]->Integral() );

    sprintf(buf,"deltaphi_%d",i);
    TH1F* deltaphiProcessX = (TH1F*) deltaphiH->Clone(buf);
    deltaphiProcessX->Sumw2();
    deltaphi.push_back(deltaphiProcessX);


    if(i < 5) {
      tree->Project(buf,"deltaPhi",selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    deltaphi[i]->Scale( expEvents[i]/deltaphi[i]->Integral() );


    sprintf(buf,"ptmax_%d",i);
    TH1F* ptmaxProcessX = (TH1F*) ptmaxH->Clone(buf);
    ptmaxProcessX->Sumw2();
    ptmax.push_back(ptmaxProcessX);

    if(i < 5) {
      tree->Project(buf,"maxPtEle",selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    ptmax[i]->Scale( expEvents[i]/ptmax[i]->Integral() );

    sprintf(buf,"ptmin_%d",i);
    TH1F* ptminProcessX = (TH1F*) ptminH->Clone(buf);
    ptminProcessX->Sumw2();
    ptmin.push_back(ptminProcessX);

    if(i < 5) {
      tree->Project(buf,"minPtEle",selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    ptmin[i]->Scale( expEvents[i]/ptmin[i]->Integral() );


  }
  

  TLegend *leg = new TLegend(0.11,0.65,0.45,0.89);
  leg->SetBorderSize(2);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(met[0],"Signal, m_{H}=160 GeV","pl");
  leg->AddEntry(met[1],"WW","f");
  leg->AddEntry(met[7],"tt","f");
  leg->AddEntry(met[6],"Z+jets","f");
  leg->AddEntry(met[3],"ZZ","f");
  leg->AddEntry(met[2],"WZ","f");
  leg->AddEntry(met[5],"W+jets","f");
  leg->AddEntry(met[4],"tW","f");

  gStyle->SetOptStat(0);

  // draw met
  TCanvas cmet;
  cmet.SetLogy();
  met[6]->SetMaximum(1000);
  met[6]->SetMinimum(0.01);
  met[6]->SetFillColor(6);
  met[6]->SetTitle("");
  met[6]->GetXaxis()->SetTitle("Missing E_{T} [GeV]");
  met[6]->GetYaxis()->SetTitle("normalized Events");
  met[6]->Draw("hist");

  met[1]->SetFillColor(5);
  met[1]->Draw("same hist");

  met[7]->SetFillColor(4);
  met[7]->Draw("same hist");

  met[5]->SetFillColor(2);
  met[5]->Draw("same hist");

  met[4]->SetFillColor(28);
  met[4]->Draw("same hist");

  met[2]->SetFillColor(7);
  met[2]->Draw("same hist");

  met[3]->SetFillColor(3);
  met[3]->Draw("same hist");
  
  met[0]->SetMarkerStyle(8);
  met[0]->Sumw2();
  met[0]->Draw("same pe1");

  leg->Draw();

  cmet.SaveAs("met.eps");
  cmet.SaveAs("met.root");


  // draw mll
  TCanvas cmll;
  cmll.SetLogy();
  mll[6]->SetMaximum(1000);
  mll[6]->SetMinimum(0.01);
  mll[6]->SetFillColor(6);
  mll[6]->SetTitle("");
  mll[6]->GetXaxis()->SetTitle("m_{ll} [GeV]");
  mll[6]->GetYaxis()->SetTitle("normalized Events");
  mll[6]->Draw("hist");

  mll[1]->SetFillColor(5);
  mll[1]->Draw("same hist");

  mll[7]->SetFillColor(4);
  mll[7]->Draw("same hist");

  mll[5]->SetFillColor(2);
  mll[5]->Draw("same hist");

  mll[4]->SetFillColor(28);
  mll[4]->Draw("same hist");

  mll[2]->SetFillColor(7);
  mll[2]->Draw("same hist");

  mll[3]->SetFillColor(3);
  mll[3]->Draw("same hist");
  
  mll[0]->SetMarkerStyle(8);
  mll[0]->Sumw2();
  mll[0]->Draw("same pe1");

  leg->Draw();

  cmll.SaveAs("mll.eps");
  cmll.SaveAs("mll.root");



  // draw deltaphi
  TCanvas cdeltaphi;
  cdeltaphi.SetLogy();
  deltaphi[6]->SetMaximum(1000);
  deltaphi[6]->SetMinimum(0.01);
  deltaphi[6]->SetFillColor(6);
  deltaphi[6]->SetTitle("");
  deltaphi[6]->GetXaxis()->SetTitle("#Delta #phi");
  deltaphi[6]->GetYaxis()->SetTitle("normalized Events");
  deltaphi[6]->Draw("hist");

  deltaphi[1]->SetFillColor(5);
  deltaphi[1]->Draw("same hist");

  deltaphi[7]->SetFillColor(4);
  deltaphi[7]->Draw("same hist");

  deltaphi[5]->SetFillColor(2);
  deltaphi[5]->Draw("same hist");

  deltaphi[4]->SetFillColor(28);
  deltaphi[4]->Draw("same hist");

  deltaphi[2]->SetFillColor(7);
  deltaphi[2]->Draw("same hist");

  deltaphi[3]->SetFillColor(3);
  deltaphi[3]->Draw("same hist");
  
  deltaphi[0]->SetMarkerStyle(8);
  deltaphi[0]->Sumw2();
  deltaphi[0]->Draw("same pe1");

  leg->Draw();

  cdeltaphi.SaveAs("deltaphi.eps");
  cdeltaphi.SaveAs("deltaphi.root");



  // draw ptmax
  TCanvas cptmax;
  cptmax.SetLogy();
  ptmax[6]->SetMaximum(200);
  ptmax[6]->SetMinimum(0.01);
  ptmax[6]->SetFillColor(6);
  ptmax[6]->SetTitle("");
  ptmax[6]->GetXaxis()->SetTitle("P^{e max}_{T} [GeV]");
  ptmax[6]->GetYaxis()->SetTitle("normalized Events");
  ptmax[6]->Draw("hist");

  ptmax[1]->SetFillColor(5);
  ptmax[1]->Draw("same hist");

  ptmax[7]->SetFillColor(4);
  ptmax[7]->Draw("same hist");

  ptmax[5]->SetFillColor(2);
  ptmax[5]->Draw("same hist");

  ptmax[4]->SetFillColor(28);
  ptmax[4]->Draw("same hist");

  ptmax[2]->SetFillColor(7);
  ptmax[2]->Draw("same hist");

  ptmax[3]->SetFillColor(3);
  ptmax[3]->Draw("same hist");
  
  ptmax[0]->SetMarkerStyle(8);
  ptmax[0]->Sumw2();
  ptmax[0]->Draw("same pe1");

  leg->Draw();

  cptmax.SaveAs("ptmax.eps");
  cptmax.SaveAs("ptmax.root");



  // draw ptmin
  TCanvas cptmin;
  cptmin.SetLogy();
  ptmin[6]->SetMaximum(200);
  ptmin[6]->SetMinimum(0.01);
  ptmin[6]->SetFillColor(6);
  ptmin[6]->SetTitle("");
  ptmin[6]->GetXaxis()->SetTitle("P^{e min}_{T} [GeV]");
  ptmin[6]->GetYaxis()->SetTitle("normalized Events");
  ptmin[6]->Draw("hist");

  ptmin[1]->SetFillColor(5);
  ptmin[1]->Draw("same hist");

  ptmin[7]->SetFillColor(4);
  ptmin[7]->Draw("same hist");

  ptmin[5]->SetFillColor(2);
  ptmin[5]->Draw("same hist");

  ptmin[4]->SetFillColor(28);
  ptmin[4]->Draw("same hist");

  ptmin[2]->SetFillColor(7);
  ptmin[2]->Draw("same hist");

  ptmin[3]->SetFillColor(3);
  ptmin[3]->Draw("same hist");
  
  ptmin[0]->SetMarkerStyle(8);
  ptmin[0]->Sumw2();
  ptmin[0]->Draw("same pe1");

  leg->Draw();

  cptmin.SaveAs("ptmin.eps");
  cptmin.SaveAs("ptmin.root");



}
