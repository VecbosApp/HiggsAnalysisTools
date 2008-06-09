// macro to draw normalized distributions
// ----
// efficiencies are hardcoded for 2 selections:
//   1) full selection
//   2) after the CJV (trigger+reco+iso+ID+CJV)
// ----
// lumi to normalize is hardcoded (100 pb-1)
// ----
// usage: 
// root -b
// .L macro/higgsPlots.cxx++
// drawKinematics("jetVeto"): draws distributions after CJV
// drawKinematics("finalSelection"): draw distributions after the full selection (but stat is poor)
// //     this uses the tight egamma electron ID (default)
// drawKinematics("jetVeto","loose"): draws distributions after CJV with loose egamma electron ID
// ---

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

std::vector<float> nEventsPresel;
std::vector<float> nEventsCJV;
std::vector<float> nEventsFinal;

void setExpectedEvents() {

  // mH = 160 GeV
  nEventsPresel.push_back(14.6); // H->WW
  nEventsPresel.push_back(72.1); // WW
  nEventsPresel.push_back(112); // WZ
  nEventsPresel.push_back(19.1);  // ZZ
  nEventsPresel.push_back(41.1);  // tW
  nEventsPresel.push_back(2628); // W+j
  nEventsPresel.push_back(2463); // Z+j
  nEventsPresel.push_back(1193); // ttbar

  nEventsCJV.push_back(4.9); // H->WW
  nEventsCJV.push_back(18.2); // WW
  nEventsCJV.push_back(3.6); // WZ
  nEventsCJV.push_back(4.2);  // ZZ
  nEventsCJV.push_back(7.1);  // tW
  nEventsCJV.push_back(10.4); // W+j
  nEventsCJV.push_back(397); // Z+j
  nEventsCJV.push_back(11.9); // ttbar

  nEventsFinal.push_back(1.90); // H->WW
  nEventsFinal.push_back(1.39); // WW
  nEventsFinal.push_back(0.17); // WZ
  nEventsFinal.push_back(0.0);  // ZZ
  nEventsFinal.push_back(0.0);  // tW
  nEventsFinal.push_back(0.20); // W+j
  nEventsFinal.push_back(0.10); // Z+j
  nEventsFinal.push_back(0.43); // ttbar

}

void drawKinematics(const char* selection) {

  setExpectedEvents();
  
  // get the expected events for each process considered
  char Selection[500];
  std::vector<float> expEvents;
  if (strcmp(selection,"Preselection")==0) {
    sprintf(Selection, "1==1");
    expEvents = nEventsPresel;
  }
  else if (strcmp(selection,"jetVeto")==0) {
    sprintf(Selection,"jetVeto");
    expEvents = nEventsCJV;
  }
  else {
    sprintf(Selection,"finalSelection");
    expEvents = nEventsFinal;
  }

  std::vector<TFile*> datasets;

  TFile *Higgs = 0, *WW_incl = 0, *WZ = 0, *ZZ_incl = 0, *tW_incl = 0, *Chowder = 0;
  Higgs   = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/HiggsH160_CMSSW_1_6_9-datasetEE.root");
  WW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/WW_incl-datasetEE.root");
  WZ      = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/WZ-datasetEE.root");
  ZZ_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/ZZ_incl-datasetEE.root");
  tW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/tW_incl-datasetEE.root");
  Chowder = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/ChowderPDElectronSkim-datasetEE.root");

  datasets.push_back(Higgs);
  datasets.push_back(WW_incl);
  datasets.push_back(WZ);
  datasets.push_back(ZZ_incl);
  datasets.push_back(tW_incl);
  datasets.push_back(Chowder);

  std::vector<TH1F*> met;
  std::vector<TH1F*> mll;
  std::vector<TH1F*> deltaphi;
  std::vector<TH1F*> ptmax;
  std::vector<TH1F*> ptmin;

  TH1F* metH = new TH1F("metH","metH",25,0,200);
  TH1F* mllH = new TH1F("mllH", "mllH",25,0,200);
  TH1F* deltaphiH = new TH1F("deltaphiH", "deltaphiH",25,0,180);
  TH1F* ptmaxH = new TH1F("ptmaxH","ptmaxH",25,0,200);
  TH1F* ptminH = new TH1F("ptminH", "ptminH",25,0,200);

  std::cout << "# of samples = " <<  expEvents.size() << std::endl;

  for(int i=0; i<(int)expEvents.size(); i++) {

    std::cout << "Expected events for sample " << i 
	      << " at step " << selection 
	      << " in 100 pb-1 = " << expEvents[i] << std::endl;

    TTree *tree;
    if(i < 5 )
      tree = (TTree*)datasets[i]->Get("T1");
    else // ALPGEN (Chowder)
      tree = (TTree*)datasets[5]->Get("T1");

    std::cout << "dataset has " << tree->GetEntries() << " entries" << std::endl;

    char buf[50];
    sprintf(buf,"met_%d",i);
    TH1F* metProcessX = (TH1F*) metH->Clone(buf);
    metProcessX->Sumw2();
    met.push_back(metProcessX);

    if(i < 5 ) {
      tree->Project(buf,"met",Selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"met",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"met",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"met",extracut);
    }

    met[i]->Scale( expEvents[i]/met[i]->Integral() );

    sprintf(buf,"mll_%d",i);
    TH1F* mllProcessX = (TH1F*) mllH->Clone(buf);
    mllProcessX->Sumw2();
    mll.push_back(mllProcessX);

    if(i < 5 ) {
      tree->Project(buf,"eleInvMass",Selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"eleInvMass",extracut);
    }
    mll[i]->Scale( expEvents[i]/mll[i]->Integral() );

    sprintf(buf,"deltaphi_%d",i);
    TH1F* deltaphiProcessX = (TH1F*) deltaphiH->Clone(buf);
    deltaphiProcessX->Sumw2();
    deltaphi.push_back(deltaphiProcessX);

    if(i < 5 ) {
      tree->Project(buf,"deltaPhi",Selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"deltaPhi",extracut);
    }
    deltaphi[i]->Scale( expEvents[i]/deltaphi[i]->Integral() );


    sprintf(buf,"ptmax_%d",i);
    TH1F* ptmaxProcessX = (TH1F*) ptmaxH->Clone(buf);
    ptmaxProcessX->Sumw2();
    ptmax.push_back(ptmaxProcessX);

    if(i < 5 ) {
      tree->Project(buf,"maxPtEle",Selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"maxPtEle",extracut);
    }
    ptmax[i]->Scale( expEvents[i]/ptmax[i]->Integral() );

    sprintf(buf,"ptmin_%d",i);
    TH1F* ptminProcessX = (TH1F*) ptminH->Clone(buf);
    ptminProcessX->Sumw2();
    ptmin.push_back(ptminProcessX);

    if(i < 5 ) {
      tree->Project(buf,"minPtEle",Selection);
    }
    else if(i==5) { // W+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=1000 && CSA07processId<2000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    else if(i==6) { // Z+jets
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=2000 && CSA07processId<3000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    else if(i==7) { // ttbar
      char extracut[200];
      sprintf(extracut,"(CSA07processId>=3000 && CSA07processId<4000 && %s)*CSA07weight",Selection);
      tree->Project(buf,"minPtEle",extracut);
    }
    ptmin[i]->Scale( expEvents[i]/ptmin[i]->Integral() );

  }
  

  TLegend *leg = new TLegend(0.11,0.65,0.45,0.89);
  leg->SetBorderSize(0);
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
  TCanvas cmet("cmet","cmet",600,600);
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
  met[7]->SetFillStyle(3004);
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
  met[0]->Draw("same pe1");

  leg->Draw();

  cmet.SaveAs("met.eps");
  cmet.SaveAs("met.root");


  // draw mll
  TCanvas cmll("cmll","cmll",600,600);
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
  mll[7]->SetFillStyle(3004);
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
  mll[0]->Draw("same pe1");

  leg->Draw();

  cmll.SaveAs("mll.eps");
  cmll.SaveAs("mll.root");



  // draw deltaphi
  TCanvas cdeltaphi("cdeltaphi","cdeltaphi",600,600);
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
  deltaphi[7]->SetFillStyle(3004);
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
  deltaphi[0]->Draw("same pe1");

  leg->Draw();

  cdeltaphi.SaveAs("deltaphi.eps");
  cdeltaphi.SaveAs("deltaphi.root");



  // draw ptmax
  TCanvas cptmax("cptmax","cptmax",600,600);
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
  ptmax[7]->SetFillStyle(3004);
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
  ptmax[0]->Draw("same pe1");

  leg->Draw();

  cptmax.SaveAs("ptmax.eps");
  cptmax.SaveAs("ptmax.root");



  // draw ptmin
  TCanvas cptmin("cptmin","cptmin",600,600);
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
  ptmin[7]->SetFillStyle(3004);
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
  ptmin[0]->Draw("same pe1");

  leg->Draw();

  cptmin.SaveAs("ptmin.eps");
  cptmin.SaveAs("ptmin.root");

}
