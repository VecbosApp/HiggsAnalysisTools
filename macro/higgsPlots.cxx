// macro to draw normalized distributions
// ----
// efficiencies are hardcoded for 2 selections:
//   1) full selection
//   2) after the CJV (trigger+reco+iso+ID+CJV)
// usage: 
// root -b
// .L macro/higgsPlots.cxx++
// drawKinematics("jetVeto"): draws distributions after CJV
// drawKinematics("finalSelection"): draw distributions after the full selection (but stat is poor)
// ---

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

std::vector<float> nEventsPresel;
std::vector<float> nEventsFinalLeptons;
std::vector<float> nEventsCJV;
std::vector<float> nEventsPreDeltaPhi;
std::vector<float> nEventsFinal;

void setExpectedEvents() {

  // events expected in 100 pb-1
  // mH = 160 GeV
  nEventsPresel.push_back(14.6); // H->WW
  nEventsPresel.push_back(72.1); // WW
  nEventsPresel.push_back(112); // WZ
  nEventsPresel.push_back(19.1);  // ZZ
  nEventsPresel.push_back(41.1);  // tW
  nEventsPresel.push_back(2628); // W+j
  nEventsPresel.push_back(2463); // Z+j (> 40 GeV)
  nEventsPresel.push_back(1193); // ttbar
  nEventsPresel.push_back(31.2); // Drell Yan (10-40 GeV)
  
  nEventsFinalLeptons.push_back(8.2); // H->WW
  nEventsFinalLeptons.push_back(26.4); // WW
  nEventsFinalLeptons.push_back(28.6); // WZ
  nEventsFinalLeptons.push_back(12.5);  // ZZ
  nEventsFinalLeptons.push_back(22.0);  // tW
  nEventsFinalLeptons.push_back(15.5); // W+j
  nEventsFinalLeptons.push_back(1538); // Z+j (> 40 GeV)
  nEventsFinalLeptons.push_back(253); // ttbar
  nEventsFinalLeptons.push_back(15.7); // Drell Yan (10-40 GeV)

  nEventsCJV.push_back(4.9); // H->WW
  nEventsCJV.push_back(18.2); // WW
  nEventsCJV.push_back(3.6); // WZ
  nEventsCJV.push_back(4.2);  // ZZ
  nEventsCJV.push_back(7.1);  // tW
  nEventsCJV.push_back(10.4); // W+j
  nEventsCJV.push_back(397); // Z+j (> 40 GeV)
  nEventsCJV.push_back(11.9); // ttbar
  nEventsCJV.push_back(7.0); // Drell Yan (10-40 GeV)

  nEventsPreDeltaPhi.push_back(2.3); // H->WW
  nEventsPreDeltaPhi.push_back(2.2); // WW
  nEventsPreDeltaPhi.push_back(0.2); // WZ
  nEventsPreDeltaPhi.push_back(0.0);  // ZZ
  nEventsPreDeltaPhi.push_back(0.03);  // tW
  nEventsPreDeltaPhi.push_back(0.62); // W+j
  nEventsPreDeltaPhi.push_back(0.11); // Z+j (> 40 GeV)
  nEventsPreDeltaPhi.push_back(0.56); // ttbar
  nEventsPreDeltaPhi.push_back(0.0); // Drell Yan (10-40 GeV)

  nEventsFinal.push_back(1.90); // H->WW
  nEventsFinal.push_back(1.39); // WW
  nEventsFinal.push_back(0.17); // WZ
  nEventsFinal.push_back(0.0);  // ZZ
  nEventsFinal.push_back(0.0);  // tW
  nEventsFinal.push_back(0.20); // W+j
  nEventsFinal.push_back(0.10); // Z+j (> 40 GeV)
  nEventsFinal.push_back(0.43); // ttbar
  nEventsFinal.push_back(0.0); // Drell Yan (10-40 GeV)

}

//! lumi in pb-1
void drawKinematics(const char* selection, float lumi=100) {

  setExpectedEvents();
  
  // get the expected events for each process considered
  char Selection[500];
  std::vector<float> expEvents;
  if (strcmp(selection,"Preselection")==0) {
    sprintf(Selection, "1==1");
    expEvents = nEventsPresel;
  }
  else if (strcmp(selection,"finalLeptons")==0) {
    sprintf(Selection, "finalLeptons");
    expEvents = nEventsFinalLeptons;
  }
  else if (strcmp(selection,"jetVeto")==0) {
    sprintf(Selection,"jetVeto");
    expEvents = nEventsCJV;
  }
  else if (strcmp(selection,"preDeltaPhi")==0) {
    sprintf(Selection,"preDeltaPhi");
    expEvents = nEventsPreDeltaPhi;
  }
  else {
    sprintf(Selection,"finalSelection");
    expEvents = nEventsFinal;
  }

  std::vector<TFile*> datasets;

  TFile *Higgs = 0, *WW_incl = 0, *WZ = 0, *ZZ_incl = 0, *tW_incl = 0, *Chowder = 0, *DY10to40 = 0;
  Higgs   = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/HiggsH160_CMSSW_1_6_9-datasetEE.root");
  WW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/WW_incl-datasetEE.root");
  WZ      = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/WZ-datasetEE.root");
  ZZ_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/ZZ_incl-datasetEE.root");
  tW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/tW_incl-datasetEE.root");
  Chowder = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/ChowderPDElectronSkim-datasetEE.root");
  DY10to40 =  TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/datasets/ForNoteSummer08_new2/H160/DrellYan_ll_10-40-datasetEE.root");

  // with deltaPhi cut at the end and the variable "beforeDeltaPhi -> all cuts with exception of dPhi"
//   Higgs   = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/HiggsH160_CMSSW_1_6_9-datasetEE.root");
//   WW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/WW_incl-datasetEE.root");
//   WZ      = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/WZ-datasetEE.root");
//   ZZ_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/ZZ_incl-datasetEE.root");
//   tW_incl = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/tW_incl-datasetEE.root");
//   Chowder = TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/ChowderPDElectronSkim-datasetEE.root");
//   DY10to40 =  TFile::Open("/cmsrm/pc18/emanuele/releases/HIGGS_RELEASES/OfflineAnalysis/HiggsAnalysisTools/tmplogs/H160/DrellYan_ll_10-40-datasetEE.root");

  datasets.push_back(Higgs); // 0
  datasets.push_back(WW_incl); // 1
  datasets.push_back(WZ); // 2 
  datasets.push_back(ZZ_incl); // 3 
  datasets.push_back(tW_incl); // 4 
  datasets.push_back(Chowder); // 5 
  datasets.push_back(DY10to40); // 6

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
    else if(i>=5 && i<=7) // ALPGEN (Chowder)
      tree = (TTree*)datasets[5]->Get("T1");
    else
      tree = (TTree*)datasets[i-2]->Get("T1");

    std::cout << "dataset has " << tree->GetEntries() << " entries" << std::endl;

    char buf[50];
    sprintf(buf,"met_%d",i);
    TH1F* metProcessX = (TH1F*) metH->Clone(buf);
    metProcessX->Sumw2();
    met.push_back(metProcessX);

    if(i < 5 || i > 7) {
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

    met[i]->Scale( expEvents[i]/met[i]->Integral() * lumi/100. );

    sprintf(buf,"mll_%d",i);
    TH1F* mllProcessX = (TH1F*) mllH->Clone(buf);
    mllProcessX->Sumw2();
    mll.push_back(mllProcessX);

    if( i < 5 || i > 7 ) {
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
    mll[i]->Scale( expEvents[i]/mll[i]->Integral() * lumi/100. );

    sprintf(buf,"deltaphi_%d",i);
    TH1F* deltaphiProcessX = (TH1F*) deltaphiH->Clone(buf);
    deltaphiProcessX->Sumw2();
    deltaphi.push_back(deltaphiProcessX);

    if( i < 5 || i > 7 ) {
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
    deltaphi[i]->Scale( expEvents[i]/deltaphi[i]->Integral() * lumi/100. );


    sprintf(buf,"ptmax_%d",i);
    TH1F* ptmaxProcessX = (TH1F*) ptmaxH->Clone(buf);
    ptmaxProcessX->Sumw2();
    ptmax.push_back(ptmaxProcessX);

    if( i < 5 || i > 7 ) {
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
    ptmax[i]->Scale( expEvents[i]/ptmax[i]->Integral() * lumi/100. );

    sprintf(buf,"ptmin_%d",i);
    TH1F* ptminProcessX = (TH1F*) ptminH->Clone(buf);
    ptminProcessX->Sumw2();
    ptmin.push_back(ptminProcessX);

    if( i < 5 || i > 7 ) {
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
    ptmin[i]->Scale( expEvents[i]/ptmin[i]->Integral() * lumi/100. );

  }
  

  TLegend *leg = new TLegend(0.11,0.65,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(met[0],"Signal, m_{H}=160 GeV","pl");
  leg->AddEntry(met[1],"WW","f");
  leg->AddEntry(met[7],"tt","f");
  leg->AddEntry(met[6],"Z+jets, m(ee)>40 GeV","f");
  leg->AddEntry(met[8], "Drell Yan, 10<m(ee)<40 GeV","f");
  leg->AddEntry(met[3],"ZZ","f");
  leg->AddEntry(met[2],"WZ","f");
  leg->AddEntry(met[5],"W+jets","f");
  leg->AddEntry(met[4],"tW","f");

  TH1F *histo0, *histo1, *histo2, *histo3, *histo4;


  gStyle->SetOptStat(0);

  // draw met
  TCanvas cmet("cmet","cmet",600,600);
  cmet.SetLogy();
  met[6]->SetMaximum(50000);
  met[6]->SetMinimum(0.01);
  met[6]->SetFillColor(6);
  met[6]->SetTitle("");
  met[6]->GetXaxis()->SetTitle("Missing E_{T} [GeV]");
  met[6]->GetYaxis()->SetTitle("normalized Events");
  met[6]->Draw("hist");

  met[8]->SetFillColor(kPink+4);
  met[8]->Draw("same hist");

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
  met[0]->Draw("same pe1");

  leg->Draw();

  TFile *fmet = new TFile("met.root","recreate");
  cmet.Write();
  histo0 = met[0];
  met[6]->Add(met[8]);
  histo1 = met[6];
  histo2 = met[7];
  histo3 = met[1];
  met[2]->Add(met[3]);
  met[2]->Add(met[4]);
  met[2]->Add(met[5]);
  histo4 = met[2];

  TCanvas cmet2("cmet2","cmet2",600,600);
  cmet2.SetLogy();

  histo0->SetName("histo0"); // H->WW
  histo1->SetName("histo1"); // Z+jets+DY_10_40
  histo2->SetName("histo2"); // ttbar
  histo3->SetName("histo3"); // WW
  histo4->SetName("histo4"); // others

  histo0->SetMarkerStyle(8);
  histo0->SetMarkerSize(1.5);
  histo1->SetFillColor(kRed+2);
  histo2->SetFillColor(kYellow-7);
  histo3->SetFillColor(kRed-4);
  histo4->SetFillColor(kBlue+3);
  histo4->SetFillStyle(3003);

  histo0->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();

  TLegend *leg2 = new TLegend(0.11,0.65,0.45,0.89);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(histo0,"Signal, m_{H}=160 GeV","pl");
  leg2->AddEntry(histo1,"Z-jets","f");
  leg2->AddEntry(histo2,"t #bar{t}-jets","f");
  leg2->AddEntry(histo3,"WW","f");
  leg2->AddEntry(histo4,"W-jets + ZZ + WZ + tW","f");

  histo1->Draw("hist");
  histo3->Draw("same hist");
  histo2->Draw("same hist");
  histo4->Draw("same hist");
  histo0->Draw("same pe1");
  leg2->Draw();
  cmet2.Write();
  fmet->Close();


  // draw mll
  TCanvas cmll("cmll","cmll",600,600);
  cmll.SetLogy();
  mll[6]->SetMaximum(50000);
  mll[6]->SetMinimum(0.01);
  mll[6]->SetFillColor(6);
  mll[6]->SetTitle("");
  mll[6]->GetXaxis()->SetTitle("m_{ll} [GeV]");
  mll[6]->GetYaxis()->SetTitle("normalized Events");
  mll[6]->Draw("hist");

  mll[8]->SetFillColor(kPink+4);
  mll[8]->Draw("same hist");

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
  mll[0]->Draw("same pe1");

  leg->Draw();

  TFile *fmll = new TFile("mll.root","recreate");
  cmll.Write();
  histo0 = mll[0];
  mll[6]->Add(mll[8]);
  histo1 = mll[6];
  histo2 = mll[7];
  histo3 = mll[1];
  mll[2]->Add(mll[3]);
  mll[2]->Add(mll[4]);
  mll[2]->Add(mll[5]);
  histo4 = mll[2];
  histo0->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();

  TCanvas cmll2("cmll2","cmll2",600,600);

  histo0->SetName("histo0"); // H->WW
  histo1->SetName("histo1"); // Z+jets+DY_10_40
  histo2->SetName("histo2"); // ttbar
  histo3->SetName("histo3"); // WW
  histo4->SetName("histo4"); // others

  histo0->SetMarkerStyle(8);
  histo0->SetMarkerSize(1.5);
  histo1->SetFillColor(kRed+2);
  histo2->SetFillColor(kYellow-7);
  histo3->SetFillColor(kRed-4);
  histo4->SetFillColor(kBlue+3);
  histo4->SetFillStyle(3003);

  cmll2.SetLogy();
  histo1->Draw("hist");
  histo3->Draw("same hist");
  histo2->Draw("same hist");
  histo4->Draw("same hist");
  histo0->Draw("same pe1");
  leg2->Draw();
  cmll2.Write();
  fmll->Close();



  // draw deltaphi
  TCanvas cdeltaphi("cdeltaphi","cdeltaphi",600,600);
  cdeltaphi.SetLogy();
  deltaphi[6]->SetMaximum(100000);
  deltaphi[6]->SetMinimum(0.1);
  deltaphi[6]->SetFillColor(6);
  deltaphi[6]->SetTitle("");
  deltaphi[6]->GetXaxis()->SetTitle("#Delta #phi");
  deltaphi[6]->GetYaxis()->SetTitle("normalized Events");
  deltaphi[6]->Draw("hist");

  deltaphi[8]->SetFillColor(kPink+4);
  deltaphi[8]->Draw("same hist");

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
  deltaphi[0]->Draw("same pe1");

  leg->Draw();

  TFile *fdeltaphi = new TFile("deltaphi.root","recreate");
  cdeltaphi.Write();
  histo0 = deltaphi[0];
  deltaphi[6]->Add(deltaphi[8]);
  histo1 = deltaphi[6];
  histo2 = deltaphi[7];
  histo3 = deltaphi[1];
  deltaphi[2]->Add(deltaphi[3]);
  deltaphi[2]->Add(deltaphi[4]);
  deltaphi[2]->Add(deltaphi[5]);
  histo4 = deltaphi[2];

  histo0->SetName("histo0"); // H->WW
  histo1->SetName("histo1"); // Z+jets+DY_10_40
  histo2->SetName("histo2"); // ttbar
  histo3->SetName("histo3"); // WW
  histo4->SetName("histo4"); // others

  histo0->SetMarkerStyle(8);
  histo0->SetMarkerSize(1.5);
  histo1->SetFillColor(kRed+2);
  histo2->SetFillColor(kYellow-7);
  histo3->SetFillColor(kRed-4);
  histo4->SetFillColor(kBlue+3);
  histo4->SetFillStyle(3003);

  histo0->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();

  TCanvas cdeltaphi2("cdeltaphi2","cdeltaphi2",600,600);
  cdeltaphi2.SetLogy();
  histo1->Draw("hist");
  histo3->Draw("same hist");
  histo2->Draw("same hist");
  histo4->Draw("same hist");
  histo0->Draw("same pe1");
  leg2->Draw();
  cdeltaphi2.Write();
  fdeltaphi->Close();



  // draw ptmax
  TCanvas cptmax("cptmax","cptmax",600,600);
  cptmax.SetLogy();
  ptmax[6]->SetMaximum(100000);
  ptmax[6]->SetMinimum(0.01);
  ptmax[6]->SetFillColor(6);
  ptmax[6]->SetTitle("");
  ptmax[6]->GetXaxis()->SetTitle("P^{e max}_{T} [GeV]");
  ptmax[6]->GetYaxis()->SetTitle("normalized Events");
  ptmax[6]->Draw("hist");

  ptmax[8]->SetFillColor(kPink+4);
  ptmax[8]->Draw("same hist");

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
  ptmax[0]->Draw("same pe1");

  leg->Draw();

  TFile *fptmax = new TFile("ptmax.root","recreate");
  cptmax.Write();
  histo0 = ptmax[0];
  ptmax[6]->Add(ptmax[8]);
  histo1 = ptmax[6];
  histo2 = ptmax[7];
  histo3 = ptmax[1];
  ptmax[2]->Add(ptmax[3]);
  ptmax[2]->Add(ptmax[4]);
  ptmax[2]->Add(ptmax[5]);
  histo4 = ptmax[2];

  histo0->SetName("histo0"); // H->WW
  histo1->SetName("histo1"); // Z+jets+DY_10_40
  histo2->SetName("histo2"); // ttbar
  histo3->SetName("histo3"); // WW
  histo4->SetName("histo4"); // others

  histo0->SetMarkerStyle(8);
  histo0->SetMarkerSize(1.5);
  histo1->SetFillColor(kRed+2);
  histo2->SetFillColor(kYellow-7);
  histo3->SetFillColor(kRed-4);
  histo4->SetFillColor(kBlue+3);
  histo4->SetFillStyle(3003);

  histo0->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();

  TCanvas cptmax2("cptmax2","cptmax2",600,600);
  cptmax2.SetLogy();
  histo1->Draw("hist");
  histo3->Draw("same hist");
  histo2->Draw("same hist");
  histo4->Draw("same hist");
  histo0->Draw("same pe1");
  leg2->Draw();
  cptmax2.Write();
  fptmax->Close();



  // draw ptmin
  TCanvas cptmin("cptmin","cptmin",600,600);
  cptmin.SetLogy();
  ptmin[6]->SetMaximum(100000);
  ptmin[6]->SetMinimum(0.01);
  ptmin[6]->SetFillColor(6);
  ptmin[6]->SetTitle("");
  ptmin[6]->GetXaxis()->SetTitle("P^{e min}_{T} [GeV]");
  ptmin[6]->GetYaxis()->SetTitle("normalized Events");
  ptmin[6]->Draw("hist");

  ptmin[8]->SetFillColor(kPink+4);
  ptmin[8]->Draw("same hist");

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
  ptmin[0]->Draw("same pe1");

  leg->Draw();

  TFile *fptmin = new TFile("ptmin.root","recreate");
  cptmin.Write();
  histo0 = ptmin[0];
  ptmin[6]->Add(ptmin[8]);
  histo1 = ptmin[6];
  histo2 = ptmin[7];
  histo3 = ptmin[1];
  ptmin[2]->Add(ptmin[3]);
  ptmin[2]->Add(ptmin[4]);
  ptmin[2]->Add(ptmin[5]);
  histo4 = ptmin[2];

  histo0->SetName("histo0"); // H->WW
  histo1->SetName("histo1"); // Z+jets+DY_10_40
  histo2->SetName("histo2"); // ttbar
  histo3->SetName("histo3"); // WW
  histo4->SetName("histo4"); // others

  histo0->SetMarkerStyle(8);
  histo0->SetMarkerSize(1.5);
  histo1->SetFillColor(kRed+2);
  histo2->SetFillColor(kYellow-7);
  histo3->SetFillColor(kRed-4);
  histo4->SetFillColor(kBlue+3);
  histo4->SetFillStyle(3003);

  histo0->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();

  TCanvas cptmin2("cptmin2","cptmin2",600,600);
  cptmin2.SetLogy();
  histo1->Draw("hist");
  histo3->Draw("same hist");
  histo2->Draw("same hist");
  histo4->Draw("same hist");
  histo0->Draw("same pe1");
  leg2->Draw();
  cptmin2.Write();
  fptmin->Close();

}


void drawSignificances() {
  
  // events/fb-1
  std::vector< std::vector<float> > events;
  std::vector<float> higgsMass;

  // mH=120
  std::vector<float> eventsH120;
  eventsH120.resize(9);
  eventsH120[0] = 3.9;   // H->WW
  eventsH120[1] = 33.6;  // WW
  eventsH120[2] = 6.9;   // ttbar
  eventsH120[3] = 24.6;  // W+jets
  eventsH120[4] = 10.3;  // Z+jets
  eventsH120[5] = 9.7;   // DY < 40 GeV
  eventsH120[6] = 2.5;   // WZ
  eventsH120[7] = 1.0;   // tW
  eventsH120[8] = 0.5;   // ZZ
  events.push_back(eventsH120);
  higgsMass.push_back(120);

  // mH=130
  std::vector<float> eventsH130;
  eventsH130.resize(9);
  eventsH130[0] = 6.2;   // H->WW
  eventsH130[1] = 34.3;  // WW
  eventsH130[2] = 6.9;   // ttbar
  eventsH130[3] = 24.6;  // W+jets
  eventsH130[4] = 10.3;  // Z+jets
  eventsH130[5] = 9.7;   // DY < 40 GeV
  eventsH130[6] = 2.5;   // WZ
  eventsH130[7] = 1.0;   // tW
  eventsH130[8] = 0.7;   // ZZ
  events.push_back(eventsH130);
  higgsMass.push_back(130);

  // mH=140
  std::vector<float> eventsH140;
  eventsH140.resize(9);
  eventsH140[0] = 10.4;   // H->WW
  eventsH140[1] = 22.4;  // WW
  eventsH140[2] = 6.4;   // ttbar
  eventsH140[3] = 13.4;  // W+jets
  eventsH140[4] = 2.0;  // Z+jets
  eventsH140[5] = 0.0;   // DY < 40 GeV
  eventsH140[6] = 2.2;   // WZ
  eventsH140[7] = 0.0;   // tW
  eventsH140[8] = 0.0;   // ZZ
  events.push_back(eventsH140);
  higgsMass.push_back(140);

  // mH=150
  std::vector<float> eventsH150;
  eventsH150.resize(9);
  eventsH150[0] = 10.0;   // H->WW
  eventsH150[1] = 13.7;  // WW
  eventsH150[2] = 4.3;   // ttbar
  eventsH150[3] = 1.0;  // W+jets
  eventsH150[4] = 1.0;  // Z+jets
  eventsH150[5] = 0.0;   // DY < 40 GeV
  eventsH150[6] = 1.7;   // WZ
  eventsH150[7] = 0.0;   // tW
  eventsH150[8] = 0.0;   // ZZ
  events.push_back(eventsH150);
  higgsMass.push_back(150);
    
  // mH=160
  std::vector<float> eventsH160;
  eventsH160.resize(9);
  eventsH160[0] = 19.0;   // H->WW
  eventsH160[1] = 13.9;  // WW
  eventsH160[2] = 4.3;   // ttbar
  eventsH160[3] = 1.0;  // W+jets
  eventsH160[4] = 1.0;  // Z+jets
  eventsH160[5] = 0.0;   // DY < 40 GeV
  eventsH160[6] = 1.7;   // WZ
  eventsH160[7] = 0.0;   // tW
  eventsH160[8] = 0.0;   // ZZ
  events.push_back(eventsH160);
  higgsMass.push_back(160);

  // mH=165
  std::vector<float> eventsH165;
  eventsH165.resize(9);
  eventsH165[0] = 21.4;   // H->WW
  eventsH165[1] = 13.0;  // WW
  eventsH165[2] = 6.0;   // ttbar
  eventsH165[3] = 1.0;  // W+jets
  eventsH165[4] = 1.8;  // Z+jets
  eventsH165[5] = 0.0;   // DY < 40 GeV
  eventsH165[6] = 1.7;   // WZ
  eventsH165[7] = 0.0;   // tW
  eventsH165[8] = 0.2;   // ZZ
  events.push_back(eventsH165);
  higgsMass.push_back(165);

  // mH=170
  std::vector<float> eventsH170;
  eventsH170.resize(9);
  eventsH170[0] = 19.8;   // H->WW
  eventsH170[1] = 14.6;  // WW
  eventsH170[2] = 7.8;   // ttbar
  eventsH170[3] = 0.0;  // W+jets
  eventsH170[4] = 1.8;  // Z+jets
  eventsH170[5] = 0.0;   // DY < 40 GeV
  eventsH170[6] = 2.9;   // WZ
  eventsH170[7] = 0.0;   // tW
  eventsH170[8] = 0.3;   // ZZ
  events.push_back(eventsH170);
  higgsMass.push_back(170);

  // mH=180
  std::vector<float> eventsH180;
  eventsH180.resize(9);
  eventsH180[0] = 14.5;   // H->WW
  eventsH180[1] = 15.0;  // WW
  eventsH180[2] = 8.2;   // ttbar
  eventsH180[3] = 1.0;  // W+jets
  eventsH180[4] = 2.8;  // Z+jets
  eventsH180[5] = 0.0;   // DY < 40 GeV
  eventsH180[6] = 3.2;   // WZ
  eventsH180[7] = 0.5;   // tW
  eventsH180[8] = 0.3;   // ZZ
  events.push_back(eventsH180);
  higgsMass.push_back(180);

  // mH=190
  std::vector<float> eventsH190;
  eventsH190.resize(9);
  eventsH190[0] = 10.1;   // H->WW
  eventsH190[1] = 15.2;  // WW
  eventsH190[2] = 8.6;   // ttbar
  eventsH190[3] = 1.0;  // W+jets
  eventsH190[4] = 2.8;  // Z+jets
  eventsH190[5] = 0.0;   // DY < 40 GeV
  eventsH190[6] = 3.2;   // WZ
  eventsH190[7] = 0.5;   // tW
  eventsH190[8] = 0.3;   // ZZ
  events.push_back(eventsH190);
  higgsMass.push_back(190);

  // mH=200
  std::vector<float> eventsH200;
  eventsH200.resize(9);
  eventsH200[0] = 8.4;   // H->WW
  eventsH200[1] = 22.2;  // WW
  eventsH200[2] = 9.9;   // ttbar
  eventsH200[3] = 3.1;  // W+jets
  eventsH200[4] = 5.0;  // Z+jets
  eventsH200[5] = 0.0;   // DY < 40 GeV
  eventsH200[6] = 4.2;   // WZ
  eventsH200[7] = 0.8;   // tW
  eventsH200[8] = 0.3;   // ZZ
  events.push_back(eventsH200);
  higgsMass.push_back(200);


  TGraph *gStatSig = new TGraph( (int)events.size() );
  TGraph *gSoverB = new TGraph( (int)events.size() );

  for(unsigned int imass=0; imass<events.size(); imass++) {
    
    std::vector<float> nEvents = events[imass];

    float nSig = nEvents[0];
    float nBkg = 0.0;

    for(unsigned int iBkg=1; iBkg<nEvents.size(); iBkg++) {
      nBkg += nEvents[iBkg];
    }

    gStatSig->SetPoint(imass, higgsMass[imass], nSig/sqrt(nBkg));
    gSoverB->SetPoint(imass, higgsMass[imass], nSig/nBkg);

    std::cout << "Higgs mass = " << higgsMass[imass] 
	      << " has S/sqrt(B) = " <<  nSig/sqrt(nBkg)
	      << " and S/B = " << nSig/nBkg << std::endl;

  }

  TCanvas csignificance;
  gStatSig->SetLineColor(34);
  gStatSig->SetLineWidth(2);
  gStatSig->SetMarkerColor(34);
  gStatSig->SetMarkerStyle(8);
  gStatSig->SetMarkerSize(2);
  gStatSig->SetMinimum(0);
  gStatSig->SetTitle("");
  gStatSig->GetXaxis()->SetTitle("Higgs mass (GeV/c^{2})");
  gStatSig->GetYaxis()->SetTitle("n_{S} / #sqrt{n_{B}}");
  gStatSig->Draw("ACP");
  csignificance.SaveAs("statSign.eps");

  TCanvas csoverb;
  gSoverB->SetLineColor(34);
  gSoverB->SetLineWidth(2);
  gSoverB->SetMarkerColor(34);
  gSoverB->SetMarkerStyle(8);
  gSoverB->SetMarkerSize(2);
  gSoverB->SetMinimum(0);
  gSoverB->SetTitle("");
  gSoverB->GetXaxis()->SetTitle("Higgs mass (GeV/c^{2})");
  gSoverB->GetYaxis()->SetTitle("n_{S} / n_{B}");
  gSoverB->Draw("ACP");
  csoverb.SaveAs("sigOverBkg.eps");

}
