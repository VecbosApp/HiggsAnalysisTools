#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include <iostream>

#define NSPECIES 5
#define NVARIABLES 8
#define NCUTS 3

void makeMCPlots(const char *finalstate, int signalFactor=1)
{
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x tdrstyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetOptTitle(0); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);

  TString suffix="";

  TString species[NSPECIES];
  species[0]="Data";
  if(signalFactor==1) species[1]="H120";
  else {
    char scaleF[10];
    sprintf(scaleF,"%dxH120",signalFactor);
    species[1]=TString(scaleF);
  }
  species[2]="others";
  species[3]="top";
  species[4]="WW";

  Color_t colors[NSPECIES];
  colors[0]=kBlack;
  colors[1]=kRed;       
  colors[2]=kAzure-1;
  colors[3]=kOrange;
  colors[4]=kViolet;

  Color_t lineColors[NSPECIES];
  lineColors[0]=kBlack;
  lineColors[1]=kRed;      
  lineColors[2]=kAzure+3;
  lineColors[3]=kOrange+3;
  lineColors[4]=kViolet+3;

  int legendOrder[NSPECIES];
  legendOrder[0]=0;
  legendOrder[1]=1;
  legendOrder[2]=4;
  legendOrder[3]=3;
  legendOrder[4]=2;

  // chiara, da sistemare
  TString files[NSPECIES];
  files[0]="results_data/merged/dataset_"+TString(finalstate)+".root";  
  files[1]="results/datasets_trees/H120_"+TString(finalstate)+".root";  
  files[2]="results/datasets_trees/others_"+TString(finalstate)+".root";
  files[3]="results/datasets_trees/top_"+TString(finalstate)+".root";
  files[4]="results/datasets_trees/WW_"+TString(finalstate)+".root";

  TString plotsDir="./HWW/"+TString(finalstate)+"/";

  TFile* fOut=new TFile("HWW_histos_"+suffix+".root","RECREATE");
  
  char icut[NCUTS][100];
  TH1F* histos[NSPECIES][NCUTS][NVARIABLES];    // 5 species, 2 cut levels, 8 variables
  
  TString variables[NVARIABLES];
  variables[0]="met";
  variables[1]="projMet";
  variables[2]="transvMass";
  variables[3]="eleInvMass";
  variables[4]="maxPtEle";
  variables[5]="minPtEle";
  variables[6]="deltaPhi";
  variables[7]="njets";

  TString units[NVARIABLES];
  units[0]="GeV";
  units[1]="GeV";
  units[2]="GeV/c^{2}";
  units[3]="GeV/c^{2}";
  units[4]="GeV/c";
  units[5]="GeV/c";
  units[6]="#circ";
  units[7]="";

  int nbins[NVARIABLES];
  nbins[0]=30;
  nbins[1]=30;
  nbins[2]=30;
  nbins[3]=30;
  nbins[4]=30;
  nbins[5]=30;
  nbins[6]=20;
  nbins[7]=7;

  float range[NVARIABLES][2]; // 8 variables, min, max
  // met
  range[0][0]=0.;
  range[0][1]=200.;
  // projected met
  range[1][0]=0.;
  range[1][1]=200.;
  // mt
  range[2][0]=0.;
  range[2][1]=250.;
  // mll
  range[3][0]=20.;
  range[3][1]=150.;
  // max pt
  range[4][0]=0.;
  range[4][1]=200.;
  // min pt
  range[5][0]=0.;
  range[5][1]=200.;
  // delta phi
  range[6][0]=0.;
  range[6][1]=180.;
  // njets
  range[7][0]=0.;
  range[7][1]=7.;

  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="MET";
  xaxisLabel[1]="projected MET";
  xaxisLabel[2]="m_{T}";
  xaxisLabel[3]="m_{ll}";
  xaxisLabel[4]="p_{T}^{max}";
  xaxisLabel[5]="p_{T}^{min}";
  xaxisLabel[6]="#Delta #phi";
  xaxisLabel[7]="n jets";

  TString binSize[NVARIABLES];

  for (int z=0;z<NVARIABLES;++z) {
    for (int j=0;j<NCUTS;++j) {
      sprintf(icut[j],"icut%d",j);
      for (int i=0;i<NSPECIES;++i) {
	histos[i][j][z]=new TH1F(variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),nbins[z],range[z][0],range[z][1]);
	if(i==0)
	  histos[i][j][z]->Sumw2();
	char binsiz[10];
	sprintf(binsiz,"%2.0f",(range[z][1]-range[z][0])/nbins[z]);
	binSize[z]=TString(binsiz);
      }
    }
  }

  TString cut[NCUTS];
  cut[0]="(finalLeptons*(eleInvMass>20))*";
  cut[1]="(jetVeto)*";
  cut[2]="(preDeltaPhi)*";

  TString intLumi="29";     
  TFile *_file[NSPECIES];
  TTree *T1[NSPECIES];
  
  TCanvas* c1= new TCanvas("test","test",800,800);
  
  for (int i=0;i<NSPECIES;++i) {
    _file[i]=TFile::Open(files[i]);
    T1[i] = (TTree*)_file[i]->Get("T1");
  }

  int nspeciesToRun=NSPECIES;

  for (int z=0;z<NVARIABLES;++z)
    {
      for (int j=0;j<NCUTS;++j)
	{
	  for (int i=0;i<nspeciesToRun;++i)
	    {
	      fOut->cd();
	      TString histoName=variables[z]+"_W_"+species[i]+"_"+TString(icut[j]);
	      std::cout << "Producing " << histoName << std::endl;
	      if (T1[i]==0)
		{
		  std::cout << "Tree not found" << std::endl;
		  return;
		}
              // T1[i]->Project(histoName,variables[z],cut[j]+"weight*"+intLumi);
              T1[i]->Project(histoName,variables[z],cut[j]+"weight");
	      std::cout << "Done " << histoName << std::endl;
	    }
	  
	  THStack histo_MC(variables[z]+"_MC",variables[z]+"_MC");
	  for (int i=2;i<nspeciesToRun;++i)
	    {
	      histos[i][j][z]->SetFillColor(colors[i]);
	      histos[i][j][z]->SetLineColor(lineColors[i]);
	      histo_MC.Add(histos[i][j][z]);
	    }
	  
	  float maximum=TMath::Max(histo_MC.GetMaximum(),histos[0][j][z]->GetMaximum());
	  histo_MC.SetMinimum(0.0001);
	  histo_MC.SetMaximum(maximum*2.);
	  histo_MC.Draw("");
	  histo_MC.GetXaxis()->SetTitle(xaxisLabel[z]+" ["+units[z]+"]");
	  histo_MC.GetYaxis()->SetTitle("Events/"+binSize[z]+" "+units[z]);

          histos[1][j][z]->Scale(signalFactor);
	  histos[1][j][z]->SetLineColor(lineColors[1]);
	  histos[1][j][z]->SetLineWidth(2.0);
	  histos[1][j][z]->Draw("SAME");
	  
	  histos[0][j][z]->SetMarkerStyle(20);
	  histos[0][j][z]->SetMarkerSize(1.3);
	  histos[0][j][z]->Draw("EPSAME");
	  //  c1->SetLogy();
	  TPaveText pt1(0.6,0.83,0.8,0.9,"NDC");
	  //	  pt1.SetTextFont(72);
	  pt1.SetTextSize(0.028);
	  pt1.SetTextAlign(12);
	  pt1.SetFillColor(0);
	  pt1.SetBorderSize(0);
	  // pt1.AddText("CMS Preliminary 2010");
	  // pt1.AddText("");
	  // pt1.AddText("");
	  // pt1.AddText("");
	  pt1.AddText("#sqrt{s}=7 TeV L_{int}="+intLumi+" pb^{-1}");
	  pt1.Draw();

	  c1->Update();
	  TLegendEntry *legge;
	  TLegend *leg;
	  leg = new TLegend(0.6,0.65,0.93,0.8,cut[j]);
	  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.025);
	  leg->SetFillColor(0);
	  //      leg->SetTitle(cut[j]);
	  for (int i=0;i<nspeciesToRun;++i)
	    {
	      if (i == 0)
		legge = leg->AddEntry(histos[legendOrder[i]][j][z],species[legendOrder[i]], "lpe");
	  else
	    legge = leg->AddEntry(histos[legendOrder[i]][j][z],species[legendOrder[i]],"f");
	    }
	  leg->Draw();
	  c1->SetLogy(0);
	  c1->Update();
	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".png");
	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".root");
	  histo_MC.SetMaximum(maximum*100);
	  c1->SetLogy(1);
	  c1->Update();
	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+"_log.png");
	}
    }
  
  fOut->Write();
  fOut->Close();
  
}
