// cvs co -d scripts UserCode/Mangano/WWAnalysis/Misc/scripts

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
#include "scripts/LatinoPlot.C"

#include <iostream>

#define NSPECIES 7
#define NVARIABLES 9
#define NCUTS 5
#define JETBINS 2

void makeDataMCPlots(const char *finalstate, float lumi, bool blindData=false, int signalFactor=1)
{
  gROOT->SetStyle("Plain");
  //  gROOT->ProcessLine(".x tdrstyle.C");
  gROOT->ProcessLine(".x scripts/LatinoStyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  //  gStyle->SetOptTitle(0); 
  
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);

  TString suffix="";

  TString species[NSPECIES];
  species[0]="Data";
  if(signalFactor==1) species[1]="H130";
  else {
    char scaleF[10];
    sprintf(scaleF,"%dxH130",signalFactor);
    species[1]=TString(scaleF);
  }
  species[2]="Wjets";
  species[3]="diboson";
  species[4]="top";
  species[5]="Zjets";
  species[6]="WW";

  TString scalefactor_datadriven[NSPECIES][JETBINS];
  scalefactor_datadriven[0][0] = "1";
  scalefactor_datadriven[1][0] = "1";
  scalefactor_datadriven[2][0] = "2.69"; // here the ee+me fake rate tree is used only for the shape. SF = (WjetsTot/ (Wjets_ee + Wjets_me))
  //scalefactor_datadriven[2][0] = "0.45"; // for ee+mm only
  scalefactor_datadriven[3][0] = "1.0"; // taken from MC
  scalefactor_datadriven[4][0] = "1.83";
  scalefactor_datadriven[5][0] = "2.0";
  scalefactor_datadriven[6][0] = "1.05";

  scalefactor_datadriven[0][1] = "1";
  scalefactor_datadriven[1][1] = "1";
  scalefactor_datadriven[2][1] = "2.69"; // here the ee+me fake rate tree is used only for the shape. SF = (WjetsTot/ (Wjets_ee + Wjets_me))
  // scalefactor_datadriven[2][1] = "0.45"; // for ee+mm only 
  scalefactor_datadriven[3][1] = "1.0"; // taken from MC
  scalefactor_datadriven[4][1] = "1.17";
  scalefactor_datadriven[5][1] = "1.0";
  scalefactor_datadriven[6][1] = "1.2";

  Color_t colors[NSPECIES];
  colors[0]=kBlack;
  colors[1]=kBlack;       
  colors[2]=kAzure-1;
  colors[3]=kOrange;
  colors[4]=kViolet;
  colors[5]=kGreen+3;
  colors[6]=kOrange+10;

  Color_t lineColors[NSPECIES];
  lineColors[0]=kBlack;
  lineColors[1]=kBlack;      
  lineColors[2]=kAzure+3;
  lineColors[3]=kOrange+3;
  lineColors[4]=kViolet+3;
  lineColors[5]=kGreen+4;
  lineColors[6]=kOrange+9;

  int legendOrder[NSPECIES];
  legendOrder[0]=0;
  legendOrder[1]=1;
  legendOrder[2]=2;
  legendOrder[3]=3;
  legendOrder[4]=4;
  legendOrder[5]=5;
  legendOrder[6]=6;

  // chiara, da sistemare
  TString files[NSPECIES];
  files[0]="results_data/datasets_trees/dataset_"+TString(finalstate)+".root";  
  files[1]="results/datasets_trees/H130_"+TString(finalstate)+".root";  
  //  files[2]="results/datasets_trees/Wjets_"+TString(finalstate)+".root";
  files[2]="results_data/datasets_trees/dataset_fake_"+TString(finalstate)+".root";
  files[3]="results/datasets_trees/others_"+TString(finalstate)+".root";
  files[4]="results/datasets_trees/top_"+TString(finalstate)+".root";
  files[5]="results/datasets_trees/Zjets_"+TString(finalstate)+".root";
  files[6]="results/datasets_trees/WW_"+TString(finalstate)+".root";

  TString plotsDir="./HWW/"+TString(finalstate)+"/";

  TFile* fOut=new TFile("HWW_histos_"+suffix+".root","RECREATE");
  
  char icut[NCUTS][100];
  TH1F* histos[NSPECIES][NCUTS][NVARIABLES];    // 5 species, 2 cut levels, 8 variables
  
  TString variables[NVARIABLES];
  variables[0]="pfMet";
  variables[1]="projMet";
  variables[2]="transvMass";
  variables[3]="eleInvMass";
  variables[4]="maxPtEle";
  variables[5]="minPtEle";
  variables[6]="deltaPhi";
  variables[7]="njets";
  variables[8]="nVtx";

  TString units[NVARIABLES];
  units[0]="GeV";
  units[1]="GeV";
  units[2]="GeV/c^{2}";
  units[3]="GeV/c^{2}";
  units[4]="GeV/c";
  units[5]="GeV/c";
  units[6]="#circ";
  units[7]="";
  units[8]="";

  int nbins[NVARIABLES];
  nbins[0]=50;
  nbins[1]=50;
  nbins[2]=50;
  nbins[3]=80;
  nbins[4]=60;
  nbins[5]=40;
  nbins[6]=50;
  nbins[7]=7;
  nbins[8]=20;


  float range[NVARIABLES][2]; // 8 variables, min, max
  // met
  range[0][0]=0.;
  range[0][1]=150.;
  // projected met
  range[1][0]=0.;
  range[1][1]=150.;
  // gamma*MR*
  range[2][0]=0.;
  range[2][1]=250.;
  // mll
  range[3][0]=0.;
  range[3][1]=200.;
  // max pt
  range[4][0]=0.;
  range[4][1]=150.;
  // min pt
  range[5][0]=0.;
  range[5][1]=100.;
  // delta phi
  range[6][0]=0.;
  range[6][1]=180.;
  // njets
  range[7][0]=0.;
  range[7][1]=7.;
  // nvtx
  range[8][0]=1.;
  range[8][1]=21.;

  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="PFMET";
  xaxisLabel[1]="projected E_{T}^{miss}";
  xaxisLabel[2]="m_{T}^{ll E_{T}^{miss}}";
  xaxisLabel[3]="m_{ll}";
  xaxisLabel[4]="p_{T}^{l,max}";
  xaxisLabel[5]="p_{T}^{l,min}";
  xaxisLabel[6]="#Delta #phi_{ll}";
  xaxisLabel[7]="n jets";
  xaxisLabel[8]="n vtx (DA)";

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
  cut[0]="(finalLeptons)*";
  cut[1]="(WWSel)*";
  cut[2]="(WWSel1j)*";
  cut[3]="(step[17] && njets==0)*"; // final 0j
  cut[4]="(step[24] && njets==1)*"; // final 1j
  //  cut[5]="(step[8] && njets==1 && projMet>30 && bTagTrackCount<2.1 && nSoftMu==0)*"; // WW 1j, relaxed MET

  char lumistr[5];
  sprintf(lumistr,"%.1f",lumi);
  TString intLumi=TString(lumistr);     
  TFile *_file[NSPECIES];
  TTree *T1[NSPECIES];
  
  if(!blindData) {
    _file[0]=TFile::Open(files[0]);
    T1[0] = (TTree*)_file[0]->Get("T1");
  } else T1[0] = 0;
  
  for (int i=1;i<NSPECIES;++i) {
    _file[i]=TFile::Open(files[i]);
    T1[i] = (TTree*)_file[i]->Get("T1");
   }

  int nspeciesToRun=NSPECIES;

  for (int z=0;z<NVARIABLES;++z)
    {
      for (int j=0;j<NCUTS;++j)
	{
          int firstSpecie = 0;
          if(blindData) firstSpecie = 1;
	  for (int i=firstSpecie;i<nspeciesToRun;++i)
	    {
	      fOut->cd();
	      TString histoName=variables[z]+"_W_"+species[i]+"_"+TString(icut[j]);
	      std::cout << "Producing " << histoName << std::endl;
	      if (T1[i]==0)
		{
		  std::cout << "Species " << i << " Tree not found" << std::endl;
		  return;
		}
              int jetbin = -1;
              if(j==1 || j==3) jetbin = 0;
              if(j==2 || j==4 || j == 5) jetbin = 1;
              if(i>0) {
                if(j>0) { // scalefactors are valid from WW level on
                  if(i!=2) T1[i]->Project(histoName,variables[z],cut[j]+TString("weight*puweight*")+scalefactor_datadriven[i][jetbin]);
                  else T1[i]->Project(histoName,variables[z],cut[j]+TString("puweight*")+scalefactor_datadriven[i][jetbin]); // W+jets uses the FR data weight only
                } else {
                  if(i!=2) T1[i]->Project(histoName,variables[z],cut[j]+TString("weight*puweight"));
                  else T1[i]->Project(histoName,variables[z],cut[j]+TString("puweight")); // W+jets uses the FR data weight only
                }
              } else T1[i]->Project(histoName,variables[z],cut[j]+"weight");
	      std::cout << "Done " << histoName << std::endl;
	    }
          
          LatinoPlot myPlot;
          myPlot.setLumi(lumi);
          myPlot.addLabel("");
          myPlot.setLabel((xaxisLabel[z]).Data());
          myPlot.setUnits((units[z]).Data());
          myPlot.setMass(130);
          myPlot.setMCHist(iHWW,     histos[1][j][z]);
          myPlot.setMCHist(iWJets,   histos[2][j][z]);
          myPlot.setMCHist(iWZ,      histos[3][j][z]);
          myPlot.setMCHist(iTop,     histos[4][j][z]);
          myPlot.setMCHist(iZJets,   histos[5][j][z]);
          myPlot.setMCHist(iWW,      histos[6][j][z]);

          myPlot.setDataHist(histos[0][j][z]);
          
	  // Draw
	  //--------------------------------------------------------------------
	  TCanvas* c1 = new TCanvas(Form("test_%d_%d_lin", z, j),
				    Form("test_%d_%d_lin", z, j));

 	  c1->SetLogy(0);

          if(j == 0 || z == 7) myPlot.Draw();
          else myPlot.Draw(2);
          //          else if(j == 1) myPlot.Draw(2);
          //          else  myPlot.Draw(4);

          c1->GetFrame()->DrawClone();

	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".lin.png");
          //	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".root");
	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".lin.eps");
          //	  c1->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".pdf");

	  TCanvas* c2 = new TCanvas(Form("test_%d_%d_log", z, j),
				    Form("test_%d_%d_log", z, j));

	  c2->SetLogy(1);

          myPlot.Draw();

          c2->GetFrame()->DrawClone();
          
	  c2->SaveAs(plotsDir+variables[z]+"MCOnly_"+TString(icut[j])+"_"+suffix+".log.png");


	}
    }
  
  fOut->Write();
  fOut->Close();
  
}
