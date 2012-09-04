#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>

enum channels { of0j, of1j, sf0j, sf1j };
using namespace std;

string getChannelSuffix(int channel) {
  if(channel==of0j) return string("of_0j");
  if(channel==of1j) return string("of_1j");
  if(channel==sf0j) return string("sf_0j");
  if(channel==sf1j) return string("sf_1j");
  return string("ERROR! Unclassified channel!");
}

string getStringChannel(int channel) {
  if(channel==of0j) return string("channel>=2 && njet==0");
  if(channel==of1j) return string("channel>=2 && njet==1");
  if(channel==sf0j) return string("channel<2 && njet==0");
  if(channel==sf1j) return string("channel<2 && njet==1");
  return string("ERROR! Unclassified channel!");
}

void cutMinimum(TH2F* h, float threshold) {
  TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
     for(int j = 1; j <= hist->GetNbinsY(); ++j) {
       float val = hist->GetBinContent(i,j);
       if(val<threshold) h->SetBinContent(i,j,threshold);
     }
  }
}

void smooth(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float count = 0.;
            float val = hist->GetBinContent(i,j);
            if (val<threshold) {
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.;
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.;
                if (i-1 != 0 && j-1 != 0)                                      val += hist->GetBinContent(i-1,j-1);
                else                                                           count -= 1.;
                if (i-1 != 0 && j+1 != hist->GetNbinsY()+1)                    val += hist->GetBinContent(i-1,j+1);
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1 && j-1 != 0)                    val += hist->GetBinContent(i+1,j-1);
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1 && j+1 != hist->GetNbinsY()+1)  val += hist->GetBinContent(i+1,j+1);
                else                                                           count -= 1.;
                val /= (9.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothSwiss(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (5.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothVertical(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_vert").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothHorizontal(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_hori").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) { 
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}


void normalize(TH2F* hist) {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double histsum = 0;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) histsum += hist->GetBinContent(i,j);
        if (histsum > 0) histsum = 1.0/histsum;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) hist->SetBinContent(i,j, hist->GetBinContent(i,j) * histsum);
    }
}

void all(int cha, float dphiMin, float dphiMax) {

  float mrMin=50;
  float mrMax=500;

  gStyle->SetOptStat(0);

  TString cut1(getStringChannel(cha).c_str());
  stringstream fitrangecut;
  fitrangecut << "mr > " << mrMin << " && mr < " << mrMax;
  fitrangecut << "&& dphillr > " << dphiMin << " && dphillr  < " << dphiMax;
  TString cut2(fitrangecut.str().c_str());
  TString cutHiggsMC = TString("(")+cut1+TString(" && ")+cut2+TString(")*effW*puW");
  TString cutMC = TString("(")+cut1+TString(" && ")+cut2+TString(")*baseW*effW*puW");
  TString cutLooseLoose = TString("(")+cut1+TString(" && ")+cut2+TString(")*fake2W");

  float mrBinSize = 10.0; // GeV
  float dphiBinSize = 5.0 * TMath::Pi() / 180.; // 5 degrees
  int yNBins = (int)(abs(dphiMax-dphiMin))/dphiBinSize;

  int xNBins = 17;
  Double_t xLowerEdges[xNBins];
  xLowerEdges[0]=mrMin;
  xLowerEdges[1]=70;
  // do constant binning up to 210 GeV (last mass point mH=180 GeV)
  //  for(Double_t x=0; x<=220; x+=mrBinSize) {
  for(int e=2; e<=15; ++e) {
    xLowerEdges[e] = 70 + ((e-1)*mrBinSize);
  }
  xLowerEdges[16]=mrMax;

  cout << "Filling the backgrounds now..." << endl;

  // input files and trees
  TFile *fileDataLooseLoose  = TFile::Open("results_data/datasets_trees/dataset_looseloose_wwbits.root");
  TFile *fileWW  = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileZjets  = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileOthers = TFile::Open("results/datasets_trees/others_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");

  TTree *treeDataLooseLoose   = (TTree*)fileDataLooseLoose->Get("latinoFitSkim");
  TTree *treeWW  = (TTree*)fileWW->Get("latinoFitSkim");
  TTree *treeZjets  = (TTree*)fileZjets->Get("latinoFitSkim");
  TTree *treeOthers = (TTree*)fileOthers->Get("latinoFitSkim");
  TTree *treeTop = (TTree*)fileTop->Get("latinoFitSkim");

  TH2F *bkg_ww     = new TH2F((string("hist2D_bkg_ww_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_top    = new TH2F((string("hist2D_bkg_top_")+getChannelSuffix(cha)).c_str(),   "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_dy     = new TH2F((string("hist2D_bkg_dy_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_wj     = new TH2F((string("hist2D_bkg_wj_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_others = new TH2F((string("hist2D_bkg_others_")+getChannelSuffix(cha)).c_str(),"",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);

  treeWW->Project(bkg_ww->GetName(),"dphillr:mr",cutMC.Data());
  treeTop->Project(bkg_top->GetName(),"dphillr:mr",cutMC.Data());
  treeZjets->Project(bkg_dy->GetName(),"dphillr:mr",cutMC.Data());
  treeOthers->Project(bkg_others->GetName(),"dphillr:mr",cutMC.Data());
  treeDataLooseLoose->Project(bkg_wj->GetName(),"dphillr:mr",cutLooseLoose.Data());

  // and now the signals
  TH2F *sig_higgs     = new TH2F((string("hist2D_sig_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  int mH[13] = {110,115,120,125,130,135,140,145,150,155,160,170,180};
  for(int i=0; i<13;i++) {
    char higgssample[1000];
    sprintf(higgssample,"results/datasets_trees/H%d_ll.root",mH[i]);
    TFile *fileHiggs  = TFile::Open(higgssample);
    TTree *treeHiggs  = (TTree*)fileHiggs->Get("latinoFitSkim");

    TH2F *sig_higgs_tmp = (TH2F*)(sig_higgs->Clone((std::string(sig_higgs->GetName())+"_tmp").c_str()));

    cout << "Filling mass mH = " << mH[i] << "..." << endl;

    // do not use the xsection to populate the 2D plane, since we need the dphi as a function of Higgs mass
    // we will normalize slices in MR later to make the conditional PDF
    treeHiggs->Project(sig_higgs_tmp->GetName(),"dphillr:mr",cutHiggsMC.Data()); 

    sig_higgs->Add(sig_higgs_tmp);

    delete sig_higgs_tmp;
  }

  // === SMOOTHING ===
  smooth(bkg_dy,     1.0);
  cutMinimum(bkg_wj, 0.0);
  smooth(bkg_wj, 1.0);

  normalize(sig_higgs);
  normalize(bkg_ww);
  normalize(bkg_top);
  normalize(bkg_dy);
  normalize(bkg_wj);
  normalize(bkg_others);

  // === SAVING ===
  TFile *fileOut = TFile::Open("hww2DShapes.root","update");
  fileOut->cd();
  sig_higgs->Write();
  bkg_ww->Write();
  bkg_top->Write();
  bkg_dy->Write();
  bkg_wj->Write();
  bkg_others->Write();

  fileOut->Close();

}

void make2DHistogram() {
  
  TFile *fileOut = TFile::Open("hww2DShapes.root","recreate");
  fileOut->Close();

  float dphiMin, dphiMax;
  dphiMin = 0.0;
  dphiMax = TMath::Pi();

  for(int i=0;i<4;++i) all(i,dphiMin,dphiMax);

}
