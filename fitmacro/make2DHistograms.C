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

#include "YieldMaker.h"
#include "FitSelection.hh"

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

void fakeAlternativeShapes(TH2F* h, TH2F *hUp, TH2F *hDn) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) { 
    for(int j = 1; j <= h->GetNbinsY(); ++j) {
      double alpha = 1.0 + 0.5*(j-0.5*h->GetNbinsY())/h->GetNbinsY();
      double beta  = 1.0/pow(alpha,1.5);
      hUp->SetBinContent(i,j, alpha*h->GetBinContent(i,j));
      hDn->SetBinContent(i,j, beta*h->GetBinContent(i,j));
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

  cout << "Filling 2D maps for channel " << cha << endl;

  FitSelection sel;
  sel.mrmin=50;
  sel.mrmax=500;
  sel.dphimin=dphiMin;
  sel.dphimax=dphiMax;

  gStyle->SetOptStat(0);

  DataYieldMaker  ymaker_data;
  YieldMaker      ymaker_hi;
  YieldMaker      ymaker_qqww;
  YieldMaker      ymaker_ggww;
  YieldMaker      ymaker_top;
  YieldMaker      ymaker_dysf;
  YieldMaker      ymaker_dyof;
  YieldMaker      ymaker_others;
  WJetsYieldMaker ymaker_wj;

  string treeFolder("latinos_tree_skim_of");

  ymaker_data   .fill(treeFolder+"/data/latino_RunA_892pbinv.root");
  ymaker_data   .fill(treeFolder+"/data/latino_RunB_4404pbinv.root");
  ymaker_data   .fill(treeFolder+"/data/latino_RunC_6807pbinv.root");
  ymaker_qqww   .fill(treeFolder+"/nominals/latino_000_WWJets2LMad.root");
  ymaker_ggww   .fill(treeFolder+"/nominals/latino_001_GluGluToWWTo4L.root");
  ymaker_top    .fill(treeFolder+"/nominals/latino_019_TTTo2L2Nu2B.root");
  ymaker_top    .fill(treeFolder+"/nominals/latino_011_TtWFullDR.root");
  ymaker_top    .fill(treeFolder+"/nominals/latino_012_TbartWFullDR.root");
  ymaker_dysf   .fill(treeFolder+"/nominals/latino_037_DY50toLLMad.root");
  ymaker_dysf   .fill(treeFolder+"/nominals/latino_036_DY10toLLMad.root");
  ymaker_dyof   .fill(treeFolder+"/nominals/latino_RunABC_DYtt_8fb.root");
  ymaker_wj     .fill(treeFolder+"/wjets/latino_RunABC_LooseLoose_skimww.root");
  ymaker_others .fill(treeFolder+"/nominals/latino_074_WZJetsMad.root");
  ymaker_others .fill(treeFolder+"/nominals/latino_075_ZZJetsMad.root");

  float mrBinSize = 10.0; // GeV
  float dphiBinSize = 5.0 * TMath::Pi() / 180.; // 5 degrees
  int yNBins = (int)(fabs(dphiMax-dphiMin))/dphiBinSize;

  int xNBins = 17;
  Double_t xLowerEdges[xNBins];
  xLowerEdges[0]=sel.mrmin;
  xLowerEdges[1]=70;
  // do constant binning up to 210 GeV (last mass point mH=180 GeV)
  //  for(Double_t x=0; x<=220; x+=mrBinSize) {
  for(int e=2; e<=15; ++e) {
    xLowerEdges[e] = 70 + ((e-1)*mrBinSize);
  }
  xLowerEdges[16]=sel.mrmax;

  cout << "Filling the backgrounds now..." << endl;

  TH2F *bkg_qqww   = new TH2F((string("hist2D_bkg_qqww_")+getChannelSuffix(cha)).c_str(),  "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_ggww   = new TH2F((string("hist2D_bkg_ggww_")+getChannelSuffix(cha)).c_str(),  "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_top    = new TH2F((string("hist2D_bkg_top_")+getChannelSuffix(cha)).c_str(),   "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_dy     = new TH2F((string("hist2D_bkg_dy_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_wj     = new TH2F((string("hist2D_bkg_wj_")+getChannelSuffix(cha)).c_str(),    "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *bkg_others = new TH2F((string("hist2D_bkg_others_")+getChannelSuffix(cha)).c_str(),"",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);
  TH2F *sig_higgs  = new TH2F((string("hist2D_sig_")+getChannelSuffix(cha)).c_str(),       "",xNBins-1,xLowerEdges,yNBins,dphiMin,dphiMax);

  ymaker_qqww.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_qqww);
  ymaker_ggww.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_ggww);
  ymaker_top.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_top);
  if(cha==of0j || cha==of1j) ymaker_dyof.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_dy);
  else ymaker_dysf.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_dy);
  ymaker_wj.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_wj);
  ymaker_others.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,bkg_others);

  // and now the signals
  int mH[13] = {110,115,120,125,130,135,140,145,150,155,160,170,180};
  for(int i=0; i<13;i++) {
    cout << "Filling mass mH = " << mH[i] << "..." << endl;
    char gghsample[1000], qqhsample[1000];
    sprintf(gghsample,"%s/nominals/latino_1%d_ggToH%dtoWWTo2LAndTau2Nu.root",treeFolder.c_str(),mH[i],mH[i]);
    sprintf(qqhsample,"%s/nominals/latino_2%d_vbfToH%dtoWWTo2LAndTau2Nu.root",treeFolder.c_str(),mH[i],mH[i]);
    ymaker_hi   .fill(gghsample); 
    ymaker_hi   .fill(qqhsample);
  }

  // do not use the xsection to populate the 2D plane, since we need the dphi as a function of Higgs mass
  // we will normalize slices in MR later to make the conditional PDF
  ymaker_hi.get2DHist(cha,sel.mrmin,sel.mrmax,sel.dphimin,sel.dphimax,sel.mtmin,sel.mtmax,sig_higgs,false); 

  // === SMOOTHING ===
  smooth(bkg_dy,     1.0);
  cutMinimum(bkg_wj, 0.0);
  smooth(bkg_wj, 1.0);

  normalize(sig_higgs);
  normalize(bkg_qqww);
  normalize(bkg_ggww);
  normalize(bkg_top);
  normalize(bkg_dy);
  normalize(bkg_wj);
  normalize(bkg_others);

  
  // === SAVING ===
  TFile *fileOut;
  if(cha==0) fileOut = TFile::Open("hww2DShapes.root","recreate");
  else fileOut = TFile::Open("hww2DShapes.root","update");
  fileOut->cd();
  sig_higgs->Write();
  bkg_qqww->Write();
  bkg_ggww->Write();
  bkg_top->Write();
  bkg_dy->Write();
  bkg_wj->Write();
  bkg_others->Write();

  // === fake alternative shapes ===
  TH2F *sig_higgs_up = (TH2F*)(sig_higgs->Clone((std::string(sig_higgs->GetName())+"_Up").c_str()));
  TH2F *sig_higgs_dn = (TH2F*)(sig_higgs->Clone((std::string(sig_higgs->GetName())+"_Dn").c_str()));
  TH2F *bkg_qqww_up = (TH2F*)(bkg_qqww->Clone((std::string(bkg_qqww->GetName())+"_Up").c_str()));
  TH2F *bkg_qqww_dn = (TH2F*)(bkg_qqww->Clone((std::string(bkg_qqww->GetName())+"_Dn").c_str()));
  TH2F *bkg_ggww_up = (TH2F*)(bkg_ggww->Clone((std::string(bkg_ggww->GetName())+"_Up").c_str()));
  TH2F *bkg_ggww_dn = (TH2F*)(bkg_ggww->Clone((std::string(bkg_ggww->GetName())+"_Dn").c_str()));
  TH2F *bkg_top_up = (TH2F*)(bkg_top->Clone((std::string(bkg_top->GetName())+"_Up").c_str()));
  TH2F *bkg_top_dn = (TH2F*)(bkg_top->Clone((std::string(bkg_top->GetName())+"_Dn").c_str()));
  TH2F *bkg_dy_up = (TH2F*)(bkg_dy->Clone((std::string(bkg_dy->GetName())+"_Up").c_str()));
  TH2F *bkg_dy_dn = (TH2F*)(bkg_dy->Clone((std::string(bkg_dy->GetName())+"_Dn").c_str()));
  TH2F *bkg_wj_up = (TH2F*)(bkg_wj->Clone((std::string(bkg_wj->GetName())+"_Up").c_str()));
  TH2F *bkg_wj_dn = (TH2F*)(bkg_wj->Clone((std::string(bkg_wj->GetName())+"_Dn").c_str()));
  TH2F *bkg_others_up = (TH2F*)(bkg_others->Clone((std::string(bkg_others->GetName())+"_Up").c_str()));
  TH2F *bkg_others_dn = (TH2F*)(bkg_others->Clone((std::string(bkg_others->GetName())+"_Dn").c_str()));

  fakeAlternativeShapes(sig_higgs,sig_higgs_up,sig_higgs_dn);
  fakeAlternativeShapes(bkg_qqww,bkg_qqww_up,bkg_qqww_dn);
  fakeAlternativeShapes(bkg_ggww,bkg_ggww_up,bkg_ggww_dn);
  fakeAlternativeShapes(bkg_top,bkg_top_up,bkg_top_dn);
  fakeAlternativeShapes(bkg_wj,bkg_wj_up,bkg_wj_dn);
  fakeAlternativeShapes(bkg_dy,bkg_dy_up,bkg_dy_dn);
  fakeAlternativeShapes(bkg_others,bkg_others_up,bkg_others_dn);

  sig_higgs_up->Write();
  bkg_qqww_up->Write();
  bkg_ggww_up->Write();
  bkg_top_up->Write();
  bkg_dy_up->Write();
  bkg_wj_up->Write();
  bkg_others_up->Write();

  sig_higgs_dn->Write();
  bkg_qqww_dn->Write();
  bkg_ggww_dn->Write();
  bkg_top_dn->Write();
  bkg_dy_dn->Write();
  bkg_wj_dn->Write();
  bkg_others_dn->Write();

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
