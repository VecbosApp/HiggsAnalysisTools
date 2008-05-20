{
// style
TStyle *tesiStyle = new TStyle("tesiStyle","");
tesiStyle->SetCanvasColor(0);
tesiStyle->SetFrameFillColor(0);
tesiStyle->SetStatColor(0);
tesiStyle->SetOptStat(0000);
tesiStyle->SetOptFit(1111);
tesiStyle->SetTitleFillColor(0);
tesiStyle->SetCanvasBorderMode(0);
tesiStyle->SetPadBorderMode(0);
tesiStyle->SetFrameBorderMode(0);
tesiStyle->cd();

int whichStudy = 0;  // 0 = high vs low pt
                     // 1 = signal vs background

// files
TFile fileSgnEBshow  ("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_sgn_showeringEB.root");
TFile fileSgnEEshow  ("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_sgn_showeringEE.root");
TFile fileSgnEBgolden("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_sgn_goldenEB.root");
TFile fileSgnEEgolden("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_sgn_goldenEE.root");
TFile fileBkgEBshow  ("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_bkg_showeringEB.root");
TFile fileBkgEEshow  ("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_bkg_showeringEE.root");
TFile fileBkgEBgolden("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_bkg_goldenEB.root");
TFile fileBkgEEgolden("/data/crovelli/HiggsAnalysis_2008/OfflineAnalysis/HiggsAnalysisTools/outHistos_bkg_goldenEE.root");

TH1F *HH_dEta[8], *HH_dPhi[8], *HH_HoE[8], *HH_S9S25[8], *HH_See[8], *HH_EoPout[8];
TH1F *HL_dEta[8], *HL_dPhi[8], *HL_HoE[8], *HL_S9S25[8], *HL_See[8], *HL_EoPout[8];

// golden, signal, EB
HH_dEta[0]   = (TH1F*)fileSgnEBgolden -> Get("HH_dEta");
HH_dPhi[0]   = (TH1F*)fileSgnEBgolden -> Get("HH_dPhi");
HH_HoE[0]    = (TH1F*)fileSgnEBgolden -> Get("HH_HoE");
HH_S9S25[0]  = (TH1F*)fileSgnEBgolden -> Get("HH_S9s25");
HH_See[0]    = (TH1F*)fileSgnEBgolden -> Get("HH_See");
HH_EoPout[0] = (TH1F*)fileSgnEBgolden -> Get("HH_EoPout");
HL_dEta[0]   = (TH1F*)fileSgnEBgolden -> Get("HL_dEta");
HL_dPhi[0]   = (TH1F*)fileSgnEBgolden -> Get("HL_dPhi");
HL_HoE[0]    = (TH1F*)fileSgnEBgolden -> Get("HL_HoE");
HL_S9S25[0]  = (TH1F*)fileSgnEBgolden -> Get("HL_S9s25");
HL_See[0]    = (TH1F*)fileSgnEBgolden -> Get("HL_See");
HL_EoPout[0] = (TH1F*)fileSgnEBgolden -> Get("HL_EoPout");

// golden, signal, EE
HH_dEta[1]   = (TH1F*)fileSgnEEgolden -> Get("HH_dEta");
HH_dPhi[1]   = (TH1F*)fileSgnEEgolden -> Get("HH_dPhi");
HH_HoE[1]    = (TH1F*)fileSgnEEgolden -> Get("HH_HoE");
HH_S9S25[1]  = (TH1F*)fileSgnEEgolden -> Get("HH_S9s25");
HH_See[1]    = (TH1F*)fileSgnEEgolden -> Get("HH_See");
HH_EoPout[1] = (TH1F*)fileSgnEEgolden -> Get("HH_EoPout");
HL_dEta[1]   = (TH1F*)fileSgnEEgolden -> Get("HL_dEta");
HL_dPhi[1]   = (TH1F*)fileSgnEEgolden -> Get("HL_dPhi");
HL_HoE[1]    = (TH1F*)fileSgnEEgolden -> Get("HL_HoE");
HL_S9S25[1]  = (TH1F*)fileSgnEEgolden -> Get("HL_S9s25");
HL_See[1]    = (TH1F*)fileSgnEEgolden -> Get("HL_See");
HL_EoPout[1] = (TH1F*)fileSgnEEgolden -> Get("HL_EoPout");

// showering, signal, EB
HH_dEta[2]   = (TH1F*)fileSgnEBshow -> Get("HH_dEta");
HH_dPhi[2]   = (TH1F*)fileSgnEBshow -> Get("HH_dPhi");
HH_HoE[2]    = (TH1F*)fileSgnEBshow -> Get("HH_HoE");
HH_S9S25[2]  = (TH1F*)fileSgnEBshow -> Get("HH_S9s25");
HH_See[2]    = (TH1F*)fileSgnEBshow -> Get("HH_See");
HH_EoPout[2] = (TH1F*)fileSgnEBshow -> Get("HH_EoPout");
HL_dEta[2]   = (TH1F*)fileSgnEBshow -> Get("HL_dEta");
HL_dPhi[2]   = (TH1F*)fileSgnEBshow -> Get("HL_dPhi");
HL_HoE[2]    = (TH1F*)fileSgnEBshow -> Get("HL_HoE");
HL_S9S25[2]  = (TH1F*)fileSgnEBshow -> Get("HL_S9s25");
HL_See[2]    = (TH1F*)fileSgnEBshow -> Get("HL_See");
HL_EoPout[2] = (TH1F*)fileSgnEBshow -> Get("HL_EoPout");

// showering, signal, EE
HH_dEta[3]   = (TH1F*)fileSgnEEshow -> Get("HH_dEta");
HH_dPhi[3]   = (TH1F*)fileSgnEEshow -> Get("HH_dPhi");
HH_HoE[3]    = (TH1F*)fileSgnEEshow -> Get("HH_HoE");
HH_S9S25[3]  = (TH1F*)fileSgnEEshow -> Get("HH_S9s25");
HH_See[3]    = (TH1F*)fileSgnEEshow -> Get("HH_See");
HH_EoPout[3] = (TH1F*)fileSgnEEshow -> Get("HH_EoPout");
HL_dEta[3]   = (TH1F*)fileSgnEEshow -> Get("HL_dEta");
HL_dPhi[3]   = (TH1F*)fileSgnEEshow -> Get("HL_dPhi");
HL_HoE[3]    = (TH1F*)fileSgnEEshow -> Get("HL_HoE");
HL_S9S25[3]  = (TH1F*)fileSgnEEshow -> Get("HL_S9s25");
HL_See[3]    = (TH1F*)fileSgnEEshow -> Get("HL_See");
HL_EoPout[3] = (TH1F*)fileSgnEEshow -> Get("HL_EoPout");

// golden, background, EB
HH_dEta[4]   = (TH1F*)fileBkgEBgolden -> Get("HH_dEta");
HH_dPhi[4]   = (TH1F*)fileBkgEBgolden -> Get("HH_dPhi");
HH_HoE[4]    = (TH1F*)fileBkgEBgolden -> Get("HH_HoE");
HH_S9S25[4]  = (TH1F*)fileBkgEBgolden -> Get("HH_S9s25");
HH_See[4]    = (TH1F*)fileBkgEBgolden -> Get("HH_See");
HH_EoPout[4] = (TH1F*)fileBkgEBgolden -> Get("HH_EoPout");
HL_dEta[4]   = (TH1F*)fileBkgEBgolden -> Get("HL_dEta");
HL_dPhi[4]   = (TH1F*)fileBkgEBgolden -> Get("HL_dPhi");
HL_HoE[4]    = (TH1F*)fileBkgEBgolden -> Get("HL_HoE");
HL_S9S25[4]  = (TH1F*)fileBkgEBgolden -> Get("HL_S9s25");
HL_See[4]    = (TH1F*)fileBkgEBgolden -> Get("HL_See");
HL_EoPout[4] = (TH1F*)fileBkgEBgolden -> Get("HL_EoPout");

// golden, background, EE
HH_dEta[5]   = (TH1F*)fileBkgEEgolden -> Get("HH_dEta");
HH_dPhi[5]   = (TH1F*)fileBkgEEgolden -> Get("HH_dPhi");
HH_HoE[5]    = (TH1F*)fileBkgEEgolden -> Get("HH_HoE");
HH_S9S25[5]  = (TH1F*)fileBkgEEgolden -> Get("HH_S9s25");
HH_See[5]    = (TH1F*)fileBkgEEgolden -> Get("HH_See");
HH_EoPout[5] = (TH1F*)fileBkgEEgolden -> Get("HH_EoPout");
HL_dEta[5]   = (TH1F*)fileBkgEEgolden -> Get("HL_dEta");
HL_dPhi[5]   = (TH1F*)fileBkgEEgolden -> Get("HL_dPhi");
HL_HoE[5]    = (TH1F*)fileBkgEEgolden -> Get("HL_HoE");
HL_S9S25[5]  = (TH1F*)fileBkgEEgolden -> Get("HL_S9s25");
HL_See[5]    = (TH1F*)fileBkgEEgolden -> Get("HL_See");
HL_EoPout[5] = (TH1F*)fileBkgEEgolden -> Get("HL_EoPout");

// showering, background, EB
HH_dEta[6]   = (TH1F*)fileBkgEBshow -> Get("HH_dEta");
HH_dPhi[6]   = (TH1F*)fileBkgEBshow -> Get("HH_dPhi");
HH_HoE[6]    = (TH1F*)fileBkgEBshow -> Get("HH_HoE");
HH_S9S25[6]  = (TH1F*)fileBkgEBshow -> Get("HH_S9s25");
HH_See[6]    = (TH1F*)fileBkgEBshow -> Get("HH_See");
HH_EoPout[6] = (TH1F*)fileBkgEBshow -> Get("HH_EoPout");
HL_dEta[6]   = (TH1F*)fileBkgEBshow -> Get("HL_dEta");
HL_dPhi[6]   = (TH1F*)fileBkgEBshow -> Get("HL_dPhi");
HL_HoE[6]    = (TH1F*)fileBkgEBshow -> Get("HL_HoE");
HL_S9S25[6]  = (TH1F*)fileBkgEBshow -> Get("HL_S9s25");
HL_See[6]    = (TH1F*)fileBkgEBshow -> Get("HL_See");
HL_EoPout[6] = (TH1F*)fileBkgEBshow -> Get("HL_EoPout");

// showering, background, EE
HH_dEta[7]   = (TH1F*)fileBkgEEshow -> Get("HH_dEta");
HH_dPhi[7]   = (TH1F*)fileBkgEEshow -> Get("HH_dPhi");
HH_HoE[7]    = (TH1F*)fileBkgEEshow -> Get("HH_HoE");
HH_S9S25[7]  = (TH1F*)fileBkgEEshow -> Get("HH_S9s25");
HH_See[7]    = (TH1F*)fileBkgEEshow -> Get("HH_See");
HH_EoPout[7] = (TH1F*)fileBkgEEshow -> Get("HH_EoPout");
HL_dEta[7]   = (TH1F*)fileBkgEEshow -> Get("HL_dEta");
HL_dPhi[7]   = (TH1F*)fileBkgEEshow -> Get("HL_dPhi");
HL_HoE[7]    = (TH1F*)fileBkgEEshow -> Get("HL_HoE");
HL_S9S25[7]  = (TH1F*)fileBkgEEshow -> Get("HL_S9s25");
HL_See[7]    = (TH1F*)fileBkgEEshow -> Get("HL_See");
HL_EoPout[7] = (TH1F*)fileBkgEEshow -> Get("HL_EoPout");


// adding HL+HH if needed - for S vs B studies
if(whichStudy){
  for(int ii=0; ii<8; ii++){
    HH_dEta[ii].  Add(HL_dEta[ii]);
    HH_dPhi[ii].  Add(HL_dPhi[ii]);
    HH_HoE[ii].   Add(HL_HoE[ii]);
    HH_S9S25[ii]. Add(HL_S9S25[ii]);
    HH_See[ii].   Add(HL_See[ii]);
    HH_EoPout[ii].Add(HL_EoPout[ii]);
  }
}

// rescaling
double scale_HH_dEta[8];
double scale_HH_dPhi[8];
double scale_HH_HoE[8];
double scale_HH_S9S25[8];
double scale_HH_See[8];
double scale_HH_EoPout[8];
double scale_HL_dEta[8];
double scale_HL_dPhi[8];
double scale_HL_HoE[8];
double scale_HL_S9S25[8];
double scale_HL_See[8];
double scale_HL_EoPout[8];

for(int ii=0; ii<8; ii++){
  scale_HH_dEta[ii]   = 1./HH_dEta[ii].Integral();
  scale_HH_dPhi[ii]   = 1./HH_dPhi[ii].Integral();
  scale_HH_HoE[ii]    = 1./HH_HoE[ii].Integral();
  scale_HH_S9S25[ii]  = 1./HH_S9S25[ii].Integral();
  scale_HH_See[ii]    = 1./HH_See[ii].Integral();
  scale_HH_EoPout[ii] = 1./HH_EoPout[ii].Integral();
  scale_HL_dEta[ii]   = 1./HL_dEta[ii].Integral();
  scale_HL_dPhi[ii]   = 1./HL_dPhi[ii].Integral();
  scale_HL_HoE[ii]    = 1./HL_HoE[ii].Integral();
  scale_HL_S9S25[ii]  = 1./HL_S9S25[ii].Integral();
  scale_HL_See[ii]    = 1./HL_See[ii].Integral();
  scale_HL_EoPout[ii] = 1./HL_EoPout[ii].Integral();
}

for(int ii=0; ii<8; ii++){
  HH_dEta[ii].  Scale(scale_HH_dEta[ii]);
  HH_dPhi[ii].  Scale(scale_HH_dPhi[ii]);
  HH_HoE[ii].   Scale(scale_HH_HoE[ii]);
  HH_S9S25[ii]. Scale(scale_HH_S9S25[ii]);
  HH_See[ii].   Scale(scale_HH_See[ii]);
  HH_EoPout[ii].Scale(scale_HH_EoPout[ii]);
  HL_dEta[ii].  Scale(scale_HL_dEta[ii]);
  HL_dPhi[ii].  Scale(scale_HL_dPhi[ii]);
  HL_HoE[ii].   Scale(scale_HL_HoE[ii]);
  HL_S9S25[ii]. Scale(scale_HL_S9S25[ii]);
  HL_See[ii].   Scale(scale_HL_See[ii]);
  HL_EoPout[ii].Scale(scale_HL_EoPout[ii]);
}

// cosmetics
if(!whichStudy){
  for(int ii=0; ii<8; ii++){
    HH_dEta[ii].  SetFillColor(1);    HH_dEta[ii].  SetFillStyle(3004);
    HH_dPhi[ii].  SetFillColor(1);    HH_dPhi[ii].  SetFillStyle(3004);
    HH_HoE[ii].   SetFillColor(1);    HH_HoE[ii].   SetFillStyle(3004);
    HH_S9S25[ii]. SetFillColor(1);    HH_S9S25[ii]. SetFillStyle(3004);
    HH_See[ii].   SetFillColor(1);    HH_See[ii].   SetFillStyle(3004);
    HH_EoPout[ii].SetFillColor(1);    HH_EoPout[ii].SetFillStyle(3004);
    HL_dEta[ii].  SetFillColor(3);    HL_dEta[ii].  SetFillStyle(3005);
    HL_dPhi[ii].  SetFillColor(3);    HL_dPhi[ii].  SetFillStyle(3005);
    HL_HoE[ii].   SetFillColor(3);    HL_HoE[ii].   SetFillStyle(3005);
    HL_S9S25[ii]. SetFillColor(3);    HL_S9S25[ii]. SetFillStyle(3005);
    HL_See[ii].   SetFillColor(3);    HL_See[ii].   SetFillStyle(3005);
    HL_EoPout[ii].SetFillColor(3);    HL_EoPout[ii].SetFillStyle(3005);
  }
}

if(whichStudy){
  for(int ii=0; ii<4; ii++){
    HH_dEta[ii].  SetFillColor(1);    HH_dEta[ii].  SetFillStyle(3004);
    HH_dPhi[ii].  SetFillColor(1);    HH_dPhi[ii].  SetFillStyle(3004);
    HH_HoE[ii].   SetFillColor(1);    HH_HoE[ii].   SetFillStyle(3004);
    HH_S9S25[ii]. SetFillColor(1);    HH_S9S25[ii]. SetFillStyle(3004);
    HH_See[ii].   SetFillColor(1);    HH_See[ii].   SetFillStyle(3004);
    HH_EoPout[ii].SetFillColor(1);    HH_EoPout[ii].SetFillStyle(3004);
  }
  for(int ii=4; ii<8; ii++){
    HH_dEta[ii].  SetFillColor(3);    HH_dEta[ii].  SetFillStyle(3005);
    HH_dPhi[ii].  SetFillColor(3);    HH_dPhi[ii].  SetFillStyle(3005);
    HH_HoE[ii].   SetFillColor(3);    HH_HoE[ii].   SetFillStyle(3005);
    HH_S9S25[ii]. SetFillColor(3);    HH_S9S25[ii]. SetFillStyle(3005);
    HH_See[ii].   SetFillColor(3);    HH_See[ii].   SetFillStyle(3005);
    HH_EoPout[ii].SetFillColor(3);    HH_EoPout[ii].SetFillStyle(3005);
  }
}

// axes titles
for(int ii=0; ii<8; ii++){
  HH_dEta[ii].  GetXaxis()->SetTitle("#Delta #eta");
  HH_dPhi[ii].  GetXaxis()->SetTitle("#Delta #phi");
  HH_HoE[ii].   GetXaxis()->SetTitle("H/E");
  HH_S9S25[ii]. GetXaxis()->SetTitle("S9/S25");
  HH_See[ii].   GetXaxis()->SetTitle("#sigma_{#eta #eta}");
  HH_EoPout[ii].GetXaxis()->SetTitle("E/P_{out}");
  HL_dEta[ii].  GetXaxis()->SetTitle("#Delta #eta");
  HL_dPhi[ii].  GetXaxis()->SetTitle("#Delta #phi");
  HL_HoE[ii].   GetXaxis()->SetTitle("H/E");
  HL_S9S25[ii]. GetXaxis()->SetTitle("S9/S25");
  HL_See[ii].   GetXaxis()->SetTitle("#sigma_{#eta #eta}");
  HL_EoPout[ii].GetXaxis()->SetTitle("E/P_{out}");
}


// ------------- plots -----------------------------------

// high vs low pt study
if(!whichStudy){

  TLegend HvsLS(0.50,0.6,0.75,0.82);
  HvsLS.AddEntry(HH_dPhi[0], "high p_{T}, signal","f");
  HvsLS.AddEntry(HL_dPhi[0], "low p_{T},  signal","f");
  HvsLS.SetFillColor(0);
  HvsLS.SetBorderSize(0.4);
  
  TLegend HvsLB(0.50,0.6,0.75,0.82);
  HvsLB.AddEntry(HH_dPhi[4], "high p_{T}, background","f");
  HvsLB.AddEntry(HL_dPhi[4], "low p_{T},  background","f");
  HvsLB.SetFillColor(0);
  HvsLB.SetBorderSize(0.4);
  
  TCanvas *csgEB = new TCanvas("csgEB", "golden, barrel, signal",1);  
  csgEB.SetLogy();
  HH_dPhi[0].Draw();
  HL_dPhi[0].Draw("same");
  HvsLS.Draw();
  csgEB->Print("EBgoldenSignal.eps");
  
  TCanvas *csgEE = new TCanvas("csgEE", "golden, endcap, signal",1);  
  csgEE.SetLogy();
  HH_dPhi[1].Draw();
  HL_dPhi[1].Draw("same");
  HvsLS.Draw();
  csgEE->Print("EEgoldenSignal.eps");
  
  TCanvas *cssEB = new TCanvas("cssEB", "showering, barrel, signal",1);  
  cssEB.SetLogy();
  HH_dPhi[2].Draw();
  HL_dPhi[2].Draw("same");
  HvsLS.Draw();
  cssEB->Print("EBshoweringSignal.eps");
  
  TCanvas *cssEE = new TCanvas("cssEE", "showering, endcap, signal",1);  
  cssEE.SetLogy();
  HH_dPhi[3].Draw();
  HL_dPhi[3].Draw("same");
  HvsLS.Draw();
  cssEE->Print("EEshoweringSignal.eps");
  
  TCanvas *cbgEB = new TCanvas("cbgEB", "golden, barrel, background",1);  
  cbgEB.SetLogy();
  HH_dPhi[4].Draw();
  HL_dPhi[4].Draw("same");
  HvsLB.Draw();
  cbgEB->Print("EBgoldenBackground.eps");
  
  TCanvas *cbgEE = new TCanvas("cbgEE", "golden, endcap, background",1);  
  cbgEE.SetLogy();
  HH_dPhi[5].Draw();
  HL_dPhi[5].Draw("same");
  HvsLB.Draw();
  cbgEE->Print("EEgoldenBackground.eps");
  
  TCanvas *cbsEB = new TCanvas("cbsEB", "showering, barrel, background",1);  
  cbsEB.SetLogy();
  HH_dPhi[6].Draw();
  HL_dPhi[6].Draw("same");
  HvsLB.Draw();
  cbsEB->Print("EBshoweringBackground.eps");
  
  TCanvas *cbsEE = new TCanvas("cbsEE", "showering, endcap, background",1);  
  cbsEE.SetLogy();
  HH_dPhi[7].Draw();
  HL_dPhi[7].Draw("same");
  HvsLB.Draw();
  cbsEE->Print("EEshoweringBackground.eps");
}

// signal vs background study
if(whichStudy){

  TLegend SvsB(0.50,0.6,0.75,0.82);
  SvsB.AddEntry(HH_S9S25[0], "signal","f");
  SvsB.AddEntry(HH_S9S25[4], "background","f");
  SvsB.SetFillColor(0);
  SvsB.SetBorderSize(0.4);
    
  TCanvas *cgEB = new TCanvas("cgEB", "golden, barrel",1);  
  cgEB.SetLogy();
  HH_S9S25[0].Draw();
  HH_S9S25[4].Draw("same");
  SvsB.Draw();
  cgEB->Print("EBgoldenSignalVsBack.eps");

  TCanvas *cgEE = new TCanvas("cgEE", "golden, endcap",1);  
  cgEE.SetLogy();
  HH_S9S25[1].Draw();
  HH_S9S25[5].Draw("same");
  SvsB.Draw();
  cgEE->Print("EEgoldenSignalVsBack.eps");

  TCanvas *csEB = new TCanvas("csEB", "showering, barrel",1);  
  csEB.SetLogy();
  HH_S9S25[2].Draw();
  HH_S9S25[6].Draw("same");
  SvsB.Draw();
  csEB->Print("EBshoweringSignalVsBack.eps");  

  TCanvas *csEE = new TCanvas("csEE", "showering, endcap",1);  
  csEE.SetLogy();
  HH_S9S25[3].Draw();
  HH_S9S25[7].Draw("same");
  SvsB.Draw();
  csEE->Print("EEshoweringSignalVsBack.eps");

}
}
