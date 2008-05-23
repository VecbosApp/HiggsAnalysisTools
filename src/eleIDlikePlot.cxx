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

// which study : 0 = 1 like; 1 = divided by pt
int whichStudy = 1;

// files
TFile fileSgn("../../Results/likelihood/outHistosSgnm160.root");
TFile fileBkg("../../Results/likelihood/outHistosBkgm160.root");

TH1F *HH_like[2], *HL_like[2];
HH_like[0] = (TH1F*)fileSgn -> Get("HH_like");
HL_like[0] = (TH1F*)fileSgn -> Get("HL_like");
HH_like[1] = (TH1F*)fileBkg -> Get("HH_like");
HL_like[1] = (TH1F*)fileBkg -> Get("HL_like");

// adding HL+HH if needed 
if(!whichStudy){
  for(int ii=0; ii<2; ii++){ HH_like[ii].Add(HL_like[ii]); }
}

// rescaling
double scale_HH[2];
double scale_HL[2];
for(int ii=0; ii<2; ii++){
  scale_HH[ii] = 1./HH_like[ii].Integral();
  scale_HL[ii] = 1./HL_like[ii].Integral();
}

for(int ii=0; ii<2; ii++){
  HH_like[ii].Scale(scale_HH[ii]);
  HL_like[ii].Scale(scale_HL[ii]);
}

// cosmetics
HH_like[0].SetFillColor(1);    HH_like[0].SetFillStyle(3004);
HL_like[0].SetFillColor(1);    HL_like[0].SetFillStyle(3004);
HH_like[1].SetFillColor(3);    HH_like[1].SetFillStyle(3005);
HL_like[1].SetFillColor(3);    HL_like[1].SetFillStyle(3005);

// axes titles
for(int ii=0; ii<2; ii++){
  HH_like[ii].GetXaxis()->SetTitle("likelihood");
  HL_like[ii].GetXaxis()->SetTitle("likelihood");
}


// ------------- plots -----------------------------------

// high vs low pt study
if(whichStudy){
  TLegend Hleg(0.50,0.6,0.75,0.82);
  Hleg.AddEntry(HH_like[0], "high p_{T}, signal","f");
  Hleg.AddEntry(HH_like[1], "high p_{T}, background","f");
  Hleg.SetFillColor(0);
  Hleg.SetBorderSize(0.4);

  TLegend Lleg(0.50,0.6,0.75,0.82);
  Lleg.AddEntry(HL_like[0], "low p_{T}, signal","f");
  Lleg.AddEntry(HL_like[1], "low p_{T}, background","f");
  Lleg.SetFillColor(0);
  Lleg.SetBorderSize(0.4);

  TCanvas *ch = new TCanvas("ch", "high p_{T}",1);  
  ch.SetLogy();
  HH_like[0].Draw();
  HH_like[1].Draw("same");
  Hleg.Draw();
  ch->Print("HighPt.png");

  TCanvas *cl = new TCanvas("cl", "low p_{T}",1);  
  cl.SetLogy();
  HL_like[0].Draw();
  HL_like[1].Draw("same");
  Lleg.Draw();
  cl->Print("LowPt.png");
}

// high and low pt together study
if(!whichStudy){
  TLegend Hleg(0.50,0.6,0.75,0.82);
  Hleg.AddEntry(HH_like[0], "signal","f");
  Hleg.AddEntry(HH_like[1], "background","f");
  Hleg.SetFillColor(0);
  Hleg.SetBorderSize(0.4);
  
  TCanvas *c = new TCanvas("c", "",1);  
  c.SetLogy();
  HH_like[0].Draw();
  HH_like[1].Draw("same");
  Hleg.Draw();
  c->Print("Like.png");
}

}
