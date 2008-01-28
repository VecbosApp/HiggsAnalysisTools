// to compile with 
// g++ TriggerEff.cpp -o triggerEff `root-config --libs --cflags --glibs` -Wno-deprecated

// input: txt inputFile with list of .root file + corresponding mass

// -------------------------------------
//
// Produce grafici con :
//         efficienza dell'OR dei trigger elettronici 
//         efficienza del trigger single electron          
//         efficienza del trigger single relaxed electron      
//         efficienza del trigger double electron         
//         efficienza del trigger double relaxed electron         
// calcolate rispetto agli eventi totali del tree ridott
// (quelli con (1)2 elettroni provenienti dalla W nell'accettanza )
// 
// Tabella con efficienza ed errore dell'OR dei 2 trigger elettronici  
//
// -------------------------------------


//! c++ includes                                                               
#include <string>
#include <stdio.h>
#include <sstream>
#include <iostream.h>
#include <unistd.h>
#include <fstream>
#include <math.h>

//! ROOT includes                                                             
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TText.h"
#include "TSelector.h"
#include "TApplication.h"
#include "TLatex.h"

int main ( int argc, char **argv)
{
  // reading input/output file names
  char inFileName[150];
  if (argc==2) { strcpy(inFileName,argv[1]); }
  else
    {
      cout << "missing argument. Insert: " << endl; 
      cout << "inputFile with ntuples name + mass" << endl; 
      return 1;
    }

  // -------------------------
  // Arrays x TGraph 
  float massErr_hlt_ele[100];
  float eff_hlt_ele[100], effErr_hlt_ele[100];
  float eff_hlt_se[100],  effErr_hlt_se[100];
  float eff_hlt_ser[100], effErr_hlt_ser[100];
  float eff_hlt_de[100],  effErr_hlt_de[100];
  float eff_hlt_der[100], effErr_hlt_der[100];

  // -------------------------
  // input file: list with ntuples and respective masses
  char Buffer[500]; 
  int IFC = 0;
  char RefFile[50][500];      
  float RefMass[50];        

  ifstream *inputFile = new ifstream(inFileName); 
  while( !(inputFile->eof()) ){  
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))){
      sscanf(Buffer,"%s %f",RefFile[IFC], &RefMass[IFC]);
      cout << RefFile[IFC] << ", " << RefMass[IFC] << endl;
      IFC++; 
    }
  }
  inputFile->close();
  delete inputFile;

  cout << endl;
  cout << endl;
  
  
  // -------------------------
  //  Loop over the reference files
  for ( int my_file = 0; my_file < IFC; my_file++ ){
    
    cout << "Working with file " << RefFile[my_file] << ", for mass = " << RefMass[my_file] << endl;

    // -------------------------
    // loading file:
    TFile *f_ntupla = TFile::Open(RefFile[my_file]);

    // -------------------------  
    // tree construction:   
    TTree *tree = (TTree*) f_ntupla->Get("T1");

    // declaration of leaves types
    int sampleOk;
    int singleElePassedTrg, singleEleRelaxPassedTrg;
    int doubleElePassedTrg, doubleEleRelaxPassedTrg;

    // declaration of branches
    TBranch *b_sampleOk;
    TBranch *b_singleElePassedTrg;
    TBranch *b_singleEleRelaxPassedTrg;
    TBranch *b_doubleElePassedTrg;
    TBranch *b_doubleEleRelaxPassedTrg;

    // setting branch addresses
    tree->SetMakeClass(1);
    tree->SetBranchAddress("sampleOk",               &sampleOk);
    tree->SetBranchAddress("singleElePassedTrg",     &singleElePassedTrg);
    tree->SetBranchAddress("singleEleRelaxPassedTrg",&singleEleRelaxPassedTrg);
    tree->SetBranchAddress("doubleElePassedTrg",     &doubleElePassedTrg);
    tree->SetBranchAddress("doubleEleRelaxPassedTrg",&doubleEleRelaxPassedTrg);

    // getting branch pointers
    b_sampleOk                = tree->GetBranch("sampleOk");
    b_singleElePassedTrg      = tree->GetBranch("singleElePassedTrg");
    b_singleEleRelaxPassedTrg = tree->GetBranch("singleEleRelaxPassedTrg");
    b_doubleElePassedTrg      = tree->GetBranch("doubleElePassedTrg");
    b_doubleEleRelaxPassedTrg = tree->GetBranch("doubleEleRelaxPassedTrg");
      
    // -------------------------  
    // counters
    float ev_total      = 0.;
    float ev_hlt_ele    = 0.;  // se || ser || de || der 
    float ev_hlt_se     = 0.;  // single ele 
    float ev_hlt_ser    = 0.;  // single ele relaxed
    float ev_hlt_de     = 0.;  // double ele 
    float ev_hlt_der    = 0.;  // double ele relaxed
      
    
    // -------------------------  
    // loop on the events
    float nEnt=tree->GetEntries();
    for (int entry=0; entry<nEnt; entry++){ 

      // if (entry%100 == 0){ cout << entry << " " << endl; }

      // charging the branches I need
      b_sampleOk                -> GetEntry(entry);
      b_singleElePassedTrg      -> GetEntry(entry);
      b_singleEleRelaxPassedTrg -> GetEntry(entry);
      b_doubleElePassedTrg      -> GetEntry(entry);
      b_doubleEleRelaxPassedTrg -> GetEntry(entry);

      // skipping events without generated electrons
      if (!sampleOk){ continue; }
      ev_total++;
  
      // to compute trigger efficiency
      if ( singleElePassedTrg       ){ ev_hlt_se++; }
      if ( singleEleRelaxPassedTrg  ){ ev_hlt_ser++; }
      if ( doubleElePassedTrg       ){ ev_hlt_de++; }
      if ( doubleEleRelaxPassedTrg  ){ ev_hlt_der++; }
      if ( singleElePassedTrg || singleEleRelaxPassedTrg || doubleElePassedTrg || doubleEleRelaxPassedTrg ){ ev_hlt_ele++; }

    } // end loop on the events
    
    
    // -------------------------        
    // efficiencies
    float this_massErr     = 0;
    float this_eff_hlt_ele    = ev_hlt_ele/ev_total;
    float this_eff_hlt_se     = ev_hlt_se/ev_total;
    float this_eff_hlt_ser    = ev_hlt_ser/ev_total;
    float this_eff_hlt_de     = ev_hlt_de/ev_total;
    float this_eff_hlt_der    = ev_hlt_der/ev_total;
    float this_effErr_hlt_ele = (pow(ev_hlt_ele*(1-this_eff_hlt_ele),0.5))/ev_hlt_ele;
    float this_effErr_hlt_se  = (pow(ev_hlt_se*(1-this_eff_hlt_se),0.5))/ev_hlt_se;
    float this_effErr_hlt_ser = (pow(ev_hlt_ser*(1-this_eff_hlt_ser),0.5))/ev_hlt_ser;
    float this_effErr_hlt_de  = (pow(ev_hlt_de*(1-this_eff_hlt_de),0.5))/ev_hlt_de;
    float this_effErr_hlt_der = (pow(ev_hlt_der*(1-this_eff_hlt_der),0.5))/ev_hlt_der;
       
    // filling the array with results      
    eff_hlt_se[my_file]      = this_eff_hlt_se;    
    eff_hlt_ser[my_file]     = this_eff_hlt_ser;    
    eff_hlt_de[my_file]      = this_eff_hlt_de;    
    eff_hlt_der[my_file]     = this_eff_hlt_der;    
    eff_hlt_ele[my_file]     = this_eff_hlt_ele;    
    effErr_hlt_ele[my_file]  = this_effErr_hlt_ele;    
    effErr_hlt_se[my_file]   = this_effErr_hlt_se;    
    effErr_hlt_ser[my_file]  = this_effErr_hlt_ser;    
    effErr_hlt_de[my_file]   = this_effErr_hlt_de;    
    effErr_hlt_der[my_file]  = this_effErr_hlt_der;    
    massErr_hlt_ele[my_file] = 0.;

    cout << "events: "  << ev_total << ", trigger ok: " << ev_hlt_ele << endl;  

  }// loop over the files

  
  // statistics
  for ( int my_file = 0; my_file<IFC; my_file++ ){
    cout << endl;
    cout << "file = " << RefFile[my_file] << endl;
    cout << "passing OR Electron HLT: " << eff_hlt_ele[my_file] << " wrt the total # of events" << endl;  
    cout << "passing SingleEl       : " << eff_hlt_se[my_file]  << " wrt the total # of events" << endl;  
    cout << "passing SingleRelaxedEl: " << eff_hlt_ser[my_file] << " wrt the total # of events" << endl;  
    cout << "passing DoubleEl       : " << eff_hlt_de[my_file]  << " wrt the total # of events" << endl;  
    cout << "passing DoubleRelaxedEl: " << eff_hlt_der[my_file] << " wrt the total # of events" << endl;  
    cout << endl;
  }

  cout << endl;
  cout << endl;
  for ( int my_file = 0; my_file<IFC; my_file++ ){
    cout << "file = " << RefFile[my_file] << ", eff HLT = " << eff_hlt_ele[my_file] << " +- " << effErr_hlt_ele[my_file] << endl;     
  }

  // TGraphs
  TGraphErrors *G_Eff_HLT_ELE = new TGraphErrors(IFC, RefMass, eff_hlt_ele, massErr_hlt_ele, effErr_hlt_ele);
  TGraphErrors *G_Eff_HLT_SE  = new TGraphErrors(IFC, RefMass, eff_hlt_se,  massErr_hlt_ele, effErr_hlt_se);
  TGraphErrors *G_Eff_HLT_SER = new TGraphErrors(IFC, RefMass, eff_hlt_ser, massErr_hlt_ele, effErr_hlt_ser);
  TGraphErrors *G_Eff_HLT_DE  = new TGraphErrors(IFC, RefMass, eff_hlt_de,  massErr_hlt_ele, effErr_hlt_de);
  TGraphErrors *G_Eff_HLT_DER = new TGraphErrors(IFC, RefMass, eff_hlt_der, massErr_hlt_ele, effErr_hlt_der);
  G_Eff_HLT_ELE->SetMarkerColor(1);  G_Eff_HLT_ELE->SetMarkerSize(1.2);  G_Eff_HLT_ELE->SetMarkerStyle(3);
  G_Eff_HLT_SE ->SetMarkerColor(2);  G_Eff_HLT_SE ->SetMarkerSize(1.);   G_Eff_HLT_SE ->SetMarkerStyle(8);
  G_Eff_HLT_DE ->SetMarkerColor(3);  G_Eff_HLT_DE ->SetMarkerSize(1.);   G_Eff_HLT_DE ->SetMarkerStyle(8);
  G_Eff_HLT_SER->SetMarkerColor(4);  G_Eff_HLT_SER->SetMarkerSize(1.);   G_Eff_HLT_SER->SetMarkerStyle(22); 
  G_Eff_HLT_DER->SetMarkerColor(6);  G_Eff_HLT_DER->SetMarkerSize(1.2);  G_Eff_HLT_DER->SetMarkerStyle(29);



  TApplication* theApp = new TApplication("App", &argc, argv);

  TStyle *tesiStyle = new TStyle("tesiStyle","");
  tesiStyle->SetCanvasColor(0);
  tesiStyle->SetFrameFillColor(0);
  tesiStyle->SetStatColor(0);
  tesiStyle->SetOptStat(0);
  tesiStyle->SetTitleFillColor(0);
  tesiStyle->SetCanvasBorderMode(0);
  tesiStyle->SetPadBorderMode(0);
  tesiStyle->SetFrameBorderMode(0);
  tesiStyle->cd();
  
  TH2F *myFuffa = new TH2F("myFuffa", "", 100, 120., 190., 100, 0., 1.);
  myFuffa -> SetTitle("Trigger Efficiency");
  myFuffa -> SetLabelSize(0.035,"X");
  myFuffa -> SetLabelSize(0.035,"Y");
  myFuffa -> GetXaxis()->SetTitle("m_{H} (GeV/c^{2})");
  myFuffa -> GetYaxis()->SetTitle("Efficiency");

  TCanvas c ("c","",116,0,206,1);
  myFuffa       -> Draw();
  G_Eff_HLT_ELE -> Draw("P");
  G_Eff_HLT_SE  -> Draw("P");
  G_Eff_HLT_DE  -> Draw("P");
  G_Eff_HLT_SER -> Draw("P");
  G_Eff_HLT_DER -> Draw("P");

  TLegend *leg =new TLegend(0.11,0.65,0.45,0.89);
  leg->SetBorderSize(2);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(G_Eff_HLT_ELE, "Single OR SingleRelaxed OR Double OR doubleRelaxed Ele","p");
  leg->AddEntry(G_Eff_HLT_SE,  "Single Ele","p");
  leg->AddEntry(G_Eff_HLT_SER, "Single Relaxed Ele","p");
  leg->AddEntry(G_Eff_HLT_DE,  "Double Ele","p");
  leg->AddEntry(G_Eff_HLT_DER, "Double Relaxed Ele","p");
  leg->Draw();

  // se si lavora con i tre fondi, scommentare questa parte 
  // TLatex l;
  // l.DrawLatex(1.,0.19,"1: ww bkg");
  // l.DrawLatex(1.,0.12,"2: wt bkg");
  // l.DrawLatex(1.,0.05,"3: t#bar{t} bkg");
  // mg->GetHistogram()->GetXaxis()->SetTitle("");

  theApp->Run(kFALSE);  
  
} // end;
  
