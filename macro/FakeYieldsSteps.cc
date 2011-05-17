#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

string fullSelCuts[25];

int UseCuts[25];

float Wj_fullSel[25];
float WjF_fullSel[25];
float data1_fullSel[25];
float data2_fullSel[25];
float data3_fullSel[25];

float Wj_eff_fullSel[25];
float WjF_eff_fullSel[25];

float Wj_finaleff_fullSel;
float WjF_finaleff_fullSel;

// xsections
float Wjets_xsec = 31314.;  // madgraph 

string sampleNames[5];

void computeYields(float lumi) {

  TChain *chains_fullSel[5];  
  for(int isample=0; isample<5; isample++) {
    char fullsel_treename[500];
    // sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_EE");
    sprintf(fullsel_treename,"FULL_SELECTION_ERRORS_EE");
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // samples name
  sampleNames[0] = "Wjets - nominal";
  sampleNames[1] = "Wjets - from fake";
  sampleNames[2] = "data, pTjet = 30";
  sampleNames[3] = "data, pTjet = 15";
  sampleNames[4] = "data, pTjet = 50";

  // mc
  chains_fullSel[0]->Add("results_nominalWjets/*Counters.root");
  chains_fullSel[1]->Add("results_wjetsFromFake/*Counters.root");
  // data
  chains_fullSel[2]->Add("results_data_nominal/*Counters.root");
  chains_fullSel[3]->Add("results_data_syst15/*Counters.root");
  chains_fullSel[4]->Add("results_data_syst50/*Counters.root");

  // initialization
  float nFullSelTot[25][5];
  for(int isample=0; isample<5; isample++) {
    for(int icut=0; icut<25; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }
  
  // full selection
  int nCutsAnaFull = 25;
  for(int isample=0; isample<5; isample++) {

    // List of branches    
    Int_t     nCutsFull;
    Float_t   nSelFull[25];   //[nCuts]
    TBranch   *b_nCutsFull;   //!
    TBranch   *b_nSelFull;    //!
    chains_fullSel[isample]->SetBranchAddress("nCuts", &nCutsFull, &b_nCutsFull);
    chains_fullSel[isample]->SetBranchAddress("nSel",  nSelFull,   &b_nSelFull);
    
    Long64_t nentriesFull = chains_fullSel[isample]->GetEntries();
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2;
      nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }
  
  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) Wj_fullSel[icut]   = lumi * Wjets_xsec * nFullSelTot[icut][0]/nFullSelTot[0][0];
    if(nFullSelTot[0][1]>0) WjF_fullSel[icut]  = lumi * Wjets_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1];
    data1_fullSel[icut] = nFullSelTot[icut][2];
    data2_fullSel[icut] = nFullSelTot[icut][3];
    data3_fullSel[icut] = nFullSelTot[icut][4];

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) Wj_eff_fullSel[icut]  = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else Wj_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][1]>0) WjF_eff_fullSel[icut] = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else WjF_eff_fullSel[icut] = 0.0;

    if(icut==0) { 
      Wj_eff_fullSel[icut]  = nFullSelTot[icut][0]/nFullSelTot[0][0];
      WjF_eff_fullSel[icut] = nFullSelTot[icut][1]/nFullSelTot[0][1];
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0) Wj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0];
  else Wj_finaleff_fullSel = 0.0;
  
  if(nFullSelTot[0][1]>0) WjF_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1];
  else WjF_finaleff_fullSel = 0.0;

  cout << "\n\nPROCESSED EVENTS:" << endl;
  for(int i=0; i<5; i++) cout << sampleNames[i] << "\t" << nFullSelTot[0][i] << endl;  
}

void setupCuts() {
  
  for(int i=0; i<25; i++) UseCuts[i] = 1;
  
  fullSelCuts[0]="event";
  fullSelCuts[1]="MCtruth";
  fullSelCuts[2]="trigger";
  fullSelCuts[3]="channel preSel.";
  fullSelCuts[4]="e/$\\mu$ ID";
  fullSelCuts[5]="e/$\\mu$ isolation";
  fullSelCuts[6]="conv. rej.";
  fullSelCuts[7]="e/$\\mu$ d0";
  fullSelCuts[8]="extra lepton veto";
  fullSelCuts[9]="$MET>20$ GeV";
  fullSelCuts[10]="$m_{ll}$";
  fullSelCuts[11]="$|m_{ll}-m_Z|>15$ GeV";
  fullSelCuts[12]="proj. MET";
  fullSelCuts[13]="MET/$p_T^{ll}$";
  fullSelCuts[14]="jet veto";
  fullSelCuts[15]="$\\mu^{soft}$ veto";
  fullSelCuts[16]="anti b-tag";
  fullSelCuts[17]="$m_{ll}2$";
  fullSelCuts[18]="sel $p_T^{max}$";
  fullSelCuts[19]="sel $p_T^{min}$";
  fullSelCuts[20]="$\\gamma^*M_{R}^{*}$";
  fullSelCuts[21]="$\\Delta \\phi$";
  fullSelCuts[22]="final";
  fullSelCuts[23]="1 jets";
  fullSelCuts[24]="$>1$ jets";
}


void printLatex(float lumi) {

  setupCuts();

  cout << " === NOW COMPUTING YIELDS ===" << endl;
  
  computeYields(lumi);
  
  
  /// ==============  print detailed breakdown  ================== ///
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;
  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{tiny}" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  textfile << "selection & data 30 & data 15 & data 50 & nominal W$(l\\nu)$+jets & fake W$(l\\nu)$+jets \\\\" << endl;
  textfile << "\\hline" << endl; 
  textfile << "\\hline" << endl;
  textfile << "\\hline" << endl;
  
  for(int icut=0; icut<25; icut++) {

    if(!UseCuts[icut]) continue;
    
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
             << data1_fullSel[icut]  << "\t&\t"
             << data2_fullSel[icut]  << "\t&\t"
             << data3_fullSel[icut]  << "\t&\t"
             << Wj_fullSel[icut]    << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << WjF_fullSel[icut]   << " (" << 100. * WjF_eff_fullSel[icut] << "\\%)" << "\t";
    textfile << "\t\\\\" << endl;
  }
    
  textfile << "\\hline" << endl;
  
  textfile << "total fullselection " << "\t&\t"
           << data1_fullSel[22]    << "\t&\t"
           << data2_fullSel[22]    << "\t&\t"
           << data3_fullSel[22]    << "\t&\t"
           << Wj_fullSel[22]       << " (" << 100. * Wj_finaleff_fullSel     << "\\%)" << "\t&\t"
           << WjF_fullSel[22]      << " (" << 100. * WjF_finaleff_fullSel    << "\\%)" << "\t";
  textfile << "\t\\\\"  << endl;
  textfile << "\\hline" << endl;
  textfile << "\t\\\\"  << endl;
    
  textfile << "\\hline"        << endl
           << "\\end{tabular}" << endl
           << "\\end{center}"  << endl
           << "\\end{tiny}"    << endl
           << "\\end{sidewaystable}" << endl;
  textfile << "\\end{document}" << endl;
  textfile.close();
}

