#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

string preSelCuts[9];
string fullSelCuts[16];

float H_preSel[9];
float Wj_preSel[9];
float ttj_preSel[9];

float H_fullSel[16];
float Wj_fullSel[16];
float ttj_fullSel[16];

float H_eff_preSel[9];
float Wj_eff_preSel[9];
float ttj_eff_preSel[9];

float H_eff_fullSel[16];
float Wj_eff_fullSel[16];
float ttj_eff_fullSel[16];

float H_finaleff_preSel;    
float Wj_finaleff_preSel;
float ttj_finaleff_preSel;

float H_finaleff_fullSel;    
float Wj_finaleff_fullSel;
float ttj_finaleff_fullSel;

float H_finaleff;    
float Wj_finaleff;
float ttj_finaleff;


// lumi of CMST3 ntuples in pb-1    
float Wj_MADGRAPH_lumi  = 206.;   // 250 * 41 / 50.;  // this is because I have run over 41 / 50 jobs (Fall08 only)
float ttj_MADGRAPH_lumi = 3010.;

// here xsec = x-sec * filter_eff (pb)
float H160_xsec = 1.31 * 1.;

void computeYields(float lumi=100.) {

  TChain *chains_preSel[3];
  TChain *chains_fullSel[3];
  for(int isample=0; isample<3; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    chains_fullSel[isample] = new TChain("FULL_SELECTION_EVENT_COUNTER_EE");
  }

  chains_preSel[0]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/signal/*Counters.root");       
  chains_preSel[1]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/WjetsMADGRAPH_Fall08/*Counters.root");   
  chains_preSel[2]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/ttjetsMADGRAPH_Fall08/*Counters.root");    

  chains_fullSel[0]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/signal/*Counters.root");  
  chains_fullSel[1]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/WjetsMADGRAPH_Fall08/*Counters.root");  
  chains_fullSel[2]->Add("/afs/cern.ch/user/c/crovelli/scratch0/HWW_2009_2X/HiggsAnalysisTools/analysisResults/newOptim_SiPreselId/ttjetsMADGRAPH_Fall08/*Counters.root");    


  float nPreSelTot[9][3];
  float nFullSelTot[16][3];

  for(int isample=0; isample<3; isample++) {
    for(int icut=0; icut<9; icut++)  { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<16; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preselections
  int nCutsAnaPre  = 0;
  int nCutsAnaFull = 0;
  for(int isample=0; isample<3; isample++) {

    cout << "Processing sample # " << isample << endl;
    
    // List of branches    
    Int_t           nCutsPre;
    Float_t         nSelPre[9];   //[nCuts]
    TBranch        *b_nCutsPre;   //!
    TBranch        *b_nSelPre;    //!
    chains_preSel[isample]->SetBranchAddress("nCuts", &nCutsPre, &b_nCutsPre);
    chains_preSel[isample]->SetBranchAddress("nSel",  nSelPre,   &b_nSelPre);

    Int_t           nCutsFull;
    Float_t         nSelFull[16];   //[nCuts]
    TBranch        *b_nCutsFull;   //!
    TBranch        *b_nSelFull;    //!
    chains_fullSel[isample]->SetBranchAddress("nCuts", &nCutsFull, &b_nCutsFull);
    chains_fullSel[isample]->SetBranchAddress("nSel",  nSelFull,   &b_nSelFull);
    
    Long64_t nentriesPre  = chains_preSel[isample]->GetEntries();
    Long64_t nentriesFull = chains_fullSel[isample]->GetEntries();
    
    // pre-selection
    for (Long64_t jentry=0; jentry<nentriesPre;jentry++) {
      Long64_t nb = chains_preSel[isample]->GetEntry(jentry);   
      nCutsAnaPre = nCutsPre;      
      for(int icut=0; icut<nCutsPre; icut++) nPreSelTot[icut][isample] += nSelPre[icut];	
    }
    
    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb = chains_fullSel[isample]->GetEntry(jentry);   
      nCutsAnaFull = nCutsFull;      
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at preselection level
  for(int icut=0; icut<nCutsAnaPre; icut++) {

    // numbers
    // N = L * x-sec * eff (eff = N_fin / N_ini)
    // weights = 1.0 assigned to MADGRAPH events
    if(nPreSelTot[0][0]>0) H_preSel[icut] = lumi * H160_xsec * nPreSelTot[icut][0]/nPreSelTot[0][0];
    Wj_preSel[icut]  = nPreSelTot[icut][1] * lumi/Wj_MADGRAPH_lumi;
    ttj_preSel[icut] = nPreSelTot[icut][2] * lumi/ttj_MADGRAPH_lumi;
    
    // efficiencies
    if(icut>0 && nPreSelTot[icut-1][0]>0) H_eff_preSel[icut]   = nPreSelTot[icut][0]/nPreSelTot[icut-1][0];
    else H_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][1]>0) Wj_eff_preSel[icut]  = nPreSelTot[icut][1]/nPreSelTot[icut-1][1];
    else Wj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][2]>0) ttj_eff_preSel[icut] = nPreSelTot[icut][2]/nPreSelTot[icut-1][2];
    else ttj_eff_preSel[icut] = 0.0;
  }

  // final efficiency at preselection level
  if(nPreSelTot[0][0]>0) H_finaleff_preSel   = nPreSelTot[nCutsAnaPre-2][0]/nPreSelTot[0][0]; 
  else H_finaleff_preSel = 0.0;
  if(nPreSelTot[0][1]>0) Wj_finaleff_preSel  = nPreSelTot[nCutsAnaPre-2][1]/nPreSelTot[0][1]; 
  else Wj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][2]>0) ttj_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][2]/nPreSelTot[0][2]; 
  else ttj_finaleff_preSel = 0.0;
  


  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) H_fullSel[icut] = lumi * H160_xsec * nFullSelTot[icut][0]/nPreSelTot[0][0];
    Wj_fullSel[icut]  = nFullSelTot[icut][1] * lumi/Wj_MADGRAPH_lumi;
    ttj_fullSel[icut] = nFullSelTot[icut][2] * lumi/ttj_MADGRAPH_lumi;

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) H_eff_fullSel[icut]   = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else H_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][1]>0) Wj_eff_fullSel[icut]  = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else Wj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][2]>0) ttj_eff_fullSel[icut] = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else ttj_eff_fullSel[icut] = 0.0;
    
    if(icut==0) { 
      H_eff_fullSel[icut]   = nFullSelTot[icut][0]/nPreSelTot[nCutsAnaPre-2][0];    
      Wj_eff_fullSel[icut]  = nFullSelTot[icut][1]/nPreSelTot[nCutsAnaPre-2][1];    
      ttj_eff_fullSel[icut] = nFullSelTot[icut][2]/nPreSelTot[nCutsAnaPre-2][2];    
    }
  }

  // final efficiency after full selections
  if(nFullSelTot[0][0]>0) H_finaleff_fullSel = nFullSelTot[nCutsAnaFull-1][0]/nFullSelTot[0][0]; 
  else H_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) Wj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-1][1]/nFullSelTot[0][1]; 
  else Wj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][2]>0) ttj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-1][2]/nFullSelTot[0][2]; 
  else ttj_finaleff_fullSel = 0.0;
  

  // final efficiency combining pre and full selections
  if(nFullSelTot[0][0]>0) H_finaleff   = nFullSelTot[nCutsAnaFull-1][0]/nPreSelTot[0][0]; 
  else H_finaleff = 0.0;
  if(nFullSelTot[0][1]>0) Wj_finaleff  = nFullSelTot[nCutsAnaFull-1][1]/nPreSelTot[0][1]; 
  else Wj_finaleff = 0.0;
  if(nFullSelTot[0][2]>0) ttj_finaleff = nFullSelTot[nCutsAnaFull-1][2]/nPreSelTot[0][2]; 
  else ttj_finaleff = 0.0;
}

void setupCuts() {
  
  preSelCuts[0]="event";
  preSelCuts[1]="MCtruth";
  preSelCuts[2]="trigger";
  preSelCuts[3]="nRecoLeptons";
  preSelCuts[4]="twoGoodRec";
  preSelCuts[5]="hardLeptonThreshold";
  preSelCuts[6]="slowLeptonThreshold";
  preSelCuts[7]="METpreselection";
  preSelCuts[8]="dileptonInvMassMin";

  fullSelCuts[0]="this channel preselected";
  fullSelCuts[1]="hardLeptonThreshold";
  fullSelCuts[2]="slowLeptonThreshold";
  fullSelCuts[3]="dileptonInvMassMin";
  fullSelCuts[4]="classDepEleId";
  fullSelCuts[5]="trackerIso";
  fullSelCuts[6]="hcalIso";
  fullSelCuts[7]="ecalIso";
  fullSelCuts[8]="jetVeto";
  fullSelCuts[9]="MET";
  fullSelCuts[10]="maxPtLepton";
  fullSelCuts[11]="minPtLepton";
  fullSelCuts[12]="dileptonInvMassMax";
  fullSelCuts[13]="detaLeptons";
  fullSelCuts[14]="deltaPhi";
}


void printLatex(float lumi) {

  setupCuts();

  computeYields(lumi);
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::trunc);
  textfile.precision(0);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;


  // W+jets events
  textfile << "\\begin{sidewaystable}[p]" << endl
	   << "\\begin{small}" << endl
	   << "\\begin{center}" << endl
	   << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl
	   << "\\hline" << endl
	   << "selection & H (160 GeV) & W$(e \\nu)$+jets & $t\\bar{t}$ \t\\\\" << endl
	   << "\\hline" << endl; 

  textfile << "\\hline"        << endl;
  textfile << "\\hline"        << endl;

  for(int icut=0; icut<9; icut++) {
    textfile << preSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << H_preSel[icut]   << " (" << 100. * H_eff_preSel[icut]   << "\\%)" << "\t&\t"
	     << Wj_preSel[icut]  << " (" << 100. * Wj_eff_preSel[icut]  << "\\%)" << "\t&\t"
	     << ttj_preSel[icut] << " (" << 100. * ttj_eff_preSel[icut] << "\\%)" 
	     << "\t\\\\" << endl;
  }

  textfile << "\\hline" << endl;

  textfile << "total preselection " << "\t&\t"
	   << H_preSel[8]   << " (" << 100. * H_finaleff_preSel  << "\\%)" << "\t&\t"
	   << Wj_preSel[8]  << " (" << 100. * Wj_finaleff_preSel << "\\%)" << "\t&\t"
	   << ttj_preSel[8] << " (" << 100. * ttj_finaleff_preSel << "\\%)" 
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;

  for(int icut=0; icut<16; icut++) {
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << H_fullSel[icut]   << " (" << 100. * H_eff_fullSel[icut]   << "\\%)" << "\t&\t"
	     << Wj_fullSel[icut]  << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << ttj_fullSel[icut] << " (" << 100. * ttj_eff_fullSel[icut] << "\\%)" 
	     << "\t\\\\" << endl;
  }

  textfile << "\\hline" << endl;

  textfile << "total fullselection " << "\t&\t"
	   << H_fullSel[14]   << " (" << 100. * H_finaleff_fullSel  << "\\%)" << "\t&\t"
	   << Wj_fullSel[14]  << " (" << 100. * Wj_finaleff_fullSel << "\\%)" << "\t&\t"
	   << ttj_fullSel[14] << " (" << 100. * ttj_finaleff_fullSel << "\\%)" 
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;

  textfile << "total " << "\t&\t"
	   << H_fullSel[14]   << " (" << 100. * H_finaleff   << "\\%)" << "\t&\t"
	   << Wj_fullSel[14]  << " (" << 100. * Wj_finaleff  << "\\%)" << "\t&\t"
	   << ttj_fullSel[14] << " (" << 100. * ttj_finaleff << "\\%)" 
	   << "\t\\\\" << endl;
  
  textfile << "\\hline" << endl
	   << "\\end{tabular}" << endl
	   << "\\end{center}" << endl
	   << "\\end{small}" << endl
	   << "\\end{sidewaystable}" << endl;

  textfile << "\\end{document}" << endl;

}
