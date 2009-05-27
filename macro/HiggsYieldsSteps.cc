#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

string preSelCuts[11];
string fullSelCuts[19];

float H_preSel[11];
float Wj_preSel[11];
float ttj_preSel[11];
float Zj_preSel[11];
float WW_preSel[11];
float ZZ4l_preSel[11];
float WZ3l_preSel[11];
float WZincl_preSel[11];

float H_fullSel[19];
float Wj_fullSel[19];
float ttj_fullSel[19];
float Zj_fullSel[19];
float WW_fullSel[19];
float ZZ4l_fullSel[19];
float WZ3l_fullSel[19];
float WZincl_fullSel[19];

float H_eff_preSel[11];
float Wj_eff_preSel[11];
float ttj_eff_preSel[11];
float Zj_eff_preSel[11];
float WW_eff_preSel[11];
float ZZ4l_eff_preSel[11];
float WZ3l_eff_preSel[11];
float WZincl_eff_preSel[11];

float H_eff_fullSel[19];
float Wj_eff_fullSel[19];
float ttj_eff_fullSel[19];
float Zj_eff_fullSel[19];
float WW_eff_fullSel[19];
float ZZ4l_eff_fullSel[19];
float WZ3l_eff_fullSel[19];
float WZincl_eff_fullSel[19];

float H_finaleff_preSel;    
float Wj_finaleff_preSel;
float ttj_finaleff_preSel;
float Zj_finaleff_preSel;
float WW_finaleff_preSel;
float ZZ4l_finaleff_preSel;
float WZ3l_finaleff_preSel;
float WZincl_finaleff_preSel;

float H_finaleff_fullSel;    
float Wj_finaleff_fullSel;
float ttj_finaleff_fullSel;
float Zj_finaleff_fullSel;
float WW_finaleff_fullSel;
float ZZ4l_finaleff_fullSel;
float WZ3l_finaleff_fullSel;
float WZincl_finaleff_fullSel;

float H_finaleff;    
float Wj_finaleff;
float ttj_finaleff;
float Zj_finaleff;
float WW_finaleff;
float ZZ4l_finaleff;
float WZ3l_finaleff;
float WZincl_finaleff;

// lumi of CMST3 ntuples in pb-1    
float Wj_MADGRAPH_lumi  = 243.336;   // 243 * 50 / 50.;  // this is because I have run over 50 / 50 jobs (Fall08 only)
float Zj_MADGRAPH_lumi  = 290.339;    
float ttj_MADGRAPH_lumi = 3381.722;

// here xsec = x-sec * filter_eff (pb)
float H160_xsec    = 1.31 * 1.;
float WW_xsec      = 74.1*1.000*0.1055*0.1055*9;
float ZZ4l_xsec    = 10.5*0.10*0.10;
float ZZ2l2nu_xsec = 10.5*0.10*0.20*2;
float WZ3l_xsec    = 32.4*0.10*0.1055*3;
float WZincl_xsec  = 32.4*1.000;

void computeYields(float lumi, const char* finalstate) {

  TChain *chains_preSel[8];
  TChain *chains_fullSel[8];
  for(int isample=0; isample<8; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }
  
  chains_preSel[0]->Add("results/H160_WW_2l/*Counters.root");       
  chains_preSel[1]->Add("results/WjetsMADGRAPH_Fall08/*Counters.root");       
  chains_preSel[2]->Add("results/ttjetsMadgraph_Fall08/*Counters.root");       
  chains_preSel[3]->Add("results/ZjetsMadgraph_Fall08/*Counters.root");       
  chains_preSel[4]->Add("results/WW_2l/*Counters.root");       
  chains_preSel[5]->Add("results/ZZ_4l/*Counters.root");       
  chains_preSel[6]->Add("results/WZ_3l/*Counters.root");       
  chains_preSel[7]->Add("results/WZ_incl/*Counters.root");       
			
  chains_fullSel[0]->Add("results/H160_WW_2l/*Counters.root");       
  chains_fullSel[1]->Add("results/WjetsMADGRAPH_Fall08/*Counters.root");       
  chains_fullSel[2]->Add("results/ttjetsMadgraph_Fall08/*Counters.root");       
  chains_fullSel[3]->Add("results/ZjetsMadgraph_Fall08/*Counters.root");       
  chains_fullSel[4]->Add("results/WW_2l/*Counters.root");       
  chains_fullSel[5]->Add("results/ZZ_4l/*Counters.root");       
  chains_fullSel[6]->Add("results/WZ_3l/*Counters.root");       
  chains_fullSel[7]->Add("results/WZ_incl/*Counters.root");       

  float nPreSelTot[11][8];
  float nFullSelTot[19][8];

  for(int isample=0; isample<8; isample++) {
    for(int icut=0; icut<11; icut++) { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<19; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preselections
  int nCutsAnaPre  = 0;
  int nCutsAnaFull = 0;
  for(int isample=0; isample<8; isample++) {

    cout << "Processing sample # " << isample << endl;
    
    // List of branches    
    Int_t           nCutsPre;
    Float_t         nSelPre[11];   //[nCuts]
    TBranch        *b_nCutsPre;   //!
    TBranch        *b_nSelPre;    //!
    chains_preSel[isample]->SetBranchAddress("nCuts", &nCutsPre, &b_nCutsPre);
    chains_preSel[isample]->SetBranchAddress("nSel",  nSelPre,   &b_nSelPre);
    
    Int_t           nCutsFull;
    Float_t         nSelFull[19];   //[nCuts]
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
      Long64_t nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      nCutsAnaFull = nCutsFull;
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at preselection level
  for(int icut=0; icut<nCutsAnaPre; icut++) {

    // numbers
    // N = L * x-sec * eff (eff = N_fin / N_ini)
    // weights = 1.0 assigned to MADGRAPH events
    if(nPreSelTot[0][0]>0) { 
      H_preSel[icut]      = lumi * H160_xsec   * nPreSelTot[icut][0]/nPreSelTot[0][0];
      WW_preSel[icut]     = lumi * WW_xsec     * nPreSelTot[icut][4]/nPreSelTot[0][4];
      ZZ4l_preSel[icut]   = lumi * ZZ4l_xsec   * nPreSelTot[icut][5]/nPreSelTot[0][5];
      WZ3l_preSel[icut]   = lumi * WZ3l_xsec   * nPreSelTot[icut][6]/nPreSelTot[0][6];
      WZincl_preSel[icut] = lumi * WZincl_xsec * nPreSelTot[icut][7]/nPreSelTot[0][7];
    }
    Wj_preSel[icut]  = nPreSelTot[icut][1] * lumi/Wj_MADGRAPH_lumi;
    ttj_preSel[icut] = nPreSelTot[icut][2] * lumi/ttj_MADGRAPH_lumi;
    Zj_preSel[icut]  = nPreSelTot[icut][3] * lumi/Zj_MADGRAPH_lumi;
    
    // efficiencies
    if(icut>0 && nPreSelTot[icut-1][0]>0) H_eff_preSel[icut]      = nPreSelTot[icut][0]/nPreSelTot[icut-1][0];
    else H_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][1]>0) Wj_eff_preSel[icut]     = nPreSelTot[icut][1]/nPreSelTot[icut-1][1];
    else Wj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][2]>0) ttj_eff_preSel[icut]    = nPreSelTot[icut][2]/nPreSelTot[icut-1][2];
    else ttj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][3]>0) Zj_eff_preSel[icut]     = nPreSelTot[icut][3]/nPreSelTot[icut-1][3];
    else Zj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][4]>0) WW_eff_preSel[icut]     = nPreSelTot[icut][4]/nPreSelTot[icut-1][4];
    else WW_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][5]>0) ZZ4l_eff_preSel[icut]   = nPreSelTot[icut][5]/nPreSelTot[icut-1][5];
    else ZZ4l_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][6]>0) WZ3l_eff_preSel[icut]   = nPreSelTot[icut][6]/nPreSelTot[icut-1][6];
    else WZ3l_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][7]>0) WZincl_eff_preSel[icut] = nPreSelTot[icut][7]/nPreSelTot[icut-1][7];
    else WZincl_eff_preSel[icut] = 0.0;
  }

  // final efficiency at preselection level
  if(nPreSelTot[0][0]>0) H_finaleff_preSel      = nPreSelTot[nCutsAnaPre-2][0]/nPreSelTot[0][0]; 
  else H_finaleff_preSel = 0.0;
  if(nPreSelTot[0][1]>0) Wj_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][1]/nPreSelTot[0][1]; 
  else Wj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][2]>0) ttj_finaleff_preSel    = nPreSelTot[nCutsAnaPre-2][2]/nPreSelTot[0][2]; 
  else ttj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][3]>0) Zj_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][3]/nPreSelTot[0][3]; 
  else Zj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][4]>0) WW_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][4]/nPreSelTot[0][4]; 
  else WW_finaleff_preSel = 0.0;
  if(nPreSelTot[0][5]>0) ZZ4l_finaleff_preSel   = nPreSelTot[nCutsAnaPre-2][5]/nPreSelTot[0][5]; 
  else ZZ4l_finaleff_preSel = 0.0;
  if(nPreSelTot[0][6]>0) WZ3l_finaleff_preSel   = nPreSelTot[nCutsAnaPre-2][6]/nPreSelTot[0][6]; 
  else WZ3l_finaleff_preSel = 0.0;
  if(nPreSelTot[0][7]>0) WZincl_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][7]/nPreSelTot[0][7]; 
  else WZincl_finaleff_preSel = 0.0;

  
  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) { 
      H_fullSel[icut]      = lumi * H160_xsec   * nFullSelTot[icut][0]/nPreSelTot[0][0];
      WW_fullSel[icut]     = lumi * WW_xsec     * nFullSelTot[icut][4]/nPreSelTot[0][4];
      ZZ4l_fullSel[icut]   = lumi * ZZ4l_xsec   * nFullSelTot[icut][5]/nPreSelTot[0][5];
      WZ3l_fullSel[icut]   = lumi * WZ3l_xsec   * nFullSelTot[icut][6]/nPreSelTot[0][6];
      WZincl_fullSel[icut] = lumi * WZincl_xsec * nFullSelTot[icut][7]/nPreSelTot[0][7];
    }
    Wj_fullSel[icut]  = nFullSelTot[icut][1] * lumi/Wj_MADGRAPH_lumi;
    ttj_fullSel[icut] = nFullSelTot[icut][2] * lumi/ttj_MADGRAPH_lumi;
    Zj_fullSel[icut]  = nFullSelTot[icut][3] * lumi/Zj_MADGRAPH_lumi;


    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) H_eff_fullSel[icut]   = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else H_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][1]>0) Wj_eff_fullSel[icut]  = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else Wj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][2]>0) ttj_eff_fullSel[icut] = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else ttj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][3]>0) Zj_eff_fullSel[icut]  = nFullSelTot[icut][3]/nFullSelTot[icut-1][3];
    else Zj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][4]>0) WW_eff_fullSel[icut]  = nFullSelTot[icut][4]/nFullSelTot[icut-1][4];
    else WW_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][5]>0) ZZ4l_eff_fullSel[icut]  = nFullSelTot[icut][5]/nFullSelTot[icut-1][5];
    else ZZ4l_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][6]>0) WZ3l_eff_fullSel[icut]  = nFullSelTot[icut][6]/nFullSelTot[icut-1][6];
    else WZ3l_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][7]>0) WZincl_eff_fullSel[icut]  = nFullSelTot[icut][7]/nFullSelTot[icut-1][7];
    else WZincl_eff_fullSel[icut] = 0.0;
    
    if(icut==0) { 
      H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nPreSelTot[nCutsAnaPre-2][0];    
      Wj_eff_fullSel[icut]     = nFullSelTot[icut][1]/nPreSelTot[nCutsAnaPre-2][1];    
      ttj_eff_fullSel[icut]    = nFullSelTot[icut][2]/nPreSelTot[nCutsAnaPre-2][2];    
      Zj_eff_fullSel[icut]     = nFullSelTot[icut][3]/nPreSelTot[nCutsAnaPre-2][3];    
      WW_eff_fullSel[icut]     = nFullSelTot[icut][4]/nPreSelTot[nCutsAnaPre-2][4];    
      ZZ4l_eff_fullSel[icut]   = nFullSelTot[icut][5]/nPreSelTot[nCutsAnaPre-2][5];    
      WZ3l_eff_fullSel[icut]   = nFullSelTot[icut][6]/nPreSelTot[nCutsAnaPre-2][6];    
      WZincl_eff_fullSel[icut] = nFullSelTot[icut][7]/nPreSelTot[nCutsAnaPre-2][7];    
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0) H_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0]; 
  else H_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) Wj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1]; 
  else Wj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][2]>0) ttj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][2]/nFullSelTot[0][2]; 
  else ttj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][3]>0) Zj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][3]/nFullSelTot[0][3]; 
  else Zj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][4]>0) WW_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][4]/nFullSelTot[0][4]; 
  else WW_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][5]>0) ZZ4l_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][5]/nFullSelTot[0][5]; 
  else ZZ4l_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][6]>0) WZ3l_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][6]/nFullSelTot[0][6]; 
  else WZ3l_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][7]>0) WZincl_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][7]/nFullSelTot[0][7]; 
  else WZincl_finaleff_fullSel = 0.0;
  
  // final efficiency combining pre and full selections
  if(nFullSelTot[0][0]>0) H_finaleff      = nFullSelTot[nCutsAnaFull-4][0]/nPreSelTot[0][0]; 
  else H_finaleff = 0.0;
  if(nFullSelTot[0][1]>0) Wj_finaleff     = nFullSelTot[nCutsAnaFull-4][1]/nPreSelTot[0][1]; 
  else Wj_finaleff = 0.0;
  if(nFullSelTot[0][2]>0) ttj_finaleff    = nFullSelTot[nCutsAnaFull-4][2]/nPreSelTot[0][2]; 
  else ttj_finaleff = 0.0;
  if(nFullSelTot[0][3]>0) Zj_finaleff     = nFullSelTot[nCutsAnaFull-4][3]/nPreSelTot[0][3]; 
  else Zj_finaleff = 0.0;
  if(nFullSelTot[0][4]>0) WW_finaleff     = nFullSelTot[nCutsAnaFull-4][4]/nPreSelTot[0][4]; 
  else WW_finaleff = 0.0;
  if(nFullSelTot[0][5]>0) ZZ4l_finaleff   = nFullSelTot[nCutsAnaFull-4][5]/nPreSelTot[0][5]; 
  else ZZ4l_finaleff = 0.0;
  if(nFullSelTot[0][6]>0) WZ3l_finaleff   = nFullSelTot[nCutsAnaFull-4][6]/nPreSelTot[0][6]; 
  else WZ3l_finaleff = 0.0;
  if(nFullSelTot[0][7]>0) WZincl_finaleff = nFullSelTot[nCutsAnaFull-4][7]/nPreSelTot[0][7]; 
  else WZincl_finaleff = 0.0;
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
  preSelCuts[9]="finalOURPreselection";
  preSelCuts[10]="preselection";

  fullSelCuts[0]="this channel preselected";
  fullSelCuts[1]="hardLeptonThreshold";
  fullSelCuts[2]="slowLeptonThreshold";
  fullSelCuts[3]="dileptonInvMassMin";
  fullSelCuts[4]="classDepEleId";
  fullSelCuts[5]="trackerIso";
  fullSelCuts[6]="hcalIso";
  fullSelCuts[7]="ecalIso";
  fullSelCuts[8]="globalIso";
  fullSelCuts[9]="MET";
  fullSelCuts[10]="maxPtLepton";
  fullSelCuts[11]="minPtLepton";
  fullSelCuts[12]="dileptonInvMassMax";
  fullSelCuts[13]="detaLeptons";
  fullSelCuts[14]="deltaPhi";
  fullSelCuts[15]="final";
  fullSelCuts[16]="0 jets";
  fullSelCuts[17]="1 jets";
  fullSelCuts[18]=">1 jets";
}


void printLatex(float lumi, const char* finalstate) {

  setupCuts();

  if(strcmp(finalstate,"EE") && strcmp(finalstate,"EM") && strcmp(finalstate,"MM")) {
    cout << "ERROR! finalstate must be one among EE/EM/MM. Exiting..." << endl;
    return;
  } else {
    cout << " === NOW COMPUTING YIELDS FOR FINAL STATE: " << finalstate << " ===" << endl;
  }

  computeYields(lumi,finalstate);
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(1);

  // W+jets events
  textfile << "\\begin{sidewaystable}[p]" << endl
	   << "\\begin{small}" << endl
	   << "\\begin{center}" << endl
	   << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}" << endl
	   << "\\hline" << endl
	   << "selection & H (160 GeV) & W$(e \\nu)$+jets & $t\\bar{t}$ & Z$(e \\nu)$+jets & WW & ZZ$->$4l & WZ$->$3l & WZ$->$incl \t\\\\" << endl
	   << "\\hline" << endl; 
  
  textfile << "\\hline"        << endl;
  textfile << "\\hline"        << endl;

  for(int icut=0; icut<9; icut++) {
    textfile << preSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << H_preSel[icut]       << " (" << 100. * H_eff_preSel[icut]      << "\\%)" << "\t&\t"
	     << Wj_preSel[icut]      << " (" << 100. * Wj_eff_preSel[icut]     << "\\%)" << "\t&\t"
	     << ttj_preSel[icut]     << " (" << 100. * ttj_eff_preSel[icut]    << "\\%)" << "\t&\t"
	     << Zj_preSel[icut]      << " (" << 100. * Zj_eff_preSel[icut]     << "\\%)" << "\t&\t"
	     << WW_preSel[icut]      << " (" << 100. * WW_eff_preSel[icut]     << "\\%)" << "\t&\t"
	     << ZZ4l_preSel[icut]    << " (" << 100. * ZZ4l_eff_preSel[icut]   << "\\%)" << "\t&\t"
	     << WZ3l_preSel[icut]    << " (" << 100. * WZ3l_eff_preSel[icut]   << "\\%)" << "\t&\t"
	     << WZincl_preSel[icut]  << " (" << 100. * WZincl_eff_preSel[icut] << "\\%)" 
	     << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total preselection " << "\t&\t"
	   << H_preSel[8]      << " (" << 100. * H_finaleff_preSel  << "\\%)"  << "\t&\t"
	   << Wj_preSel[8]     << " (" << 100. * Wj_finaleff_preSel << "\\%)"  << "\t&\t"
	   << ttj_preSel[8]    << " (" << 100. * ttj_finaleff_preSel << "\\%)" << "\t&\t"
	   << Zj_preSel[8]     << " (" << 100. * Zj_finaleff_preSel << "\\%)"  << "\t&\t"
	   << WW_preSel[8]     << " (" << 100. * WW_finaleff_preSel << "\\%)"  << "\t&\t"
	   << ZZ4l_preSel[8]   << " (" << 100. * ZZ4l_finaleff_preSel << "\\%)"  << "\t&\t"
	   << WZ3l_preSel[8]   << " (" << 100. * WZ3l_finaleff_preSel << "\\%)"  << "\t&\t"
	   << WZincl_preSel[8] << " (" << 100. * WZincl_finaleff_preSel << "\\%)"  
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;
  
  for(int icut=0; icut<16; icut++) {
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << H_fullSel[icut]      << " (" << 100. * H_eff_fullSel[icut]   << "\\%)" << "\t&\t"
	     << Wj_fullSel[icut]     << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << ttj_fullSel[icut]    << " (" << 100. * ttj_eff_fullSel[icut] << "\\%)" << "\t&\t"
	     << Zj_fullSel[icut]     << " (" << 100. * Zj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << WW_fullSel[icut]     << " (" << 100. * WW_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << ZZ4l_fullSel[icut]   << " (" << 100. * ZZ4l_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << WZ3l_fullSel[icut]   << " (" << 100. * WZ3l_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << WZincl_fullSel[icut] << " (" << 100. * WZincl_eff_fullSel[icut]  << "\\%)" 
	     << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;

  textfile << "total fullselection " << "\t&\t"
	   << H_fullSel[15]      << " (" << 100. * H_finaleff_fullSel  << "\\%)"     << "\t&\t"
	   << Wj_fullSel[15]     << " (" << 100. * Wj_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << ttj_fullSel[15]    << " (" << 100. * ttj_finaleff_fullSel << "\\%)"    << "\t&\t"
	   << Zj_fullSel[15]     << " (" << 100. * Zj_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << WW_fullSel[15]     << " (" << 100. * WW_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << ZZ4l_fullSel[15]   << " (" << 100. * ZZ4l_finaleff_fullSel << "\\%)"   << "\t&\t"
	   << WZ3l_fullSel[15]   << " (" << 100. * WZ3l_finaleff_fullSel << "\\%)"   << "\t&\t"
	   << WZincl_fullSel[15] << " (" << 100. * WZincl_finaleff_fullSel << "\\%)" 
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;

  textfile << "total " << "\t&\t"
	   << H_fullSel[15]      << " (" << 100. * H_finaleff   << "\\%)"    << "\t&\t"
	   << Wj_fullSel[15]     << " (" << 100. * Wj_finaleff  << "\\%)"    << "\t&\t"
	   << ttj_fullSel[15]    << " (" << 100. * ttj_finaleff << "\\%)"    << "\t&\t"
	   << Zj_fullSel[15]     << " (" << 100. * Zj_finaleff << "\\%)"     << "\t&\t"
	   << WW_fullSel[15]     << " (" << 100. * WW_finaleff << "\\%)"     << "\t&\t"
	   << ZZ4l_fullSel[15]   << " (" << 100. * ZZ4l_finaleff << "\\%)"   << "\t&\t"
	   << WZ3l_fullSel[15]   << " (" << 100. * WZ3l_finaleff << "\\%)"   << "\t&\t"
	   << WZincl_fullSel[15] << " (" << 100. * WZincl_finaleff << "\\%)" 
	   << "\t\\\\" << endl;

  textfile << "0 jets bin " << "\t&\t"
	   << H_fullSel[16]      << "\t&\t"
	   << Wj_fullSel[16]     << "\t&\t"
	   << ttj_fullSel[16]    << "\t&\t"
	   << Zj_fullSel[16]     << "\t&\t"
	   << WW_fullSel[16]     << "\t&\t"
	   << ZZ4l_fullSel[16]   << "\t&\t"
	   << WZ3l_fullSel[16]   << "\t&\t"
	   << WZincl_fullSel[16] << "\t\\\\" << endl;
	   
  textfile << "1 jets bin " << "\t&\t"
	   << H_fullSel[17]      << "\t&\t"
	   << Wj_fullSel[17]     << "\t&\t"
	   << ttj_fullSel[17]    << "\t&\t"
	   << Zj_fullSel[17]     << "\t&\t"
	   << WW_fullSel[17]     << "\t&\t"
	   << ZZ4l_fullSel[17]   << "\t&\t"
	   << WZ3l_fullSel[17]   << "\t&\t"
	   << WZincl_fullSel[17] << "\t\\\\" << endl;

  textfile << "$>1$ jets bin " << "\t&\t"
	   << H_fullSel[18]      << "\t&\t"
	   << Wj_fullSel[18]     << "\t&\t"
	   << ttj_fullSel[18]    << "\t&\t"
	   << Zj_fullSel[18]     << "\t&\t"
	   << WW_fullSel[18]     << "\t&\t"
	   << ZZ4l_fullSel[18]   << "\t&\t"
	   << WZ3l_fullSel[18]   << "\t&\t"
	   << WZincl_fullSel[18] << "\t\\\\" << endl;
	   
  
  textfile << "\\hline" << endl
	   << "\\end{tabular}" << endl
	   << "\\end{center}" << endl
	   << "\\end{small}" << endl
           << "\\caption{Breakdown of signal and backgrounds events in "
           << lumi << " $pb^{-1}$ for " << finalstate << " final state.} " << endl 
	   << "\\end{sidewaystable}" << endl;

}

void printLatex(float lumi) {
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::trunc);
  textfile.precision(1);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;

  textfile.close();

  printLatex(lumi, "EE");
  printLatex(lumi, "MM");
  printLatex(lumi, "EM");

  textfile.open(namefile, ios_base::app);
  textfile << "\\end{document}" << endl;
  textfile.close();

}

