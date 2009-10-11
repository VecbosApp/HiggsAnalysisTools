#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

string preSelCuts[10];
string fullSelCuts[20];

float Wenu_preSel[10];
float ttbar_preSel[10];
float Zee_preSel[10];
float QCD_preSel[10];

float Wenu_fullSel[20];
float ttbar_fullSel[20];
float Zee_fullSel[20];
float QCD_fullSel[20];

float Wenu_eff_preSel[10];
float ttbar_eff_preSel[10];
float Zee_eff_preSel[10];
float QCD_eff_preSel[10];

float Wenu_eff_fullSel[20];
float ttbar_eff_fullSel[20];
float Zee_eff_fullSel[20];
float QCD_eff_fullSel[20];

float Wenu_finaleff_preSel;
float ttbar_finaleff_preSel;
float Zee_finaleff_preSel;
float WW_finaleff_preSel;
float QCD_finaleff_preSel;

float Wenu_finaleff_fullSel;
float ttbar_finaleff_fullSel;
float Zee_finaleff_fullSel;
float QCD_finaleff_fullSel;

float Wenu_finaleff;
float ttbar_finaleff;
float Zee_finaleff;
float QCD_finaleff;

float lumiWeight = 10.0; // the luminosity used to evaluate the event weights

// here weight = cross-section * lumi (10pb-1) * prescale / initial number of events
// http://ceballos.web.cern.ch/ceballos/hwwlnln/hww_SDs_Oct09.txt
float Wenu_weight  = 0.81550;
float ttbar_weight = 0.70788;
float Zee_weight   = 0.72474;

float QCD_em20to30_weight = 0.95129;
float QCD_em30to80_weight = 1.22521;
float QCD_em80to170_weight = 0.49742;
float QCD_bctoe20to30_weight = 0.80542;
float QCD_bctoe30to80_weight = 1.17929;
float QCD_bctoe80to170_weight = 0.21964;

void computeYields(float lumi, const char* finalstate) {

  TChain *chains_preSel[9];
  TChain *chains_fullSel[9];
  for(int isample=0; isample<9; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }
  
  chains_preSel[0]->Add("results/Wenu/*Counters.root");       
  chains_preSel[1]->Add("results/TTbarEle/*Counters.root");       
  chains_preSel[2]->Add("results/Zee/*Counters.root");       
  chains_preSel[3]->Add("results/QCD_EMEnriched_Pt20to30/*Counters.root");       
  chains_preSel[4]->Add("results/QCD_EMEnriched_Pt30to80/*Counters.root");       
  chains_preSel[5]->Add("results/QCD_EMEnriched_Pt80to170/*Counters.root");       
  chains_preSel[6]->Add("results/QCD_BCtoE_Pt20to30/*Counters.root");       
  chains_preSel[7]->Add("results/QCD_BCtoE_Pt30to80/*Counters.root");       
  chains_preSel[8]->Add("results/QCD_BCtoE_Pt80to170/*Counters.root");       

  chains_fullSel[0]->Add("results/Wenu/*Counters.root");       
  chains_fullSel[1]->Add("results/TTbarEle/*Counters.root");       
  chains_fullSel[2]->Add("results/Zee/*Counters.root");       
  chains_fullSel[3]->Add("results/QCD_EMEnriched_Pt20to30/*Counters.root");       
  chains_fullSel[4]->Add("results/QCD_EMEnriched_Pt30to80/*Counters.root");       
  chains_fullSel[5]->Add("results/QCD_EMEnriched_Pt80to170/*Counters.root");       
  chains_fullSel[6]->Add("results/QCD_BCtoE_Pt20to30/*Counters.root");       
  chains_fullSel[7]->Add("results/QCD_BCtoE_Pt30to80/*Counters.root");       
  chains_fullSel[8]->Add("results/QCD_BCtoE_Pt80to170/*Counters.root");       
			
  float nPreSelTot[10][9];
  float nFullSelTot[20][9];

  for(int isample=0; isample<9; isample++) {
    for(int icut=0; icut<10; icut++) { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<20; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preselections
  int nCutsAnaPre  = 10;
  int nCutsAnaFull = 20;
  for(int isample=0; isample<9; isample++) {

    cout << "Processing sample # " << isample << endl;
    
    // List of branches    
    Int_t           nCutsPre;
    Float_t         nSelPre[10];   //[nCuts]
    TBranch        *b_nCutsPre;   //!
    TBranch        *b_nSelPre;    //!
    chains_preSel[isample]->SetBranchAddress("nCuts", &nCutsPre, &b_nCutsPre);
    chains_preSel[isample]->SetBranchAddress("nSel",  nSelPre,   &b_nSelPre);
    
    Int_t           nCutsFull;
    Float_t         nSelFull[20];   //[nCuts]
    TBranch        *b_nCutsFull;   //!
    TBranch        *b_nSelFull;    //!
    chains_fullSel[isample]->SetBranchAddress("nCuts", &nCutsFull, &b_nCutsFull);
    chains_fullSel[isample]->SetBranchAddress("nSel",  nSelFull,   &b_nSelFull);
    
    Long64_t nentriesPre  = chains_preSel[isample]->GetEntries();
    Long64_t nentriesFull = chains_fullSel[isample]->GetEntries();
    
    // pre-selection
    for (Long64_t jentry=0; jentry<nentriesPre;jentry++) {
      Long64_t nb = chains_preSel[isample]->GetEntry(jentry);   
      //      nCutsAnaPre = nCutsPre;      
      for(int icut=0; icut<nCutsPre; icut++) nPreSelTot[icut][isample] += nSelPre[icut];	
    }
    
    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      // nCutsAnaFull = nCutsFull;
      for(int icut=0; icut<nCutsAnaFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at preselection level
  for(int icut=0; icut<nCutsAnaPre; icut++) {

    // numbers
    // N = Ni * weight * lumi / 10 (weights computed for 10pb-1)
    Wenu_preSel[icut]   = lumi / lumiWeight * Wenu_weight  * nPreSelTot[icut][0];
    ttbar_preSel[icut]  = lumi / lumiWeight * ttbar_weight * nPreSelTot[icut][1];
    Zee_preSel[icut]    = lumi / lumiWeight * Zee_weight   * nPreSelTot[icut][2];
    
    float qcd_tmp=0.;
    qcd_tmp += lumi / lumiWeight * QCD_em20to30_weight     * nPreSelTot[icut][3];
    qcd_tmp += lumi / lumiWeight * QCD_em30to80_weight     * nPreSelTot[icut][4];
    qcd_tmp += lumi / lumiWeight * QCD_em80to170_weight    * nPreSelTot[icut][5];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe20to30_weight  * nPreSelTot[icut][6];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe30to80_weight  * nPreSelTot[icut][7];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe80to170_weight * nPreSelTot[icut][8];
    QCD_preSel[icut] = qcd_tmp;
    
    
    // efficiencies
    if(icut>0 && nPreSelTot[icut-1][0]>0) Wenu_eff_preSel[icut]   = nPreSelTot[icut][0]/nPreSelTot[icut-1][0];
    else Wenu_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][1]>0) ttbar_eff_preSel[icut]  = nPreSelTot[icut][1]/nPreSelTot[icut-1][1];
    else ttbar_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][2]>0) Zee_eff_preSel[icut]    = nPreSelTot[icut][2]/nPreSelTot[icut-1][2];
    else Zee_eff_preSel[icut] = 0.0;
    if(icut>0 && QCD_preSel[icut-1]>0) QCD_eff_preSel[icut]       = QCD_preSel[icut]/QCD_preSel[icut-1];
    else QCD_eff_preSel[icut] = 0.0;
  }

  // final efficiency at preselection level
  if(nPreSelTot[0][0]>0) Wenu_finaleff_preSel   = nPreSelTot[nCutsAnaPre-2][0]/nPreSelTot[0][0]; 
  else Wenu_finaleff_preSel = 0.0;
  if(nPreSelTot[0][1]>0) ttbar_finaleff_preSel  = nPreSelTot[nCutsAnaPre-2][1]/nPreSelTot[0][1]; 
  else ttbar_finaleff_preSel = 0.0;
  if(nPreSelTot[0][2]>0) Zee_finaleff_preSel    = nPreSelTot[nCutsAnaPre-2][2]/nPreSelTot[0][2]; 
  else Zee_finaleff_preSel = 0.0;
  if(QCD_preSel[0]>0) QCD_finaleff_preSel     = QCD_preSel[nCutsAnaPre-2]/QCD_preSel[0];
  else QCD_finaleff_preSel = 0.0;
  
  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    Wenu_fullSel[icut]   = lumi / lumiWeight * Wenu_weight   * nFullSelTot[icut][0];
    ttbar_fullSel[icut]  = lumi / lumiWeight * ttbar_weight  * nFullSelTot[icut][1];
    Zee_fullSel[icut]    = lumi / lumiWeight * Zee_weight    * nFullSelTot[icut][2];

    float qcd_tmp=0.;
    qcd_tmp += lumi / lumiWeight * QCD_em20to30_weight     * nFullSelTot[icut][3];
    qcd_tmp += lumi / lumiWeight * QCD_em30to80_weight     * nFullSelTot[icut][4];
    qcd_tmp += lumi / lumiWeight * QCD_em80to170_weight    * nFullSelTot[icut][5];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe20to30_weight  * nFullSelTot[icut][6];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe30to80_weight  * nFullSelTot[icut][7];
    qcd_tmp += lumi / lumiWeight * QCD_bctoe80to170_weight * nFullSelTot[icut][8];
    QCD_fullSel[icut] = qcd_tmp;

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) Wenu_eff_fullSel[icut]  = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else Wenu_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][1]>0) ttbar_eff_fullSel[icut] = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else ttbar_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][2]>0) Zee_eff_fullSel[icut]   = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else Zee_eff_fullSel[icut] = 0.0;
    if(icut>0 && QCD_fullSel[icut-1]>0) QCD_eff_fullSel[icut]      = QCD_fullSel[icut]/QCD_fullSel[icut-1];
    else QCD_eff_fullSel[icut] = 0.0;

    if(icut==0) { 
      Wenu_eff_fullSel[icut]   = nFullSelTot[icut][0]/nPreSelTot[nCutsAnaPre-2][0];    
      ttbar_eff_fullSel[icut]  = nFullSelTot[icut][1]/nPreSelTot[nCutsAnaPre-2][1];    
      Zee_eff_fullSel[icut]    = nFullSelTot[icut][2]/nPreSelTot[nCutsAnaPre-2][2];    
      QCD_eff_fullSel[icut]    = QCD_fullSel[icut]/QCD_preSel[nCutsAnaPre-2];    
    }

  }

  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0) Wenu_finaleff_fullSel  = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0]; 
  else Wenu_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) ttbar_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1]; 
  else ttbar_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][2]>0) Zee_finaleff_fullSel   = nFullSelTot[nCutsAnaFull-4][2]/nFullSelTot[0][2]; 
  else Zee_finaleff_fullSel = 0.0;
  if(QCD_fullSel[0]>0) QCD_finaleff_fullSel      = QCD_fullSel[nCutsAnaPre-4]/QCD_fullSel[0];
  else QCD_finaleff_fullSel = 0.0;

  // final efficiency combining pre and full selections
  if(nFullSelTot[0][0]>0) Wenu_finaleff   = nFullSelTot[nCutsAnaFull-4][0]/nPreSelTot[0][0]; 
  else Wenu_finaleff = 0.0;
  if(nFullSelTot[0][1]>0) ttbar_finaleff  = nFullSelTot[nCutsAnaFull-4][1]/nPreSelTot[0][1]; 
  else ttbar_finaleff = 0.0;
  if(nFullSelTot[0][2]>0) Zee_finaleff    = nFullSelTot[nCutsAnaFull-4][2]/nPreSelTot[0][2]; 
  else Zee_finaleff = 0.0;
  if(QCD_fullSel[0]>0) QCD_finaleff        = QCD_fullSel[nCutsAnaFull-4]/QCD_preSel[0]; 
  else QCD_finaleff = 0.0;


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

  fullSelCuts[0]="this channel preselected";
  fullSelCuts[1]="hardLeptonThreshold";
  fullSelCuts[2]="slowLeptonThreshold";
  fullSelCuts[3]="dileptonInvMassMin";
  fullSelCuts[4]="classDepEleId";
  fullSelCuts[5]="trackerIso";
  fullSelCuts[6]="hcalIso";
  fullSelCuts[7]="ecalIso";
  fullSelCuts[8]="globalIso";
  fullSelCuts[9]="lepton d0";
  fullSelCuts[10]="MET";
  fullSelCuts[11]="maxPtLepton";
  fullSelCuts[12]="minPtLepton";
  fullSelCuts[13]="dileptonInvMassMax";
  fullSelCuts[14]="detaLeptons";
  fullSelCuts[15]="deltaPhi";
  fullSelCuts[16]="final";
  fullSelCuts[17]="0 jets";
  fullSelCuts[18]="1 jets";
  fullSelCuts[19]=">1 jets";
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
	   << "\\begin{tabular}{|c|c|c|c|c|}" << endl
	   << "\\hline" << endl
	   << "selection & W$(e \\nu)$ & $t\\bar{t}$ & Z$(e^+e^-)$ & QCD \t\\\\" << endl
	   << "\\hline" << endl; 
  
  textfile << "\\hline"        << endl;
  textfile << "\\hline"        << endl;

  for(int icut=0; icut<9; icut++) {
    textfile << preSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << Wenu_preSel[icut]  << " (" << 100. * Wenu_eff_preSel[icut]     << "\\%)" << "\t&\t"
	     << ttbar_preSel[icut] << " (" << 100. * ttbar_eff_preSel[icut]    << "\\%)" << "\t&\t"
	     << Zee_preSel[icut]   << " (" << 100. * Zee_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << QCD_preSel[icut]   << " (" << 100. * QCD_eff_preSel[icut]   << "\\%)" 
	     << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total preselection " << "\t&\t"
	   << Wenu_preSel[8]  << " (" << 100. * Wenu_finaleff_preSel << "\\%)"  << "\t&\t"
	   << ttbar_preSel[8] << " (" << 100. * ttbar_finaleff_preSel << "\\%)" << "\t&\t"
	   << Zee_preSel[8]   << " (" << 100. * Zee_finaleff_preSel << "\\%)"  << "\t&\t"
	   << QCD_preSel[8]   << " (" << 100. * QCD_finaleff_preSel << "\\%)"
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;
  
  for(int icut=0; icut<17; icut++) {
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << Wenu_fullSel[icut]  << " (" << 100. * Wenu_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << ttbar_fullSel[icut] << " (" << 100. * ttbar_eff_fullSel[icut] << "\\%)" << "\t&\t"
	     << Zee_fullSel[icut]   << " (" << 100. * Zee_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << QCD_fullSel[icut]   << " (" << 100. * QCD_eff_fullSel[icut]  << "\\%)"
	     << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;

  textfile << "total fullselection " << "\t&\t"
	   << Wenu_fullSel[16]  << " (" << 100. * Wenu_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << ttbar_fullSel[16] << " (" << 100. * ttbar_finaleff_fullSel << "\\%)"    << "\t&\t"
	   << Zee_fullSel[16]   << " (" << 100. * Zee_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << QCD_fullSel[16]   << " (" << 100. * QCD_finaleff_fullSel << "\\%)"
	   << "\t\\\\" << endl;

  textfile << "\\hline" << endl;

  textfile << "total " << "\t&\t"
	   << Wenu_fullSel[16]  << " (" << 100. * Wenu_finaleff  << "\\%)"    << "\t&\t"
	   << ttbar_fullSel[16] << " (" << 100. * ttbar_finaleff << "\\%)"    << "\t&\t"
	   << Zee_fullSel[16]   << " (" << 100. * Zee_finaleff << "\\%)"     << "\t&\t"
	   << QCD_fullSel[16]   << " (" << 100. * QCD_finaleff << "\\%)"
	   << "\t\\\\" << endl;

  textfile << "0 jets bin " << "\t&\t"
	   << Wenu_fullSel[17]  << "\t&\t"
	   << ttbar_fullSel[17] << "\t&\t"
	   << Zee_fullSel[17]   << "\t&\t"
	   << QCD_fullSel[17]
	   << "\t\\\\" << endl;
	   
  textfile << "1 jets bin " << "\t&\t"
	   << Wenu_fullSel[18]  << "\t&\t"
	   << ttbar_fullSel[18] << "\t&\t"
	   << Zee_fullSel[18]   << "\t&\t"
	   << QCD_fullSel[18]
	   << "\t\\\\" << endl;

  textfile << "$>1$ jets bin " << "\t&\t"
	   << Wenu_fullSel[19]  << "\t&\t"
	   << ttbar_fullSel[19] << "\t&\t"
	   << Zee_fullSel[19]   << "\t&\t"
	   << QCD_fullSel[19]
           << "\t\\\\" << endl;
	   
  
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

