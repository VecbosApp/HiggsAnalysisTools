#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "massDependentCuts.cc"

using namespace std;

enum { ee=0, mm=1, em=2, me=3, ll=4 };

float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void HiggsYields(int mH, int njets, bool showData);
void printLatex(bool showData);

void printLatex(bool showData) {

  ofstream textfile;
  textfile.open("yields.tex", ios_base::trunc);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{fullpage}" << endl;
  textfile << "\\begin{document}" << endl << endl;
  textfile << "\\section{Yields for 0 jet}" << endl;

  textfile << "\\subsection{Yields at WW selection level}" << std::endl;
  std::cout << "Evaluating yields at WW selection level" << std::endl;
  HiggsYields(-1, 0, showData);

  textfile.open("yields.tex", ios_base::app);
  textfile << "\\subsection{Yields at full H selection level}" << std::endl;

  int mH[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  for(int i=0; i<17;i++) {
    std::cout << "mH = " << mH[i] << "\t0 jet" << std::endl;
    HiggsYields(mH[i], 0, showData);
  }
  textfile.close();

  ofstream textfile2;
  textfile2.open("yields.tex", ios_base::app);
  textfile2 << "\\cleardoublepage" << endl << endl;
  textfile2 << "\\section{Yields for 1 jet}" << endl;
  textfile2 << "\\subsection{Yields at WW selection level}" << std::endl;
  std::cout << "Evaluating yields at WW selection level" << std::endl;
  HiggsYields(-1, 1, showData);

  textfile2.open("yields.tex", ios_base::app);
  textfile2 << "\\subsection{Yields at full H selection level}" << std::endl;
  for(int i=0; i<17;i++) {
    std::cout << "mH = " << mH[i] << "\t1 jet" << std::endl;
    HiggsYields(mH[i], 1, showData);
  }
  textfile2.close();

  ofstream textfile3;
  textfile3.open("yields.tex", ios_base::app);
  textfile3 << "\\end{document}" << endl << endl;

}

void HiggsYields(int mH, int njets, bool showData) {

  std::vector<std::vector<double> > yields;
  std::vector<std::vector<double> > yields_err;

  std::vector<float> yields_bkgtot, yields_bkgtot_err;
  std::vector<float> yields_data;
  for(int icha=0; icha<5; icha++) {
    yields_bkgtot.push_back(0.);
    yields_bkgtot_err.push_back(0.);
    yields_data.push_back(0.);
  }

  std::vector<std::string> sampleName;
  sampleName.push_back("Z/$\\gamma^*$");
  sampleName.push_back("$t \\bar{t}$");
  sampleName.push_back("single t");
  sampleName.push_back("W+jets");
  sampleName.push_back("WZ+ZZ");
  sampleName.push_back("qq$\\to$ WW");
  sampleName.push_back("gg$\\to$ WW");
  char mass[10];
  sprintf(mass,"H%d",mH);
  if(mH!=-1) sampleName.push_back(std::string(mass));
  else  sampleName.push_back("H130"); // take one mass as example

  // open files
  TFile *fileZj = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileTTbar = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileSingleTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWj = TFile::Open("results/datasets_trees/Wjets_ll.root");
  TFile *fileOthers = TFile::Open("results/datasets_trees/others_ll.root");
  TFile *fileqqWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileggWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileData = TFile::Open("results_data/datasets_trees/dataset_ll.root");

  char signalFile[200];
  if(mH!=-1) sprintf(signalFile, "results/datasets_trees/H%d_ll.root", mH);
  else sprintf(signalFile, "results/datasets_trees/H130_ll.root"); // take one mass as example
  TFile *fileH = TFile::Open(signalFile);

  // get trees
  TTree* treeZj = (TTree*)fileZj->Get("T1");
  TTree* treeTTbar = (TTree*)fileTTbar->Get("T1");
  TTree* treeSingleTop = (TTree*)fileSingleTop->Get("T1");
  TTree* treeWj = (TTree*)fileWj->Get("T1");
  TTree* treeOthers = (TTree*)fileOthers->Get("T1");
  TTree* treeqqWW = (TTree*)fileqqWW->Get("T1");
  TTree* treeggWW = (TTree*)fileggWW->Get("T1");
  TTree* treeH = (TTree*)fileH->Get("T1");
  TTree* treeData = (TTree*)fileData->Get("T1");

  std::vector<TTree*> trees;
  trees.push_back(treeZj); // 0
  trees.push_back(treeTTbar); // 1
  trees.push_back(treeSingleTop); // 2
  trees.push_back(treeWj); // 3
  trees.push_back(treeOthers); // 4
  trees.push_back(treeqqWW); // 5
  trees.push_back(treeggWW); // 6
  trees.push_back(treeH); // 7

  TH1F *histo = new TH1F("histo","histo",100,0,180);
  
  char njcut[30];
  sprintf(njcut, "njets==%d", njets);
  char wwselcut[30];
  if(njets==0) sprintf(wwselcut,"WWSel");
  else if(njets==1) sprintf(wwselcut,"WWSel1j");
  else {
    std::cout << "Jet bin must be 0/1" << std::endl;
    return;
  }

  TString higgsMassDependentCut;
  if(mH==-1) higgsMassDependentCut = TString("1"); // i.e. give the yields at WW level
  else higgsMassDependentCut = higgsCuts(mH,true);

  for(int isample=0; isample<(int)trees.size(); isample++) {
    TTree *tree = trees[isample];

    std::vector<TString> cutChannel;
    TString HCut_ee = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==0)*weight*puweight");
    TString HCut_mm = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==1)*weight*puweight");
    TString HCut_em = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==2)*weight*puweight");
    TString HCut_me = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==3)*weight*puweight");
    TString HCut_all = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(")*weight*puweight");
    
    cutChannel.push_back(HCut_ee);
    cutChannel.push_back(HCut_mm);
    cutChannel.push_back(HCut_em);
    cutChannel.push_back(HCut_me);
    cutChannel.push_back(HCut_all);

    std::vector<double> sampleYield, sampleYield_err;
    for(int icha=0; icha<5; icha++) {

      // some ROOT file contains more than 1 bkg
      // split top / ttbar
      TString addCut=TString("");
      if(isample==1) addCut = TString("*(process==11)");
      if(isample==2) addCut = TString("*(process>=8 && process<=10)");
      // split qqWW and ggWW
      if(isample==5) addCut = TString("*(process==13)");
      if(isample==6) addCut = TString("*(process==14)");

      TString TheFinalCut = cutChannel[icha]+addCut;

      //      cout << "final cut = " << TheFinalCut.Data() << endl;

      tree->Project("histo","deltaPhi",TheFinalCut);
      double yield = histo->Integral();
      double yield_err = yieldErrPoisson(yield,histo->GetEntries());
      sampleYield.push_back(yield);
      sampleYield_err.push_back(yield_err);


      // sum the backgrounds (signal is the last)
      if(isample != (int)trees.size()-1) {
        yields_bkgtot[icha] += yield;
        yields_bkgtot_err[icha] += yield_err * yield_err;
      }
    }

    yields.push_back(sampleYield);
    yields_err.push_back(sampleYield_err);

  }

  // quadrature sum of the tot bkg error
  for(int isample=0; isample<(int)trees.size(); isample++) {
    yields_bkgtot_err[isample] = sqrt(yields_bkgtot_err[isample]);
  }

  // data counting
  for(int icha=0; icha<5; icha++) {
    std::vector<TString> cutChannel;
    TString HCut_ee = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==0)");
    TString HCut_mm = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==1)");
    TString HCut_em = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==2)");
    TString HCut_me = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(" && finalstate==3)");
    TString HCut_all = TString("(")+TString(wwselcut)+TString(" && ")+higgsMassDependentCut+TString(")");
    
    cutChannel.push_back(HCut_ee);
    cutChannel.push_back(HCut_mm);
    cutChannel.push_back(HCut_em);
    cutChannel.push_back(HCut_me);
    cutChannel.push_back(HCut_all);
    
    treeData->Project("histo","deltaPhi",cutChannel[icha]);
    double yield = histo->Integral();
    yields_data[icha] = yield;
  }

  ofstream textfile;
  textfile.open("yields.tex", ios_base::app);
  textfile.setf(ios::fixed,ios::floatfield);
  textfile.precision(1);
  textfile << "\\begin{table}[p]" << endl
           << "\\begin{small}" << endl
           << "\\begin{center}" << endl;
  if(showData) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  else textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;

  // header
  for(int isample=0; isample<(int)trees.size(); isample++) {
    if(isample==0) textfile << "\t & " << sampleName[isample] << " & ";
    else if(isample==(int)trees.size()-1 && !showData) textfile << sampleName[isample] << " \\\\ " << std::endl;
    else if(isample==(int)trees.size()-1 && showData) textfile << sampleName[isample] << " & data \\\\ " << std::endl;
    else if(isample==(int)trees.size()-2) textfile << sampleName[isample] << " & bkg. tot. & ";
    else textfile << sampleName[isample] << " & ";
  }

  textfile << "\\hline" << std::endl;

  std::vector<TString> chanName;
  chanName.push_back("ee");
  chanName.push_back("$\\mu\\mu$");
  chanName.push_back("e$\\mu$");
  chanName.push_back("$\\mu$e");
  chanName.push_back("all");

  // yields
  for(int icha=0; icha<5; icha++) {
    for(int isample=0; isample<(int)trees.size(); isample++) {
      std::vector<double> sampleYiled = yields[isample];
      std::vector<double> sampleYiled_err = yields_err[isample];

      double val = sampleYiled[icha];
      double err = sampleYiled_err[icha];
      
      if(isample==0) textfile << chanName[icha] << " & " << val << " $\\pm$  " << err << " & ";
      else if(isample==(int)trees.size()-1 && !showData) textfile << val << " $\\pm$  " << err << " \\\\ " << std::endl;
      else if(isample==(int)trees.size()-1 && showData) textfile << val << " $\\pm$  " << err << " & " << yields_data[icha] << " \\\\ " << std::endl;
      else if(isample==(int)trees.size()-2) textfile << val << " $\\pm$  " << err << " & " << yields_bkgtot[icha] << " $\\pm$  " << yields_bkgtot_err[icha] << " & ";
      else textfile << val << " $\\pm$  " << err << " & ";
    }
    if(icha==3) textfile << "\\hline" << std::endl;
  }

  // trailer
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl;
  if(mH==-1) textfile << "\\caption{WW selection level, " << njets << " jet.}" << std::endl;
  else textfile << "\\caption{Higgs $m_H$ = " << mH << " GeV/c$^2$, " << njets << " jet.}" << std::endl;
  textfile << "\\end{center}" << endl
           << "\\end{small}" << endl
           << "\\end{table}" << endl;

  delete histo;

}

float yieldErrPoisson(float nEst1, float n1, float nEst2, float n2, float nEst3, float n3, float nEst4, float n4, float nEst5, float n5, float nEst6, float n6) {

  float sum=0;
  if(n1>0) sum += pow(nEst1,2)/n1;
  if(n2>0) sum += pow(nEst2,2)/n2;
  if(n3>0) sum += pow(nEst3,2)/n3;
  if(n4>0) sum += pow(nEst4,2)/n4;
  if(n5>0) sum += pow(nEst5,2)/n5;
  if(n6>0) sum += pow(nEst6,2)/n6;
  
  return sqrt(sum);
}
