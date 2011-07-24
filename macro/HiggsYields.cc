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
void HiggsYields(int mH, int njets);
void printLatex();

void printLatex() {

  ofstream textfile;
  textfile.open("yields.tex", ios_base::trunc);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl << endl;
  textfile << "\\section{Yields for 0 jet}" << endl;

  int mH[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  for(int i=0; i<17;i++) {
    std::cout << "mH = " << mH[i] << "\t0 jet" << std::endl;
    HiggsYields(mH[i], 0);
  }

  textfile << endl << endl;
  textfile << "\\cleardoublepage" << endl << endl;
  textfile << "\\section{Yields for 1 jet}" << endl;
  for(int i=0; i<17;i++) {
    std::cout << "mH = " << mH[i] << "\t1 jet" << std::endl;
    HiggsYields(mH[i], 1);
  }

  textfile << "\\end{document}" << endl << endl;

}

void HiggsYields(int mH, int njets) {

  std::vector<std::vector<double> > yields;
  std::vector<std::vector<double> > yields_err;

  std::vector<std::string> sampleName;
  sampleName.push_back("Z/$\\gamma^*$");
  sampleName.push_back("$t \\bar{t}$");
  sampleName.push_back("single t");
  sampleName.push_back("W+jets");
  sampleName.push_back("WZ+ZZ");
  sampleName.push_back("gg$\\to$ WW");
  sampleName.push_back("qq$\\to$ WW");
  char mass[10];
  sprintf(mass,"H%d",mH);
  sampleName.push_back(std::string(mass));


  // open files
  TFile *fileZj = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileTTbar = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileSingleTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWj = TFile::Open("results/datasets_trees/Wjets_ll.root");
  TFile *fileOthers = TFile::Open("results/datasets_trees/others_ll.root");
  TFile *fileqqWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileggWW = TFile::Open("results/datasets_trees/WW_ll.root");

  char signalFile[200];
  sprintf(signalFile, "results/datasets_trees/H%d_ll.root", mH);
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

  for(int isample=0; isample<(int)trees.size(); isample++) {
    TTree *tree = trees[isample];

    std::vector<TString> cutChannel;
    TString HCut_ee = TString("(")+TString(wwselcut)+TString(" && ")+higgsCuts(mH,true)+TString(" && finalstate==0)*weight*puweight");
    TString HCut_mm = TString("(")+TString(wwselcut)+TString(" && ")+higgsCuts(mH,true)+TString(" && finalstate==1)*weight*puweight");
    TString HCut_em = TString("(")+TString(wwselcut)+TString(" && ")+higgsCuts(mH,true)+TString(" && finalstate==2)*weight*puweight");
    TString HCut_me = TString("(")+TString(wwselcut)+TString(" && ")+higgsCuts(mH,true)+TString(" && finalstate==3)*weight*puweight");
    TString HCut_all = TString("(")+TString(wwselcut)+TString(" && ")+higgsCuts(mH,true)+TString(")*weight*puweight");
    
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
    }

    yields.push_back(sampleYield);
    yields_err.push_back(sampleYield_err);

  }
  
  ofstream textfile;
  textfile.open("yields.tex", ios_base::app);
  textfile.precision(2);
  textfile << "\\begin{table}[p]" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;

  // header
  for(int isample=0; isample<(int)trees.size(); isample++) {
    if(isample==0) textfile << "\t & " << sampleName[isample] << " & ";
    else if(isample==(int)trees.size()-1) textfile << sampleName[isample] << " \\\\ " << std::endl;
    else textfile << sampleName[isample] << " & ";
  }

  textfile << "\\hline" << std::endl;

  std::vector<TString> chanName;
  chanName.push_back("ee");
  chanName.push_back("$\\mu\\mu$");
  chanName.push_back("e$\\mu$");
  chanName.push_back("$\\mu$ e");
  chanName.push_back("all");

  // yields
  for(int icha=0; icha<5; icha++) {
    for(int isample=0; isample<(int)trees.size(); isample++) {
      std::vector<double> sampleYiled = yields[isample];
      std::vector<double> sampleYiled_err = yields_err[isample];

      double val = sampleYiled[icha];
      double err = sampleYiled_err[icha];
      
      if(isample==0) textfile << chanName[icha] << " & " << val << " $\\pm$  " << err << " & ";
      else if(isample==(int)trees.size()-1) textfile << val << " $\\pm$  " << err << " \\\\ " << std::endl;
      else textfile << val << " $\\pm$  " << err << " & ";
    }
    if(icha==3) textfile << "\\hline" << std::endl;
  }

  // trailer
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\caption{Higgs $m_H$ = " << mH << " GeV/c$^2$, " << njets << " jet.}" << std::endl
           << "\\end{center}" << endl
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
