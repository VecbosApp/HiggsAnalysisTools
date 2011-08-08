#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include "massDependentCuts.cc"

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
TString kinematicCut(int ibin);
void estimateWjets(int njets);

enum { ee=0, mm=1, em=2, me=3 };
enum { barrellowpt=0, barrelhighpt=1, endcaplowpt=2, endcaphighpt=3 };

void printLatex() {

  ofstream textfile;
  textfile.open("Wjets_yields.tex", ios_base::trunc);
  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{fullpage}" << endl;
  textfile << "\\begin{document}" << endl << endl;
  std::cout << "now estimating 0 jet bin " << std::endl;
  estimateWjets(0);

  std::cout << "now estimating 1 jet bin " << std::endl;
  estimateWjets(1);

  ofstream textfile2;
  textfile2.open("Wjets_yields.tex", ios_base::app);
  textfile2 << "\\end{document}" << endl << endl;

}

void estimateWjets(int njets) {
  
  char njetscut[20];
  sprintf(njetscut,"(njets==%d)",njets);

  float yield_WWSel[5][4]; // [bin][icha] bin = 5 => total
  float yield_WWSel_staterr[5][4]; // [bin][icha]

  // for now only ee is done
  TFile *fileEE = TFile::Open("results_data/datasets_trees_skim/dataset_fakes_ee.root");
  TFile *fileMM = TFile::Open("results_data/datasets_trees_skim/dataset_fakes_ee.root");
  TFile *fileEM = TFile::Open("results_data/datasets_trees_skim/dataset_fakes_ee.root");
  TFile *fileME = TFile::Open("results_data/datasets_trees_skim/dataset_fakes_ee.root");

  TTree *treeEE = (TTree*)fileEE->Get("T1");
  TTree *treeMM = (TTree*)fileMM->Get("T1");
  TTree *treeEM = (TTree*)fileEM->Get("T1");
  TTree *treeME = (TTree*)fileME->Get("T1");

  std::vector<TTree*> trees;
  trees.push_back(treeEE);
  trees.push_back(treeMM);
  trees.push_back(treeEM);
  trees.push_back(treeME);

  TH1F *histo = new TH1F("histo","",100,0.,180.);

  // yields at WW level
  for(int icha=0; icha<4; icha++) {
    for(int ibin=0; ibin<4; ibin++) {
      
      TString fpCut = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP");
      TString fpCutStatErr = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP");

      trees[icha]->Project("histo","dphill",fpCut);
      yield_WWSel[ibin][icha] = histo->Integral();

      trees[icha]->Project("histo","dphill",fpCutStatErr);
      yield_WWSel_staterr[ibin][icha] = sqrt(histo->Integral());

    }
  }  

  for(int icha=0; icha<4; icha++) {
    for(int ibin=0; ibin<4; ibin++) {
      if(ibin==0) { 
        yield_WWSel[4][icha] = yield_WWSel[ibin][icha];
        yield_WWSel_staterr[4][icha] = pow(yield_WWSel_staterr[ibin][icha],2);
      }
      else {
        yield_WWSel[4][icha] += yield_WWSel[ibin][icha];
        yield_WWSel_staterr[4][icha] += pow(yield_WWSel_staterr[ibin][icha],2);}
    }
  }

  for(int icha=0; icha<4; icha++) yield_WWSel_staterr[4][icha] = sqrt(yield_WWSel_staterr[4][icha]);

  std::vector<TString> labels;
  labels.push_back(TString("Barrel, 10 $\\leq p_{T} < $ 20 GeV"));
  labels.push_back(TString("Barrel, $p_{T} \\geq $ 20 GeV"));
  labels.push_back(TString("Endcap, 10 $\\leq p_{T} < $ 20 GeV"));
  labels.push_back(TString("Endcap, $p_{T} \\geq $ 20 GeV"));
  labels.push_back(TString("Total"));

  ofstream textfile;
  textfile.open("Wjets_yields.tex", ios_base::app);
  textfile.setf(ios::fixed,ios::floatfield);
  textfile.precision(1);

  // ee and me
  textfile << "\\begin{table}[p]" << endl;
  textfile << "\\begin{small}" << endl;
  textfile << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  textfile << "Fake Lepton Bin \t & \t Electron + Fake Electron \t & \t Muon + Fake Electron \\\\" << std::endl;
  textfile << "\\hline" << endl;
  for(int ibin=0; ibin<5; ibin++) {
    if(ibin==4) textfile << "\\hline" << endl;
    textfile << labels[ibin] << " & "
             << yield_WWSel[ibin][ee] << " $ \\pm $ " << yield_WWSel_staterr[ibin][ee] << " (stat.) & "
             << yield_WWSel[ibin][me] << " $ \\pm $ " << yield_WWSel_staterr[ibin][me] << " (stat.) \\\\ "
             << std::endl;
  }
  textfile << "\\hline" << endl;
  textfile << "\\end{tabular}" << endl;
  textfile << "\\caption{W+jets yields in the ee and $\\mu$e channels at WW selection level in the " 
           << njets << " jet bin.}" << std::endl;
  textfile << "\\end{center}" << endl;
  textfile << "\\end{small}" << endl;
  textfile << "\\end{table}" << endl;


  // mm and em
  textfile << "\\begin{table}[p]" << endl;
  textfile << "\\begin{small}" << endl;
  textfile << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  textfile << "Fake Lepton Bin \t & \t Electron + Fake Muon \t & \t Muon + Fake Muon \\\\" << std::endl;
  textfile << "\\hline" << endl;
  for(int ibin=0; ibin<5; ibin++) {
    if(ibin==4) textfile << "\\hline" << endl;
    textfile << labels[ibin] << " & "
             << yield_WWSel[ibin][em] << " $ \\pm $ " << yield_WWSel_staterr[ibin][em] << " (stat.) & "
             << yield_WWSel[ibin][mm] << " $ \\pm $ " << yield_WWSel_staterr[ibin][mm] << " (stat.) \\\\ "
             << std::endl;
  }
  textfile << "\\hline" << endl;
  textfile << "\\end{tabular}" << endl;
  textfile << "\\caption{W+jets yields in the $\\mu$e and $\\mu\\mu$ channels at WW selection level in the " 
           << njets << " jet bin.}" << std::endl;
  textfile << "\\end{center}" << endl;
  textfile << "\\end{small}" << endl;
  textfile << "\\end{table}" << endl;


  char nameFileTable[100];
  sprintf(nameFileTable, "WJetsYieldsData_ForTable_%dj.txt",njets);
  ofstream tablefile;
  tablefile.open(nameFileTable, ios_base::trunc);
  tablefile.precision(2);

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    TString higgsMassDependentCut = higgsCuts(mass,true);
    
    // yields at WW level
    for(int icha=0; icha<4; icha++) {
        TString fpCut = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP");
        TString fpCutStatErr = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP");
        
        trees[icha]->Project("histo","dphill",fpCut);
        yield_WWSel[4][icha] = histo->Integral();
        
        trees[icha]->Project("histo","dphill",fpCutStatErr);
        yield_WWSel_staterr[4][icha] = sqrt(histo->Integral());
    }

    // total for all the channels
    float yield_tot = 0.;
    float yield_tot_err = 0.;
    for(int icha=0; icha<4; icha++) {
      yield_tot += yield_WWSel[4][icha];
      yield_tot_err += pow(yield_WWSel_staterr[4][icha],2);
    }
    yield_tot_err = sqrt(yield_tot_err);

    // summary table for limits
    if (i==0) { 
      tablefile << "# " << njets << " jets bin data" << endl;
      tablefile << "# \t\t mumu \t\t mue \t\t emu \t\t ee \t\t ll" << endl;
    }
    tablefile << mass 
              << "\t\t" << yield_WWSel[4][1] << " +/- " <<  yield_WWSel_staterr[4][1] 
              << "\t\t" << yield_WWSel[4][3] << " +/- " <<  yield_WWSel_staterr[4][3] 
              << "\t\t" << yield_WWSel[4][2] << " +/- " <<  yield_WWSel_staterr[4][2] 
              << "\t\t" << yield_WWSel[4][0] << " +/- " <<  yield_WWSel_staterr[4][0] 
              << "\t\t" << yield_tot << " +/- " <<  yield_tot_err  
              << std::endl;

  }

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
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

TString kinematicCut(int ibin) {
  if(ibin==barrellowpt) return TString("(abs(eta2)<=1.479 && pt2<20)");
  if(ibin==barrelhighpt) return TString("(abs(eta2)<=1.479 && pt2>=20)");
  if(ibin==endcaplowpt) return TString("(abs(eta2)>1.479 && pt2<20)");
  if(ibin==endcaphighpt) return TString("(abs(eta2)>1.479 && pt2>=20)");
  return TString("wrong bin chosen. Check method kinematicCut(int ibin)");
}
