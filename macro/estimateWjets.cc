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
  char wwLevelCut[20];
  if(njets==0) sprintf(wwLevelCut,"WWSel");
  if(njets==1) sprintf(wwLevelCut,"WWSel1j");

  float yield_WWSel[5][4]; // [bin][icha] bin = 5 => total
  float yield_WWSel_staterr[5][4]; // [bin][icha]

  // for now only ee and me are done
  TFile *fileEE = TFile::Open("results_data/datasets_trees_fake/dataset_fake_ee.root");
  TFile *fileMM = TFile::Open("results_data/datasets_trees_fake/dataset_fake_mm.root");
  TFile *fileEM = TFile::Open("results_data/datasets_trees_fake/dataset_fake_em.root");
  TFile *fileME = TFile::Open("results_data/datasets_trees_fake/dataset_fake_me.root");

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

      TString fpCut, fpCutStatErr;
      if(icha==0 || icha==3) { // chiara has the hlt bit and WWSel to be applied
        fpCut = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP*hlt*") + TString(wwLevelCut);
        fpCutStatErr = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP*hlt*") + TString(wwLevelCut);
      } else { // for mm, em HLT already applied
        fpCut = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP");
        fpCutStatErr = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP");
      }
      // for MC closure test
      // TString fpCut = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*1.545*weightFP*baseW*") + TString(wwLevelCut);
      // TString fpCutStatErr = TString("(") + kinematicCut(ibin) + TString(" && ") + TString(njetscut) + TString(")") + TString("*1.545*weightStatFP*baseW*") + TString(wwLevelCut);
      
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

  // total for all the channels
  float yield_WWSel_tot = 0.;
  float yield_WWSel_tot_err = 0.;
  for(int icha=0; icha<4; icha++) {
    yield_WWSel_tot += yield_WWSel[4][icha];
    yield_WWSel_tot_err += pow(yield_WWSel_staterr[4][icha],2);
  }
  yield_WWSel_tot_err = sqrt(yield_WWSel_tot_err);

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
  textfile << "\\caption{W+jets yields in the e$\\mu$ and $\\mu\\mu$ channels at WW selection level in the " 
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
      TString fpCut, fpCutStatErr;
      if(icha==0 || icha==3) { // chiara has the hlt bit
        fpCut = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP*hlt*") + TString(wwLevelCut);
        fpCutStatErr = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP*hlt*") + TString(wwLevelCut);
      } else { // alicia has already applied WWSel + HLT
        fpCut = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightFP");
        fpCutStatErr = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*weightStatFP");
      }
      // for MC closure test
      // TString fpCut = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*1.545*weightFP*baseW*") + TString(wwLevelCut);
      // TString fpCutStatErr = TString("(") + higgsMassDependentCut + TString(" && ") + TString(njetscut) + TString(")") + TString("*1.545*weightStatFP*baseW*") + TString(wwLevelCut);
      
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


  /// ===> compute !b-veto efficiency on W+jets 
  // this is useful for W+jets contamination in the top-tagged control sample
  // use a channel independent estimate
  // for now only ee is done
  TFile *fileEEAll = TFile::Open("results_data/datasets_trees_fake/dataset_fake_ee.root");
  TFile *fileMMAll = TFile::Open("results_data/datasets_trees_fake/dataset_fake_ee.root");
  TFile *fileEMAll = TFile::Open("results_data/datasets_trees_fake/dataset_fake_ee.root");
  TFile *fileMEAll = TFile::Open("results_data/datasets_trees_fake/dataset_fake_me.root");

  TTree *treeEEAll = (TTree*)fileEEAll->Get("T1");
  TTree *treeMMAll = (TTree*)fileMMAll->Get("T1");
  TTree *treeEMAll = (TTree*)fileEMAll->Get("T1");
  TTree *treeMEAll = (TTree*)fileMEAll->Get("T1");

  trees.clear();
  trees.push_back(treeEEAll);
  trees.push_back(treeMMAll);
  trees.push_back(treeEMAll);
  trees.push_back(treeMEAll);


  // yields at WW level: tagged events
  float yield_WWSel_denom[4]; // [icha]
  float yield_WWSel_denom_staterr[4]; // [icha]

  float yield_WWSel_denom_tot = 0.;
  float yield_WWSel_denom_tot_staterr = 0.;

  float yield_WWSel_1jnum_tot = 0.;
  float yield_WWSel_1jnum_tot_staterr = 0.;

  TH1F *histo2 = new TH1F("histo2","",100,0.,180.); // for the denominator

  for(int icha=0; icha<4; icha++) {
    // step[9] is all before b-veto and njet cut 

    TString fpCut = TString("(") + TString("step[9] && ") + TString(njetscut) + TString(")") + TString("*weightFP*hlt");            
    TString fpCutStatErr = TString("(") + TString("step[9] && ") + TString(njetscut) + TString(")") + TString("*weightStatFP*hlt"); 
    // for MC closure test
    // TString fpCut = TString("(") + TString("step[9] && ") + TString(njetscut) + TString(")") + TString("*1.545*weightFP*baseW");               
    // TString fpCutStatErr = TString("(") + TString("step[9] && ") + TString(njetscut) + TString(")") + TString("*1.545*weightStatFP*baseW");    
    
    trees[icha]->Project("histo2","dphill",fpCut);
    yield_WWSel_denom[icha] = histo2->Integral();
    yield_WWSel_denom_tot += yield_WWSel_denom[icha];
    
    trees[icha]->Project("histo2","dphill",fpCutStatErr);
    yield_WWSel_denom_staterr[icha] = sqrt(histo2->Integral());
    yield_WWSel_denom_tot_staterr += pow(yield_WWSel_denom_staterr[icha],2);

    // as the tagged 1 jet bin
    // for the 1 jet numerator: one event enter the top control region if the leading jet is btagged && all the rest are not 

    fpCut = TString("(") + TString("step[9] && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && ") + TString(njetscut) + TString(")") + TString("*weightFP*hlt");
    fpCutStatErr = TString("(") + TString("step[9] && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && ") + TString(njetscut) + TString(")") + TString("*weightStatFP*hlt");
    // for MC closure test
    // fpCut = TString("(") + TString("step[9] && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && ") + TString(njetscut) + TString(")") + TString("*1.545*weightFP*baseW");
    // fpCutStatErr = TString("(") + TString("step[9] && leadingJetBTagTrackCount>2.1 && subleadingJetBTagTrackCount<=2.1 && ") + TString(njetscut) + TString(")") + TString("*1.545*weightStatFP*baseW");

    trees[icha]->Project("histo2","dphill",fpCut);
    yield_WWSel_1jnum_tot += histo2->Integral();
    
    trees[icha]->Project("histo2","dphill",fpCutStatErr);
    yield_WWSel_1jnum_tot_staterr += pow(histo2->Integral(),2);
  }
  
  yield_WWSel_denom_tot_staterr = sqrt(yield_WWSel_denom_tot_staterr);
  yield_WWSel_1jnum_tot_staterr = sqrt(yield_WWSel_1jnum_tot_staterr);

  float numerator = 0;
  if(njets==0) numerator = yield_WWSel_denom_tot - yield_WWSel_tot; // events that are btagged
  else numerator = yield_WWSel_1jnum_tot_staterr;

  float btageff = numerator / yield_WWSel_denom_tot;
  float btageff_err = sqrt(btageff*(1-btageff)/yield_WWSel_denom_tot);
  std::cout << "On the TF sample: btagging (IP + soft mu) efficiency = " << btageff << " +/- " << btageff_err << std::endl;

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
