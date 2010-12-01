#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

// parameters to configure
int exampleHiggsMass = 200;
TString castordir_data("/cmsrm/pc21_2/emanuele/data/Higgs3.9.X/Data_HiggsRev_V8/");
TString castordir_mc("/cmsrm/pc23_2/crovelli/data/Higgs3.9.X/mc_higgsReview_v8/");

using namespace std;

enum { ee=0, mm=1, em=2 };

string preSelCuts[10];
string fullSelCuts[24];

int UsePreSelCuts[3][10];
int UseCuts[3][24];

float H_preSel[11];
float Wj_preSel[11];
float ttj_preSel[11];
float SingleTop_preSel[11];
float Zj_preSel[11];
float WW_preSel[11];
float ZZ_preSel[11];
float WZ_preSel[11];
float Wgamma_preSel[11];
float QCDem_preSel[11];
float QCDbc_preSel[11];
float Photj_preSel[11];
float QCDmu_preSel[11];
float ttbar_preSel[11];
float data_preSel[11];

float H_fullSel[24];
float Wj_fullSel[24];
float ttj_fullSel[24];
float SingleTop_fullSel[24];
float Zj_fullSel[24];
float WW_fullSel[24];
float ZZ_fullSel[24];
float WZ_fullSel[24];
float Wgamma_fullSel[24];
float QCDem_fullSel[24];
float QCDbc_fullSel[24];
float Photj_fullSel[24];
float QCDmu_fullSel[24];
float ttbar_fullSel[24];
float data_fullSel[24];

float H_eff_preSel[11];
float Wj_eff_preSel[11];
float ttj_eff_preSel[11];
float SingleTop_eff_preSel[11];
float Zj_eff_preSel[11];
float WW_eff_preSel[11];
float ZZ_eff_preSel[11];
float WZ_eff_preSel[11];
float Wgamma_eff_preSel[11];
float QCDem_eff_preSel[11];
float QCDbc_eff_preSel[11];
float Photj_eff_preSel[11];
float QCDmu_eff_preSel[11];
float ttbar_eff_preSel[11];

float H_eff_fullSel[24];
float Wj_eff_fullSel[24];
float ttj_eff_fullSel[24];
float SingleTop_eff_fullSel[24];
float Zj_eff_fullSel[24];
float WW_eff_fullSel[24];
float ZZ_eff_fullSel[24];
float WZ_eff_fullSel[24];
float Wgamma_eff_fullSel[24];
float QCDem_eff_fullSel[24];
float QCDbc_eff_fullSel[24];
float Photj_eff_fullSel[24];
float QCDmu_eff_fullSel[24];
float ttbar_eff_fullSel[24];

float H_finaleff_preSel;    
float Wj_finaleff_preSel;
float ttj_finaleff_preSel;
float SingleTop_finaleff_preSel;
float Zj_finaleff_preSel;
float WW_finaleff_preSel;
float ZZ_finaleff_preSel;
float WZ_finaleff_preSel;
float Wgamma_finaleff_preSel;
float QCDem_finaleff_preSel;
float QCDbc_finaleff_preSel;
float Photj_finaleff_preSel;
float QCDmu_finaleff_preSel;
float ttbar_finaleff_preSel;

float H_finaleff_fullSel;    
float Wj_finaleff_fullSel;
float ttj_finaleff_fullSel;
float SingleTop_finaleff_fullSel;
float Zj_finaleff_fullSel;
float WW_finaleff_fullSel;
float ZZ_finaleff_fullSel;
float WZ_finaleff_fullSel;
float Wgamma_finaleff_fullSel;
float QCDem_finaleff_fullSel;
float QCDbc_finaleff_fullSel;
float Photj_finaleff_fullSel;
float QCDmu_finaleff_fullSel;
float ttbar_finaleff_fullSel;

float H_finaleff;    
float Wj_finaleff;
float ttj_finaleff;
float SingleTop_finaleff;
float Zj_finaleff;
float WW_finaleff;
float ZZ_finaleff;
float WZ_finaleff;
float Wgamma_finaleff;
float QCDem_finaleff;
float QCDbc_finaleff;
float Photj_finaleff;
float QCDmu_finaleff;
float ttbar_finaleff;

float H_final[3];
float Wj_final[3];
float ttj_final[3];
float SingleTop_final[3];
float Zj_final[3];
float WW_final[3];
float ZZ_final[3];
float WZ_final[3];
float Wgamma_final[3];
float data_final[3];

// xsections
std::map<int,float> Higgs_xsec_masses;

float Wgamma_xsec = 41.76;
float Wlnu_xsec = 31314./3. * 0.742; // NLO * filtereff (BR W->lnu included);
float ZjetsLoMass_xsec = 2659./3.;  // production page...
float ZjetsHiMass_xsec = 4998./3.; // MCFM mll > 20 GeV
float TTjets_xsec = 157.5;
float ggWW_xsec = 0.1538; // gg->WW->4l
float WW_xsec = 4.50347; // WW_2l2nu
float WZ_xsec = 0.599442; // WZ_3l
//float ZZ_xsec = 0.25252; // ZZ_2l2nu
float ZZ_xsec = 5.9; // inclusive
float TTbar_xsec = 165;

// here xsec = x-sec * filter_eff (pb)
float QCD_EMenriched_Pt20to30_xsec = 235.5E+06 * 0.0073;
float QCD_EMenriched_Pt30to80_xsec = 59.3E+06 * 0.059;
float QCD_EMenriched_Pt80to170_xsec = 0.906E+06 * 0.148;
float QCD_BCtoE_Pt20to30_xsec = 235.5E+06 * 0.00046;
float QCD_BCtoE_Pt30to80_xsec = 59.3E+06 * 0.00234;
float QCD_BCtoE_Pt80to170_xsec = 0.906E+06 * 0.0104;

// xsec (pb)
float SingleTopS_xsec = 4.21 * (0.1080*3);
float SingleTopT_xsec = 64.6 * (0.1080*3);
float SingleTopTW_xsec = 10.6 * (0.1080*3);

// xsec (pb)
float PhotonJet_Pt0to15_xsec = 84.46E+06;
float PhotonJet_Pt15to20_xsec = 114700;
float PhotonJet_Pt20to30_xsec = 57180;
float PhotonJet_Pt30to50_xsec = 16520;
float PhotonJet_Pt50to80_xsec = 2723;
float PhotonJet_Pt80to120_xsec = 446.2;
float PhotonJet_Pt120to170_xsec = 84.43;
float PhotonJet_Pt170to300_xsec = 22.55;
float PhotonJet_Pt300to500_xsec = 1.545;
float PhotonJet_Pt500toInf_xsec = 0.0923;

float InclusiveMu15_xsec = 48.44*0.00176*1.0E+09; // ppMuX

string sampleNames[21];

void computeYields(float lumi, const char* finalstate, int mass=0) {

  TChain *chains_preSel[21];
  TChain *chains_fullSel[21];

  for(int isample=0; isample<21; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // signal
  sampleNames[0] = "Higgs";
  // backgrounds
  sampleNames[1] = "TTbar+jets";
  sampleNames[2] = "Z(ee) 10-20 GeV";
  sampleNames[3] = "Z(mumu) 10-20 GeV";
  sampleNames[4] = "Z(tautau) 10-20 GeV ";
  sampleNames[5] = "Z(ee) >20 GeV";
  sampleNames[6] = "Z(mumu) >20 GeV";
  sampleNames[7] = "Z(tautau) >20 GeV ";
  sampleNames[8] = "WW";
  sampleNames[9] = "gg->WW";
  sampleNames[10] = "ZZ";
  sampleNames[11] = "WZ";
  sampleNames[12] = "Wgamma";
  sampleNames[13] = "SingleTop_sChannel";
  sampleNames[14] = "SingleTop_tChannel";
  sampleNames[15] = "SingleTop_tWChannel";
  sampleNames[16] = "Wenu";
  sampleNames[17] = "Wmunu";
  sampleNames[18] = "Wtaunu";
  sampleNames[19] = "data Sep13Rereco";
  sampleNames[20] = "data 11To40 /pb";

  float Higgs_xsec;
  Higgs_xsec_masses.insert(std::make_pair(120,0.247143));
  Higgs_xsec_masses.insert(std::make_pair(130,0.452859));
  Higgs_xsec_masses.insert(std::make_pair(140,0.64926));
  Higgs_xsec_masses.insert(std::make_pair(150,0.787871));
  Higgs_xsec_masses.insert(std::make_pair(160,0.897043));
  Higgs_xsec_masses.insert(std::make_pair(170,0.808914));
  Higgs_xsec_masses.insert(std::make_pair(200,0.422487));
  Higgs_xsec_masses.insert(std::make_pair(300,0.181931));
  Higgs_xsec_masses.insert(std::make_pair(400,0.125106));

  if(mass==0) { // use the default mass to print the cut-by cut table: the one pointed by results/ dir
    Higgs_xsec = Higgs_xsec_masses[exampleHiggsMass] * 4./9.; // 4/9 because we are considering only the samples containing e-mu combinations.

    char dir_mc[1000];
    sprintf(dir_mc,"%s/OptimMH%d/mc_higgsReview_v8/OptimMH%d/",castordir_mc.Data(),exampleHiggsMass,exampleHiggsMass);
    char HiggsSample[500];
    sprintf(HiggsSample,"HiggsWW/GluGluToHToWWTo2L2Nu_M-%d/*Counters.root",exampleHiggsMass);

    // signal
    chains_preSel[0]->Add(TString(dir_mc)+TString(HiggsSample));       
    // backgrounds
    chains_preSel[1]->Add(TString(dir_mc)+TString("/TTbar/TTJets_TuneD6T/*Counters.root"));       
    chains_preSel[2]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[3]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[4]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[5]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_preSel[6]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_preSel[7]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*Counters.root"));
    chains_preSel[8]->Add(TString(dir_mc)+TString("/DiBosons/WWTo2L2Nu_TuneZ2/*Counters.root"));
    chains_preSel[9]->Add(TString(dir_mc)+TString("/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*Counters.root"));
    chains_preSel[10]->Add(TString(dir_mc)+TString("/DiBosons/ZZtoAnything_TuneZ2/*Counters.root"));   
    chains_preSel[11]->Add(TString(dir_mc)+TString("/DiBosons/WZTo3LNu_TuneZ2/*Counters.root"));
    //    chains_preSel[12]->Add(TString(dir_mc)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_preSel[13]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[14]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[15]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[16]->Add(TString(dir_mc)+TString("/WPYTHIA/WToENu_TuneZ2/*Counters.root"));
    chains_preSel[17]->Add(TString(dir_mc)+TString("/WPYTHIA/WToMuNu_TuneZ2/*Counters.root"));
    chains_preSel[18]->Add(TString(dir_mc)+TString("/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));

    // data
    char dir_data[1000];
    sprintf(dir_data,"%s/OptimMH%d/Data_HiggsRev_V8/OptimMH%d/Data7TeV/",castordir_data.Data(),exampleHiggsMass,exampleHiggsMass);
    chains_preSel[19]->Add(TString(dir_data)+TString("/dataset_eg_Sep3rdReReco/*Counters.root"));
    chains_preSel[20]->Add(TString(dir_data)+TString("/PDElectron_11pbTo40pb/*Counters.root"));

    // signal
    chains_fullSel[0]->Add(TString(dir_mc)+TString(HiggsSample));       
    // backgrounds
    chains_fullSel[1]->Add(TString(dir_mc)+TString("/TTbar/TTJets_TuneD6T/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir_mc)+TString("/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*Counters.root"));
    chains_fullSel[8]->Add(TString(dir_mc)+TString("/DiBosons/WWTo2L2Nu_TuneZ2/*Counters.root"));
    chains_fullSel[9]->Add(TString(dir_mc)+TString("/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*Counters.root"));
    chains_fullSel[10]->Add(TString(dir_mc)+TString("/DiBosons/ZZtoAnything_TuneZ2/*Counters.root"));   
    chains_fullSel[11]->Add(TString(dir_mc)+TString("/DiBosons/WZTo3LNu_TuneZ2/*Counters.root"));
    //    chains_fullSel[12]->Add(TString(dir_mc)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_fullSel[13]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[14]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[15]->Add(TString(dir_mc)+TString("/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[16]->Add(TString(dir_mc)+TString("/WPYTHIA/WToENu_TuneZ2/*Counters.root"));
    chains_fullSel[17]->Add(TString(dir_mc)+TString("/WPYTHIA/WToMuNu_TuneZ2/*Counters.root"));
    chains_fullSel[18]->Add(TString(dir_mc)+TString("/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));

    // data
    chains_fullSel[19]->Add(TString(dir_data)+TString("/dataset_eg_Sep3rdReReco/*Counters.root"));
    chains_fullSel[20]->Add(TString(dir_data)+TString("/PDElectron_11pbTo40pb/*Counters.root"));

  } else {

    Higgs_xsec = Higgs_xsec_masses[mass] * 4./9.; // 4/9 because we are considering only the samples containing e-mu combinations.

    char dir[1000];
    sprintf(dir,"%s/OptimMH%d/mc_higgsReview_v8/OptimMH%d/",castordir_mc.Data(),mass,mass);
    char higgsDir[100];
    sprintf(higgsDir,"HiggsWW/GluGluToHToWWTo2L2Nu_M-%d/",mass);

    // signal
    chains_preSel[0]->Add(TString(dir)+TString(higgsDir)+TString("*Counters.root"));
    // backgrounds
    chains_preSel[1]->Add(TString(dir)+TString("/TTbar/TTJets_TuneD6T/*Counters.root"));       
    chains_preSel[2]->Add(TString(dir)+TString("/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[3]->Add(TString(dir)+TString("/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[4]->Add(TString(dir)+TString("/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*Counters.root"));       
    chains_preSel[5]->Add(TString(dir)+TString("/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_preSel[6]->Add(TString(dir)+TString("/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_preSel[7]->Add(TString(dir)+TString("/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*Counters.root"));
    chains_preSel[8]->Add(TString(dir)+TString("/DiBosons/WWTo2L2Nu_TuneZ2/*Counters.root"));
    chains_preSel[9]->Add(TString(dir)+TString("/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*Counters.root"));
    chains_preSel[10]->Add(TString(dir)+TString("/DiBosons/ZZtoAnything_TuneZ2/*Counters.root"));   
    chains_preSel[11]->Add(TString(dir)+TString("/DiBosons/WZTo3LNu_TuneZ2/*Counters.root"));
    //    chains_preSel[12]->Add(TString(dir)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_preSel[13]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[14]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[15]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_preSel[16]->Add(TString(dir)+TString("/WPYTHIA/WToENu_TuneZ2/*Counters.root"));
    chains_preSel[17]->Add(TString(dir)+TString("/WPYTHIA/WToMuNu_TuneZ2/*Counters.root"));
    chains_preSel[18]->Add(TString(dir)+TString("/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));


    // data
    char dir_data[1000];
    sprintf(dir_data,"%s/OptimMH%d/Data_HiggsRev_V8/OptimMH%d/Data7TeV/",castordir_data.Data(),mass,mass);
    chains_preSel[19]->Add(TString(dir_data)+TString("/dataset_eg_Sep3rdReReco/*Counters.root"));
    chains_preSel[20]->Add(TString(dir_data)+TString("/PDElectron_11pbTo40pb/*Counters.root"));

    // signal
    chains_fullSel[0]->Add(TString(dir)+TString(higgsDir)+TString("*Counters.root"));
    // backgrounds
    chains_fullSel[1]->Add(TString(dir)+TString("/TTbar/TTJets_TuneD6T/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir)+TString("/ZPYTHIA/DYToEE_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir)+TString("/ZPYTHIA/DYToMuMu_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir)+TString("/ZPYTHIA/DYToTauTau_M-10To20_TuneZ2/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir)+TString("/ZPYTHIA/DYToEE_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir)+TString("/ZPYTHIA/DYToMuMu_M-20_CT10_TuneZ2_PU/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir)+TString("/ZPYTHIA/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/*Counters.root"));
    chains_fullSel[8]->Add(TString(dir)+TString("/DiBosons/WWTo2L2Nu_TuneZ2/*Counters.root"));
    chains_fullSel[9]->Add(TString(dir)+TString("/DiBosons/GluGluToWWTo4L_TuneZ2_PU/*Counters.root"));
    chains_fullSel[10]->Add(TString(dir)+TString("/DiBosons/ZZtoAnything_TuneZ2/*Counters.root"));   
    chains_fullSel[11]->Add(TString(dir)+TString("/DiBosons/WZTo3LNu_TuneZ2/*Counters.root"));
    //    chains_fullSel[12]->Add(TString(dir)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_fullSel[13]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[14]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[15]->Add(TString(dir)+TString("/SingleTop/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[16]->Add(TString(dir)+TString("/WPYTHIA/WToENu_TuneZ2/*Counters.root"));
    chains_fullSel[17]->Add(TString(dir)+TString("/WPYTHIA/WToMuNu_TuneZ2/*Counters.root"));
    chains_fullSel[18]->Add(TString(dir)+TString("/WPYTHIA/WToTauNu_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));

    // data 
    chains_fullSel[19]->Add(TString(dir_data)+TString("/dataset_eg_Sep3rdReReco/*Counters.root"));
    chains_fullSel[20]->Add(TString(dir_data)+TString("/PDElectron_11pbTo40pb/*Counters.root"));

  }

  float nPreSelTot[10][21];
  float nFullSelTot[24][21];

  for(int isample=0; isample<21; isample++) {
    for(int icut=0; icut<10; icut++) { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<28; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preSelections
  int nCutsAnaPre  = 10;
  int nCutsAnaFull = 24;
  for(int isample=0; isample<21; isample++) {

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
      Long64_t nb;
      nb = chains_preSel[isample]->GetEntry(jentry);   
      // nCutsAnaPre = nCutsPre;      
      for(int icut=0; icut<nCutsPre; icut++) nPreSelTot[icut][isample] += nSelPre[icut];	
    }
    
    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2;
      nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      // nCutsAnaFull = nCutsFull;
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at preSelection level
  for(int icut=0; icut<nCutsAnaPre; icut++) {

    // numbers
    // N = L * x-sec * eff (eff = N_fin / N_ini)
    if(nPreSelTot[0][0]>0) { 
      H_preSel[icut]      = lumi * Higgs_xsec  * nPreSelTot[icut][0]/nPreSelTot[0][0];
      // Wj_preSel[icut]     = lumi * Wjets_xsec  * nPreSelTot[icut][1]/nPreSelTot[0][1];
      float W_tmp=0.;
      W_tmp += lumi * Wlnu_xsec  * nPreSelTot[icut][16]/nPreSelTot[0][16]; // W->enu
      W_tmp += lumi * Wlnu_xsec  * nPreSelTot[icut][17]/nPreSelTot[0][17]; // W->munu
      W_tmp += lumi * Wlnu_xsec  * nPreSelTot[icut][18]/nPreSelTot[0][18]; // W->taunu
      Wj_preSel[icut]     = W_tmp;
      ttj_preSel[icut]    = lumi * TTjets_xsec * nPreSelTot[icut][1]/nPreSelTot[0][1];
      float Z_tmp=0.;
      Z_tmp += lumi * ZjetsLoMass_xsec  * nPreSelTot[icut][2]/nPreSelTot[0][2]; // Z->ee
      Z_tmp += lumi * ZjetsLoMass_xsec  * nPreSelTot[icut][3]/nPreSelTot[0][3]; // Z->mumu
      Z_tmp += lumi * ZjetsLoMass_xsec  * nPreSelTot[icut][4]/nPreSelTot[0][4]; // Z->tautau
      Z_tmp += lumi * ZjetsHiMass_xsec  * nPreSelTot[icut][5]/nPreSelTot[0][5]; // Z->ee
      Z_tmp += lumi * ZjetsHiMass_xsec  * nPreSelTot[icut][6]/nPreSelTot[0][6]; // Z->mumu
      Z_tmp += lumi * ZjetsHiMass_xsec  * nPreSelTot[icut][7]/nPreSelTot[0][7]; // Z->tautau
      Zj_preSel[icut]     = Z_tmp;
      float WW_tmp=0.;
      WW_tmp += lumi * WW_xsec     * nPreSelTot[icut][8]/nPreSelTot[0][8];
      WW_tmp += lumi * ggWW_xsec     * nPreSelTot[icut][9]/nPreSelTot[0][9];
      WW_preSel[icut]     = WW_tmp;
      ZZ_preSel[icut]     = lumi * ZZ_xsec     * nPreSelTot[icut][10]/nPreSelTot[0][10];
      WZ_preSel[icut]     = lumi * WZ_xsec     * nPreSelTot[icut][11]/nPreSelTot[0][11];
      Wgamma_preSel[icut] = lumi * Wgamma_xsec * nPreSelTot[icut][12]/nPreSelTot[0][12];
      
      float singletop_tmp=0.;
      singletop_tmp += lumi * SingleTopS_xsec * nPreSelTot[icut][13]/nPreSelTot[0][13];
      singletop_tmp += lumi * SingleTopT_xsec * nPreSelTot[icut][14]/nPreSelTot[0][14];
      singletop_tmp += lumi * SingleTopTW_xsec * nPreSelTot[icut][15]/nPreSelTot[0][15];
      SingleTop_preSel[icut] = singletop_tmp;
      
      // data yields
      data_preSel[icut] = nPreSelTot[icut][19] + nPreSelTot[icut][20];

    }

    
    // efficiencies
    if(icut>0 && nPreSelTot[icut-1][0]>0) H_eff_preSel[icut]      = nPreSelTot[icut][0]/nPreSelTot[icut-1][0];
    else H_eff_preSel[icut] = 0.0;
    if(icut>0 && Wj_preSel[icut-1]>0) Wj_eff_preSel[icut]     = Wj_preSel[icut]/Wj_preSel[icut-1];
    else Wj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][1]>0) ttj_eff_preSel[icut]    = nPreSelTot[icut][1]/nPreSelTot[icut-1][1];
    else ttj_eff_preSel[icut] = 0.0;
    if(icut>0 && Zj_preSel[icut-1]>0) Zj_eff_preSel[icut]     = Zj_preSel[icut]/Zj_preSel[icut-1];
    else Zj_eff_preSel[icut] = 0.0;
    if(icut>0 && WW_preSel[icut-1]>0) WW_eff_preSel[icut]     = WW_preSel[icut]/WW_preSel[icut-1];
    else WW_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][5]>0) ZZ_eff_preSel[icut]     = nPreSelTot[icut][10]/nPreSelTot[icut-1][10];
    else ZZ_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][6]>0) WZ_eff_preSel[icut]     = nPreSelTot[icut][11]/nPreSelTot[icut-1][11];
    else WZ_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][7]>0) Wgamma_eff_preSel[icut] = nPreSelTot[icut][12]/nPreSelTot[icut-1][12];
    else Wgamma_eff_preSel[icut] = 0.0;
    if(icut>0 && SingleTop_preSel[icut-1]>0) SingleTop_eff_preSel[icut] = SingleTop_preSel[icut]/SingleTop_preSel[icut-1];
    else SingleTop_eff_preSel[icut] = 0.0;

  }

  // final efficiency at preSelection level
  if(nPreSelTot[0][0]>0) H_finaleff_preSel      = nPreSelTot[nCutsAnaPre-2][0]/nPreSelTot[0][0];
  else H_finaleff_preSel = 0.0;
  if(Wj_preSel[0]>0) Wj_finaleff_preSel     = Wj_preSel[nCutsAnaPre-2]/Wj_preSel[0];
  else Wj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][1]>0) ttj_finaleff_preSel    = nPreSelTot[nCutsAnaPre-2][1]/nPreSelTot[0][1];
  else ttj_finaleff_preSel = 0.0;
  if(Zj_preSel[0]>0) Zj_finaleff_preSel     = Zj_preSel[nCutsAnaPre-2]/Zj_preSel[0];
  else Zj_finaleff_preSel = 0.0;
  if(WW_preSel[0]>0) WW_finaleff_preSel     = WW_preSel[nCutsAnaPre-2]/WW_preSel[0];
  else WW_finaleff_preSel = 0.0;
  if(nPreSelTot[0][5]>0) ZZ_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][10]/nPreSelTot[0][10];
  else ZZ_finaleff_preSel = 0.0;
  if(nPreSelTot[0][6]>0) WZ_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][11]/nPreSelTot[0][11];
  else WZ_finaleff_preSel = 0.0;
  if(nPreSelTot[0][7]>0) Wgamma_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][12]/nPreSelTot[0][12];
  else Wgamma_finaleff_preSel = 0.0;
  if(SingleTop_preSel[0]>0) SingleTop_finaleff_preSel = SingleTop_preSel[nCutsAnaPre-2]/SingleTop_preSel[0];
  else SingleTop_finaleff_preSel = 0.0;

  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) { 
      H_fullSel[icut]      = lumi * Higgs_xsec  * nFullSelTot[icut][0]/nPreSelTot[0][0];
      float W_tmp=0.;
      W_tmp += lumi * Wlnu_xsec * nFullSelTot[icut][16]/nPreSelTot[0][16]; // W->enu
      W_tmp += lumi * Wlnu_xsec * nFullSelTot[icut][17]/nPreSelTot[0][17]; // W->munu
      W_tmp += lumi * Wlnu_xsec * nFullSelTot[icut][18]/nPreSelTot[0][18]; // W->taunu
      Wj_fullSel[icut]     = W_tmp;
      ttj_fullSel[icut]    = lumi * TTjets_xsec * nFullSelTot[icut][1]/nPreSelTot[0][1];
      float Z_tmp=0.;
      Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][2]/nPreSelTot[0][2]; // Z->ee
      Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][3]/nPreSelTot[0][3]; // Z->mumu
      Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][4]/nPreSelTot[0][4]; // Z->tautau
      Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][5]/nPreSelTot[0][5]; // Z->ee
      Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][6]/nPreSelTot[0][6]; // Z->mumu
      Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][7]/nPreSelTot[0][7]; // Z->tautau
      Zj_fullSel[icut]     = Z_tmp;
      float WW_tmp=0.;
      WW_tmp += lumi * WW_xsec     * nFullSelTot[icut][8]/nPreSelTot[0][8];
      WW_tmp += lumi * ggWW_xsec   * nFullSelTot[icut][9]/nPreSelTot[0][9];
      WW_fullSel[icut]     = WW_tmp;
      ZZ_fullSel[icut]     = lumi * ZZ_xsec     * nFullSelTot[icut][10]/nPreSelTot[0][10];
      WZ_fullSel[icut]     = lumi * WZ_xsec     * nFullSelTot[icut][11]/nPreSelTot[0][11];
      Wgamma_fullSel[icut] = lumi * Wgamma_xsec * nFullSelTot[icut][12]/nPreSelTot[0][12];

      float singletop_tmp=0.;
      singletop_tmp += lumi * SingleTopS_xsec * nFullSelTot[icut][13]/nPreSelTot[0][13];
      singletop_tmp += lumi * SingleTopT_xsec * nFullSelTot[icut][14]/nPreSelTot[0][14];
      singletop_tmp += lumi * SingleTopTW_xsec * nFullSelTot[icut][15]/nPreSelTot[0][15];
      SingleTop_fullSel[icut] = singletop_tmp;

      // data
      data_fullSel[icut] = nFullSelTot[icut][19] + nFullSelTot[icut][20];

    }

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else H_eff_fullSel[icut] = 0.0;
    if(icut>0 && Wj_fullSel[icut-1]>0) Wj_eff_fullSel[icut]     = Wj_fullSel[icut]/Wj_fullSel[icut-1];
    else Wj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][1]>0) ttj_eff_fullSel[icut]    = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else ttj_eff_fullSel[icut] = 0.0;
    if(icut>0 && Zj_fullSel[icut-1]>0) Zj_eff_fullSel[icut]     = Zj_fullSel[icut]/Zj_fullSel[icut-1];
    else Zj_eff_fullSel[icut] = 0.0;
    if(icut>0 && WW_fullSel[icut-1]>0) WW_eff_fullSel[icut]     = WW_fullSel[icut]/WW_fullSel[icut-1];
    else WW_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][5]>0) ZZ_eff_fullSel[icut]     = nFullSelTot[icut][10]/nFullSelTot[icut-1][10];
    else ZZ_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][6]>0) WZ_eff_fullSel[icut]     = nFullSelTot[icut][11]/nFullSelTot[icut-1][11];
    else WZ_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][7]>0) Wgamma_eff_fullSel[icut] = nFullSelTot[icut][12]/nFullSelTot[icut-1][12];
    else Wgamma_eff_fullSel[icut] = 0.0;
    if(icut>0 && SingleTop_fullSel[icut-1]>0) SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[icut-1];
    else SingleTop_eff_fullSel[icut] = 0.0;

    if(icut==0) { 
      H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nPreSelTot[nCutsAnaPre-2][0];
      Wj_eff_fullSel[icut]     = Wj_fullSel[icut]/Wj_preSel[nCutsAnaPre-2];
      ttj_eff_fullSel[icut]    = nFullSelTot[icut][1]/nPreSelTot[nCutsAnaPre-2][1];
      Zj_eff_fullSel[icut]     = Zj_fullSel[icut]/Zj_preSel[nCutsAnaPre-2];
      WW_eff_fullSel[icut]     = WW_fullSel[icut]/WW_fullSel[nCutsAnaPre-2];
      ZZ_eff_fullSel[icut]     = nFullSelTot[icut][10]/nPreSelTot[nCutsAnaPre-2][10];
      WZ_eff_fullSel[icut]     = nFullSelTot[icut][11]/nPreSelTot[nCutsAnaPre-2][11];
      Wgamma_eff_fullSel[icut] = nFullSelTot[icut][12]/nPreSelTot[nCutsAnaPre-2][12];
      SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_preSel[nCutsAnaPre-2];
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0) H_finaleff_fullSel      = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0];
  else H_finaleff_fullSel = 0.0;
  if(Wj_fullSel[0]>0) Wj_finaleff_fullSel = Wj_fullSel[nCutsAnaFull-4]/Wj_fullSel[0];
  else Wj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) ttj_finaleff_fullSel    = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1];
  else ttj_finaleff_fullSel = 0.0;
  if(Zj_fullSel[0]>0) Zj_finaleff_fullSel = Zj_fullSel[nCutsAnaFull-4]/Zj_fullSel[0];
  else Zj_finaleff_fullSel = 0.0;
  if(WW_fullSel[0]>0) WW_finaleff_fullSel     = WW_fullSel[nCutsAnaFull-4]/WW_fullSel[0];
  else WW_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][5]>0) ZZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][10]/nFullSelTot[0][10];
  else ZZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][6]>0) WZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][11]/nFullSelTot[0][11];
  else WZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][7]>0) Wgamma_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][12]/nFullSelTot[0][12];
  else Wgamma_finaleff_fullSel = 0.0;
  if(SingleTop_fullSel[0]>0) SingleTop_finaleff_fullSel = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_fullSel[0];
  else SingleTop_finaleff_fullSel = 0.0;

  // final efficiency combining pre and full selections
  if(nPreSelTot[0][0]>0) H_finaleff      = nFullSelTot[nCutsAnaFull-4][0]/nPreSelTot[0][0];
  else H_finaleff = 0.0;
  if(Wj_preSel[0]>0) Wj_finaleff = Wj_fullSel[nCutsAnaFull-4]/Wj_preSel[0];
  else Wj_finaleff = 0.0;
  if(nPreSelTot[0][1]>0) ttj_finaleff    = nFullSelTot[nCutsAnaFull-4][1]/nPreSelTot[0][1];
  else ttj_finaleff = 0.0;
  if(Zj_preSel[0]>0) Zj_finaleff = Zj_fullSel[nCutsAnaFull-4]/Zj_preSel[0];
  else Zj_finaleff = 0.0;
  if(WW_preSel[0]>0) WW_finaleff     = WW_fullSel[nCutsAnaFull-4]/WW_preSel[0];
  else WW_finaleff = 0.0;
  if(nPreSelTot[0][5]>0) ZZ_finaleff     = nFullSelTot[nCutsAnaFull-4][10]/nPreSelTot[0][10];
  else ZZ_finaleff = 0.0;
  if(nPreSelTot[0][6]>0) WZ_finaleff     = nFullSelTot[nCutsAnaFull-4][11]/nPreSelTot[0][11];
  else WZ_finaleff = 0.0;
  if(nPreSelTot[0][7]>0) Wgamma_finaleff = nFullSelTot[nCutsAnaFull-4][12]/nPreSelTot[0][12];
  else Wgamma_finaleff = 0.0;
  if(SingleTop_preSel[0]>0) SingleTop_finaleff = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_preSel[0];
  else SingleTop_finaleff = 0.0;

  cout << "\n\nPROCESSED EVENTS:" << endl;
  for(int i=0; i<21; i++) {
    cout << sampleNames[i] << "\t" << nPreSelTot[0][i] << endl;
  }

}

void setupCuts() {
  
  for(int ichan=0; ichan<3; ichan++) {
    for(int i=0; i<10; i++) {
      UsePreSelCuts[ichan][i] = 1;
    }
  }

  for(int ichan=0; ichan<3; ichan++) {
    for(int i=0; i<24; i++) {
      UseCuts[ichan][i] = 1;
    }
  }

  // unusued ee cuts
  UsePreSelCuts[ee][1] = UsePreSelCuts[ee][7] = UsePreSelCuts[ee][8] = 0; 
  UseCuts[ee][5] = UseCuts[ee][6] = UseCuts[ee][7] = UseCuts[ee][8] = UseCuts[ee][21] = UseCuts[ee][11] = UseCuts[ee][15] = 0; // muon iso

  // unusued mm cuts
  UsePreSelCuts[mm][1] = UsePreSelCuts[mm][7] = UsePreSelCuts[mm][8] = 0; 
  UseCuts[mm][4] = UseCuts[mm][10] = UseCuts[mm][5] = UseCuts[mm][6] = UseCuts[mm][7] = UseCuts[mm][21] = UseCuts[mm][11] = UseCuts[mm][15] = 0; // ele ID and conv. rej. and separate mu iso

  // unusued em cuts
  UsePreSelCuts[em][1] = UsePreSelCuts[em][7] = UsePreSelCuts[em][8] = UseCuts[em][13] = UseCuts[em][21] = UseCuts[em][11] = UseCuts[em][15] = 0;
  UseCuts[em][5] = UseCuts[em][6] = UseCuts[em][7] = 0; // separate mu iso 
  
  preSelCuts[0]="event";
  preSelCuts[1]="MCtruth";
  preSelCuts[2]="trigger";
  preSelCuts[3]="$\\geq 2$ leptons";
  preSelCuts[4]="acceptance";
  preSelCuts[5]="$p_T^{max}$";
  preSelCuts[6]="$p_T^{min}$";
  preSelCuts[7]="MET preSelection";
  preSelCuts[8]="$m_{ll}>12$ GeV";
  preSelCuts[9]="final preSel.";

  fullSelCuts[0]="channel preSel.";
  fullSelCuts[1]="sel $p_T^{max}$";
  fullSelCuts[2]="sel $p_T^{min}$";
  fullSelCuts[3]="e/$\\mu$ d0";
  fullSelCuts[4]="e isolation";
  fullSelCuts[5]="$\\mu$ tracker Iso";
  fullSelCuts[6]="$\\mu$ HCAL Iso";
  fullSelCuts[7]="$\\mu$ ECAL Iso";
  fullSelCuts[8]="$\\mu$ isolation";
  fullSelCuts[9]="e/$\\mu$ ID";
  fullSelCuts[10]="conv. rej.";
  fullSelCuts[11]="$MET>20$ GeV";
  fullSelCuts[12]="$m_{ll}$";
  fullSelCuts[13]="$|m_{ll}-m_Z|>15$ GeV";
  fullSelCuts[14]="proj. MET";
  fullSelCuts[15]="MET/$p_T^{ll}$";
  fullSelCuts[16]="jet veto";
  fullSelCuts[17]="anti b-tag";
  fullSelCuts[18]="$\\mu^{soft}$ veto";
  fullSelCuts[19]="extra lepton veto";
  fullSelCuts[20]="$\\Delta \\phi$";
  fullSelCuts[21]="final";
  fullSelCuts[22]="1 jets";
  fullSelCuts[23]="$>1$ jets";

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


  /// ==============  print detailed breakdown  ================== ///
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{tiny}" << endl
           << "\\begin{center}" << endl;
  // if QCD is considered
  //   if(!strcmp(finalstate,"EE")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  //   if(!strcmp(finalstate,"EM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  //   if(!strcmp(finalstate,"MM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  // if QCD is considered
  //   if(!strcmp(finalstate,"EE")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  //   if(!strcmp(finalstate,"EM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  //   if(!strcmp(finalstate,"MM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) \t\\\\" << endl;
  textfile << "selection & data & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ \t\\\\" << endl;
  textfile << "\\hline" << endl; 
    
  textfile << "\\hline"        << endl;
  textfile << "\\hline"        << endl;
    
  for(int icut=0; icut<9; icut++) {

    if(!strcmp(finalstate,"EE")) {
      if(!UsePreSelCuts[ee][icut]) continue;
    }
    if(!strcmp(finalstate,"MM")) {
      if(!UsePreSelCuts[mm][icut]) continue;
    }
    if(!strcmp(finalstate,"EM")) {
      if(!UsePreSelCuts[em][icut]) continue;
    }

    textfile << preSelCuts[icut] << "\t&\t";
      
    textfile << fixed
             << data_preSel[icut]    << "\t&\t"
             << H_preSel[icut]       << " (" << 100. * H_eff_preSel[icut]      << "\\%)" << "\t&\t"
             << Wj_preSel[icut]      << " (" << 100. * Wj_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << ttj_preSel[icut]     << " (" << 100. * ttj_eff_preSel[icut]    << "\\%)" << "\t&\t"
             << SingleTop_preSel[icut]  << " (" << 100. * SingleTop_eff_preSel[icut]    << "\\%)" << "\t&\t"
             << Zj_preSel[icut]      << " (" << 100. * Zj_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << WW_preSel[icut]      << " (" << 100. * WW_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << ZZ_preSel[icut]      << " (" << 100. * ZZ_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << WZ_preSel[icut]      << " (" << 100. * WZ_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << Wgamma_preSel[icut]  << " (" << 100. * Wgamma_eff_preSel[icut] << "\\%)" << "\t";
    //     if(!strcmp(finalstate,"MM")) textfile << QCDmu_preSel[icut]   << " (" << 100. * QCDmu_eff_preSel[icut]  << "\\%)";
    //     if(!strcmp(finalstate,"EM")) textfile << QCDmu_preSel[icut]   << " (" << 100. * QCDmu_eff_preSel[icut]  << "\\%)" << "\t&\t";
    //     if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    //       textfile << QCDem_preSel[icut]   << " (" << 100. * QCDem_eff_preSel[icut]  << "\\%)" << "\t&\t"
    //                << QCDbc_preSel[icut]   << " (" << 100. * QCDbc_eff_preSel[icut]  << "\\%)" << "\t&\t"
    //                << Photj_preSel[icut]   << " (" << 100. * Photj_eff_preSel[icut]  << "\\%)";
    //     }
    textfile << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total preSelection " << "\t&\t"
           << data_preSel[8]   << "\t&\t"
           << H_preSel[8]      << " (" << 100. * H_finaleff_preSel  << "\\%)"  << "\t&\t"
           << Wj_preSel[8]     << " (" << 100. * Wj_finaleff_preSel << "\\%)"  << "\t&\t"
           << ttj_preSel[8]    << " (" << 100. * ttj_finaleff_preSel << "\\%)" << "\t&\t"
           << SingleTop_preSel[8]    << " (" << 100. * SingleTop_finaleff_preSel << "\\%)" << "\t&\t"
           << Zj_preSel[8]     << " (" << 100. * Zj_finaleff_preSel << "\\%)"  << "\t&\t"
           << WW_preSel[8]     << " (" << 100. * WW_finaleff_preSel << "\\%)"  << "\t&\t"
           << ZZ_preSel[8]     << " (" << 100. * ZZ_finaleff_preSel << "\\%)"  << "\t&\t"
           << WZ_preSel[8]     << " (" << 100. * WZ_finaleff_preSel << "\\%)"  << "\t&\t"
           << Wgamma_preSel[8] << " (" << 100. * Wgamma_finaleff_preSel << "\\%)"  << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_preSel[8]   << " (" << 100. * QCDmu_eff_preSel[8]  << "\\%)";
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_preSel[8]   << " (" << 100. * QCDmu_eff_preSel[8]  << "\\%)" << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_preSel[8]  << " (" << 100. * QCDem_finaleff_preSel << "\\%)"  << "\t&\t"
  //              << QCDbc_preSel[8]  << " (" << 100. * QCDbc_finaleff_preSel << "\\%)"  << "\t&\t"
  //              << Photj_preSel[8]  << " (" << 100. * Photj_finaleff_preSel << "\\%)";
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl;
    
  for(int icut=0; icut<22; icut++) {

    if(!strcmp(finalstate,"EE")) {
      if(!UseCuts[ee][icut]) continue;
    }
    if(!strcmp(finalstate,"MM")) {
      if(!UseCuts[mm][icut]) continue;
    }
    if(!strcmp(finalstate,"EM")) {
      if(!UseCuts[em][icut]) continue;
    }

    textfile << fullSelCuts[icut] << "\t&\t";
      
    textfile << fixed
             << data_fullSel[icut]   << "\t&\t"
             << H_fullSel[icut]      << " (" << 100. * H_eff_fullSel[icut]   << "\\%)" << "\t&\t"
             << Wj_fullSel[icut]     << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << ttj_fullSel[icut]    << " (" << 100. * ttj_eff_fullSel[icut] << "\\%)" << "\t&\t"
             << SingleTop_fullSel[icut]    << " (" << 100. * SingleTop_eff_fullSel[icut] << "\\%)" << "\t&\t"
             << Zj_fullSel[icut]     << " (" << 100. * Zj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << WW_fullSel[icut]     << " (" << 100. * WW_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << ZZ_fullSel[icut]     << " (" << 100. * ZZ_eff_fullSel[icut]     << "\\%)" << "\t&\t"
             << WZ_fullSel[icut]     << " (" << 100. * WZ_eff_fullSel[icut]     << "\\%)" << "\t&\t"
             << Wgamma_fullSel[icut] << " (" << 100. * Wgamma_eff_fullSel[icut] << "\\%)" << "\t";
    //     if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[icut]   << " (" << 100. * QCDmu_eff_fullSel[icut]  << "\\%)";
    //     if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[icut]   << " (" << 100. * QCDmu_eff_fullSel[icut]  << "\\%)" << "\t&\t";
    //     if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    //       textfile << QCDem_fullSel[icut]  << " (" << 100. * QCDem_eff_fullSel[icut]  << "\\%)" << "\t&\t"
    //                << QCDbc_fullSel[icut]  << " (" << 100. * QCDbc_eff_fullSel[icut]  << "\\%)" << "\t&\t"
    //                << Photj_fullSel[icut]  << " (" << 100. * Photj_eff_fullSel[icut]  << "\\%)";
    //     }
    textfile << "\t\\\\" << endl;
  }
    
  textfile << "\\hline" << endl;

  textfile << "total fullselection " << "\t&\t"
           << data_fullSel[21]   << "\t&\t"
           << H_fullSel[21]      << " (" << 100. * H_finaleff_fullSel  << "\\%)"     << "\t&\t"
           << Wj_fullSel[21]     << " (" << 100. * Wj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ttj_fullSel[21]    << " (" << 100. * ttj_finaleff_fullSel << "\\%)"    << "\t&\t"
           << SingleTop_fullSel[21]    << " (" << 100. * SingleTop_finaleff_fullSel << "\\%)"    << "\t&\t"
           << Zj_fullSel[21]     << " (" << 100. * Zj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WW_fullSel[21]     << " (" << 100. * WW_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ZZ_fullSel[21]     << " (" << 100. * ZZ_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WZ_fullSel[21]     << " (" << 100. * WZ_finaleff_fullSel << "\\%)"   << "\t&\t"
           << Wgamma_fullSel[21] << " (" << 100. * Wgamma_finaleff_fullSel << "\\%)" << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)";
  //   if(!strcmp(finalstate,"EM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)" << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[21]  << " (" << 100. * QCDem_finaleff_fullSel << "\\%)" << "\t&\t"
  //              << QCDbc_fullSel[21]  << " (" << 100. * QCDbc_finaleff_fullSel << "\\%)" << "\t&\t"
  //              << Photj_fullSel[21]  << " (" << 100. * Photj_finaleff_fullSel << "\\%)";
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl;
    
  textfile << "total " << "\t&\t"
           << data_fullSel[21]   << "\t&\t"
           << H_fullSel[21]      << " (" << 100. * H_finaleff   << "\\%)"    << "\t&\t"
           << Wj_fullSel[21]     << " (" << 100. * Wj_finaleff  << "\\%)"    << "\t&\t"
           << ttj_fullSel[21]    << " (" << 100. * ttj_finaleff << "\\%)"    << "\t&\t"
           << SingleTop_fullSel[21]    << " (" << 100. * SingleTop_finaleff << "\\%)"    << "\t&\t"
           << Zj_fullSel[21]     << " (" << 100. * Zj_finaleff << "\\%)"     << "\t&\t"
           << WW_fullSel[21]     << " (" << 100. * WW_finaleff << "\\%)"     << "\t&\t"
           << ZZ_fullSel[21]     << " (" << 100. * ZZ_finaleff << "\\%)"     << "\t&\t"
           << WZ_fullSel[21]     << " (" << 100. * WZ_finaleff << "\\%)"     << "\t&\t"
           << Wgamma_fullSel[21] << " (" << 100. * Wgamma_finaleff << "\\%)" << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[21] << " (" << 100. * QCDmu_finaleff << "\\%)";
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[21] << " (" << 100. * QCDmu_finaleff << "\\%)" << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[21]  << " (" << 100. * QCDem_finaleff << "\\%)" << "\t&\t"
  //              << QCDbc_fullSel[21]  << " (" << 100. * QCDbc_finaleff << "\\%)" << "\t&\t"
  //              << Photj_fullSel[21]  << " (" << 100. * Photj_finaleff << "\\%)";
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "0 jets bin " << "\t&\t"
           << data_fullSel[21]   << "\t&\t"
           << H_fullSel[21]      << "\t&\t"
           << Wj_fullSel[21]     << "\t&\t"
           << ttj_fullSel[21]    << "\t&\t"
           << SingleTop_fullSel[21]    << "\t&\t"
           << Zj_fullSel[21]     << "\t&\t"
           << WW_fullSel[21]     << "\t&\t"
           << ZZ_fullSel[21]     << "\t&\t"
           << WZ_fullSel[21]     << "\t&\t"
           << Wgamma_fullSel[21] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[21];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[21] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[21]  << "\t&\t"
  //              << QCDbc_fullSel[21]  << "\t&\t"
  //              << Photj_fullSel[21];
  //   }
  textfile << "\t\\\\" << endl;

  textfile << "1 jets bin " << "\t&\t"
           << data_fullSel[22]   << "\t&\t"
           << H_fullSel[22]      << "\t&\t"
           << Wj_fullSel[22]     << "\t&\t"
           << ttj_fullSel[22]    << "\t&\t"
           << SingleTop_fullSel[22]    << "\t&\t"
           << Zj_fullSel[22]     << "\t&\t"
           << WW_fullSel[22]     << "\t&\t"
           << ZZ_fullSel[22]   << "\t&\t"
           << WZ_fullSel[22]   << "\t&\t"
           << Wgamma_fullSel[22] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[22];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[22] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[22]  << "\t&\t"
  //              << QCDbc_fullSel[22]  << "\t&\t"
  //              << Photj_fullSel[22];
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "$>1$ jets bin " << "\t&\t"
           << data_fullSel[23]   << "\t&\t"
           << H_fullSel[23]      << "\t&\t"
           << Wj_fullSel[23]     << "\t&\t"
           << ttj_fullSel[23]    << "\t&\t"
           << SingleTop_fullSel[23]    << "\t&\t"
           << Zj_fullSel[23]     << "\t&\t"
           << WW_fullSel[23]     << "\t&\t"
           << ZZ_fullSel[23]   << "\t&\t"
           << WZ_fullSel[23]   << "\t&\t"
           << Wgamma_fullSel[23] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[23];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[23] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[23]  << "\t&\t"
  //              << QCDbc_fullSel[23]  << "\t&\t"
  //              << Photj_fullSel[23];
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\end{center}" << endl
           << "\\end{tiny}" << endl
           << "\\caption{Breakdown of signal and backgrounds events in "
           << lumi << " $pb^{-1}$ for " << finalstate << " final state.} " << endl 
           << "\\end{sidewaystable}" << endl;

  // assign the final yields
  int channel=-1;
  if(!strcmp(finalstate,"MM")) channel = mm;
  if(!strcmp(finalstate,"EM")) channel = em;
  if(!strcmp(finalstate,"EE")) channel = ee;

  data_final[channel] = data_fullSel[21];
  H_final[channel] = H_fullSel[21];
  Wj_final[channel] = Wj_fullSel[21];
  ttj_final[channel] = ttj_fullSel[21];
  SingleTop_final[channel] = SingleTop_fullSel[21];
  Zj_final[channel] = Zj_fullSel[21];
  WW_final[channel] = WW_fullSel[21];
  ZZ_final[channel] = ZZ_fullSel[21];
  WZ_final[channel] = WZ_fullSel[21];
  Wgamma_final[channel] = Wgamma_fullSel[21];

}

void printShortBkgSummary(float lumi) {

  /// ==============  print short summary  ================== ///

  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  std::string channelName[3];
  channelName[mm] = "$\\mu\\mu$";
  channelName[ee] = "$ee$";
  channelName[em] = "$e\\mu$";

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;

  textfile << "final state & H & WW & WZ & ZZ & Z+jets \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<3; ++ichan) {
    textfile << channelName[ichan] << "\t&\t";
    textfile << fixed 
             << H_final[ichan] << "\t&\t"
             << WW_final[ichan] << "\t&\t"
             << WZ_final[ichan] << "\t&\t"
             << ZZ_final[ichan] << "\t&\t"
             << Zj_final[ichan] << "\t\\\\"
             << endl;
  }
  textfile << "ll" << "\t&\t";
  textfile << fixed
           << H_final[mm] + H_final[ee] + H_final[em] << "\t&\t"
           << WW_final[mm] + WW_final[ee] + WW_final[em] << "\t&\t"
           << WZ_final[mm] + WZ_final[ee] + WZ_final[em] << "\t&\t"
           << ZZ_final[mm] + ZZ_final[ee] + ZZ_final[em] << "\t&\t"
           << Zj_final[mm] + Zj_final[ee] + Zj_final[em] << "\t\\\\"
           << endl;
  textfile << "\\hline" << endl; 

  textfile << "final state & W+jets & W+$\\gamma$ & $t\\bar t$ & single $t$ & data \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<3; ++ichan) {
    textfile << channelName[ichan] << "\t&\t";
    textfile << fixed 
             << Wj_final[ichan] << "\t&\t"
             << Wgamma_final[ichan] << "\t&\t"
             << ttj_final[ichan] << "\t&\t"
             << SingleTop_final[ichan] << "\t&\t"
             << data_final[ichan] << "\t\\\\"
             << endl;
  }
  textfile << "ll" << "\t&\t";
  textfile << fixed
           << Wj_final[mm] + Wj_final[ee] + Wj_final[em] << "\t&\t"
           << Wgamma_final[mm] + Wgamma_final[ee] + Wgamma_final[em] << "\t&\t"
           << ttj_final[mm] + ttj_final[ee] + ttj_final[em] << "\t&\t"
           << SingleTop_final[mm] + SingleTop_final[ee] + SingleTop_final[em] << "\t&\t"
           << data_final[mm] + data_final[ee] + data_final[em] << "\t\\\\"
           << endl;
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\end{center}" << endl
           << "\\caption{Expected backgrounds events in "
           << lumi << " $pb^{-1}$.} " << endl 
           << "\\end{sidewaystable}" << endl;

}

void printLatex(float lumi) {
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::trunc);
  textfile.precision(2);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;

  textfile.close();

  printLatex(lumi, "EE");
  printLatex(lumi, "MM");
  printLatex(lumi, "EM");

  printShortBkgSummary(lumi);

  textfile.open(namefile, ios_base::app);
  textfile << "\\end{document}" << endl;
  textfile.close();

}

void printSuperSummary(float lumi, int massset) {

  int masses[3];
  if(massset==0) {
    masses[0] = 120;
    masses[1] = 130;
    masses[2] = 140;
  } else if(massset==1) {
    masses[0] = 150;
    masses[1] = 160;
    masses[2] = 170;
  } else if(massset==2) {
    masses[0] = 200;
    masses[1] = 300;
    masses[2] = 400;
  }

  char namefile[200];
  sprintf(namefile,"yieldsSummary_byCut.tex");
  ofstream textfile;
  if(massset==0) textfile.open(namefile, ios_base::trunc);
  else textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  if(massset==0) {
    textfile << "\\documentclass{article}" << endl;
    textfile << "\\setlength\\textheight{9.8in}" << endl;
    textfile << "\\usepackage{rotating}" << endl;
    textfile << "\\begin{document}" << endl;
  }

  std::string channelName[3];
  channelName[mm] = "$\\mu\\mu$";
  channelName[ee] = "$ee$";
  channelName[em] = "$e\\mu$";
  
  textfile << "\\begin{table}[p]" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline\\hline" << endl;

  char title[1000];
  sprintf(title,"final state & \\multicolumn{3}{|c|} {$m_H=%d$ GeV} & \\multicolumn{3}{|c|} {$m_H=%d$ GeV} & \\multicolumn{3}{|c|} {$m_H=%d$ GeV} \t\\\\",masses[0],masses[1],masses[2]);
  textfile << title << endl;
  textfile << "\\hline" << endl; 
  textfile << " & signal & background & data & signal & background & data & signal & background & data \\\\" << endl;
  textfile << "\\hline" << endl;

  for(int ichan=0; ichan<3; ++ichan) {

    textfile << channelName[ichan] << "\t&\t";

    for(int i=0; i<3; i++) {
      int mass = masses[i];
      cout << "======> Now analyzing mass = " << mass << "<======" << endl;
      if(ichan==ee) computeYields(lumi,"EE",mass);
      if(ichan==mm) computeYields(lumi,"MM",mass);
      if(ichan==em) computeYields(lumi,"EM",mass);

      H_final[ichan] = H_fullSel[21];
      Wj_final[ichan] = Wj_fullSel[21];
      ttj_final[ichan] = ttj_fullSel[21];
      SingleTop_final[ichan] = SingleTop_fullSel[21];
      Zj_final[ichan] = Zj_fullSel[21];
      WW_final[ichan] = WW_fullSel[21];
      ZZ_final[ichan] = ZZ_fullSel[21];
      WZ_final[ichan] = WZ_fullSel[21];
      Wgamma_final[ichan] = Wgamma_fullSel[21];
      data_final[ichan] = data_fullSel[21];

      textfile << fixed 
               << H_final[ichan] << "\t&\t"
               << WW_final[ichan] + 
        WZ_final[ichan] +
        ZZ_final[ichan] +
        Zj_final[ichan] +
        Wj_final[ichan] +
        ttj_final[ichan] +
        SingleTop_final[ichan] << "\t&\t"
               << data_final[ichan];
      if(i<2) textfile << "\t&\t";
      else textfile << "\t\\\\" << endl;
      cout << "$$$$$$$> Done with mass = " << mass << "<$$$$$$$" << endl;
    }
    cout << "done with ichannel = " << ichan << endl;
    textfile << "\\hline" << endl;
  }
    
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\end{center}" << endl
           << "\\caption{Breakdown of signal and backgrounds events for an integrated luminosity of " << lumi << "pb$^-1$.} "
           << "\\end{table}" << endl;

  if(massset==2) {
    textfile << "\\end{document}" << endl;
    textfile.close();
  }

}

