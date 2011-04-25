#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

// parameters to configure
int exampleHiggsMass = -1;
bool runStandalone = true;

using namespace std;

enum { ee=0, mm=1, em=2, me=3 };

string fullSelCuts[24];

int UseCuts[4][24];

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

float H_final[4];
float Wj_final[4];
float ttj_final[4];
float SingleTop_final[4];
float Zj_final[4];
float WW_final[4];
float ZZ_final[4];
float WZ_final[4];
float Wgamma_final[4];
float data_final[4];

// xsections
std::map<int,float> Higgs_xsec_masses;

float Wgamma_xsec      = 41.76;      // not used
// float Wlnu_xsec        = 31314.;     // madgraph 
float Wlnu_xsec        = 31314./3. * 0.742;  // pythia
float ZjetsLoMass_xsec = 2659./3.;           // pythia <20
float ZjetsHiMass_xsec = 4998./3.;           // pythia >20
float TTjets_xsec      = 157.5;
float ggWW_xsec        = 0.1538; // gg->WW->4l
float WW_xsec          = 4.50347; // WW_2l2nu
float WZ_xsec          = 0.599442; // WZ_3l
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

string sampleNames[20];

void computeYields(float lumi, const char* finalstate, int mass=0) {

  TChain *chains_fullSel[20];

  for(int isample=0; isample<20; isample++) {
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // signal
  sampleNames[0]  = "Higgs";
  // backgrounds
  sampleNames[1]  = "TTbar+jets";
  sampleNames[2]  = "Z(ee) 10-20 GeV";
  sampleNames[3]  = "Z(mumu) 10-20 GeV";
  sampleNames[4]  = "Z(tautau) 10-20 GeV ";
  sampleNames[5]  = "Z(ee) >20 GeV";
  sampleNames[6]  = "Z(mumu) >20 GeV";
  sampleNames[7]  = "Z(tautau) >20 GeV ";
  sampleNames[8]  = "WW";
  sampleNames[9]  = "gg->WW";
  sampleNames[10] = "ZZ";
  sampleNames[11] = "WZ";
  sampleNames[12] = "Wgamma";
  sampleNames[13] = "SingleTop_sChannel";
  sampleNames[14] = "SingleTop_tChannel";
  sampleNames[15] = "SingleTop_tWChannel";
  sampleNames[16] = "WToENu";
  sampleNames[17] = "WToMuNu";
  sampleNames[18] = "WToTauNu";
  sampleNames[19] = "data";

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

    // mc
    char dir_mc[1000];                  
    sprintf(dir_mc,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/MC2011_WP90_V2/OptimMH%d/Spring11/",exampleHiggsMass);
    char HiggsSample[500];
    sprintf(HiggsSample,"GluGluToHToWWTo2L2Nu_M-%d_7TeV-powheg-pythia6/*Counters.root",exampleHiggsMass);

    // data
    char dir_data[1000]; 
    sprintf(dir_data,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/Data2011_WP90_V2/OptimMH%d/Data7TeV",exampleHiggsMass);

    // signal
    chains_fullSel[0]->Add(TString(dir_mc)+TString(HiggsSample));       
    // backgrounds
    chains_fullSel[1]->Add(TString(dir_mc)+TString("/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir_mc)+TString("/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir_mc)+TString("/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir_mc)+TString("/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir_mc)+TString("/DYToEE_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir_mc)+TString("/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir_mc)+TString("/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[8]->Add(TString(dir_mc)+TString("/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[9]->Add(TString(dir_mc)+TString("/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*Counters.root"));
    chains_fullSel[10]->Add(TString(dir_mc)+TString("/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));
    chains_fullSel[11]->Add(TString(dir_mc)+TString("/WZTo3LNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    // chains_fullSel[12]->Add(TString(dir_mc)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_fullSel[13]->Add(TString(dir_mc)+TString("TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[14]->Add(TString(dir_mc)+TString("TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[15]->Add(TString(dir_mc)+TString("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[16]->Add(TString(dir_mc)+TString("WToENu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[17]->Add(TString(dir_mc)+TString("WToMuNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[18]->Add(TString(dir_mc)+TString("WToTauNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    // data
    if (finalstate=="EE") chains_fullSel[19]->Add(TString(dir_data)+TString("/DoubleElectron/*Counters.root"));
    if (finalstate=="MM") chains_fullSel[19]->Add(TString(dir_data)+TString("/DoubleMu/*Counters.root"));
    if (finalstate=="EM" || finalstate=="ME") chains_fullSel[19]->Add(TString(dir_data)+TString("/MuEG/*Counters.root"));

  } else {

    Higgs_xsec = Higgs_xsec_masses[mass] * 4./9.; // 4/9 because we are considering only the samples containing e-mu combinations.

    // mc
    char dir[1000];                  
    sprintf(dir,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/MC2011_WP90_V2/OptimMH%d/Spring11/",mass);
    char HiggsDir[500];
    sprintf(HiggsDir,"GluGluToHToWWTo2L2Nu_M-%d_7TeV-powheg-pythia6/*Counters.root",mass);

    // data
    char dir_data[1000]; 
    if (finalstate=="EE") sprintf(dir_data,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/Data2011_WP90_V2/OptimMH%d/Data7TeV/DoubleElectron",mass);
    if (finalstate=="MM") sprintf(dir_data,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/Data2011_WP90_V2/OptimMH%d/Data7TeV/DoubleMu",mass);
    if (finalstate=="EM" || finalstate=="ME") sprintf(dir_data,"/cmsrm/pc23_2/emanuele/data/Higgs4.1.X/Data2011_WP90_V2/OptimMH%d/Data7TeV/MuEG",mass);

    // signal
    chains_fullSel[0]->Add(TString(dir)+TString(HiggsDir)+TString("*Counters.root"));
    // backgrounds
    chains_fullSel[1]->Add(TString(dir)+TString("/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir)+TString("/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir)+TString("/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir)+TString("/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir)+TString("/DYToEE_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir)+TString("/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir)+TString("/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*Counters.root"));       
    chains_fullSel[8]->Add(TString(dir)+TString("/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[9]->Add(TString(dir)+TString("/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*Counters.root"));
    chains_fullSel[10]->Add(TString(dir)+TString("/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*Counters.root"));
    chains_fullSel[11]->Add(TString(dir)+TString("/WZTo3LNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    // chains_fullSel[12]->Add(TString(dir)+TString("/DiBosons/WgammaXXXX/*Counters.root"));
    chains_fullSel[13]->Add(TString(dir)+TString("TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[14]->Add(TString(dir)+TString("TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[15]->Add(TString(dir)+TString("TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root"));
    chains_fullSel[16]->Add(TString(dir)+TString("WToENu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[17]->Add(TString(dir)+TString("WToMuNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    chains_fullSel[18]->Add(TString(dir)+TString("WToTauNu_TuneZ2_7TeV-pythia6/*Counters.root"));
    // data
    if (finalstate=="EE") chains_fullSel[19]->Add(TString(dir_data)+TString("/DoubleElectron/*Counters.root"));
    if (finalstate=="MM") chains_fullSel[19]->Add(TString(dir_data)+TString("/DoubleMu/*Counters.root"));
    if (finalstate=="EM" || finalstate=="ME") chains_fullSel[19]->Add(TString(dir_data)+TString("/MuEG/*Counters.root"));
  }

  float nFullSelTot[24][20];

  for(int isample=0; isample<20; isample++) {
    for(int icut=0; icut<24; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }
  
  // full selection
  int nCutsAnaFull = 24;
  for(int isample=0; isample<20; isample++) {

    // List of branches    
    Int_t           nCutsFull;
    Float_t         nSelFull[24];   //[nCuts]
    TBranch        *b_nCutsFull;   //!
    TBranch        *b_nSelFull;    //!
    chains_fullSel[isample]->SetBranchAddress("nCuts", &nCutsFull, &b_nCutsFull);
    chains_fullSel[isample]->SetBranchAddress("nSel",  nSelFull,   &b_nSelFull);
    
    Long64_t nentriesFull = chains_fullSel[isample]->GetEntries();

    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2;
      nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) { 
      
      H_fullSel[icut]   = lumi * Higgs_xsec  * nFullSelTot[icut][0]/nFullSelTot[0][0];

      ttj_fullSel[icut] = lumi * TTjets_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1];

      float Z_tmp=0.;
      Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][2]/nFullSelTot[0][2];
      Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][3]/nFullSelTot[0][3];
      // Z_tmp += lumi * ZjetsLoMass_xsec  * nFullSelTot[icut][4]/nFullSelTot[0][4];
      Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][5]/nFullSelTot[0][5];
      Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][6]/nFullSelTot[0][6];
      // Z_tmp += lumi * ZjetsHiMass_xsec  * nFullSelTot[icut][7]/nFullSelTot[0][7];
      Zj_fullSel[icut]  = Z_tmp;

      float WW_tmp=0.;
      WW_tmp += lumi * WW_xsec   * nFullSelTot[icut][8]/nFullSelTot[0][8];
      WW_tmp += lumi * ggWW_xsec * nFullSelTot[icut][9]/nFullSelTot[0][9];
      WW_fullSel[icut] = WW_tmp;

      ZZ_fullSel[icut]     = lumi * ZZ_xsec     * nFullSelTot[icut][10]/nFullSelTot[0][10];
      WZ_fullSel[icut]     = lumi * WZ_xsec     * nFullSelTot[icut][11]/nFullSelTot[0][11];
      Wgamma_fullSel[icut] = lumi * Wgamma_xsec * nFullSelTot[icut][12]/nFullSelTot[0][12];

      float singletop_tmp=0.;
      singletop_tmp += lumi * SingleTopS_xsec  * nFullSelTot[icut][13]/nFullSelTot[0][13];
      singletop_tmp += lumi * SingleTopT_xsec  * nFullSelTot[icut][14]/nFullSelTot[0][14];
      singletop_tmp += lumi * SingleTopTW_xsec * nFullSelTot[icut][15]/nFullSelTot[0][15];
      SingleTop_fullSel[icut] = singletop_tmp;

      // Wj_fullSel[icut]  = lumi * Wlnu_xsec * nFullSelTot[icut][16]/nFullSelTot[0][16]; 
      float W_tmp=0.;
      W_tmp += lumi * Wlnu_xsec  * nFullSelTot[icut][16]/nFullSelTot[0][16]; // W->enu                                
      W_tmp += lumi * Wlnu_xsec  * nFullSelTot[icut][17]/nFullSelTot[0][17]; // W->munu                               
      // W_tmp += lumi * Wlnu_xsec  * nFullSelTot[icut][18]/nFullSelTot[0][18]; // W->taunu                              
      Wj_fullSel[icut]  = W_tmp;

      // data
      data_fullSel[icut] = nFullSelTot[icut][19];
    }

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0)  H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else H_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][1]>0)  ttj_eff_fullSel[icut]    = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else ttj_eff_fullSel[icut] = 0.0;

    if(icut>0 && Zj_fullSel[icut-1]>0)      Zj_eff_fullSel[icut]     = Zj_fullSel[icut]/Zj_fullSel[icut-1];
    else Zj_eff_fullSel[icut] = 0.0;

    if(icut>0 && WW_fullSel[icut-1]>0)      WW_eff_fullSel[icut]     = WW_fullSel[icut]/WW_fullSel[icut-1];
    else WW_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][10]>0) ZZ_eff_fullSel[icut]     = nFullSelTot[icut][10]/nFullSelTot[icut-1][10];
    else ZZ_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][11]>0) WZ_eff_fullSel[icut]     = nFullSelTot[icut][11]/nFullSelTot[icut-1][11];
    else WZ_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][12]>0) Wgamma_eff_fullSel[icut] = nFullSelTot[icut][12]/nFullSelTot[icut-1][12];
    else Wgamma_eff_fullSel[icut] = 0.0;

    if(icut>0 && SingleTop_fullSel[icut-1]>0) SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[icut-1];
    else SingleTop_eff_fullSel[icut] = 0.0;

    if(icut>0 && Wj_fullSel[icut-1]>0)       Wj_eff_fullSel[icut]     = Wj_fullSel[icut]/Wj_fullSel[icut-1];
    else Wj_eff_fullSel[icut] = 0.0;

    
    if(icut==0) { 
      H_eff_fullSel[icut]         = nFullSelTot[icut][0]/nFullSelTot[0][0];
      ttj_eff_fullSel[icut]       = nFullSelTot[icut][1]/nFullSelTot[0][1];
      Zj_eff_fullSel[icut]        = Zj_fullSel[icut]/Zj_fullSel[0];
      WW_eff_fullSel[icut]        = WW_fullSel[icut]/WW_fullSel[0];
      ZZ_eff_fullSel[icut]        = nFullSelTot[icut][10]/nFullSelTot[0][10];
      WZ_eff_fullSel[icut]        = nFullSelTot[icut][11]/nFullSelTot[0][11];
      Wgamma_eff_fullSel[icut]    = nFullSelTot[icut][12]/nFullSelTot[0][12];
      SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[0];
      Wj_eff_fullSel[icut]        = Wj_fullSel[icut]/Wj_fullSel[0];
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0)  H_finaleff_fullSel  = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0];
  else H_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) ttj_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1];
  else ttj_finaleff_fullSel = 0.0;
  if(Zj_fullSel[0]>0) Zj_finaleff_fullSel      = Zj_fullSel[nCutsAnaFull-4]/Zj_fullSel[0];
  else Zj_finaleff_fullSel = 0.0;
  if(WW_fullSel[0]>0) WW_finaleff_fullSel     = WW_fullSel[nCutsAnaFull-4]/WW_fullSel[0];
  else WW_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][10]>0) ZZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][10]/nFullSelTot[0][10];
  else ZZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][11]>0) WZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][11]/nFullSelTot[0][11];
  else WZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][12]>0) Wgamma_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][12]/nFullSelTot[0][12];
  else Wgamma_finaleff_fullSel = 0.0;
  if(SingleTop_fullSel[0]>0) SingleTop_finaleff_fullSel = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_fullSel[0];
  else SingleTop_finaleff_fullSel = 0.0;
  if(Wj_fullSel[0]>0) Wj_finaleff_fullSel      = Wj_fullSel[nCutsAnaFull-4]/Wj_fullSel[0];
  else Wj_finaleff_fullSel = 0.0;
  
  // final efficiency combining pre and full selections
  if(nFullSelTot[0][0]>0) H_finaleff       = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0];
  else H_finaleff = 0.0;
  if(nFullSelTot[0][1]>0) ttj_finaleff     = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1];
  else ttj_finaleff = 0.0;
  if(Zj_fullSel[0]>0)     Zj_finaleff      = Zj_fullSel[nCutsAnaFull-4]/Zj_fullSel[0];
  else Zj_finaleff = 0.0;
  if(WW_fullSel[0]>0)     WW_finaleff      = WW_fullSel[nCutsAnaFull-4]/WW_fullSel[0];
  else WW_finaleff = 0.0;
  if(nFullSelTot[0][10]>0) ZZ_finaleff     = nFullSelTot[nCutsAnaFull-4][10]/nFullSelTot[0][10];
  else ZZ_finaleff = 0.0;
  if(nFullSelTot[0][11]>0) WZ_finaleff     = nFullSelTot[nCutsAnaFull-4][11]/nFullSelTot[0][11];
  else WZ_finaleff = 0.0;
  if(nFullSelTot[0][12]>0) Wgamma_finaleff = nFullSelTot[nCutsAnaFull-4][12]/nFullSelTot[0][12];
  else Wgamma_finaleff = 0.0;
  if(SingleTop_fullSel[0]>0) SingleTop_finaleff = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_fullSel[0];
  else SingleTop_finaleff = 0.0;
  if(Wj_fullSel[0]>0)     Wj_finaleff      = Wj_fullSel[nCutsAnaFull-4]/Wj_fullSel[0];
  else Wj_finaleff = 0.0;

  cout << "\n\nPROCESSED EVENTS:" << endl;
  for(int i=0; i<20; i++) {
    cout << sampleNames[i] << "\t" << nFullSelTot[0][i] << endl;
  }

}

void setupCuts() {
  
  for(int ichan=0; ichan<4; ichan++) {
    for(int i=0; i<24; i++) {
      UseCuts[ichan][i] = 1;
    }
  }

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
  fullSelCuts[20]="$\\Delta \\phi$";
  fullSelCuts[21]="final";
  fullSelCuts[22]="1 jets";
  fullSelCuts[23]="$>1$ jets";
}


void printLatex(float lumi, const char* finalstate) {

  setupCuts();

  if(strcmp(finalstate,"EE") && strcmp(finalstate,"EM") && strcmp(finalstate,"MM") && strcmp(finalstate,"ME")) {
    cout << "ERROR! finalstate must be one among EE/EM/MM/ME. Exiting..." << endl;
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
  //   if(!strcmp(finalstate,"ME")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
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
    
  for(int icut=0; icut<24; icut++) {

    if(!strcmp(finalstate,"EE")) {
      if(!UseCuts[ee][icut]) continue;
    }
    if(!strcmp(finalstate,"MM")) {
      if(!UseCuts[mm][icut]) continue;
    }
    if(!strcmp(finalstate,"EM")) {
      if(!UseCuts[em][icut]) continue;
    }
    if(!strcmp(finalstate,"ME")) {
      if(!UseCuts[me][icut]) continue;
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
  //     textfile << QCDem_fullSel[23]  << "\t&\t"
  //              << QCDbc_fullSel[23]  << "\t&\t"
  //              << Photj_fullSel[23B];
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
           << "\\caption{Higgs $m_H$ = " << exampleHiggsMass << " GeV/c$^2$. Breakdown of signal and backgrounds events in "
           << lumi << " $pb^{-1}$ for " << finalstate << " final state.} " << endl 
           << "\\end{sidewaystable}" << endl;

  // assign the final yields
  int channel=-1;
  if(!strcmp(finalstate,"MM")) channel = mm;
  if(!strcmp(finalstate,"EM")) channel = em;
  if(!strcmp(finalstate,"EE")) channel = ee;
  if(!strcmp(finalstate,"ME")) channel = me;

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

  std::string channelName[4];
  channelName[mm] = "$\\mu\\mu$";
  channelName[ee] = "$ee$";
  channelName[em] = "$e\\mu$";
  channelName[me] = "$\\mu e$";

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;

  textfile << "final state & H & WW & WZ & ZZ & Z+jets \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<4; ++ichan) {
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
           << H_final[mm]  + H_final[ee]  + H_final[em]  + H_final[me]  << "\t&\t"
           << WW_final[mm] + WW_final[ee] + WW_final[em] + WW_final[me] << "\t&\t"
           << WZ_final[mm] + WZ_final[ee] + WZ_final[em] + WZ_final[me] << "\t&\t"
           << ZZ_final[mm] + ZZ_final[ee] + ZZ_final[em] + ZZ_final[me] << "\t&\t"
           << Zj_final[mm] + Zj_final[ee] + Zj_final[em] + Zj_final[me] << "\t\\\\"
           << endl;
  textfile << "\\hline" << endl; 

  textfile << "final state & W+jets & W+$\\gamma$ & $t\\bar t$ & single $t$ & data \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<4; ++ichan) {
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
           << Wj_final[mm]     + Wj_final[ee]     + Wj_final[em]     + Wj_final[me]      << "\t&\t"
           << Wgamma_final[mm] + Wgamma_final[ee] + Wgamma_final[em] + Wgamma_final[me]  << "\t&\t"
           << ttj_final[mm]    + ttj_final[ee]    + ttj_final[em]    + ttj_final[me]     << "\t&\t"
           << SingleTop_final[mm] + SingleTop_final[ee] + SingleTop_final[em] + SingleTop_final[me] << "\t&\t"
           << data_final[mm] + data_final[ee] + data_final[em]  + data_final[me]  << "\t\\\\"
           << endl;
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\end{center}" << endl
           << "\\caption{Higgs $m_H$ = " << exampleHiggsMass << " GeV/c$^2$. Expected backgrounds events in "
           << lumi << " $pb^{-1}$.} " << endl 
           << "\\end{sidewaystable}" << endl;

}

void printLatex(float lumi, int HiggsMassOptim) {

  exampleHiggsMass = HiggsMassOptim;
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::trunc);
  textfile.precision(2);

  if(runStandalone) {
    textfile << "\\documentclass{article}" << endl;
    textfile << "\\setlength\\textheight{9.8in}" << endl;
    textfile << "\\usepackage{rotating}" << endl;
    textfile << "\\begin{document}" << endl;
  }

  textfile.close();

  printLatex(lumi, "EE");
  printLatex(lumi, "MM");
  printLatex(lumi, "EM");
  printLatex(lumi, "ME");

  printShortBkgSummary(lumi);

  textfile.open(namefile, ios_base::app);
  textfile << "\\clearpage" << endl;
  textfile.close();  

  if(runStandalone) {
    textfile.open(namefile, ios_base::app);
    textfile << "\\end{document}" << endl;
    textfile.close();
  }
}

void printSuperSummary(float lumi, int massset) {

  int masses[3];
  if(massset==0) {
    masses[0] = 120;
    masses[1] = 130;
    masses[2] = 140;
  } else if(massset==1) {
    masses[0] = 160;
    masses[1] = 120;
    masses[2] = 130;
  } 

  char namefile[200];
  sprintf(namefile,"yieldsSummary_byCut.tex");
  ofstream textfile;
  if(massset==0) textfile.open(namefile, ios_base::trunc);
  else textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  if(massset==0) {
    if(runStandalone) {
      textfile << "\\documentclass{article}" << endl;
      textfile << "\\setlength\\textheight{9.8in}" << endl;
      textfile << "\\usepackage{rotating}" << endl;
      textfile << "\\begin{document}" << endl;
    }
  }

  std::string channelName[4];
  channelName[mm] = "$\\mu\\mu$";
  channelName[ee] = "$ee$";
  channelName[em] = "$e\\mu$";
  channelName[me] = "$\\mu e$";
  
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

  for(int ichan=0; ichan<4; ++ichan) {

    textfile << channelName[ichan] << "\t&\t";

    for(int i=0; i<3; i++) {
      int mass = masses[i];
      cout << "======> Now analyzing mass = " << mass << "<======" << endl;
      exampleHiggsMass = mass;
      /*
      if(ichan==ee) computeYields(lumi,"EE",mass);
      if(ichan==mm) computeYields(lumi,"MM",mass);
      if(ichan==em) computeYields(lumi,"EM",mass);
      if(ichan==me) computeYields(lumi,"ME",mass);
      */
      if(ichan==ee) computeYields(lumi,"EE");
      if(ichan==mm) computeYields(lumi,"MM");
      if(ichan==em) computeYields(lumi,"EM");
      if(ichan==me) computeYields(lumi,"ME");

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
        // WZ_final[ichan] +
        // ZZ_final[ichan] +
        Zj_final[ichan] +
        Wj_final[ichan] +
        ttj_final[ichan] +
        SingleTop_final[ichan] 
	       << "\t&\t"
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

  if(massset==1) {
    textfile << "\\end{document}" << endl;
    textfile.close();
  }

}

void setRunStandalone(bool what) { runStandalone = what; }
