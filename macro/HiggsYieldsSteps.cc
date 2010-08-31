#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

enum { ee=0, mm=1, em=2 };

string preSelCuts[10];
string fullSelCuts[22];

int UsePreSelCuts[3][10];
int UseCuts[3][22];

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

float H_fullSel[22];
float Wj_fullSel[22];
float ttj_fullSel[22];
float SingleTop_fullSel[22];
float Zj_fullSel[22];
float WW_fullSel[22];
float ZZ_fullSel[22];
float WZ_fullSel[22];
float Wgamma_fullSel[22];
float QCDem_fullSel[22];
float QCDbc_fullSel[22];
float Photj_fullSel[22];
float QCDmu_fullSel[22];
float ttbar_fullSel[22];

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

float H_eff_fullSel[22];
float Wj_eff_fullSel[22];
float ttj_eff_fullSel[22];
float SingleTop_eff_fullSel[22];
float Zj_eff_fullSel[22];
float WW_eff_fullSel[22];
float ZZ_eff_fullSel[22];
float WZ_eff_fullSel[22];
float Wgamma_eff_fullSel[22];
float QCDem_eff_fullSel[22];
float QCDbc_eff_fullSel[22];
float Photj_eff_fullSel[22];
float QCDmu_eff_fullSel[22];
float ttbar_eff_fullSel[22];

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

// xsections
float H120_xsec = 0.247143;
float H130_xsec = 0.452859;
float H140_xsec = 0.64926;
float H150_xsec = 0.787871;
float H155_xsec = 0.842093;
float H160_xsec = 0.897043;
float H165_xsec = 0.867591;
float H170_xsec = 0.808914;
float H175_xsec = 0.751628;
float H180_xsec = 0.685617;
float H190_xsec = 0.503611;
float H200_xsec = 0; // samples not yet produced
float H210_xsec = 0;
float H220_xsec = 0;
float H230_xsec = 0;
float H240_xsec = 0;
float H250_xsec = 0;
float H275_xsec = 0;
float H300_xsec = 0;
float H350_xsec = 0;
float H400_xsec = 0;
float H450_xsec = 0;
float H500_xsec = 0;
float H550_xsec = 0;
float H600_xsec = 0;
float Wgamma_xsec = 41.76;
float Wjets_xsec = 31314;
float Zjets_xsec = 3048;
float TTjets_xsec = 157.5;
float WW_xsec = 4.50347; // WW_2l2nu
float WZ_xsec = 0.599442; // WZ_3l
float ZZ_xsec = 0.25252; // ZZ_2l2nu
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
float SingleTopTW_xsec = 10.6;

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

float Higgs_xsec = H160_xsec;

string sampleNames[29];

void computeYields(float lumi, const char* finalstate) {

  TChain *chains_preSel[29];
  TChain *chains_fullSel[29];
  for(int isample=0; isample<29; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // signal
  sampleNames[0] = "Higgs";
  // backgrounds
  sampleNames[1] = "W+jets";
  sampleNames[2] = "TTbar+jets";
  sampleNames[3] = "Z+jets";
  sampleNames[4] = "WW";
  sampleNames[5] = "ZZ";
  sampleNames[6] = "WZ";
  sampleNames[7] = "Wgamma";
  sampleNames[8] = "SingleTop_sChannel";
  sampleNames[9] = "SingleTop_tChannel";
  sampleNames[10] = "SingleTop_tWChannel";
  sampleNames[11] = "QCD_EMEnriched_Pt20to30";
  sampleNames[12] = "QCD_EMEnriched_Pt30to80";
  sampleNames[13] = "QCD_EMEnriched_Pt80to170";
  sampleNames[14] = "QCD_BCtoE_Pt20to30";
  sampleNames[15] = "QCD_BCtoE_Pt30to80";
  sampleNames[16] = "QCD_BCtoE_Pt80to170";
  sampleNames[17] = "PhotonJet_Pt0to15";
  sampleNames[18] = "PhotonJet_Pt15to20";
  sampleNames[19] = "PhotonJet_Pt20to30";
  sampleNames[20] = "PhotonJet_Pt30to50";
  sampleNames[21] = "PhotonJet_Pt50to80";
  sampleNames[22] = "PhotonJet_Pt80to120";
  sampleNames[23] = "PhotonJet_Pt120to170";
  sampleNames[24] = "PhotonJet_Pt170to300";
  sampleNames[25] = "PhotonJet_Pt300to500";
  sampleNames[26] = "PhotonJet_Pt500toInf";
  sampleNames[27] = "InclusiveMu15";
  // backgrounds (alternative to some of the above ones)
  sampleNames[28] = "TTbar";

  // signal
  chains_preSel[0]->Add("results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*Counters.root");       
  // backgrounds
  chains_preSel[1]->Add("results/WJetsMADGRAPH/WJets-madgraph/*Counters.root");       
  chains_preSel[2]->Add("results/TTbar/TTbarJets-madgraph/*Counters.root");       
  chains_preSel[3]->Add("results/ZJetsMADGRAPH/ZJets-madgraph/*Counters.root");       
  chains_preSel[4]->Add("results/DiBosons/WW_2l_7TeV/*Counters.root");   
  chains_preSel[5]->Add("results/DiBosons/ZZ_2l2nu/*Counters.root");   
  chains_preSel[6]->Add("results/DiBosons/WZ_3l_7TeV/*Counters.root");
  chains_preSel[7]->Add("results/DiBosons/Wgamma/*Counters.root");
  chains_preSel[8]->Add("results/SingleTop/SingleTop_sChannel-madgraph/*Counters.root");
  chains_preSel[9]->Add("results/SingleTop/SingleTop_tChannel-madgraph/*Counters.root");
  chains_preSel[10]->Add("results/SingleTop/SingleTop_tWChannel-madgraph/*Counters.root");
  chains_preSel[11]->Add("results/QCD_EMEnriched_Pt20to30/*Counters.root");
  chains_preSel[12]->Add("results/QCD_EMEnriched_Pt30to80/*Counters.root");
  chains_preSel[13]->Add("results/QCD_EMEnriched_Pt80to170/*Counters.root");
  chains_preSel[14]->Add("results/QCD_BCtoE_Pt20to30/*Counters.root");
  chains_preSel[15]->Add("results/QCD_BCtoE_Pt30to80/*Counters.root");
  chains_preSel[16]->Add("results/QCD_BCtoE_Pt80to170/*Counters.root");
  chains_preSel[17]->Add("results/PhotonJet_Pt0to15/*Counters.root");
  chains_preSel[18]->Add("results/PhotonJet_Pt15to20/*Counters.root");
  chains_preSel[19]->Add("results/PhotonJet_Pt20to30/*Counters.root");
  chains_preSel[20]->Add("results/PhotonJet_Pt30to50/*Counters.root");
  chains_preSel[21]->Add("results/PhotonJet_Pt50to80/*Counters.root");
  chains_preSel[22]->Add("results/PhotonJet_Pt80to120/*Counters.root");
  chains_preSel[23]->Add("results/PhotonJet_Pt120to170/*Counters.root");
  chains_preSel[24]->Add("results/PhotonJet_Pt170to300/*Counters.root");
  chains_preSel[25]->Add("results/PhotonJet_Pt300to500/*Counters.root");
  chains_preSel[26]->Add("results/PhotonJet_Pt500toInf/*Counters.root");
  chains_preSel[27]->Add("results/InclusiveMu15/*Counters.root");
  // backgrounds (alternative to some of the above ones)
  chains_preSel[28]->Add("results/TTbar/*Counters.root");       

  chains_fullSel[0]->Add("results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*Counters.root");       
  // backgrounds
  chains_fullSel[1]->Add("results/WJetsMADGRAPH/WJets-madgraph/*Counters.root");       
  chains_fullSel[2]->Add("results/TTbar/TTbarJets-madgraph/*Counters.root");       
  chains_fullSel[3]->Add("results/ZJetsMADGRAPH/ZJets-madgraph/*Counters.root");       
  chains_fullSel[4]->Add("results/DiBosons/WW_2l_7TeV/*Counters.root");   
  chains_fullSel[5]->Add("results/DiBosons/ZZ_2l2nu/*Counters.root");   
  chains_fullSel[6]->Add("results/DiBosons/WZ_3l_7TeV/*Counters.root");
  chains_fullSel[7]->Add("results/DiBosons/Wgamma/*Counters.root");
  chains_fullSel[8]->Add("results/SingleTop/SingleTop_sChannel-madgraph/*Counters.root");
  chains_fullSel[9]->Add("results/SingleTop/SingleTop_tChannel-madgraph/*Counters.root");
  chains_fullSel[10]->Add("results/SingleTop/SingleTop_tWChannel-madgraph/*Counters.root");
  chains_fullSel[11]->Add("results/QCD_EMEnriched_Pt20to30/*Counters.root");
  chains_fullSel[12]->Add("results/QCD_EMEnriched_Pt30to80/*Counters.root");
  chains_fullSel[13]->Add("results/QCD_EMEnriched_Pt80to170/*Counters.root");
  chains_fullSel[14]->Add("results/QCD_BCtoE_Pt20to30/*Counters.root");
  chains_fullSel[15]->Add("results/QCD_BCtoE_Pt30to80/*Counters.root");
  chains_fullSel[16]->Add("results/QCD_BCtoE_Pt80to170/*Counters.root");
  chains_fullSel[17]->Add("results/PhotonJet_Pt0to15/*Counters.root");
  chains_fullSel[18]->Add("results/PhotonJet_Pt15to20/*Counters.root");
  chains_fullSel[19]->Add("results/PhotonJet_Pt20to30/*Counters.root");
  chains_fullSel[20]->Add("results/PhotonJet_Pt30to50/*Counters.root");
  chains_fullSel[21]->Add("results/PhotonJet_Pt50to80/*Counters.root");
  chains_fullSel[22]->Add("results/PhotonJet_Pt80to120/*Counters.root");
  chains_fullSel[23]->Add("results/PhotonJet_Pt120to170/*Counters.root");
  chains_fullSel[24]->Add("results/PhotonJet_Pt170to300/*Counters.root");
  chains_fullSel[25]->Add("results/PhotonJet_Pt300to500/*Counters.root");
  chains_fullSel[26]->Add("results/PhotonJet_Pt500toInf/*Counters.root");
  chains_fullSel[27]->Add("results/InclusiveMu15/*Counters.root");
  // backgrounds (alternative to some of the above ones)
  chains_fullSel[28]->Add("results/TTbar/*Counters.root");       


  float nPreSelTot[10][29];
  float nFullSelTot[22][29];

  for(int isample=0; isample<29; isample++) {
    for(int icut=0; icut<10; icut++) { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<28; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preselections
  int nCutsAnaPre  = 10;
  int nCutsAnaFull = 22;
  for(int isample=0; isample<29; isample++) {

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
      // nCutsAnaPre = nCutsPre;      
      for(int icut=0; icut<nCutsPre; icut++) nPreSelTot[icut][isample] += nSelPre[icut];	
    }
    
    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      // nCutsAnaFull = nCutsFull;
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }

  // eff at preselection level
  for(int icut=0; icut<nCutsAnaPre; icut++) {

    // numbers
    // N = L * x-sec * eff (eff = N_fin / N_ini)
    if(nPreSelTot[0][0]>0) { 
      H_preSel[icut]      = lumi * Higgs_xsec  * nPreSelTot[icut][0]/nPreSelTot[0][0];
      Wj_preSel[icut]     = lumi * Wjets_xsec  * nPreSelTot[icut][1]/nPreSelTot[0][1];
      ttj_preSel[icut]    = lumi * TTjets_xsec * nPreSelTot[icut][2]/nPreSelTot[0][2];
      Zj_preSel[icut]     = lumi * Zjets_xsec  * nPreSelTot[icut][3]/nPreSelTot[0][3];
      WW_preSel[icut]     = lumi * WW_xsec     * nPreSelTot[icut][4]/nPreSelTot[0][4];
      ZZ_preSel[icut]     = lumi * ZZ_xsec     * nPreSelTot[icut][5]/nPreSelTot[0][5];
      WZ_preSel[icut]     = lumi * WZ_xsec     * nPreSelTot[icut][6]/nPreSelTot[0][6];
      Wgamma_preSel[icut] = lumi * Wgamma_xsec * nPreSelTot[icut][7]/nPreSelTot[0][7];
      
      float singletop_tmp=0.;
      singletop_tmp += lumi * SingleTopS_xsec * nPreSelTot[icut][8]/nPreSelTot[0][8];
      singletop_tmp += lumi * SingleTopT_xsec * nPreSelTot[icut][9]/nPreSelTot[0][9];
      singletop_tmp += lumi * SingleTopTW_xsec * nPreSelTot[icut][10]/nPreSelTot[0][10];
      SingleTop_preSel[icut] = singletop_tmp;

      float qcd_em_tmp=0.;
      qcd_em_tmp += lumi * QCD_EMenriched_Pt20to30_xsec * nPreSelTot[icut][11]/nPreSelTot[0][11];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt30to80_xsec * nPreSelTot[icut][12]/nPreSelTot[0][12];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt80to170_xsec * nPreSelTot[icut][13]/nPreSelTot[0][13];
      QCDem_preSel[icut] = qcd_em_tmp;
      
      float qcd_bc_tmp=0.;
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt20to30_xsec * nPreSelTot[icut][14]/nPreSelTot[0][14];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt30to80_xsec * nPreSelTot[icut][15]/nPreSelTot[0][15];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt80to170_xsec * nPreSelTot[icut][16]/nPreSelTot[0][16];
      QCDbc_preSel[icut] = qcd_bc_tmp;

      float photj_tmp=0.;
      photj_tmp += lumi * PhotonJet_Pt0to15_xsec * nPreSelTot[icut][17]/nPreSelTot[0][17];
      photj_tmp += lumi * PhotonJet_Pt15to20_xsec * nPreSelTot[icut][18]/nPreSelTot[0][18];
      photj_tmp += lumi * PhotonJet_Pt20to30_xsec * nPreSelTot[icut][19]/nPreSelTot[0][19];
      photj_tmp += lumi * PhotonJet_Pt30to50_xsec * nPreSelTot[icut][20]/nPreSelTot[0][20];
      photj_tmp += lumi * PhotonJet_Pt50to80_xsec * nPreSelTot[icut][21]/nPreSelTot[0][21];
      photj_tmp += lumi * PhotonJet_Pt80to120_xsec * nPreSelTot[icut][22]/nPreSelTot[0][22];
      photj_tmp += lumi * PhotonJet_Pt120to170_xsec * nPreSelTot[icut][23]/nPreSelTot[0][23];
      photj_tmp += lumi * PhotonJet_Pt170to300_xsec * nPreSelTot[icut][24]/nPreSelTot[0][24];
      photj_tmp += lumi * PhotonJet_Pt300to500_xsec * nPreSelTot[icut][25]/nPreSelTot[0][25];
      photj_tmp += lumi * PhotonJet_Pt500toInf_xsec * nPreSelTot[icut][26]/nPreSelTot[0][26];
      Photj_preSel[icut] = photj_tmp;

      QCDmu_preSel[icut] = lumi * InclusiveMu15_xsec * nPreSelTot[icut][27]/nPreSelTot[0][27];

      ttbar_preSel[icut] = lumi * TTbar_xsec * nPreSelTot[icut][28]/nPreSelTot[0][28];

    }

    
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
    if(icut>0 && nPreSelTot[icut-1][5]>0) ZZ_eff_preSel[icut]     = nPreSelTot[icut][5]/nPreSelTot[icut-1][5];
    else ZZ_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][6]>0) WZ_eff_preSel[icut]     = nPreSelTot[icut][6]/nPreSelTot[icut-1][6];
    else WZ_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][7]>0) Wgamma_eff_preSel[icut] = nPreSelTot[icut][7]/nPreSelTot[icut-1][7];
    else Wgamma_eff_preSel[icut] = 0.0;
    if(icut>0 && SingleTop_preSel[icut-1]>0) SingleTop_eff_preSel[icut] = SingleTop_preSel[icut]/SingleTop_preSel[icut-1];
    else SingleTop_eff_preSel[icut] = 0.0;
    if(icut>0 && QCDem_preSel[icut-1]>0) QCDem_eff_preSel[icut] = QCDem_preSel[icut]/QCDem_preSel[icut-1];
    else QCDem_eff_preSel[icut] = 0.0;
    if(icut>0 && QCDbc_preSel[icut-1]>0) QCDbc_eff_preSel[icut] = QCDbc_preSel[icut]/QCDbc_preSel[icut-1];
    else QCDbc_eff_preSel[icut] = 0.0;
    if(icut>0 && Photj_preSel[icut-1]>0) Photj_eff_preSel[icut] = Photj_preSel[icut]/Photj_preSel[icut-1];
    else Photj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][24]>0) QCDmu_eff_preSel[icut] = nPreSelTot[icut][27]/nPreSelTot[icut-1][27];
    else QCDmu_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][25]>0) ttbar_eff_preSel[icut] = nPreSelTot[icut][28]/nPreSelTot[icut-1][28];
    else ttbar_eff_preSel[icut] = 0.0;

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
  if(nPreSelTot[0][5]>0) ZZ_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][5]/nPreSelTot[0][5];
  else ZZ_finaleff_preSel = 0.0;
  if(nPreSelTot[0][6]>0) WZ_finaleff_preSel     = nPreSelTot[nCutsAnaPre-2][6]/nPreSelTot[0][6];
  else WZ_finaleff_preSel = 0.0;
  if(nPreSelTot[0][7]>0) Wgamma_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][7]/nPreSelTot[0][7];
  else Wgamma_finaleff_preSel = 0.0;
  if(SingleTop_preSel[0]>0) SingleTop_finaleff_preSel = SingleTop_preSel[nCutsAnaPre-2]/SingleTop_preSel[0];
  else SingleTop_finaleff_preSel = 0.0;
  if(QCDem_preSel[0]>0) QCDem_finaleff_preSel = QCDem_preSel[nCutsAnaPre-2]/QCDem_preSel[0];
  else QCDem_finaleff_preSel = 0.0;
  if(QCDbc_preSel[0]>0) QCDbc_finaleff_preSel = QCDbc_preSel[nCutsAnaPre-2]/QCDbc_preSel[0];
  else QCDbc_finaleff_preSel = 0.0;
  if(Photj_preSel[0]>0) Photj_finaleff_preSel = Photj_preSel[nCutsAnaPre-2]/Photj_preSel[0];
  else Photj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][24]>0) QCDmu_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][27]/nPreSelTot[0][27];
  else QCDmu_finaleff_preSel = 0.0;
  if(nPreSelTot[0][25]>0) ttbar_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][28]/nPreSelTot[0][28];
  else ttbar_finaleff_preSel = 0.0;

  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    // numbers
    if(nFullSelTot[0][0]>0) { 
      H_fullSel[icut]      = lumi * Higgs_xsec  * nFullSelTot[icut][0]/nPreSelTot[0][0];
      Wj_fullSel[icut]     = lumi * Wjets_xsec  * nFullSelTot[icut][1]/nPreSelTot[0][1];
      ttj_fullSel[icut]    = lumi * TTjets_xsec * nFullSelTot[icut][2]/nPreSelTot[0][2];
      Zj_fullSel[icut]     = lumi * Zjets_xsec  * nFullSelTot[icut][3]/nPreSelTot[0][3];
      WW_fullSel[icut]     = lumi * WW_xsec     * nFullSelTot[icut][4]/nPreSelTot[0][4];
      ZZ_fullSel[icut]     = lumi * ZZ_xsec     * nFullSelTot[icut][5]/nPreSelTot[0][5];
      WZ_fullSel[icut]     = lumi * WZ_xsec     * nFullSelTot[icut][6]/nPreSelTot[0][6];
      Wgamma_fullSel[icut] = lumi * Wgamma_xsec * nFullSelTot[icut][7]/nPreSelTot[0][7];

      float singletop_tmp=0.;
      singletop_tmp += lumi * SingleTopS_xsec * nFullSelTot[icut][8]/nPreSelTot[0][8];
      singletop_tmp += lumi * SingleTopT_xsec * nFullSelTot[icut][9]/nPreSelTot[0][9];
      singletop_tmp += lumi * SingleTopTW_xsec * nFullSelTot[icut][10]/nPreSelTot[0][10];
      SingleTop_fullSel[icut] = singletop_tmp;

      float qcd_em_tmp=0.;
      qcd_em_tmp += lumi * QCD_EMenriched_Pt20to30_xsec * nFullSelTot[icut][11]/nPreSelTot[0][11];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt30to80_xsec * nFullSelTot[icut][12]/nPreSelTot[0][12];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt80to170_xsec * nFullSelTot[icut][13]/nPreSelTot[0][13];
      QCDem_fullSel[icut] = qcd_em_tmp;
      
      float qcd_bc_tmp=0.;
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt20to30_xsec * nFullSelTot[icut][14]/nPreSelTot[0][14];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt30to80_xsec * nFullSelTot[icut][15]/nPreSelTot[0][15];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt80to170_xsec * nFullSelTot[icut][16]/nPreSelTot[0][16];
      QCDbc_fullSel[icut] = qcd_bc_tmp;

      float photj_tmp=0.;
      photj_tmp += lumi * PhotonJet_Pt0to15_xsec * nFullSelTot[icut][17]/nPreSelTot[0][17];
      photj_tmp += lumi * PhotonJet_Pt15to20_xsec * nFullSelTot[icut][18]/nPreSelTot[0][18];
      photj_tmp += lumi * PhotonJet_Pt20to30_xsec * nFullSelTot[icut][19]/nPreSelTot[0][19];
      photj_tmp += lumi * PhotonJet_Pt30to50_xsec * nFullSelTot[icut][20]/nPreSelTot[0][20];
      photj_tmp += lumi * PhotonJet_Pt50to80_xsec * nFullSelTot[icut][21]/nPreSelTot[0][21];
      photj_tmp += lumi * PhotonJet_Pt80to120_xsec * nFullSelTot[icut][22]/nPreSelTot[0][22];
      photj_tmp += lumi * PhotonJet_Pt120to170_xsec * nFullSelTot[icut][23]/nPreSelTot[0][23];
      photj_tmp += lumi * PhotonJet_Pt170to300_xsec * nFullSelTot[icut][24]/nPreSelTot[0][24];
      photj_tmp += lumi * PhotonJet_Pt300to500_xsec * nFullSelTot[icut][25]/nPreSelTot[0][25];
      photj_tmp += lumi * PhotonJet_Pt500toInf_xsec * nFullSelTot[icut][26]/nPreSelTot[0][26];
      Photj_fullSel[icut] = photj_tmp;

      QCDmu_fullSel[icut] = lumi * InclusiveMu15_xsec * nFullSelTot[icut][27]/nPreSelTot[0][27];

      ttbar_fullSel[icut] = lumi * TTbar_xsec * nFullSelTot[icut][28]/nPreSelTot[0][28];
    }

    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][0]>0) H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nFullSelTot[icut-1][0];
    else H_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][1]>0) Wj_eff_fullSel[icut]     = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else Wj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][2]>0) ttj_eff_fullSel[icut]    = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else ttj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][3]>0) Zj_eff_fullSel[icut]     = nFullSelTot[icut][3]/nFullSelTot[icut-1][3];
    else Zj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][4]>0) WW_eff_fullSel[icut]     = nFullSelTot[icut][4]/nFullSelTot[icut-1][4];
    else WW_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][5]>0) ZZ_eff_fullSel[icut]     = nFullSelTot[icut][5]/nFullSelTot[icut-1][5];
    else ZZ_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][6]>0) WZ_eff_fullSel[icut]     = nFullSelTot[icut][6]/nFullSelTot[icut-1][6];
    else WZ_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][7]>0) Wgamma_eff_fullSel[icut] = nFullSelTot[icut][7]/nFullSelTot[icut-1][7];
    else Wgamma_eff_fullSel[icut] = 0.0;
    if(icut>0 && SingleTop_fullSel[icut-1]>0) SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[icut-1];
    else SingleTop_eff_fullSel[icut] = 0.0;
    if(icut>0 && QCDem_fullSel[icut-1]>0) QCDem_eff_fullSel[icut] = QCDem_fullSel[icut]/QCDem_fullSel[icut-1];
    else QCDem_eff_fullSel[icut] = 0.0;
    if(icut>0 && QCDbc_fullSel[icut-1]>0) QCDbc_eff_fullSel[icut] = QCDbc_fullSel[icut]/QCDbc_fullSel[icut-1];
    else QCDbc_eff_fullSel[icut] = 0.0;
    if(icut>0 && Photj_fullSel[icut-1]>0) Photj_eff_fullSel[icut] = Photj_fullSel[icut]/Photj_fullSel[icut-1];
    else Photj_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][24]>0) QCDmu_eff_fullSel[icut] = nFullSelTot[icut][24]/nFullSelTot[icut-1][24];
    else QCDmu_eff_fullSel[icut] = 0.0;
    if(icut>0 && nFullSelTot[icut-1][25]>0) ttbar_eff_fullSel[icut] = nFullSelTot[icut][25]/nFullSelTot[icut-1][25];
    else ttbar_eff_fullSel[icut] = 0.0;
    
    if(icut==0) { 
      H_eff_fullSel[icut]      = nFullSelTot[icut][0]/nPreSelTot[nCutsAnaPre-2][0];
      Wj_eff_fullSel[icut]     = nFullSelTot[icut][1]/nPreSelTot[nCutsAnaPre-2][1];
      ttj_eff_fullSel[icut]    = nFullSelTot[icut][2]/nPreSelTot[nCutsAnaPre-2][2];
      Zj_eff_fullSel[icut]     = nFullSelTot[icut][3]/nPreSelTot[nCutsAnaPre-2][3];
      WW_eff_fullSel[icut]     = nFullSelTot[icut][4]/nPreSelTot[nCutsAnaPre-2][4];
      ZZ_eff_fullSel[icut]     = nFullSelTot[icut][5]/nPreSelTot[nCutsAnaPre-2][5];
      WZ_eff_fullSel[icut]     = nFullSelTot[icut][6]/nPreSelTot[nCutsAnaPre-2][6];
      Wgamma_eff_fullSel[icut] = nFullSelTot[icut][7]/nPreSelTot[nCutsAnaPre-2][7];
      SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_preSel[nCutsAnaPre-2];
      QCDem_eff_fullSel[icut] = QCDem_fullSel[icut]/QCDem_preSel[nCutsAnaPre-2];
      QCDbc_eff_fullSel[icut] = QCDbc_fullSel[icut]/QCDbc_preSel[nCutsAnaPre-2];
      Photj_eff_fullSel[icut] = Photj_fullSel[icut]/Photj_preSel[nCutsAnaPre-2];
      QCDmu_eff_fullSel[icut] = nFullSelTot[icut][24]/nPreSelTot[nCutsAnaPre-2][24];
      ttbar_eff_fullSel[icut] = nFullSelTot[icut][25]/nPreSelTot[nCutsAnaPre-2][25];
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][0]>0) H_finaleff_fullSel      = nFullSelTot[nCutsAnaFull-4][0]/nFullSelTot[0][0];
  else H_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][1]>0) Wj_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][1]/nFullSelTot[0][1];
  else Wj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][2]>0) ttj_finaleff_fullSel    = nFullSelTot[nCutsAnaFull-4][2]/nFullSelTot[0][2];
  else ttj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][3]>0) Zj_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][3]/nFullSelTot[0][3];
  else Zj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][4]>0) WW_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][4]/nFullSelTot[0][4];
  else WW_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][5]>0) ZZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][5]/nFullSelTot[0][5];
  else ZZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][6]>0) WZ_finaleff_fullSel     = nFullSelTot[nCutsAnaFull-4][6]/nFullSelTot[0][6];
  else WZ_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][7]>0) Wgamma_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][7]/nFullSelTot[0][7];
  else Wgamma_finaleff_fullSel = 0.0;
  if(SingleTop_fullSel[0]>0) SingleTop_finaleff_fullSel = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_fullSel[0];
  else SingleTop_finaleff_fullSel = 0.0;
  if(QCDem_fullSel[0]>0) QCDem_finaleff_fullSel = QCDem_fullSel[nCutsAnaFull-4]/QCDem_fullSel[0];
  else QCDem_finaleff_fullSel = 0.0;
  if(QCDbc_fullSel[0]>0) QCDbc_finaleff_fullSel = QCDbc_fullSel[nCutsAnaFull-4]/QCDbc_fullSel[0];
  else QCDbc_finaleff_fullSel = 0.0;
  if(Photj_fullSel[0]>0) Photj_finaleff_fullSel = Photj_fullSel[nCutsAnaFull-4]/Photj_fullSel[0];
  else Photj_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][24]>0) QCDmu_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][24]/nFullSelTot[0][24];
  else QCDmu_finaleff_fullSel = 0.0;
  if(nFullSelTot[0][25]>0) ttbar_finaleff_fullSel = nFullSelTot[nCutsAnaFull-4][25]/nFullSelTot[0][25];
  else ttbar_finaleff_fullSel = 0.0;

  // final efficiency combining pre and full selections
  if(nPreSelTot[0][0]>0) H_finaleff      = nFullSelTot[nCutsAnaFull-4][0]/nPreSelTot[0][0];
  else H_finaleff = 0.0;
  if(nPreSelTot[0][1]>0) Wj_finaleff     = nFullSelTot[nCutsAnaFull-4][1]/nPreSelTot[0][1];
  else Wj_finaleff = 0.0;
  if(nPreSelTot[0][2]>0) ttj_finaleff    = nFullSelTot[nCutsAnaFull-4][2]/nPreSelTot[0][2];
  else ttj_finaleff = 0.0;
  if(nPreSelTot[0][3]>0) Zj_finaleff     = nFullSelTot[nCutsAnaFull-4][3]/nPreSelTot[0][3];
  else Zj_finaleff = 0.0;
  if(nPreSelTot[0][4]>0) WW_finaleff     = nFullSelTot[nCutsAnaFull-4][4]/nPreSelTot[0][4];
  else WW_finaleff = 0.0;
  if(nPreSelTot[0][5]>0) ZZ_finaleff     = nFullSelTot[nCutsAnaFull-4][5]/nPreSelTot[0][5];
  else ZZ_finaleff = 0.0;
  if(nPreSelTot[0][6]>0) WZ_finaleff     = nFullSelTot[nCutsAnaFull-4][6]/nPreSelTot[0][6];
  else WZ_finaleff = 0.0;
  if(nPreSelTot[0][7]>0) Wgamma_finaleff = nFullSelTot[nCutsAnaFull-4][7]/nPreSelTot[0][7];
  else Wgamma_finaleff = 0.0;
  if(SingleTop_preSel[0]>0) SingleTop_finaleff = SingleTop_fullSel[nCutsAnaFull-4]/SingleTop_preSel[0];
  else SingleTop_finaleff = 0.0;
  if(QCDem_preSel[0]>0) QCDem_finaleff = QCDem_fullSel[nCutsAnaFull-4]/QCDem_preSel[0];
  else QCDem_finaleff = 0.0;
  if(QCDbc_preSel[0]>0) QCDbc_finaleff = QCDbc_fullSel[nCutsAnaFull-4]/QCDbc_preSel[0];
  else QCDbc_finaleff = 0.0;
  if(Photj_preSel[0]>0) Photj_finaleff = Photj_fullSel[nCutsAnaFull-4]/Photj_preSel[0];
  else Photj_finaleff = 0.0;
  if(nPreSelTot[0][24]>0) QCDmu_finaleff = nFullSelTot[nCutsAnaFull-4][24]/nPreSelTot[0][24];
  else QCDmu_finaleff = 0.0;
  if(nPreSelTot[0][25]>0) ttbar_finaleff = nFullSelTot[nCutsAnaFull-4][25]/nPreSelTot[0][25];
  else ttbar_finaleff = 0.0;

  cout << "\n\nPROCESSED EVENTS:" << endl;
  for(int i=0; i<29; i++) {
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
    for(int i=0; i<22; i++) {
      UseCuts[ichan][i] = 1;
    }
  }

  // unusued ee cuts
  UsePreSelCuts[ee][1] = UsePreSelCuts[ee][7] = UsePreSelCuts[ee][8] = 0; 
  UseCuts[ee][5] = UseCuts[ee][6] = UseCuts[ee][7] = UseCuts[ee][8] = UseCuts[ee][19] = 0; // muon iso

  // unusued mm cuts
  UsePreSelCuts[mm][1] = UsePreSelCuts[mm][7] = UsePreSelCuts[mm][8] = 0; 
  UseCuts[mm][4] = UseCuts[mm][10] = UseCuts[mm][5] = UseCuts[mm][6] = UseCuts[mm][7] = UseCuts[mm][19] = 0; // ele ID and conv. rej. and separate mu iso

  // unusued em cuts
  UsePreSelCuts[em][1] = UsePreSelCuts[em][7] = UsePreSelCuts[em][8] = UseCuts[em][13] = UseCuts[em][19] = 0;
  UseCuts[em][5] = UseCuts[em][6] = UseCuts[em][7] = 0; // separate mu iso 
  
  preSelCuts[0]="event";
  preSelCuts[1]="MCtruth";
  preSelCuts[2]="trigger";
  preSelCuts[3]="$\\geq 2$ leptons";
  preSelCuts[4]="acceptance";
  preSelCuts[5]="$p_T^{max}>20$ GeV";
  preSelCuts[6]="$p_T^{min}>20$ GeV";
  preSelCuts[7]="MET preselection";
  preSelCuts[8]="$m_{ll}>12$ GeV";
  preSelCuts[9]="final presel.";

  fullSelCuts[0]="channel presel.";
  fullSelCuts[1]="$p_T^{max}>20$ GeV";
  fullSelCuts[2]="$p_T^{min}>20$ GeV";
  fullSelCuts[3]="e/$\\mu$ d0";
  fullSelCuts[4]="e isolation";
  fullSelCuts[5]="$\\mu$ tracker Iso";
  fullSelCuts[6]="$\\mu$ HCAL Iso";
  fullSelCuts[7]="$\\mu$ ECAL Iso";
  fullSelCuts[8]="$\\mu$ isolation";
  fullSelCuts[9]="e/$\\mu$ ID";
  fullSelCuts[10]="conv. rej.";
  fullSelCuts[11]="$MET>20$ GeV";
  fullSelCuts[12]="$m_{ll}>12$ GeV";
  fullSelCuts[13]="$|m_{ll}-m_Z|>15$ GeV";
  fullSelCuts[14]="tight (p)MET";
  fullSelCuts[15]="jet veto";
  fullSelCuts[16]="$\\mu^{soft}$ veto";
  fullSelCuts[17]="extra lepton veto";
  fullSelCuts[18]="$\\Delta \\phi$";
  fullSelCuts[19]="final";
  fullSelCuts[20]="1 jets";
  fullSelCuts[21]="$>1$ jets";

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
  textfile.precision(1);

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{tiny}" << endl
           << "\\begin{center}" << endl;
  // if QCD is considered
  //   if(!strcmp(finalstate,"EE")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  //   if(!strcmp(finalstate,"EM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  //   if(!strcmp(finalstate,"MM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  // if QCD is considered
  //   if(!strcmp(finalstate,"EE")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  //   if(!strcmp(finalstate,"EM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  //   if(!strcmp(finalstate,"MM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) \t\\\\" << endl;
  textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & single t & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ \t\\\\" << endl;
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
  
  textfile << "total preselection " << "\t&\t"
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
    
  for(int icut=0; icut<20; icut++) {

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
           << H_fullSel[19]      << " (" << 100. * H_finaleff_fullSel  << "\\%)"     << "\t&\t"
           << Wj_fullSel[19]     << " (" << 100. * Wj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ttj_fullSel[19]    << " (" << 100. * ttj_finaleff_fullSel << "\\%)"    << "\t&\t"
           << SingleTop_fullSel[19]    << " (" << 100. * SingleTop_finaleff_fullSel << "\\%)"    << "\t&\t"
           << Zj_fullSel[19]     << " (" << 100. * Zj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WW_fullSel[19]     << " (" << 100. * WW_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ZZ_fullSel[19]     << " (" << 100. * ZZ_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WZ_fullSel[19]     << " (" << 100. * WZ_finaleff_fullSel << "\\%)"   << "\t&\t"
           << Wgamma_fullSel[19] << " (" << 100. * Wgamma_finaleff_fullSel << "\\%)" << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)";
  //   if(!strcmp(finalstate,"EM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)" << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[19]  << " (" << 100. * QCDem_finaleff_fullSel << "\\%)" << "\t&\t"
  //              << QCDbc_fullSel[19]  << " (" << 100. * QCDbc_finaleff_fullSel << "\\%)" << "\t&\t"
  //              << Photj_fullSel[19]  << " (" << 100. * Photj_finaleff_fullSel << "\\%)";
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl;
    
  textfile << "total " << "\t&\t"
           << H_fullSel[19]      << " (" << 100. * H_finaleff   << "\\%)"    << "\t&\t"
           << Wj_fullSel[19]     << " (" << 100. * Wj_finaleff  << "\\%)"    << "\t&\t"
           << ttj_fullSel[19]    << " (" << 100. * ttj_finaleff << "\\%)"    << "\t&\t"
           << SingleTop_fullSel[19]    << " (" << 100. * SingleTop_finaleff << "\\%)"    << "\t&\t"
           << Zj_fullSel[19]     << " (" << 100. * Zj_finaleff << "\\%)"     << "\t&\t"
           << WW_fullSel[19]     << " (" << 100. * WW_finaleff << "\\%)"     << "\t&\t"
           << ZZ_fullSel[19]     << " (" << 100. * ZZ_finaleff << "\\%)"     << "\t&\t"
           << WZ_fullSel[19]     << " (" << 100. * WZ_finaleff << "\\%)"     << "\t&\t"
           << Wgamma_fullSel[19] << " (" << 100. * Wgamma_finaleff << "\\%)" << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[19] << " (" << 100. * QCDmu_finaleff << "\\%)";
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[19] << " (" << 100. * QCDmu_finaleff << "\\%)" << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[19]  << " (" << 100. * QCDem_finaleff << "\\%)" << "\t&\t"
  //              << QCDbc_fullSel[19]  << " (" << 100. * QCDbc_finaleff << "\\%)" << "\t&\t"
  //              << Photj_fullSel[19]  << " (" << 100. * Photj_finaleff << "\\%)";
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "0 jets bin " << "\t&\t"
           << H_fullSel[19]      << "\t&\t"
           << Wj_fullSel[19]     << "\t&\t"
           << ttj_fullSel[19]    << "\t&\t"
           << SingleTop_fullSel[19]    << "\t&\t"
           << Zj_fullSel[19]     << "\t&\t"
           << WW_fullSel[19]     << "\t&\t"
           << ZZ_fullSel[19]     << "\t&\t"
           << WZ_fullSel[19]     << "\t&\t"
           << Wgamma_fullSel[19] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[19];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[19] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[19]  << "\t&\t"
  //              << QCDbc_fullSel[19]  << "\t&\t"
  //              << Photj_fullSel[19];
  //   }
  textfile << "\t\\\\" << endl;

  textfile << "1 jets bin " << "\t&\t"
           << H_fullSel[20]      << "\t&\t"
           << Wj_fullSel[20]     << "\t&\t"
           << ttj_fullSel[20]    << "\t&\t"
           << SingleTop_fullSel[20]    << "\t&\t"
           << Zj_fullSel[20]     << "\t&\t"
           << WW_fullSel[20]     << "\t&\t"
           << ZZ_fullSel[20]   << "\t&\t"
           << WZ_fullSel[20]   << "\t&\t"
           << Wgamma_fullSel[20] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[20];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[20] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[20]  << "\t&\t"
  //              << QCDbc_fullSel[20]  << "\t&\t"
  //              << Photj_fullSel[20];
  //   }
  textfile << "\t\\\\" << endl;
    
  textfile << "$>1$ jets bin " << "\t&\t"
           << H_fullSel[21]      << "\t&\t"
           << Wj_fullSel[21]     << "\t&\t"
           << ttj_fullSel[21]    << "\t&\t"
           << SingleTop_fullSel[21]    << "\t&\t"
           << Zj_fullSel[21]     << "\t&\t"
           << WW_fullSel[21]     << "\t&\t"
           << ZZ_fullSel[21]   << "\t&\t"
           << WZ_fullSel[21]   << "\t&\t"
           << Wgamma_fullSel[21] << "\t";
  //   if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[21];
  //   if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[21] << "\t&\t";
  //   if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
  //     textfile << QCDem_fullSel[21]  << "\t&\t"
  //              << QCDbc_fullSel[21]  << "\t&\t"
  //              << Photj_fullSel[21];
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

  H_final[channel] = H_fullSel[19];
  Wj_final[channel] = Wj_fullSel[19];
  ttj_final[channel] = ttj_fullSel[19];
  SingleTop_final[channel] = SingleTop_fullSel[19];
  Zj_final[channel] = Zj_fullSel[19];
  WW_final[channel] = WW_fullSel[19];
  ZZ_final[channel] = ZZ_fullSel[19];
  WZ_final[channel] = WZ_fullSel[19];
  Wgamma_final[channel] = Wgamma_fullSel[19];

}

void printShortBkgSummary(float lumi) {

  /// ==============  print short summary  ================== ///

  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(1);

  std::string channelName[3];
  channelName[mm] = "$\\mu\\mu$";
  channelName[ee] = "$ee$";
  channelName[em] = "$e\\mu$";

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;

  textfile << "final state & WW & WZ & ZZ & Z+jets \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<3; ++ichan) {
    textfile << channelName[ichan] << "\t&\t";
    textfile << fixed 
             << WW_final[ichan] << "\t&\t"
             << WZ_final[ichan] << "\t&\t"
             << ZZ_final[ichan] << "\t&\t"
             << Zj_final[ichan] << "\t\\\\"
             << endl;
  }
  textfile << "ll" << "\t&\t";
  textfile << fixed
           << WW_final[mm] + WW_final[ee] + WW_final[em] << "\t&\t"
           << WZ_final[mm] + WZ_final[ee] + WZ_final[em] << "\t&\t"
           << ZZ_final[mm] + ZZ_final[ee] + ZZ_final[em] << "\t&\t"
           << Zj_final[mm] + Zj_final[ee] + Zj_final[em] << "\t\\\\"
           << endl;
  textfile << "\\hline" << endl; 

  textfile << "final state & W+jets & W+$\\gamma$ & $t\\bar t$ & single $t$ \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  for(int ichan=0; ichan<3; ++ichan) {
    textfile << channelName[ichan] << "\t&\t";
    textfile << fixed 
             << Wj_final[ichan] << "\t&\t"
             << Wgamma_final[ichan] << "\t&\t"
             << ttj_final[ichan] << "\t&\t"
             << SingleTop_final[ichan] << "\t\\\\"
             << endl;
  }
  textfile << "ll" << "\t&\t";
  textfile << fixed
           << Wj_final[mm] + Wj_final[ee] + Wj_final[em] << "\t&\t"
           << Wgamma_final[mm] + Wgamma_final[ee] + Wgamma_final[em] << "\t&\t"
           << ttj_final[mm] + ttj_final[ee] + ttj_final[em] << "\t&\t"
           << SingleTop_final[mm] + SingleTop_final[ee] + SingleTop_final[em] << "\t\\\\"
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
  textfile.precision(1);

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

