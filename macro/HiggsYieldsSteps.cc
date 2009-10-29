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

float H_preSel[11];
float Wj_preSel[11];
float ttj_preSel[11];
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

float H_fullSel[19];
float Wj_fullSel[19];
float ttj_fullSel[19];
float Zj_fullSel[19];
float WW_fullSel[19];
float ZZ_fullSel[19];
float WZ_fullSel[19];
float Wgamma_fullSel[19];
float QCDem_fullSel[19];
float QCDbc_fullSel[19];
float Photj_fullSel[19];
float QCDmu_fullSel[19];
float ttbar_fullSel[19];

float H_eff_preSel[11];
float Wj_eff_preSel[11];
float ttj_eff_preSel[11];
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

float H_eff_fullSel[19];
float Wj_eff_fullSel[19];
float ttj_eff_fullSel[19];
float Zj_eff_fullSel[19];
float WW_eff_fullSel[19];
float ZZ_eff_fullSel[19];
float WZ_eff_fullSel[19];
float Wgamma_eff_fullSel[19];
float QCDem_eff_fullSel[19];
float QCDbc_eff_fullSel[19];
float Photj_eff_fullSel[19];
float QCDmu_eff_fullSel[19];
float ttbar_eff_fullSel[19];

float H_finaleff_preSel;    
float Wj_finaleff_preSel;
float ttj_finaleff_preSel;
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

// xsections
float H120_xsec = 0.429291;
float H130_xsec = 0.798372;
float H140_xsec = 1.16031;
float H150_xsec = 1.42712;
float H155_xsec = 1.53503;
float H160_xsec = 1.64532;
float H165_xsec = 1.60112;
float H170_xsec = 1.50192;
float H175_xsec = 1.40384;
float H180_xsec = 1.28808;
float H190_xsec = 0.956866;
float H200_xsec = 0.811859;
float H210_xsec = 0.720998;
float H220_xsec = 0.652178;
float H230_xsec = 0.596385;
float H240_xsec = 0.549726;
float H250_xsec = 0.510142;
float H275_xsec = 0.434763;
float H300_xsec = 0.385833;
float H350_xsec = 0.382384;
float H400_xsec = 0.294458;
float H450_xsec = 0.194491;
float H500_xsec = 0.126372;
float H550_xsec = 0.082537;
float H600_xsec = 0.054517;
float Wgamma_xsec = 11960;
float Wjets_xsec = 46050; // NLO estimate
float Zjets_xsec = 7164; // NLO estimate
float TTjets_xsec = 415; // NLO estimate
float WW_xsec = 44.8;
float WZ_xsec = 17.4;
float ZZ_xsec = 7.1;
float TTbar_xsec = 375.0;
float QCD_EMenriched_Pt20to30_xsec = 0.4000*0.00800*1.0E+09;
float QCD_EMenriched_Pt30to80_xsec = 0.1000*0.04700*1.0E+09;
float QCD_EMenriched_Pt80to170_xsec = 0.0019*0.15000*1.0E+09;
float QCD_BCtoE_Pt20to30_xsec = 0.4000*0.00048*1.0E+09;
float QCD_BCtoE_Pt30to80_xsec = 0.1000*0.00240*1.0E+09;
float QCD_BCtoE_Pt80to170_xsec = 0.0019*0.01200*1.0E+09;
float PhotonJet_Pt0to15_xsec = 100.5*1.0E+06;
float PhotonJet_Pt15to20_xsec = 168.0*1.0E+03;
float PhotonJet_Pt20to30_xsec = 84.87*1.0E+03;
float PhotonJet_Pt30to50_xsec = 26.32*1.0E+03;
float PhotonJet_Pt50to80_xsec = 4.589*1.0E+03;
float PhotonJet_Pt80to120_xsec = 0.7864*1.0E+03;
float PhotonJet_Pt120to170_xsec = 0.1648*1.0E+03;
float PhotonJet_Pt170to300_xsec = 0.04596*1.0E+03;
float PhotonJet_Pt300to500_xsec = 3.708;
float PhotonJet_Pt500toInf_xsec = 0.3285;
float InclusiveMu15_xsec = 0.5091*0.0002881*1.0E+09;

float Higgs_xsec = H165_xsec;

void computeYields(float lumi, const char* finalstate) {

  TChain *chains_preSel[26];
  TChain *chains_fullSel[26];
  for(int isample=0; isample<26; isample++) {
    chains_preSel[isample]  = new TChain("PRESELECTION_EVENT_COUNTER");
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // signal
  chains_preSel[0]->Add("results/H165/*Counters.root");       
  // backgrounds
  chains_preSel[1]->Add("results/WjetsMadgraph/*Counters.root");       
  chains_preSel[2]->Add("results/TTbarJetsMadgraph/*Counters.root");       
  chains_preSel[3]->Add("results/ZjetsMadgraph/*Counters.root");       
  chains_preSel[4]->Add("results/WW/*Counters.root");   
  chains_preSel[5]->Add("results/ZZ/*Counters.root");   
  chains_preSel[6]->Add("results/WZ/*Counters.root");
  chains_preSel[7]->Add("results/Wgamma/*Counters.root");
  chains_preSel[8]->Add("results/QCD_EMEnriched_Pt20to30Ele10/*Counters.root");
  chains_preSel[9]->Add("results/QCD_EMEnriched_Pt30to80Ele10/*Counters.root");
  chains_preSel[10]->Add("results/QCD_EMEnriched_Pt80to170Ele10/*Counters.root");
  chains_preSel[11]->Add("results/QCD_BCtoE_Pt20to30Ele10/*Counters.root");
  chains_preSel[12]->Add("results/QCD_BCtoE_Pt30to80Ele10/*Counters.root");
  chains_preSel[13]->Add("results/QCD_BCtoE_Pt80to170Ele10/*Counters.root");
  chains_preSel[14]->Add("results/PhotonJet_Pt0to15/*Counters.root");
  chains_preSel[15]->Add("results/PhotonJet_Pt15to20/*Counters.root");
  chains_preSel[16]->Add("results/PhotonJet_Pt20to30/*Counters.root");
  chains_preSel[17]->Add("results/PhotonJet_Pt30to50/*Counters.root");
  chains_preSel[18]->Add("results/PhotonJet_Pt50to80/*Counters.root");
  chains_preSel[19]->Add("results/PhotonJet_Pt80to120/*Counters.root");
  chains_preSel[20]->Add("results/PhotonJet_Pt120to170/*Counters.root");
  chains_preSel[21]->Add("results/PhotonJet_Pt170to300/*Counters.root");
  chains_preSel[22]->Add("results/PhotonJet_Pt300to500/*Counters.root");
  chains_preSel[23]->Add("results/PhotonJet_Pt500toInf/*Counters.root");
  chains_preSel[24]->Add("results/InclusiveMu15/*Counters.root");
  // backgrounds (alternative to some of the above ones)
  chains_preSel[25]->Add("results/TTbar/*Counters.root");       

  // signal
  chains_fullSel[0]->Add("results/H165/*Counters.root");       
  // backgrounds
  chains_fullSel[1]->Add("results/WjetsMadgraph/*Counters.root");       
  chains_fullSel[2]->Add("results/TTbarJetsMadgraph/*Counters.root");       
  chains_fullSel[3]->Add("results/ZjetsMadgraph/*Counters.root");       
  chains_fullSel[4]->Add("results/WW/*Counters.root");   
  chains_fullSel[5]->Add("results/ZZ/*Counters.root");   
  chains_fullSel[6]->Add("results/WZ/*Counters.root");
  chains_fullSel[7]->Add("results/Wgamma/*Counters.root");
  chains_fullSel[8]->Add("results/QCD_EMEnriched_Pt20to30Ele10/*Counters.root");
  chains_fullSel[9]->Add("results/QCD_EMEnriched_Pt30to80Ele10/*Counters.root");
  chains_fullSel[10]->Add("results/QCD_EMEnriched_Pt80to170Ele10/*Counters.root");
  chains_fullSel[11]->Add("results/QCD_BCtoE_Pt20to30Ele10/*Counters.root");
  chains_fullSel[12]->Add("results/QCD_BCtoE_Pt30to80Ele10/*Counters.root");
  chains_fullSel[13]->Add("results/QCD_BCtoE_Pt80to170Ele10/*Counters.root");
  chains_fullSel[14]->Add("results/PhotonJet_Pt0to15/*Counters.root");
  chains_fullSel[15]->Add("results/PhotonJet_Pt15to20/*Counters.root");
  chains_fullSel[16]->Add("results/PhotonJet_Pt20to30/*Counters.root");
  chains_fullSel[17]->Add("results/PhotonJet_Pt30to50/*Counters.root");
  chains_fullSel[18]->Add("results/PhotonJet_Pt50to80/*Counters.root");
  chains_fullSel[19]->Add("results/PhotonJet_Pt80to120/*Counters.root");
  chains_fullSel[20]->Add("results/PhotonJet_Pt120to170/*Counters.root");
  chains_fullSel[21]->Add("results/PhotonJet_Pt170to300/*Counters.root");
  chains_fullSel[22]->Add("results/PhotonJet_Pt300to500/*Counters.root");
  chains_fullSel[23]->Add("results/PhotonJet_Pt500toInf/*Counters.root");
  chains_fullSel[24]->Add("results/InclusiveMu15/*Counters.root");
  // backgrounds (alternative to some of the above ones)
  chains_fullSel[25]->Add("results/TTbar/*Counters.root");       

  float nPreSelTot[10][26];
  float nFullSelTot[20][26];

  for(int isample=0; isample<26; isample++) {
    for(int icut=0; icut<10; icut++) { nPreSelTot[icut][isample]  = 0.0; }
    for(int icut=0; icut<20; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }

  // preselections
  int nCutsAnaPre  = 10;
  int nCutsAnaFull = 20;
  for(int isample=0; isample<26; isample++) {

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

      float qcd_em_tmp=0.;
      qcd_em_tmp += lumi * QCD_EMenriched_Pt20to30_xsec * nPreSelTot[icut][8]/nPreSelTot[0][8];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt30to80_xsec * nPreSelTot[icut][9]/nPreSelTot[0][9];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt80to170_xsec * nPreSelTot[icut][10]/nPreSelTot[0][10];
      QCDem_preSel[icut] = qcd_em_tmp;
      
      float qcd_bc_tmp=0.;
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt20to30_xsec * nPreSelTot[icut][11]/nPreSelTot[0][11];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt30to80_xsec * nPreSelTot[icut][12]/nPreSelTot[0][12];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt80to170_xsec * nPreSelTot[icut][13]/nPreSelTot[0][13];
      QCDbc_preSel[icut] = qcd_bc_tmp;

      float photj_tmp=0.;
      photj_tmp += lumi * PhotonJet_Pt0to15_xsec * nPreSelTot[icut][14]/nPreSelTot[0][14];
      photj_tmp += lumi * PhotonJet_Pt15to20_xsec * nPreSelTot[icut][15]/nPreSelTot[0][15];
      photj_tmp += lumi * PhotonJet_Pt20to30_xsec * nPreSelTot[icut][16]/nPreSelTot[0][16];
      photj_tmp += lumi * PhotonJet_Pt30to50_xsec * nPreSelTot[icut][17]/nPreSelTot[0][17];
      photj_tmp += lumi * PhotonJet_Pt50to80_xsec * nPreSelTot[icut][18]/nPreSelTot[0][18];
      photj_tmp += lumi * PhotonJet_Pt80to120_xsec * nPreSelTot[icut][19]/nPreSelTot[0][19];
      photj_tmp += lumi * PhotonJet_Pt120to170_xsec * nPreSelTot[icut][20]/nPreSelTot[0][20];
      photj_tmp += lumi * PhotonJet_Pt170to300_xsec * nPreSelTot[icut][21]/nPreSelTot[0][21];
      photj_tmp += lumi * PhotonJet_Pt300to500_xsec * nPreSelTot[icut][22]/nPreSelTot[0][22];
      photj_tmp += lumi * PhotonJet_Pt500toInf_xsec * nPreSelTot[icut][23]/nPreSelTot[0][23];
      Photj_preSel[icut] = photj_tmp;

      QCDmu_preSel[icut] = lumi * InclusiveMu15_xsec * nPreSelTot[icut][24]/nPreSelTot[0][24];

      ttbar_preSel[icut] = lumi * TTbar_xsec * nPreSelTot[icut][25]/nPreSelTot[0][25];

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
    if(icut>0 && QCDem_preSel[icut-1]>0) QCDem_eff_preSel[icut] = QCDem_preSel[icut]/QCDem_preSel[icut-1];
    else QCDem_eff_preSel[icut] = 0.0;
    if(icut>0 && QCDbc_preSel[icut-1]>0) QCDbc_eff_preSel[icut] = QCDbc_preSel[icut]/QCDbc_preSel[icut-1];
    else QCDbc_eff_preSel[icut] = 0.0;
    if(icut>0 && Photj_preSel[icut-1]>0) Photj_eff_preSel[icut] = Photj_preSel[icut]/Photj_preSel[icut-1];
    else Photj_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][24]>0) QCDmu_eff_preSel[icut] = nPreSelTot[icut][24]/nPreSelTot[icut-1][24];
    else QCDmu_eff_preSel[icut] = 0.0;
    if(icut>0 && nPreSelTot[icut-1][25]>0) ttbar_eff_preSel[icut] = nPreSelTot[icut][25]/nPreSelTot[icut-1][25];
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
  if(QCDem_preSel[0]>0) QCDem_finaleff_preSel = QCDem_preSel[nCutsAnaPre-2]/QCDem_preSel[0];
  else QCDem_finaleff_preSel = 0.0;
  if(QCDbc_preSel[0]>0) QCDbc_finaleff_preSel = QCDbc_preSel[nCutsAnaPre-2]/QCDbc_preSel[0];
  else QCDbc_finaleff_preSel = 0.0;
  if(Photj_preSel[0]>0) Photj_finaleff_preSel = Photj_preSel[nCutsAnaPre-2]/Photj_preSel[0];
  else Photj_finaleff_preSel = 0.0;
  if(nPreSelTot[0][24]>0) QCDmu_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][24]/nPreSelTot[0][24];
  else QCDmu_finaleff_preSel = 0.0;
  if(nPreSelTot[0][25]>0) ttbar_finaleff_preSel = nPreSelTot[nCutsAnaPre-2][25]/nPreSelTot[0][25];
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

      float qcd_em_tmp=0.;
      qcd_em_tmp += lumi * QCD_EMenriched_Pt20to30_xsec * nFullSelTot[icut][8]/nPreSelTot[0][8];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt30to80_xsec * nFullSelTot[icut][9]/nPreSelTot[0][9];
      qcd_em_tmp += lumi * QCD_EMenriched_Pt80to170_xsec * nFullSelTot[icut][10]/nPreSelTot[0][10];
      QCDem_fullSel[icut] = qcd_em_tmp;
      
      float qcd_bc_tmp=0.;
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt20to30_xsec * nFullSelTot[icut][11]/nPreSelTot[0][11];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt30to80_xsec * nFullSelTot[icut][12]/nPreSelTot[0][12];
      qcd_bc_tmp += lumi * QCD_BCtoE_Pt80to170_xsec * nFullSelTot[icut][13]/nPreSelTot[0][13];
      QCDbc_fullSel[icut] = qcd_bc_tmp;

      float photj_tmp=0.;
      photj_tmp += lumi * PhotonJet_Pt0to15_xsec * nFullSelTot[icut][14]/nPreSelTot[0][14];
      photj_tmp += lumi * PhotonJet_Pt15to20_xsec * nFullSelTot[icut][15]/nPreSelTot[0][15];
      photj_tmp += lumi * PhotonJet_Pt20to30_xsec * nFullSelTot[icut][16]/nPreSelTot[0][16];
      photj_tmp += lumi * PhotonJet_Pt30to50_xsec * nFullSelTot[icut][17]/nPreSelTot[0][17];
      photj_tmp += lumi * PhotonJet_Pt50to80_xsec * nFullSelTot[icut][18]/nPreSelTot[0][18];
      photj_tmp += lumi * PhotonJet_Pt80to120_xsec * nFullSelTot[icut][19]/nPreSelTot[0][19];
      photj_tmp += lumi * PhotonJet_Pt120to170_xsec * nFullSelTot[icut][20]/nPreSelTot[0][20];
      photj_tmp += lumi * PhotonJet_Pt170to300_xsec * nFullSelTot[icut][21]/nPreSelTot[0][21];
      photj_tmp += lumi * PhotonJet_Pt300to500_xsec * nFullSelTot[icut][22]/nPreSelTot[0][22];
      photj_tmp += lumi * PhotonJet_Pt500toInf_xsec * nFullSelTot[icut][23]/nPreSelTot[0][23];
      Photj_fullSel[icut] = photj_tmp;

      QCDmu_fullSel[icut] = lumi * InclusiveMu15_xsec * nFullSelTot[icut][24]/nPreSelTot[0][24];

      ttbar_fullSel[icut] = lumi * TTbar_xsec * nFullSelTot[icut][25]/nPreSelTot[0][25];
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

  textfile << "\\begin{sidewaystable}[p]" << endl
           << "\\begin{tiny}" << endl
           << "\\begin{center}" << endl;
  if(!strcmp(finalstate,"EE")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  if(!strcmp(finalstate,"EM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  if(!strcmp(finalstate,"MM")) textfile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  if(!strcmp(finalstate,"EE")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  if(!strcmp(finalstate,"EM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) & QCD(e.m.) & QCD(b,c) & $\\gamma$+jets \t\\\\" << endl;
  if(!strcmp(finalstate,"MM")) textfile << "selection & H(WW) & W$(l\\nu)$+jets & $t\\bar{t}$ & Z$(ll)$+jets & WW & ZZ & WZ & W$\\gamma$ & QCD($\\mu$) \t\\\\" << endl;
  textfile << "\\hline" << endl; 
    
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
             << ZZ_preSel[icut]      << " (" << 100. * ZZ_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << WZ_preSel[icut]      << " (" << 100. * WZ_eff_preSel[icut]     << "\\%)" << "\t&\t"
             << Wgamma_preSel[icut]  << " (" << 100. * Wgamma_eff_preSel[icut] << "\\%)" << "\t&\t";
    if(!strcmp(finalstate,"MM")) textfile << QCDmu_preSel[icut]   << " (" << 100. * QCDmu_eff_preSel[icut]  << "\\%)";
    if(!strcmp(finalstate,"EM")) textfile << QCDmu_preSel[icut]   << " (" << 100. * QCDmu_eff_preSel[icut]  << "\\%)" << "\t&\t";
    if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
      textfile << QCDem_preSel[icut]   << " (" << 100. * QCDem_eff_preSel[icut]  << "\\%)" << "\t&\t"
               << QCDbc_preSel[icut]   << " (" << 100. * QCDbc_eff_preSel[icut]  << "\\%)" << "\t&\t"
               << Photj_preSel[icut]   << " (" << 100. * Photj_eff_preSel[icut]  << "\\%)";
    }
    textfile << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total preselection " << "\t&\t"
           << H_preSel[8]      << " (" << 100. * H_finaleff_preSel  << "\\%)"  << "\t&\t"
           << Wj_preSel[8]     << " (" << 100. * Wj_finaleff_preSel << "\\%)"  << "\t&\t"
           << ttj_preSel[8]    << " (" << 100. * ttj_finaleff_preSel << "\\%)" << "\t&\t"
           << Zj_preSel[8]     << " (" << 100. * Zj_finaleff_preSel << "\\%)"  << "\t&\t"
           << WW_preSel[8]     << " (" << 100. * WW_finaleff_preSel << "\\%)"  << "\t&\t"
           << ZZ_preSel[8]     << " (" << 100. * ZZ_finaleff_preSel << "\\%)"  << "\t&\t"
           << WZ_preSel[8]     << " (" << 100. * WZ_finaleff_preSel << "\\%)"  << "\t&\t"
           << Wgamma_preSel[8] << " (" << 100. * Wgamma_finaleff_preSel << "\\%)"  << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << QCDmu_preSel[8]   << " (" << 100. * QCDmu_eff_preSel[8]  << "\\%)";
  if(!strcmp(finalstate,"EM")) textfile << QCDmu_preSel[8]   << " (" << 100. * QCDmu_eff_preSel[8]  << "\\%)" << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_preSel[8]  << " (" << 100. * QCDem_finaleff_preSel << "\\%)"  << "\t&\t"
             << QCDbc_preSel[8]  << " (" << 100. * QCDbc_finaleff_preSel << "\\%)"  << "\t&\t"
             << Photj_preSel[8]  << " (" << 100. * Photj_finaleff_preSel << "\\%)";
  }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl;
    
  for(int icut=0; icut<17; icut++) {
    textfile << fullSelCuts[icut] << "\t&\t";
      
    textfile << fixed
             << H_fullSel[icut]      << " (" << 100. * H_eff_fullSel[icut]   << "\\%)" << "\t&\t"
             << Wj_fullSel[icut]     << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << ttj_fullSel[icut]    << " (" << 100. * ttj_eff_fullSel[icut] << "\\%)" << "\t&\t"
             << Zj_fullSel[icut]     << " (" << 100. * Zj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << WW_fullSel[icut]     << " (" << 100. * WW_eff_fullSel[icut]  << "\\%)" << "\t&\t"
             << ZZ_fullSel[icut]     << " (" << 100. * ZZ_eff_fullSel[icut]     << "\\%)" << "\t&\t"
             << WZ_fullSel[icut]     << " (" << 100. * WZ_eff_fullSel[icut]     << "\\%)" << "\t&\t"
             << Wgamma_fullSel[icut] << " (" << 100. * Wgamma_eff_fullSel[icut] << "\\%)" << "\t&\t";
    if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[icut]   << " (" << 100. * QCDmu_eff_fullSel[icut]  << "\\%)";
    if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[icut]   << " (" << 100. * QCDmu_eff_fullSel[icut]  << "\\%)" << "\t&\t";
    if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
      textfile << QCDem_fullSel[icut]  << " (" << 100. * QCDem_eff_fullSel[icut]  << "\\%)" << "\t&\t"
               << QCDbc_fullSel[icut]  << " (" << 100. * QCDbc_eff_fullSel[icut]  << "\\%)" << "\t&\t"
               << Photj_fullSel[icut]  << " (" << 100. * Photj_eff_fullSel[icut]  << "\\%)";
    }
    textfile << "\t\\\\" << endl;
  }
    
  textfile << "\\hline" << endl;

  textfile << "total fullselection " << "\t&\t"
           << H_fullSel[16]      << " (" << 100. * H_finaleff_fullSel  << "\\%)"     << "\t&\t"
           << Wj_fullSel[16]     << " (" << 100. * Wj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ttj_fullSel[16]    << " (" << 100. * ttj_finaleff_fullSel << "\\%)"    << "\t&\t"
           << Zj_fullSel[16]     << " (" << 100. * Zj_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WW_fullSel[16]     << " (" << 100. * WW_finaleff_fullSel << "\\%)"     << "\t&\t"
           << ZZ_fullSel[16]     << " (" << 100. * ZZ_finaleff_fullSel << "\\%)"     << "\t&\t"
           << WZ_fullSel[16]     << " (" << 100. * WZ_finaleff_fullSel << "\\%)"   << "\t&\t"
           << Wgamma_fullSel[16] << " (" << 100. * Wgamma_finaleff_fullSel << "\\%)" << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)";
  if(!strcmp(finalstate,"EM")) textfile << " ( " << 100 * QCDmu_finaleff_fullSel << "\\%)" << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_fullSel[16]  << " (" << 100. * QCDem_finaleff_fullSel << "\\%)" << "\t&\t"
             << QCDbc_fullSel[16]  << " (" << 100. * QCDbc_finaleff_fullSel << "\\%)" << "\t&\t"
             << Photj_fullSel[16]  << " (" << 100. * Photj_finaleff_fullSel << "\\%)";
  }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl;
    
  textfile << "total " << "\t&\t"
           << H_fullSel[16]      << " (" << 100. * H_finaleff   << "\\%)"    << "\t&\t"
           << Wj_fullSel[16]     << " (" << 100. * Wj_finaleff  << "\\%)"    << "\t&\t"
           << ttj_fullSel[16]    << " (" << 100. * ttj_finaleff << "\\%)"    << "\t&\t"
           << Zj_fullSel[16]     << " (" << 100. * Zj_finaleff << "\\%)"     << "\t&\t"
           << WW_fullSel[16]     << " (" << 100. * WW_finaleff << "\\%)"     << "\t&\t"
           << ZZ_fullSel[16]     << " (" << 100. * ZZ_finaleff << "\\%)"     << "\t&\t"
           << WZ_fullSel[16]     << " (" << 100. * WZ_finaleff << "\\%)"     << "\t&\t"
           << Wgamma_fullSel[16] << " (" << 100. * Wgamma_finaleff << "\\%)" << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[16] << " (" << 100. * QCDmu_finaleff << "\\%)";
  if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[16] << " (" << 100. * QCDmu_finaleff << "\\%)" << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_fullSel[16]  << " (" << 100. * QCDem_finaleff << "\\%)" << "\t&\t"
             << QCDbc_fullSel[16]  << " (" << 100. * QCDbc_finaleff << "\\%)" << "\t&\t"
             << Photj_fullSel[16]  << " (" << 100. * Photj_finaleff << "\\%)";
  }
  textfile << "\t\\\\" << endl;
    
  textfile << "0 jets bin " << "\t&\t"
           << H_fullSel[17]      << "\t&\t"
           << Wj_fullSel[17]     << "\t&\t"
           << ttj_fullSel[17]    << "\t&\t"
           << Zj_fullSel[17]     << "\t&\t"
           << WW_fullSel[17]     << "\t&\t"
           << ZZ_fullSel[17]     << "\t&\t"
           << WZ_fullSel[17]     << "\t&\t"
           << Wgamma_fullSel[17] << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[17];
  if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[17] << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_fullSel[17]  << "\t&\t"
             << QCDbc_fullSel[17]  << "\t&\t"
             << Photj_fullSel[17];
  }
  textfile << "\t\\\\" << endl;

  textfile << "1 jets bin " << "\t&\t"
           << H_fullSel[18]      << "\t&\t"
           << Wj_fullSel[18]     << "\t&\t"
           << ttj_fullSel[18]    << "\t&\t"
           << Zj_fullSel[18]     << "\t&\t"
           << WW_fullSel[18]     << "\t&\t"
           << ZZ_fullSel[18]   << "\t&\t"
           << WZ_fullSel[18]   << "\t&\t"
           << Wgamma_fullSel[18] << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[18];
  if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[18] << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_fullSel[18]  << "\t&\t"
             << QCDbc_fullSel[18]  << "\t&\t"
             << Photj_fullSel[18];
  }
  textfile << "\t\\\\" << endl;
    
  textfile << "$>1$ jets bin " << "\t&\t"
           << H_fullSel[19]      << "\t&\t"
           << Wj_fullSel[19]     << "\t&\t"
           << ttj_fullSel[19]    << "\t&\t"
           << Zj_fullSel[19]     << "\t&\t"
           << WW_fullSel[19]     << "\t&\t"
           << ZZ_fullSel[19]   << "\t&\t"
           << WZ_fullSel[19]   << "\t&\t"
           << Wgamma_fullSel[19] << "\t&\t";
  if(!strcmp(finalstate,"MM")) textfile << QCDmu_fullSel[19];
  if(!strcmp(finalstate,"EM")) textfile << QCDmu_fullSel[19] << "\t&\t";
  if(!strcmp(finalstate,"EE") || !strcmp(finalstate,"EM")) {
    textfile << QCDem_fullSel[19]  << "\t&\t"
             << QCDbc_fullSel[19]  << "\t&\t"
             << Photj_fullSel[19];
  }
  textfile << "\t\\\\" << endl;
    
  textfile << "\\hline" << endl
           << "\\end{tabular}" << endl
           << "\\end{center}" << endl
           << "\\end{tiny}" << endl
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
  textfile << "\\setlength\\textheight{9.8in}" << endl;
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

