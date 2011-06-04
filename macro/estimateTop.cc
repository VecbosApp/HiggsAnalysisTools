#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass, const char* channel);

void estimateTop() {

  // constants
  float eff_2b = 0.438;
  float eff_2b_err = 0.069; 
  float eff_1b_softmu = 0.053; // MC

  // scale factors for the backgrounds
  float WWDataOverMC = 85./64.; // estimation 150/pb
  float DYDataOverMC = 22.0/11.4; // estimation 190/pb 
  float WjDataTot[4] = {4.49, 6.65, 7.53, 10.16};

  TFile *fileData = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileDY = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileWj = TFile::Open("results/datasets_trees/Wjets_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeTop= (TTree*)fileTop->Get("T1");
  TTree *treeWW= (TTree*)fileWW->Get("T1");
  TTree *treeDY= (TTree*)fileDY->Get("T1");
  TTree *treeWj= (TTree*)fileWj->Get("T1");

  // these used  for the channel-split estimates (not now)
  TH1F *topHEE = new TH1F("topHEE","",50,0,180);
  TH1F *topHMM = new TH1F("topHMM","",50,0,180);
  TH1F *topHEM = new TH1F("topHEM","",50,0,180);
  TH1F *topHME = new TH1F("topHME","",50,0,180);
  TH1F *btagHDataEE = new TH1F("btagHDataEE","",50,0,180);
  TH1F *btagHDataMM = new TH1F("btagHDataMM","",50,0,180);
  TH1F *btagHDataEM = new TH1F("btagHDataEM","",50,0,180);
  TH1F *btagHDataME = new TH1F("btagHDataME","",50,0,180);

  // thses for the ll estimate
  TH1F *btagHDataLL = new TH1F("btagHDataLL","",50,0,180);

  // backgrounds in the tagged region
  TH1F *btagWWHLL = new TH1F("btagWWHLL","",50,0,180);
  TH1F *btagDYHLL = new TH1F("btagDYHLL","",50,0,180);
  TH1F *CJVWjHEE = new TH1F("CJVWjHEE","",50,0,180);
  TH1F *CJVWjHMM = new TH1F("CJVWjHMM","",50,0,180);
  TH1F *CJVWjHEM = new TH1F("CJVWjHEM","",50,0,180);
  TH1F *CJVWjHME = new TH1F("CJVWjHME","",50,0,180);
  TH1F *btagWjHEE = new TH1F("btagWjHEE","",50,0,180);
  TH1F *btagWjHMM = new TH1F("btagWjHMM","",50,0,180);
  TH1F *btagWjHEM = new TH1F("btagWjHEM","",50,0,180);
  TH1F *btagWjHME = new TH1F("btagWjHME","",50,0,180);

  treeTop->Project("topHLL","deltaPhi","WWSel*weight*puweight");
  treeData->Project("btagHDataLL","deltaPhi","jetVeto && bTagTrackCount>2.1");

  treeTop->Project("topHEE","deltaPhi","(WWSel && finalstate==0)*weight*puweight");
  treeTop->Project("topHMM","deltaPhi","(WWSel && finalstate==1)*weight*puweight");
  treeTop->Project("topHEM","deltaPhi","(WWSel && finalstate==2)*weight*puweight");
  treeTop->Project("topHME","deltaPhi","(WWSel && finalstate==3)*weight*puweight");
  treeData->Project("btagHDataEE","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==0");
  treeData->Project("btagHDataMM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==1");
  treeData->Project("btagHDataEM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==2");
  treeData->Project("btagHDataME","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==3");

  treeWW->Project("btagWWHLL","deltaPhi","(jetVeto && bTagTrackCount>2.1)*weight*puweight");
  treeDY->Project("btagDYHLL","deltaPhi","(jetVeto && bTagTrackCount>2.1)*weight*puweight");

  treeWj->Project("btagWjHEE","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==0");
  treeWj->Project("btagWjHMM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==1");
  treeWj->Project("btagWjHEM","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==2");
  treeWj->Project("btagWjHME","deltaPhi","jetVeto && bTagTrackCount>2.1 && finalstate==3");
  treeWj->Project("CJVWjHEE","deltaPhi","bTagTrackCount>2.1 && finalstate==0");
  treeWj->Project("CJVWjHMM","deltaPhi","bTagTrackCount>2.1 && finalstate==1");
  treeWj->Project("CJVWjHEM","deltaPhi","bTagTrackCount>2.1 && finalstate==2");
  treeWj->Project("CJVWjHME","deltaPhi","bTagTrackCount>2.1 && finalstate==3");
  
  // backgrounds in the tagged region (inclusive)
  float effBtagWj[4] = { btagWjHEE->Integral()/CJVWjHEE->Integral(), 
                         btagWjHMM->Integral()/CJVWjHMM->Integral(),
                         btagWjHEM->Integral()/CJVWjHEM->Integral(),
                         btagWjHME->Integral()/CJVWjHME->Integral() };

  float WW_tot = btagWWHLL->Integral() * WWDataOverMC;
  float WW_tot_err = yieldErrPoisson(WW_tot,btagWWHLL->Integral());
  float DY_tot = btagDYHLL->Integral() * DYDataOverMC;
  float DY_tot_err = yieldErrPoisson(DY_tot,btagDYHLL->Integral());
  float Wjets_tot = WjDataTot[ee] * effBtagWj[ee] + WjDataTot[mm] * effBtagWj[mm] +
    WjDataTot[em] * effBtagWj[em] + WjDataTot[me] * effBtagWj[me];
  float Wjets_tot_err = 0.5 * Wjets_tot; // approximation

  float tagBkg_tot = WW_tot + DY_tot + Wjets_tot;
  float tagBkg_tot_err = quadrSum(WW_tot_err,DY_tot_err,Wjets_tot_err);

  std::cout << "--- Background estimations: ---" << std::endl;
  std::cout << "WW = " << WW_tot << " +/-" << WW_tot_err << std::endl;
  std::cout << "DY = " << DY_tot << " +/-" << DY_tot_err << std::endl;
  std::cout << "Wjets = " << Wjets_tot << " +/-" << Wjets_tot_err << std::endl;
  std::cout << "Tot background to tagged events = " << tagBkg_tot << " +/- " << tagBkg_tot_err << std::endl;
  std::cout << "---- end backgrounds ---\n\n" << std::endl;

  ///// TOP ESTIMATION ///////
  float nTopMC[4];
  float nTopMC_err[4];
  
  nTopMC[ee] = topHEE->Integral();
  nTopMC_err[ee] = yieldErrPoisson(nTopMC[ee],topHEE->GetEntries());
  nTopMC[mm] = topHMM->Integral();
  nTopMC_err[mm] = yieldErrPoisson(nTopMC[mm],topHMM->GetEntries());
  nTopMC[em] = topHEM->Integral();
  nTopMC_err[em] = yieldErrPoisson(nTopMC[em],topHEM->GetEntries());
  nTopMC[me] = topHME->Integral();
  nTopMC_err[me] = yieldErrPoisson(nTopMC[me],topHME->GetEntries());
  
  float nTopMC_tot = nTopMC[ee] + nTopMC[mm] + nTopMC[em] + nTopMC[me];
  float nTopMC_tot_err = quadrSum(nTopMC_err[ee],nTopMC_err[mm],nTopMC_err[em],nTopMC_err[me]);

  /// TOP CHANNEL FRACTIONS FROM MC ///
  // they are used in case of low stat to us.e the summed WW estimation in the 4 sub-channels
  std::cout << "Using channel fractions from MC!" << std::endl;
  float frac[4];
  for(int icha=0; icha<4; icha++) {
    frac[icha] = nTopMC[icha]/nTopMC_tot;
    std::cout << "Fraction of Top in channel " << icha << " = " << frac[icha] << std::endl;
  }

  // top estimation from data (0-jet bin method)
  float nTopData[4];
  float nTopData_err[4];
  float eff_2b_softmu = 1 - pow(1-eff_1b_softmu,2);

  // EE
  float nBTagTag_data = btagHDataEE->Integral();
  float nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  float nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  float nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  float nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EE = " << nBTagTag_data << std::endl;

  nTopData[ee] = nTopSoftMuVeto_data;
  nTopData_err[ee] = nTopSoftMuVeto_data_err;

  // MM
  nBTagTag_data = btagHDataMM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "MM = " << nBTagTag_data << std::endl;

  nTopData[mm] = nTopSoftMuVeto_data;
  nTopData_err[mm] = nTopSoftMuVeto_data_err;

  // EM
  nBTagTag_data = btagHDataEM->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "EM = " << nBTagTag_data << std::endl;

  nTopData[em] = nTopSoftMuVeto_data;
  nTopData_err[em] = nTopSoftMuVeto_data_err;

  // ME
  nBTagTag_data = btagHDataME->Integral();
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  std::cout << "Tagged events in data: " << "ME = " << nBTagTag_data << std::endl;

  float nBTagTag_data_tot =  btagHDataEE->Integral() +  btagHDataMM->Integral() +  btagHDataEM->Integral() +  btagHDataME->Integral();

  nTopData[me] = nTopSoftMuVeto_data;
  nTopData_err[me] = nTopSoftMuVeto_data_err;


  // here summing up the separate ones
  //   float nTopData_tot = nTopData[ee] + nTopData[mm] + nTopData[em] + nTopData[me];
  //   float nTopData_tot_err = quadrSum(nTopData_err[ee],nTopData_err[mm],nTopData_err[em],nTopData_err[me]);

  // LL
  nBTagTag_data = btagHDataLL->Integral() - tagBkg_tot;
  std::cout << "number of tagged events after bkg subtraction = " << nBTagTag_data << std::endl;
  nTopBTagVeto_data = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).first;
  nTopBTagVeto_data_err = (nVeto(nBTagTag_data, eff_2b, eff_2b_err)).second; 
  nTopSoftMuVeto_data = nTopBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  nTopSoftMuVeto_data_err = nTopBTagVeto_data_err * (1-eff_2b_softmu); 

  float nTopData_tot = nTopSoftMuVeto_data;
  float nTopData_tot_err = nTopSoftMuVeto_data_err;

  std::cout << "Number of Top events from data at W+W- level: " << std::endl
    // commented because the bkg estimation is done only for the cumulative sample
//             << "*\tEE = " << nTopData[ee] << " +/- " << nTopData_err[ee] << std::endl
//             << "*\tMM = " << nTopData[mm] << " +/- " << nTopData_err[mm] << std::endl
//             << "*\tEM = " << nTopData[em] << " +/- " << nTopData_err[em] << std::endl
//             << "*\tME = " << nTopData[me] << " +/- " << nTopData_err[me] << std::endl
            << "*\tTOT = " << nTopData_tot << " +/- " << nTopData_tot_err << std::endl;

  std::cout << "Number of Top events from MC at W+W- level: " << std::endl
//             << "*\tEE = " << nTopMC[ee] << " +/- " << nTopMC_err[ee] << std::endl
//             << "*\tMM = " << nTopMC[mm] << " +/- " << nTopMC_err[mm] << std::endl
//             << "*\tEM = " << nTopMC[em] << " +/- " << nTopMC_err[em] << std::endl
//             << "*\tME = " << nTopMC[me] << " +/- " << nTopMC_err[me] << std::endl
            << "*\tTOT = " << nTopMC_tot << " +/- " << nTopMC_tot_err << std::endl;

  ofstream textfile;
  textfile.open("TopYieldsData.txt", ios_base::trunc);
  textfile.precision(2);

  ofstream tablefile1;
  tablefile1.open("TopYieldsData_ForTable_0j.txt", ios_base::trunc);
  tablefile1.precision(2);

  ofstream tablefile2;
  tablefile2.open("TopYieldsData_ForTable_1j.txt", ios_base::trunc);
  tablefile2.precision(2);

  ofstream tablefile3;
  tablefile3.open("TopYieldsMC_ForTable_0j.txt", ios_base::trunc);
  tablefile3.precision(2);

  ofstream tablefile4;
  tablefile4.open("TopYieldsMC_ForTable_1j.txt", ios_base::trunc);
  tablefile4.precision(2);

  // these are for limits
  ofstream cardfile[4][2]; //[cha][jetbin]
  for(int icha=0; icha<4; icha++) {
    for(int j=0; j<2; j++) {
      char fileName[2];
      if(icha==ee) sprintf(fileName,"TopCard_ee_%dj.txt",j);
      if(icha==mm) sprintf(fileName,"TopCard_mm_%dj.txt",j);
      if(icha==em) sprintf(fileName,"TopCard_em_%dj.txt",j);
      if(icha==me) sprintf(fileName,"TopCard_me_%dj.txt",j);
      cardfile[icha][j].open(fileName, ios_base::trunc);
      cardfile[icha][j].precision(2);
    }
  }
    

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    countEvents(mass,"EE");
    countEvents(mass,"MM");
    countEvents(mass,"EM");
    countEvents(mass,"ME");

    float nTopData_HiggsSel_0j[4], nTopData_HiggsSel_1j[4];
    float nTopData_HiggsSel_0j_err[4], nTopData_HiggsSel_1j_err[4];

    float nTopMC_HiggsSel_0j[4],  nTopMC_HiggsSel_1j[4];
    float nTopMC_HiggsSel_0j_err[4], nTopMC_HiggsSel_1j_err[4];

    for(int icha=0;icha<4;icha++) {
      float eff_0j = (nEv_endWW[icha]==0) ? 0. : nEv_end0j[icha]/nEv_endWW[icha];
      float eff_0j_err = (nEv_endWW[icha]==0) ? 0. : sqrt(eff_0j*(1-eff_0j)/nEv_endWW[icha]);
      float effErrRel = (eff_0j==0) ? 0. : eff_0j_err/eff_0j;
    
      // this is the correct esztimation for when we have sufficient stat
      // nTopData_HiggsSel_0j[icha] = nTopData[icha] * eff_0j;
      // float topDataErrRel = (nTopData[icha]==0) ? 0. : nTopData_err[icha]/nTopData[icha];
      // nTopData_HiggsSel_0j_err[icha] = nTopData_HiggsSel_0j[icha] * quadrSum(topDataErrRel,effErrRel);
      nTopData_HiggsSel_0j[icha] = nTopData_tot * frac[icha] * eff_0j;
      nTopData_HiggsSel_0j_err[icha] = nTopData_HiggsSel_0j[icha] * quadrSum(nTopData_tot_err/nTopData_tot,eff_0j_err/eff_0j);

      nTopMC_HiggsSel_0j[icha] = nTopMC[icha] * eff_0j;
      float topMCErrRel = (nTopMC[icha]==0) ? 0. : nTopMC_err[icha]/nTopMC[icha];
      nTopMC_HiggsSel_0j_err[icha] = nTopMC_HiggsSel_0j[icha] * quadrSum(topMCErrRel,effErrRel);

      float eff_1j = nEv_end1j[icha]/nEv_endWW[icha];
      float eff_1j_err = sqrt(eff_1j*(1-eff_1j)/nEv_endWW[icha]);
      effErrRel = (eff_1j==0) ? 0. : eff_1j_err/eff_1j;

      // this is the correct esztimation for when we have sufficient stat
      // nTopData_HiggsSel_1j[icha] = nTopData[icha] * eff_1j;
      // topDataErrRel = (nTopData[icha]==0) ? 0. : nTopData_err[icha]/nTopData[icha];
      // nTopData_HiggsSel_1j_err[icha] = nTopData_HiggsSel_1j[icha] * quadrSum(topDataErrRel,effErrRel);
      nTopData_HiggsSel_1j[icha] = nTopData_tot * frac[icha] * eff_1j;
      nTopData_HiggsSel_1j_err[icha] =  nTopData_HiggsSel_1j[icha] * quadrSum(nTopData_tot_err/nTopData_tot,eff_1j_err/eff_1j);      

      nTopMC_HiggsSel_1j[icha] = nTopMC[icha] * eff_1j;
      topMCErrRel = (nTopMC[icha]==0) ? 0. : nTopMC_err[icha]/nTopMC[icha];
      nTopMC_HiggsSel_1j_err[icha] = nTopMC_HiggsSel_1j[icha] * quadrSum(topMCErrRel,effErrRel);

      char channelName[2];
      if(icha==ee) sprintf(channelName,"EE");
      if(icha==mm) sprintf(channelName,"MM");
      if(icha==em) sprintf(channelName,"EM");
      if(icha==me) sprintf(channelName,"ME");
      
      // for Giovanni
      float alpha_0j = frac[icha] * (1-eff_2b)/eff_2b * (1-eff_2b_softmu) * eff_0j;
      float alpha_0j_err = alpha_0j * quadrSum(eff_2b_err/pow(eff_2b,2), eff_0j_err/eff_0j); // assuming no err on softmu

      float alpha_1j = frac[icha] * (1-eff_2b)/eff_2b * (1-eff_2b_softmu) * eff_1j;
      float alpha_1j_err = alpha_1j * quadrSum(eff_2b_err/pow(eff_2b,2), eff_1j_err/eff_1j); // assuming no err on softmu

      cardfile[icha][0] << mass 
                        << "\t" << nBTagTag_data_tot << "\t" << alpha_0j
                        << "\t" <<  alpha_0j_err 
                        << std::endl;

      cardfile[icha][1] << mass 
                        << "\t" << nBTagTag_data_tot << "\t" << alpha_1j
                        << "\t" <<  alpha_1j_err 
                        << std::endl;

      ///////////////

      textfile << channelName << ": Higgs Mass = " << mass 
               << "\tdata 0 jet = " << nTopData_HiggsSel_0j[icha] << " +/- " << nTopData_HiggsSel_0j_err[icha] 
               << "\tdata 1 jet = " << nTopData_HiggsSel_1j[icha] << " +/- " << nTopData_HiggsSel_1j_err[icha] 
               << "\tMC 0 jet = " << nTopMC_HiggsSel_0j[icha] << " +/- " << nTopMC_HiggsSel_0j_err[icha] 
               << "\tMC 1 jet = " << nTopMC_HiggsSel_1j[icha] << " +/- " << nTopMC_HiggsSel_1j_err[icha]
               << std::endl;
    }

    // summary table for limits
    if (i==0) { 
      tablefile1 << "# zero jets bin data" << endl;
      tablefile1 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile1 << mass 
	       << " " << "\t\t" << nTopData_HiggsSel_0j[1] << " +/- " <<  nTopData_HiggsSel_0j_err[1] 
	       << " " << "\t\t" << nTopData_HiggsSel_0j[3] << " +/- " <<  nTopData_HiggsSel_0j_err[3] 
	       << " " << "\t\t" << nTopData_HiggsSel_0j[2] << " +/- " <<  nTopData_HiggsSel_0j_err[2] 
	       << " " << "\t\t" << nTopData_HiggsSel_0j[0] << " +/- " <<  nTopData_HiggsSel_0j_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile2 << "#one jet bin data" << endl;
      tablefile2 << "#\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile2 << mass 
	       << " " << "\t\t" << nTopData_HiggsSel_1j[1] << " +/- " <<  nTopData_HiggsSel_1j_err[1] 
	       << " " << "\t\t" << nTopData_HiggsSel_1j[3] << " +/- " <<  nTopData_HiggsSel_1j_err[3] 
	       << " " << "\t\t" << nTopData_HiggsSel_1j[2] << " +/- " <<  nTopData_HiggsSel_1j_err[2] 
	       << " " << "\t\t" << nTopData_HiggsSel_1j[0] << " +/- " <<  nTopData_HiggsSel_1j_err[0] 
	       << std::endl;

    if (i==0) { 
      tablefile3 << "# zero jets bin MC" << endl;
      tablefile3 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile3 << mass 
	       << " " << "\t\t" << nTopMC_HiggsSel_0j[1] << " +/- " <<  nTopMC_HiggsSel_0j_err[1] 
	       << " " << "\t\t" << nTopMC_HiggsSel_0j[3] << " +/- " <<  nTopMC_HiggsSel_0j_err[3] 
	       << " " << "\t\t" << nTopMC_HiggsSel_0j[2] << " +/- " <<  nTopMC_HiggsSel_0j_err[2] 
	       << " " << "\t\t" << nTopMC_HiggsSel_0j[0] << " +/- " <<  nTopMC_HiggsSel_0j_err[0] 
	       << std::endl;
    
    if (i==0) { 
      tablefile4 << "#one jet bin MC" << endl;
      tablefile4 << "#\t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile4 << mass 
	       << " " << "\t\t" << nTopMC_HiggsSel_1j[1] << " +/- " <<  nTopMC_HiggsSel_1j_err[1] 
	       << " " << "\t\t" << nTopMC_HiggsSel_1j[3] << " +/- " <<  nTopMC_HiggsSel_1j_err[3] 
	       << " " << "\t\t" << nTopMC_HiggsSel_1j[2] << " +/- " <<  nTopMC_HiggsSel_1j_err[2] 
	       << " " << "\t\t" << nTopMC_HiggsSel_1j[0] << " +/- " <<  nTopMC_HiggsSel_1j_err[0] 
	       << std::endl;

    float nTopData_HiggsSel_0j_Tot = nTopData_HiggsSel_0j[ee] + nTopData_HiggsSel_0j[mm] + nTopData_HiggsSel_0j[em] + nTopData_HiggsSel_0j[me];
    float nTopData_HiggsSel_0j_Tot_err = quadrSum(nTopData_HiggsSel_0j_err[ee],nTopData_HiggsSel_0j_err[mm],nTopData_HiggsSel_0j_err[em],nTopData_HiggsSel_0j_err[me]);

    float nTopMC_HiggsSel_0j_Tot = nTopMC_HiggsSel_0j[ee] + nTopMC_HiggsSel_0j[mm] + nTopMC_HiggsSel_0j[em] + nTopMC_HiggsSel_0j[me];
    float nTopMC_HiggsSel_0j_Tot_err = quadrSum(nTopMC_HiggsSel_0j_err[ee],nTopMC_HiggsSel_0j_err[mm],nTopMC_HiggsSel_0j_err[em],nTopMC_HiggsSel_0j_err[me]);

    float nTopData_HiggsSel_1j_Tot = nTopData_HiggsSel_1j[ee] + nTopData_HiggsSel_1j[mm] + nTopData_HiggsSel_1j[em] + nTopData_HiggsSel_1j[me];
    float nTopData_HiggsSel_1j_Tot_err = quadrSum(nTopData_HiggsSel_1j_err[ee],nTopData_HiggsSel_1j_err[mm],nTopData_HiggsSel_1j_err[em],nTopData_HiggsSel_1j_err[me]);

    float nTopMC_HiggsSel_1j_Tot = nTopMC_HiggsSel_1j[ee] + nTopMC_HiggsSel_1j[mm] + nTopMC_HiggsSel_1j[em] + nTopMC_HiggsSel_1j[me];
    float nTopMC_HiggsSel_1j_Tot_err = quadrSum(nTopMC_HiggsSel_1j_err[ee],nTopMC_HiggsSel_1j_err[mm],nTopMC_HiggsSel_1j_err[em],nTopMC_HiggsSel_1j_err[me]);
    
    textfile.precision(2);
    textfile << "\t===>> TOTAL: Higgs Mass = " << mass 
             << "\tdata 0 jet = " << nTopData_HiggsSel_0j_Tot << " +/- " << nTopData_HiggsSel_0j_Tot_err 
             << "\tdata 1 jet = " << nTopData_HiggsSel_1j_Tot << " +/- " << nTopData_HiggsSel_1j_Tot_err
             << "\tMC 0 jet = " << nTopMC_HiggsSel_0j_Tot << " +/- " << nTopMC_HiggsSel_0j_Tot_err 
             << "\tMC 1 jet = " << nTopMC_HiggsSel_1j_Tot << " +/- " << nTopMC_HiggsSel_1j_Tot_err
             << std::endl;

  }

  std::cout << "Full top yields in data in:  TopYieldsData.txt " << std::endl;

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr) {
  float val = ntag * (1-eff2b) / eff2b;
  float err = ntag * eff2berr / pow(eff2b,2);
  return std::make_pair(val,err);
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

void countEvents(int mass, const char *channel) {

  // taking the EE or ME trees for the wanted mass
  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_%s",channel);  
  TChain *theChain = new TChain(nametree);

  // assume that the final selection efficiency is the same for all the top samples and use average of it
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V14/OptimMH%d/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V14/OptimMH%d/Spring11_V2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V14/OptimMH%d/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root",mass);  
  theChain->Add(file_mc);
  sprintf(file_mc,"/cmsrm/pc21_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V14/OptimMH%d/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root",mass);  
  theChain->Add(file_mc);
  //  cout << "reading tree " << nametree << " from file " << file_mc << endl;    
  
  int theCha=-1;
  if(TString(channel).Contains("EE")) theCha=ee;
  if(TString(channel).Contains("MM")) theCha=mm;
  if(TString(channel).Contains("EM")) theCha=em;
  if(TString(channel).Contains("ME")) theCha=me;

  // number of events at the wanted step of the selection
  nEv_endWW[theCha] = 0.0;
  nEv_end0j[theCha] = 0.0;
  nEv_end1j[theCha] = 0.0;

  // reading the tree
  Int_t    nCuts;
  Float_t  nSel[25];   //[nCuts]                                      
  TBranch  *b_nCuts;   
  TBranch  *b_nSel;    

  theChain->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
  theChain->SetBranchAddress("nSel", nSel, &b_nSel);
  Long64_t nentries = theChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = theChain->GetEntry(jentry);   
    nbytes    += nb;
    nEv_endWW[theCha] += nSel[16];
    nEv_end0j[theCha] += nSel[22];
    nEv_end1j[theCha] += nSel[23];
  }
}

