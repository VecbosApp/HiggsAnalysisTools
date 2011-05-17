#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

// numbers filled from counters
float nEv_endWW;
float nEv_end0j;
float nEv_end1j;

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass);

void estimateWW() {

  // constants
  float eff_2b = 0.4983;
  float eff_2b_err = 0.03; 
  float eff_2b_softmu = 0.1677; // MC 

  float eff_softmu_Z = 0.867;

  float Rmm = 0.187398;
  float Rmm_err = 0.00355062;
  float kmm = 0.592288;
  float kmm_err = 0.0139356;
  float Ree = 0.152089;
  float Ree_err = 0.00338336;
  float kee = 0.422092;
  float kee_err = 0.00874687;


  TFile *fileData = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileWW = TFile::Open("results/datasets_trees/WW_ll.root");
  TFile *fileTop = TFile::Open("results/datasets_trees/top_ll.root");
  TFile *fileWjets = TFile::Open("results/datasets_trees/Wjets_ll.root");
  TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");
  TFile *fileDiBosons = TFile::Open("results/datasets_trees/others_ll.root");

  TTree *treeData = (TTree*)fileData->Get("T1");
  TTree *treeWW = (TTree*)fileWW->Get("T1");
  TTree *treeTop= (TTree*)fileTop->Get("T1");
  TTree *treeWjets = (TTree*)fileWjets->Get("T1");
  TTree *treeZjets = (TTree*)fileZjets->Get("T1");
  TTree *treeDiBosons = (TTree*)fileDiBosons->Get("T1");

  std::vector<TFile*> bakgrounds;
  bakgrounds.push_back(fileTop);
  bakgrounds.push_back(fileWjets);
  bakgrounds.push_back(fileZjets);
  bakgrounds.push_back(fileDiBosons);

  std::vector<float> backgroundsUnc; // dummy values so far
  backgroundsUnc.push_back(0.5);
  backgroundsUnc.push_back(0.5);
  backgroundsUnc.push_back(0.1);
  backgroundsUnc.push_back(0.1);
  
  TH1F *dataH = new TH1F("dataH","",50,0,180);
  TH1F *WWCH = new TH1F("WWCH","",50,0,180); // control region
  TH1F *WWH = new TH1F("WWH","",50,0,180); // all
  TH1F *topH = new TH1F("topH","",50,0,180);
  TH1F *btagHData = new TH1F("btagHData","",50,0,180);
  TH1F *WjetsEEH = new TH1F("WjetsEEH","",50,0,180);
  TH1F *WjetsMEH = new TH1F("WjetsMEH","",50,0,180);
  TH1F *WjetsMMH = new TH1F("WjetsMMH","",50,0,180);
  TH1F *WjetsEMH = new TH1F("WjetsEMH","",50,0,180);
  TH1F *ZeejetsH = new TH1F("ZeejetsH","",50,0,180);
  TH1F *ZmmjetsH = new TH1F("ZmmjetsH","",50,0,180);
  TH1F *ZemjetsH = new TH1F("ZemjetsH","",50,0,180);
  TH1F *neeInH = new TH1F("neeInH","",50,0,180);
  TH1F *nmmInH = new TH1F("nmmInH","",50,0,180);
  TH1F *nemInH = new TH1F("nemInH","",50,0,180);
  TH1F *DiBosonsH = new TH1F("DiBosonsH","",50,0,180);

  treeData->Project("dataH","deltaPhi","eleInvMass>100 && WWSel");
  treeWW->Project("WWCH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeWW->Project("WWH","deltaPhi","WWSel*weight*puweight");
  treeTop->Project("topH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");
  treeData->Project("btagHData","deltaPhi","eleInvMass>100 && jetVeto && bTagTrackCount>2.1");
  treeWjets->Project("WjetsEEH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==0)*weight*puweight");
  treeWjets->Project("WjetsMEH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==3)*weight*puweight");
  treeWjets->Project("WjetsMMH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==1)*weight*puweight");
  treeWjets->Project("WjetsEMH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==2)*weight*puweight");
  treeZjets->Project("ZeejetsH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==0)*weight*puweight");
  treeZjets->Project("ZmmjetsH","deltaPhi","(eleInvMass>100 && WWSel && finalstate==1)*weight*puweight");
  treeZjets->Project("ZemjetsH","deltaPhi","(eleInvMass>100 && WWSel && (finalstate==2 || finalstate==3))*weight*puweight");
  treeData->Project("neeInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nmmInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nemInH","deltaPhi","eleInvMass>100 && eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && (finalstate==2 || finalstate==3)"); // missing softmu... not available in red trees... hopefully small contrib here
  treeDiBosons->Project("DiBosonsH","deltaPhi","(eleInvMass>100 && WWSel)*weight*puweight");

  ///// TOP ESTIMATION ///////
  float nTopOut = topH->Integral();
  float nTopOut_err = yieldErrPoisson(nTopOut,topH->GetEntries());

  // top estimation from data (0-jet bin method)
  float nBTagTagOut_data = btagHData->Integral();
  float nTopOutBTagVeto_data = (nVeto(nBTagTagOut_data, eff_2b, eff_2b_err)).first;
  float nTopOutBTagVeto_data_err = (nVeto(nBTagTagOut_data, eff_2b, eff_2b_err)).second; 
  float nTopOutSoftMuVeto_data = nTopOutBTagVeto_data * (1-eff_2b_softmu); // efficiency of passing the soft muon veto (both the b's).
  float nTopOutSoftMuVeto_data_err = nTopOutBTagVeto_data_err * (1-eff_2b_softmu); 


  std::cout << "TOP ESTIMATION..." << std::endl;
  std::cout << "Using eff_2b = " << eff_2b << " +/- " << eff_2b_err << std::endl;
  std::cout << "Using eff_2b_softmu = " << eff_2b_softmu << std::endl;
  std::cout << "Tagged events = " << nBTagTagOut_data
            << "   Top out from data after btag veto = " << nTopOutBTagVeto_data << " +/- " << nTopOutBTagVeto_data_err 
            << "   Top out from data after soft mu veto =  " << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;

  std::cout << "Top from MC = " << nTopOut << " +/- " << nTopOut_err << "\tTop from data = " << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;
  std::cout << "END TOP ESTIMATION." << std::endl;
  ///// END TOP /////////


  /// W + JETS ESTIMATION ////
  // MC estimation
  float nWjetsEEOut = WjetsEEH->Integral();
  float nWjetsEEOut_err = yieldErrPoisson(nWjetsEEOut,WjetsEEH->GetEntries());
  float nWjetsMEOut = WjetsMEH->Integral();
  float nWjetsMEOut_err = yieldErrPoisson(nWjetsMEOut,WjetsMEH->GetEntries());
  float nWjetsMMOut = WjetsMMH->Integral();
  float nWjetsMMOut_err = yieldErrPoisson(nWjetsMMOut,WjetsMMH->GetEntries());
  float nWjetsEMOut = WjetsEMH->Integral();
  float nWjetsEMOut_err = yieldErrPoisson(nWjetsEMOut,WjetsEMH->GetEntries());

  float nWjetsOut = nWjetsEEOut + nWjetsMEOut + nWjetsMMOut + nWjetsEMOut;
  float nWjetsOut_err = quadrSum(nWjetsEEOut_err,nWjetsMEOut_err,nWjetsMMOut_err,nWjetsEMOut_err);

  // data estimation (ee)
  float nWjetsEEOutData = 0.612778;
  float nWjetsEEOutData_err = 0.113244;
  float nWjetsMEOutData = 0.887411;
  float nWjetsMEOutData_err = 0.13546;

  // still MC fort mumu
  float nWjetsOutData = nWjetsEEOutData + nWjetsMEOutData + nWjetsMMOut + nWjetsEMOut;
  float nWjetsOutData_err = quadrSum(nWjetsEEOutData_err,nWjetsMEOutData_err,nWjetsMMOut_err,nWjetsEMOut_err);
  /////////////////////////////
  
  // DY estimation /////
  std::cout << "DY ESTIMATION..." << std::endl;
  float nZeejetsOut = ZeejetsH->Integral();
  float nZeejetsOut_err = yieldErrPoisson(nZeejetsOut,ZeejetsH->GetEntries());
  float nZmmjetsOut = ZmmjetsH->Integral();
  float nZmmjetsOut_err = yieldErrPoisson(nZmmjetsOut,ZmmjetsH->GetEntries());
  float nZemjetsOut = ZemjetsH->Integral();
  float nZemjetsOut_err = yieldErrPoisson(nZemjetsOut,ZemjetsH->GetEntries());

  float neeIn = neeInH->Integral() * eff_softmu_Z;
  float nmmIn = nmmInH->Integral() * eff_softmu_Z;
  float nemIn = nemInH->Integral() * eff_softmu_Z;


  float neeExp = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).first;
  float neeExp_err = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).second;
  float nmmExp = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).first;
  float nmmExp_err = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).second;

  // ee and mm from data and em from MC
  float nemExp = nZemjetsOut;
  float nemExp_err = yieldErrPoisson(nZemjetsOut,ZemjetsH->GetEntries());

  std::cout << "neeIn = " << neeIn << "\tnmmIn = " << nmmIn << "\tnemIn = " << nemIn << std::endl;
  std::cout << "nEE MC = " << nZeejetsOut << "\tData = " << neeExp << " +/- " << neeExp_err << std::endl;
  std::cout << "nMM MC = " << nZmmjetsOut << "\tData = " << nmmExp << " +/- " << nmmExp_err << std::endl;
  std::cout << "nEM MC = " << nemExp << " +/- " << nemExp_err << std::endl; 

  std::cout << "END DY ESTIMATION." << std::endl;
  ////////// END DY ///////////

  float nDiBosonsOut = DiBosonsH->Integral();
  float nDiBosonsOut_err = yieldErrPoisson(nDiBosonsOut,DiBosonsH->GetEntries());
  std::cout << "--> Dibosons from MC = " << nDiBosonsOut << " +/- " << nDiBosonsOut_err << std::endl;
  
  // sum of the backgrounds ///
  // data estimation (where possible)
  float DYTot = neeExp + nmmExp + nemExp;
  float DYTot_err = quadrSum(neeExp_err,nmmExp_err,nemExp_err);
  float bkgTot = nWjetsOutData + nTopOutSoftMuVeto_data + DYTot + nDiBosonsOut;
  float bkgTot_err = quadrSum(nWjetsOutData_err,nTopOutSoftMuVeto_data_err,DYTot_err,nDiBosonsOut_err);

  // MC estimation
  float DYTotMC = nZeejetsOut + nZmmjetsOut + nZemjetsOut;
  float DYTotMC_err = quadrSum(nZeejetsOut_err,nZmmjetsOut_err,nZemjetsOut_err);
  float bkgTotMC = nWjetsOut + nTopOut + DYTotMC + nDiBosonsOut;
  float bkgTotMC_err = quadrSum(nWjetsOut_err,nTopOut_err,DYTotMC_err,nDiBosonsOut_err);

  ///// SUMMARY BKG ////
  std::cout << "---->  BACKGROUND SUMMARY  <-------" << std::endl;
  std::cout << "bkg\t\tMC\t\t\tdata" << std::endl;   
  std::cout.precision(3);
  std::cout << "W+jets =\t" << nWjetsOut << " +/- " << nWjetsOut_err << "\t\t" << nWjetsOutData << " +/- " << nWjetsOutData_err << std::endl;
  std::cout << "top =\t\t" << nTopOut << " +/- " << nTopOut_err << "\t\t" << nTopOutSoftMuVeto_data << " +/-" << nTopOutSoftMuVeto_data_err << std::endl;
  std::cout << "DY =\t\t" << DYTotMC << " +/- " << DYTotMC_err << "\t\t\t" << DYTot << " +/- " << DYTot_err << std::endl;
  std::cout << "WZ,ZZ =\t\t" << nDiBosonsOut << " +/- " << nDiBosonsOut_err << "\tn.a." << std::endl;
  std::cout << "TOTAL:\n\t\t" << bkgTotMC << " +/- " << bkgTotMC_err << "\t\t" << bkgTot << " +/-" << bkgTot_err << std::endl;
  std::cout << "-----------------------------------\n\n\n" << std::endl;
  /////////////////////

  // WW 
  float nDataOut = dataH->Integral();
  // since we extrapolate Out -> In, assign the stat error to Out
  float nDataOut_err = yieldErrPoisson(nDataOut,dataH->GetEntries());
  float nWWOutData = nDataOut - bkgTot;
  float nWWOutData_err = quadrSum(nDataOut_err,bkgTot_err);

  // MC estimation of WW
  float nWWOutMC = WWCH->Integral();
  float nWWOutMC_err = yieldErrPoisson(nWWOutMC,WWCH->GetEntries());
  std::cout << "WW CONTROL REGION:" << std::endl;
  std::cout << "MC = " << nWWOutMC << " +/- " << nWWOutMC_err 
            << "\t\tData = " << nWWOutData << " +/- " << nWWOutData_err << std::endl;


  float dataOmc = nWWOutData/nWWOutMC;
  float dataOmc_err = dataOmc * quadrSum(nWWOutData_err/nWWOutData,nWWOutMC_err/nWWOutMC);

  std::cout << "Scale factor data / MC = " << dataOmc << " +/- " << dataOmc_err << std::endl;

  // Ratio Signal region / Control region
  float nWWMC = WWH->Integral();
  float nWWMC_err = yieldErrPoisson(nWWMC,WWH->GetEntries());
  float RSC = nWWMC/nWWOutMC;
  float RSC_err = RSC * quadrSum(nWWMC_err/nWWMC,nWWOutMC_err/nWWOutMC);

  // The WW for the wole region estimated from data at the WW sel level
  float nWWData_WWSel = RSC * nWWOutData;
  float nWWData_WWSel_err = nWWData_WWSel * quadrSum(RSC_err/RSC, nWWOutData_err/nWWOutData);

  std::cout << "===> WW ESTIMATION AT THE WW SELECTION LEVEL <===" << std::endl;
  std::cout << "MC = " << nWWMC << " +/- " << nWWMC_err << std::endl;
  std::cout << "data = " << nWWData_WWSel << " +/- " << nWWData_WWSel_err << std::endl;
  std::cout << "=================================================" << std::endl; 

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    countEvents(mass);

    float eff_0j = nEv_end0j/nEv_endWW;
    float eff_0j_err = sqrt(eff_0j*(1-eff_0j)/nEv_endWW);
    
    float nWWData_HiggsSel = nWWData_WWSel * eff_0j;
    float nWWData_HiggsSel_err = nWWData_HiggsSel * quadrSum(nWWData_WWSel_err/nWWData_WWSel,eff_0j_err/eff_0j);

    float nWWMC_HiggsSel = nWWMC * eff_0j;
    float nWWMC_HiggsSel_err = nWWMC_HiggsSel * quadrSum(nWWMC_err/nWWMC,eff_0j_err/eff_0j);

    cout << "Higgs Mass = " << mass 
         << "\tMC = " << nWWMC_HiggsSel << " +/- " << nWWMC_HiggsSel_err
         << "\tdata = " << nWWData_HiggsSel << " +/- " << nWWData_HiggsSel_err << std::endl;
  }

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
}

std::pair<float,float> nVeto(float ntag, float eff2b, float eff2berr) {
  float val = ntag * (1-eff2b) / eff2b;
  float err = ntag * eff2berr / pow(eff2b,2);
  return std::make_pair(val,err);
}

std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK) {
  float val = R * (nDYin - 0.5 * nemu * K);
  float err = sqrt(pow(sigmaR,2) * pow(nDYin - 0.5 * nemu * K , 2) + 
                   1./4. * pow(nemu,2) * R*R * sigmaK*sigmaK + 
                   nDYin * R*R +
                   nemu * 1./4. * K*K * R*R); 
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

void countEvents(int mass) {

  // taking the EE or ME trees for the wanted mass
  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EE");  
  // sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_ME");
  TChain *theChain = new TChain(nametree);
  
  char file_mc[1000];
  sprintf(file_mc,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH%d/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*Counters.root",mass);  
  theChain->Add(file_mc);
  //  cout << "reading tree " << nametree << " from file " << file_mc << endl;    
  
  // number of events at the wanted step of the selection
  nEv_endWW = 0.0;
  nEv_end0j = 0.0;
  nEv_end1j = 0.0;

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
    nEv_endWW += nSel[16];
    nEv_end0j += nSel[22];
    nEv_end1j += nSel[23];
  }
}

