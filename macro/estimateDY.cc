#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>

enum { ee=0, mm=1, em=2, me=3 };

// numbers filled from counters
float nEv_endWW[4];
float nEv_end0j[4];
float nEv_end1j[4];

float quadrSum(float x1, float x2, float x3=0, float x4=0, float x5=0, float x6=0, float x7=0, float x8=0);
std::pair<float,float> nDYout(float nDYin, float nemu, float R, float sigmaR, float K, float sigmaK);
float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);
void countEvents(int mass, const char* channel);

void estimateDY() {

  // constants
  float eff_softmu_Z = 0.867;
  //  float eff_softmu_Z = 1;
  float Rmm = 0.187398;
  float Rmm_err = 0.00355062;
  float kmm = 0.592288;
  float kmm_err = 0.0139356;
  float Ree = 0.152089;
  float Ree_err = 0.00338336;
  float kee = 0.422092;
  float kee_err = 0.00874687;


  TFile *fileData  = TFile::Open("results_data/merged/dataset_ll.root");
  TFile *fileZjets = TFile::Open("results/datasets_trees/Zjets_ll.root");

  TTree *treeData  = (TTree*)fileData->Get("T1");
  TTree *treeZjets = (TTree*)fileZjets->Get("T1");

  TH1F *ZeejetsH = new TH1F("ZeejetsH","",50,0,180);
  TH1F *ZmmjetsH = new TH1F("ZmmjetsH","",50,0,180);
  TH1F *ZemjetsH = new TH1F("ZemjetsH","",50,0,180);
  TH1F *ZmejetsH = new TH1F("ZmejetsH","",50,0,180);

  TH1F *neeInH = new TH1F("neeInH","",50,0,180);
  TH1F *nmmInH = new TH1F("nmmInH","",50,0,180);
  TH1F *nemInH = new TH1F("nemInH","",50,0,180);

  treeZjets->Project("ZeejetsH","deltaPhi","(WWSel && finalstate==0)*weight*puweight");
  treeZjets->Project("ZmmjetsH","deltaPhi","(WWSel && finalstate==1)*weight*puweight");
  treeZjets->Project("ZemjetsH","deltaPhi","(WWSel && (finalstate==2))*weight*puweight");
  treeZjets->Project("ZmejetsH","deltaPhi","(WWSel && (finalstate==3))*weight*puweight");
  
  treeData->Project("neeInH","deltaPhi","eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==0"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nmmInH","deltaPhi","eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && finalstate==1"); // missing softmu... not available in red trees... hopefully small contrib here
  treeData->Project("nemInH","deltaPhi","eleInvMass>12 && finalLeptons && pfMet>30 && projMet>35 && njets==0 && bTagTrackCount<2.1 && abs(eleInvMass-91.1876)<15 && (finalstate==2 || finalstate==3)"); // missing softmu... not available in red trees... hopefully small contrib here

  // DY estimation /////
  std::cout << "DY ESTIMATION..." << std::endl;
  float nZeejetsMC = ZeejetsH->Integral();
  float nZeejetsMC_err = yieldErrPoisson(nZeejetsMC,ZeejetsH->GetEntries());
  float nZmmjetsMC = ZmmjetsH->Integral();
  float nZmmjetsMC_err = yieldErrPoisson(nZmmjetsMC,ZmmjetsH->GetEntries());
  float nZemjetsMC = ZemjetsH->Integral();
  float nZemjetsMC_err = yieldErrPoisson(nZemjetsMC,ZemjetsH->GetEntries());
  float nZmejetsMC = ZmejetsH->Integral();
  float nZmejetsMC_err = yieldErrPoisson(nZmejetsMC,ZmejetsH->GetEntries());
  
  float nDYMC[4], nDYMC_err[4];
  nDYMC[ee] = ZeejetsH->Integral();
  nDYMC_err[ee] = yieldErrPoisson(nDYMC[ee],ZeejetsH->GetEntries());
  nDYMC[mm] = ZmmjetsH->Integral();
  nDYMC_err[mm] = yieldErrPoisson(nDYMC[mm],ZmmjetsH->GetEntries());
  nDYMC[em] = ZemjetsH->Integral();
  nDYMC_err[em] = yieldErrPoisson(nDYMC[em],ZemjetsH->GetEntries());
  nDYMC[me] = ZmejetsH->Integral();
  nDYMC_err[me] = yieldErrPoisson(nDYMC[me],ZmejetsH->GetEntries());

  float neeIn = neeInH->Integral() * eff_softmu_Z;
  float nmmIn = nmmInH->Integral() * eff_softmu_Z;
  float nemIn = (nemInH->Integral()) * eff_softmu_Z;

  // ee and mm from data
  float neeExp = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).first;
  float neeExp_err = (nDYout(neeIn, nemIn, Ree, Ree_err, kee, kee_err)).second;
  float nmmExp = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).first;
  float nmmExp_err = (nDYout(nmmIn, nemIn, Rmm, Rmm_err, kmm, kmm_err)).second;

  // and em from MC
  float nemExp = nZemjetsMC;
  float nemExp_err = yieldErrPoisson(nZemjetsMC,ZemjetsH->GetEntries());
  float nmeExp = nZmejetsMC;
  float nmeExp_err = yieldErrPoisson(nZmejetsMC,ZmejetsH->GetEntries());

  float nemmeExp = nZemjetsMC+nZmejetsMC;
  float nemmeExp_err = quadrSum(nZemjetsMC_err,nZmejetsMC_err);

  std::cout << "Number of DY->ll events at W+W- level: " << std::endl;
  std::cout << "neeIn = "  << ( 150./126.) * neeIn << "\tnmmIn = "  << (150./126.) * nmmIn << "\tnemIn = " << (150./126.) * nemIn << std::endl;
  std::cout << "nEE MC = " << ( 150./126.) * nZeejetsMC  << " +/- " << nZeejetsMC_err
            << "\tData = " << ( 150./126.) * neeExp      << " +/- " << neeExp_err << std::endl;
  std::cout << "nMM MC = " << ( 150./126.) * nZmmjetsMC  << " +/- " << nZmmjetsMC_err
            << "\tData = " << ( 150./126.) * nmmExp      << " +/- " << nmmExp_err << std::endl;
  std::cout << "nEM+nME MC = " << (150./126.) * nemmeExp << " +/- " << nemmeExp_err << std::endl; 

  std::cout << "END DY ESTIMATION." << std::endl;
  ////////// END DY ///////////

  float nDYData[4];
  float nDYData_err[4];
  nDYData[ee] = neeExp;
  nDYData_err[ee] = neeExp_err;
  nDYData[mm] = nmmExp;
  nDYData_err[mm] = nmmExp_err;

  // these are MC
  nDYData[em] = nemExp;
  nDYData_err[em] = nemExp_err;
  nDYData[me] = nmeExp;
  nDYData_err[me] = nmeExp_err;

  
  ofstream tablefile1;
  tablefile1.open("DYYieldsData_ForTable_0j.txt", ios_base::app);
  tablefile1.precision(2);

  ofstream tablefile3;
  tablefile3.open("DYYieldsMC_ForTable_0j.txt", ios_base::app);
  tablefile3.precision(2);

  int masses[17] = {120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600};
  // -------------------------------------------------------------------
  // now considering all masses to estimate the number of events at the end of the HWW selection
  for (int i=0; i<17; i++) {
    
    int mass = masses[i];

    countEvents(mass,"EE");
    countEvents(mass,"MM");
    countEvents(mass,"EM");
    countEvents(mass,"ME");

    float nDYData_HiggsSel[4];
    float nDYData_HiggsSel_err[4];

    float nDYMC_HiggsSel[4];
    float nDYMC_HiggsSel_err[4];

    for(int icha=0;icha<4;icha++) {
      float eff_0j = (nEv_endWW[icha]==0) ? 0. : nEv_end0j[icha]/nEv_endWW[icha];
      float eff_0j_err = (nEv_endWW[icha]==0) ? 0. : sqrt(eff_0j*(1-eff_0j)/nEv_endWW[icha]);
      float effErrRel = (eff_0j==0) ? 0. : eff_0j_err/eff_0j;
    
      nDYData_HiggsSel[icha] = nDYData[icha] * eff_0j;
      float DYDataErrRel = (nDYData[icha]==0) ? 0. : nDYData_err[icha]/nDYData[icha];
      nDYData_HiggsSel_err[icha] = nDYData_HiggsSel[icha] * quadrSum(DYDataErrRel,effErrRel);

      nDYMC_HiggsSel[icha] = nDYMC[icha] * eff_0j;
      float DYMCErrRel = (nDYMC[icha]==0) ? 0. : nDYMC_err[icha]/nDYMC[icha];
      nDYMC_HiggsSel_err[icha] = nDYMC_HiggsSel[icha] * quadrSum(DYMCErrRel,effErrRel);
    }
  
    cout.precision(2);
    cout << "Zero jets bin" << endl;
    cout << "Higgs Mass = " << mass << endl;
    cout << "\tee: MC = "   << 150./126. * nDYMC_HiggsSel[ee]   << " +/- " << nDYMC_HiggsSel_err[ee]
	 << "\tee: data = " << 150./126. * nDYData_HiggsSel[ee] << " +/- " << nDYData_HiggsSel_err[ee] << std::endl;
    cout << "\tmm: MC = "   << 150./126. * nDYMC_HiggsSel[mm]   << " +/- " << nDYMC_HiggsSel_err[mm]
	 << "\tmm: data = " << 150./126. * nDYData_HiggsSel[mm] << " +/- " << nDYData_HiggsSel_err[mm] << std::endl;
    cout << "\tem: MC = "   << 150./126. * nDYMC_HiggsSel[em]   << " +/- " << nDYMC_HiggsSel_err[em]
	 << "\tem: data = " << 150./126. * nDYData_HiggsSel[em] << " +/- " << nDYData_HiggsSel_err[em] << std::endl;
    cout << "\tme: MC = "   << 150./126. * nDYMC_HiggsSel[me]   << " +/- " << nDYMC_HiggsSel_err[me]
	 << "\tme: data = " << 150./126. * nDYData_HiggsSel[me] << " +/- " << nDYData_HiggsSel_err[me] << std::endl;
    cout << endl;


    // summary table for limits                                                                                                             
    if (i==0) {
      tablefile1 << "# zero jets bin data" << endl;
      tablefile1 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile1 << mass
               << " " << "\t" << (150./126.) * nDYData_HiggsSel[1] << " +/- " <<  nDYData_HiggsSel_err[1]
               << " " << "\t" << (150./126.) * nDYData_HiggsSel[3] << " +/- " <<  nDYData_HiggsSel_err[3]
               << " " << "\t" << (150./126.) * nDYData_HiggsSel[2] << " +/- " <<  nDYData_HiggsSel_err[2]
               << " " << "\t" << (150./126.) * nDYData_HiggsSel[0] << " +/- " <<  nDYData_HiggsSel_err[0]
               << std::endl;

    if (i==0) {
      tablefile3 << "# zero jets bin MC" << endl;
      tablefile3 << "# \t mumu \t mue \t emu \t ee" << endl;
    }
    tablefile3 << mass
               << " " << "\t" << (150./126.) * nDYMC_HiggsSel[1] << " +/- " <<  nDYMC_HiggsSel_err[1]
               << " " << "\t" << (150./126.) * nDYMC_HiggsSel[3] << " +/- " <<  nDYMC_HiggsSel_err[3]
               << " " << "\t" << (150./126.) * nDYMC_HiggsSel[2] << " +/- " <<  nDYMC_HiggsSel_err[2]
               << " " << "\t" << (150./126.) * nDYMC_HiggsSel[0] << " +/- " <<  nDYMC_HiggsSel_err[0]
               << std::endl;

    float nDYData_HiggsSel_Tot = nDYData_HiggsSel[ee] + nDYData_HiggsSel[mm] + nDYData_HiggsSel[em] + nDYData_HiggsSel[me];
    float nDYData_HiggsSel_Tot_err = quadrSum(nDYData_HiggsSel_err[ee],nDYData_HiggsSel_err[mm],nDYData_HiggsSel_err[em],nDYData_HiggsSel_err[me]);

    float nDYMC_HiggsSel_Tot = nDYMC_HiggsSel[ee] + nDYMC_HiggsSel[mm] + nDYMC_HiggsSel[em] + nDYMC_HiggsSel[me];
    float nDYMC_HiggsSel_Tot_err = quadrSum(nDYMC_HiggsSel_err[ee],nDYMC_HiggsSel_err[mm],nDYMC_HiggsSel_err[em],nDYMC_HiggsSel_err[me]);
    
    cout.precision(2);
    cout << "Higgs Mass = " << mass 
         << "\tMC = "   << 150./126. * nDYMC_HiggsSel_Tot   << " +/- " << nDYMC_HiggsSel_Tot_err
         << "\tdata = " << 150./126. * nDYData_HiggsSel_Tot << " +/- " << nDYData_HiggsSel_Tot_err << std::endl;
  }

}

float quadrSum(float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8) {
  return sqrt(x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8);
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

void countEvents(int mass, const char *channel) {

  // taking the EE or ME trees for the wanted mass
  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_%s",channel);  
  TChain *theChain = new TChain(nametree);
  
  char sampleName1[100];
  char sampleName2[100];
  // some approximation in the following...
  if(TString(channel).Contains("EE")) sprintf(sampleName1,"DYToEE_M-20_TuneZ2_7TeV-pythia6");
  if(TString(channel).Contains("MM")) sprintf(sampleName1,"DYToMuMu_M-20_TuneZ2_7TeV-pythia6");
  if(TString(channel).Contains("EM") || TString(channel).Contains("ME")) sprintf(sampleName1,"DYToTauTau_M-20_TuneZ2_7TeV-pythia6");

  if(TString(channel).Contains("EE")) sprintf(sampleName2,"DYToEE_M-10To20_TuneZ2_7TeV-pythia6");
  if(TString(channel).Contains("MM")) sprintf(sampleName2,"DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6");
  if(TString(channel).Contains("EM") || TString(channel).Contains("ME")) sprintf(sampleName2,"DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6");

  char file_mc1[1000];
  sprintf(file_mc1,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH%d/Spring11_V2/%s/*Counters.root",mass,sampleName1);  

  char file_mc2[1000];
  sprintf(file_mc2,"/cmsrm/pc24_2/emanuele/data/Higgs4.1.X/MC2011_LHLoose_V13/OptimMH%d/Spring11_V2/%s/*Counters.root",mass,sampleName2);  

  theChain->Add(file_mc1);
  theChain->Add(file_mc2);
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

