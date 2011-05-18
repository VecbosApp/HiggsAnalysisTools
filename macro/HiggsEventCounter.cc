#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#define NSAMPLES 17

using namespace std;

double weight(double ngen, double xsec, double filtereff, double lumi = 1);

void countEvents(int mH=160) {

  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EE");

  cout << "nametree = " << nametree << endl;

  TChain *chains[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    chains[isample] = new TChain(nametree);
  }

  chains[0]->Add("results/Spring11_V2/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root");

  chains[1]->Add("results/Spring11_V2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root");
  chains[2]->Add("results/Spring11_V2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root");
  chains[3]->Add("results/Spring11_V2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6/*Counters.root");

  chains[4]->Add("results/Spring11_V2/DYToEE_M-20_TuneZ2_7TeV-pythia6/*Counters.root");
  chains[5]->Add("results/Spring11_V2/DYToMuMu_M-20_TuneZ2_7TeV-pythia6/*Counters.root");
  chains[6]->Add("results/Spring11_V2/DYToTauTau_M-20_TuneZ2_7TeV-pythia6/*Counters.root");

  chains[7]->Add("results/Spring11_V2/ToBLNu_TuneZ2_s-channel_7TeV-madgraph/*Counters.root");
  chains[8]->Add("results/Spring11_V2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph/*Counters.root");
  chains[9]->Add("results/Spring11_V2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/*Counters.root");

  chains[10]->Add("results/Spring11_V2/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root");

  TString hSample("results/Spring11_V2/GluGluToHToWWTo2L2Nu_M-");
  char mass[5];
  sprintf(mass,"%d",mH);
  hSample += TString(mass)+TString("_7TeV-powheg-pythia6/*Counters.root");

  chains[11]->Add(hSample.Data());

  chains[12]->Add("results/Spring11_V2/PhotonVJets_7TeV-madgraph/*Counters.root");
  //  chains[13]->Add("results/Spring11_V2/WWTo2L2Nu_TuneZ2_7TeV-pythia6/*Counters.root"); // PYTHIA sample
  chains[13]->Add("results/Spring11_V2/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola-WWFilter/*Counters.root"); // MADGRAPH sample
  chains[14]->Add("results/Spring11_V2/GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6/*Counters.root");
  chains[15]->Add("results/Spring11_V2/WZTo3LNu_TuneZ2_7TeV-pythia6/*Counters.root");
  chains[16]->Add("results/Spring11_V2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/*Counters.root");

  cout << "chains added. " << endl;

  std::vector<TString> sampleName;

  sampleName.push_back("results/merged/Wjets_ee.root");

  sampleName.push_back("results/merged/Zee_Lo_ee.root");
  sampleName.push_back("results/merged/Zmm_Lo_ee.root");
  sampleName.push_back("results/merged/Ztautau_Lo_ee.root");

  sampleName.push_back("results/merged/Zee_Hi_ee.root");
  sampleName.push_back("results/merged/Zmm_Hi_ee.root");
  sampleName.push_back("results/merged/Ztautau_Hi_ee.root");

  sampleName.push_back("results/merged/SingleTop_sChannel_ee.root");
  sampleName.push_back("esults/merged/SingleTop_tChannel_ee.root");
  sampleName.push_back("results/merged/SingleTop_tWChannel_ee.root");

  sampleName.push_back("results/merged/TTbar_ee.root");

  sampleName.push_back(TString("results/merged/H")+TString(mass)+TString("_ee.root"));

  sampleName.push_back("results/merged/Wgamma_ee.root");
  sampleName.push_back("results/merged/WW_ee.root");
  sampleName.push_back("results/merged/ggWW_ee.root");
  sampleName.push_back("results/merged/WZ_ee.root");
  sampleName.push_back("results/merged/ZZ_ee.root");

  std::map<int,float> Higgs_xsec_masses;
  // samples are emu only
  Higgs_xsec_masses.insert(std::make_pair(120,0.249642*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(130,0.452090*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(140,0.641773*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(150,0.770471*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(160,0.866443*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(170,0.782962*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(180,0.659328*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(190,0.486486*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(200,0.408305*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(210,0.358465*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(220,0.321398*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(230,0.290454*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(250,0.243724*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(300,0.175652*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(350,0.160052*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(400,0.124330*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(450,0.078433*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(500,0.048702*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(550,0.030364*4./9.));
  Higgs_xsec_masses.insert(std::make_pair(600,0.019184*4./9.));

  std::vector<float> sampleXsec;
  sampleXsec.push_back(31314.);
  sampleXsec.push_back(3457./3.);
  sampleXsec.push_back(3457./3.);
  sampleXsec.push_back(3457./3.);
  sampleXsec.push_back(4998./3.);
  sampleXsec.push_back(4998./3.);
  sampleXsec.push_back(4998./3.);
  sampleXsec.push_back(4.21 * (0.1080*3));
  sampleXsec.push_back(64.6 * (0.1080*3));
  sampleXsec.push_back(10.6);
  sampleXsec.push_back(157.5);
  sampleXsec.push_back(Higgs_xsec_masses[mH]);
  sampleXsec.push_back(165.);
  sampleXsec.push_back(4.50347);
  sampleXsec.push_back(0.1538);
  sampleXsec.push_back(0.599442);
  sampleXsec.push_back(7.67); // 5.9*1.3 is to consider the ratio (m_LL>12/m_LL>40) https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopDileptonRefAnalysis2010Pass6

  std::cout << "For mH = " << mH << " GeV found xsec x BR = " << Higgs_xsec_masses[mH] << " pb " << std::endl;

  float nEv[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<NSAMPLES; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;

    Int_t           nCuts;
    Float_t         nSel[25];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    // loop over files (>1 if VecBos in batch, splitted in many jobs)
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;

      nEv[isample] += nSel[0];
    }
  }

  if (sampleXsec.size() != sampleName.size() ) cout << "nasty error! check sizes..." << endl;

  std::ofstream weightsFile;
  weightsFile.open("weightTrees.sh");
  weightsFile << "#! /bin/sh\n\n" << std::endl;
  weightsFile << "lumiEE=$1" << std::endl;
  weightsFile << "lumiMM=$2" << std::endl;
  weightsFile << "lumiEM=$3" << std::endl;
  
  weightsFile << "echo \"===> THIS WEIGHTING IS FOR MH = " << mH << " <====\" "<< std::endl;

  // now write ee
  weightsFile << "echo \"Adding weights for ee datasets for \" $lumiEE \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    weightsFile << "addWeights(\"" << sampleName[isample].Data() << "\", " << w << "*$lumiEE, 0);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;
  
  // now write mm
  weightsFile << "echo \"Adding weights for mm datasets for \" $lumiMM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int isample=0; isample<NSAMPLES; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMM = sampleName[isample].ReplaceAll("_ee","_mm");
    weightsFile << "addWeights(\"" << sampleNameMM.Data() << "\", " << w << "*$lumiMM, 1);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;


  // now write em
  weightsFile << "echo \"Adding weights for em datasets for \" $lumiEM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEM = sampleName[isample].ReplaceAll("_mm","_em");
    weightsFile << "addWeights(\"" << sampleNameEM.Data() << "\", " << w << "*$lumiEM, 2);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;


  // now write me
  weightsFile << "echo \"Adding weights for me datasets for \" $lumiEM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEM = sampleName[isample].ReplaceAll("_em","_me");
    weightsFile << "addWeights(\"" << sampleNameEM.Data() << "\", " << w << "*$lumiEM, 3);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;

  weightsFile << "echo \"done weighting.\"" << std::endl;

}

double weight(double ngen, double xsec, double filtereff, double lumi) {

  if(ngen==0) return 0;
  return xsec * filtereff * lumi / ngen;

}

