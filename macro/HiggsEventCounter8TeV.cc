#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#define NSAMPLES 11

using namespace std;

double weight(double ngen, double xsec, double filtereff, double lumi = 1);

void countEvents() {

  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EE");

  cout << "nametree = " << nametree << endl;

  // signals
  std::vector<std::vector<TChain*> > signalChains;
  int mH[27] = {110,115,120,125,130,135,140,145,150,155,160,170,180,190,200,250,300,350,400,450,500,550,600,700,800,900,1000};
  for(int imass=0; imass<27;imass++) {
    char mass[5];
    sprintf(mass,"%d",mH[imass]);

    std::vector<TChain*> massChain;

    TChain *sigChains[2];
    for(int i=0; i<2; i++) sigChains[i] = new TChain(nametree);

    TString hSample("results/Summer12_V14_52X/GluGluToHToWWTo2LAndTau2Nu_M-");
    hSample += TString(mass)+TString("_8TeV-powheg-pythia6/*Counters.root");
    sigChains[0]->Add(hSample.Data());
    
    hSample = TString("results/Summer12_V14_52X/VBF_HToWWTo2LAndTau2Nu_M-");
    hSample += TString(mass)+TString("_8TeV-powheg-pythia6/*Counters.root");
    sigChains[1]->Add(hSample.Data());
    
    for(int i=0; i<2; i++) massChain.push_back(sigChains[i]);
    signalChains.push_back(massChain);
  }
  // backgrounds
  TChain *chains[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    chains[isample] = new TChain(nametree);
  }

  // nominal sample first, then the systematics ones
  chains[0]->Add("results/Summer12_V14_52X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/*Counters.root");
  chains[1]->Add("results/Summer12_V14_52X/DYJetsToLL_M-10To50filter_8TeV-madgraph/*Counters.root");
  chains[2]->Add("results/Summer12_V14_52X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/*Counters.root");

  chains[3]->Add("results/Summer12_V14_52X/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/*Counters.root");
  chains[4]->Add("results/Summer12_V14_52X/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/*Counters.root");
  chains[5]->Add("results/Summer12_V14_52X/TTJets_TuneZ2star_8TeV-madgraph-tauola/*Counters.root"); 

  chains[6]->Add("results/Summer12_V14_52X/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/*Counters.root"); // nominal WW MADGRAPH sample
  chains[7]->Add("results/Summer12_V14_52X/GluGluToWWTo4L_TuneZ2star_8TeV-gg2ww-pythia6/*Counters.root");
  chains[8]->Add("results/Summer12_V14_52X/WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola/*Counters.root");
  chains[9]->Add("results/Summer12_V14_52X/ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola/*Counters.root");

  // samples for systematics
  chains[10]->Add("results/Summer12_V14_52X/WWTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola/*Counters.root"); // PYTHIA sample 
  cout << "chains added. " << endl;

  std::vector<std::vector<TString> > signalSampleName;
  for(int imass=0; imass<27;imass++) {
    char mass[5];
    sprintf(mass,"%d",mH[imass]);
    
    std::vector<TString> massSampleName;
    massSampleName.push_back(TString("results/merged/ggH")+TString(mass)+TString("2LAndTau2Nu_ee.root"));
    massSampleName.push_back(TString("results/merged/qqH")+TString(mass)+TString("2LAndTau2Nu_ee.root"));

    signalSampleName.push_back(massSampleName);
  }

  std::vector<TString> sampleName;
  sampleName.push_back("results/merged/Wjets_ee.root"); // 0
  sampleName.push_back("results/merged/Zll_Lo_ee.root"); // 1
  sampleName.push_back("results/merged/Zll_Hi_ee.root"); // 2
  sampleName.push_back("results/merged/SingleT_tWChannel_ee.root"); // 3
  sampleName.push_back("results/merged/SingleTbar_tWChannel_ee.root"); // 4
  sampleName.push_back("results/merged/TTbar_ee.root"); // 5
  sampleName.push_back("results/merged/WW_ee.root"); // 6
  sampleName.push_back("results/merged/ggWW_ee.root"); // 7
  sampleName.push_back("results/merged/WZ_ee.root"); // 8
  sampleName.push_back("results/merged/ZZ_ee.root"); // 9
  sampleName.push_back("results/merged/WW_pythia_ee.root"); // 10

  std::map<int,float> ggHiggs_xsec;
  // sigma_effective for gg->H and qqH = sigma * BR(H->WW) * BR(W->lnu)^2 (pb)
  ggHiggs_xsec.insert(std::make_pair(110,0.100387));
  ggHiggs_xsec.insert(std::make_pair(115,0.165009));
  ggHiggs_xsec.insert(std::make_pair(120,0.249642));
  ggHiggs_xsec.insert(std::make_pair(125,0.347151));
  ggHiggs_xsec.insert(std::make_pair(130,0.452090));
  ggHiggs_xsec.insert(std::make_pair(135,0.553354));
  ggHiggs_xsec.insert(std::make_pair(140,0.641773));
  ggHiggs_xsec.insert(std::make_pair(145,0.713397));
  ggHiggs_xsec.insert(std::make_pair(150,0.770471));
  ggHiggs_xsec.insert(std::make_pair(155,0.818479));
  ggHiggs_xsec.insert(std::make_pair(160,0.866443));
  ggHiggs_xsec.insert(std::make_pair(170,0.782962));
  ggHiggs_xsec.insert(std::make_pair(180,0.659328));
  ggHiggs_xsec.insert(std::make_pair(190,0.486486));
  ggHiggs_xsec.insert(std::make_pair(200,0.408305));
  ggHiggs_xsec.insert(std::make_pair(250,0.243724));
  ggHiggs_xsec.insert(std::make_pair(300,0.175652));
  ggHiggs_xsec.insert(std::make_pair(350,0.160052));
  ggHiggs_xsec.insert(std::make_pair(400,0.124330));
  ggHiggs_xsec.insert(std::make_pair(450,0.078433));
  ggHiggs_xsec.insert(std::make_pair(500,0.048702));
  ggHiggs_xsec.insert(std::make_pair(550,0.030364));
  ggHiggs_xsec.insert(std::make_pair(600,0.019184));
  ggHiggs_xsec.insert(std::make_pair(700,0.007995));
  ggHiggs_xsec.insert(std::make_pair(800,0.003532));
  ggHiggs_xsec.insert(std::make_pair(900,0.001637));
  ggHiggs_xsec.insert(std::make_pair(1000,0.000789));

  std::map<int,float> qqHiggs_xsec;
  // sigma_effective for gg->H and qqH = sigma * BR(H->WW) * BR(W->lnu)^2 (pb)
  qqHiggs_xsec.insert(std::make_pair(110,0.007073));
  qqHiggs_xsec.insert(std::make_pair(115,0.012126));
  qqHiggs_xsec.insert(std::make_pair(120,0.019057));
  qqHiggs_xsec.insert(std::make_pair(125,0.027450));
  qqHiggs_xsec.insert(std::make_pair(130,0.036942));
  qqHiggs_xsec.insert(std::make_pair(135,0.046536));
  qqHiggs_xsec.insert(std::make_pair(140,0.055686));
  qqHiggs_xsec.insert(std::make_pair(145,0.713397));
  qqHiggs_xsec.insert(std::make_pair(150,0.070561));
  qqHiggs_xsec.insert(std::make_pair(155,0.818479));
  qqHiggs_xsec.insert(std::make_pair(160,0.083858));
  qqHiggs_xsec.insert(std::make_pair(170,0.082794));
  qqHiggs_xsec.insert(std::make_pair(180,0.073183));
  qqHiggs_xsec.insert(std::make_pair(190,0.057139));
  qqHiggs_xsec.insert(std::make_pair(200,0.049566));
  qqHiggs_xsec.insert(std::make_pair(250,0.031672));
  qqHiggs_xsec.insert(std::make_pair(300,0.021873));
  qqHiggs_xsec.insert(std::make_pair(350,0.015113));
  qqHiggs_xsec.insert(std::make_pair(400,0.009885));
  qqHiggs_xsec.insert(std::make_pair(450,0.007143));
  qqHiggs_xsec.insert(std::make_pair(500,0.005439));
  qqHiggs_xsec.insert(std::make_pair(550,0.004249));
  qqHiggs_xsec.insert(std::make_pair(600,0.003374));
  qqHiggs_xsec.insert(std::make_pair(700,0.002205));
  qqHiggs_xsec.insert(std::make_pair(800,0.001478));
  qqHiggs_xsec.insert(std::make_pair(900,0.001016));
  qqHiggs_xsec.insert(std::make_pair(1000,0.000724));

  std::map<int,float> Higgs_8TevOver7TeV_ratio;
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(110,1.27394)); // using the same as 115 GeV
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(115,1.27394));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(120,1.27742));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(125,1.28090));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(130,1.28437));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(135,1.28775));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(140,1.29112));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(145,1.29439));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(150,1.29766));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(155,1.29766)); // using the same as 150 GeV
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(160,1.30403));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(170,1.31031));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(180,1.31641));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(190,1.32237));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(200,1.32822));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(250,1.35530));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(300,1.38238));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(350,1.40661));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(400,1.43085));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(450,1.45471));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(500,1.47857));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(550,1.50211));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(600,1.52566));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(700,1.52566)); // for mH>600 GeV, using the 600 GeV value
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(800,1.52566));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(900,1.52566));
  Higgs_8TevOver7TeV_ratio.insert(std::make_pair(10000,1.52566));

  // these for the 2LAndTauNu samples
  std::vector<std::vector<double> > signalXSec;
  for(int imass=0; imass<27;imass++) {
    std::vector<double> massXsec;
    double ggHiggs_xsec_2landtau2nu = ggHiggs_xsec[mH[imass]];
    double qqHiggs_xsec_2landtau2nu = qqHiggs_xsec[mH[imass]];
    double scaleto8TeV = Higgs_8TevOver7TeV_ratio[mH[imass]];

    massXsec.push_back(ggHiggs_xsec_2landtau2nu*scaleto8TeV);
    massXsec.push_back(qqHiggs_xsec_2landtau2nu*scaleto8TeV);

    signalXSec.push_back(massXsec);
  }

  std::vector<float> sampleXsec;
  sampleXsec.push_back(37509.); // 0
  sampleXsec.push_back(860.5013); // 1
  sampleXsec.push_back(3532.8149); // 2
  sampleXsec.push_back(11.1773); // 3
  sampleXsec.push_back(11.1773); // 4
  sampleXsec.push_back(225.1967); // 5
  sampleXsec.push_back(5.8123); // 6
  sampleXsec.push_back(0.182852); // 7
  sampleXsec.push_back(0.7346); // 8
  sampleXsec.push_back(0.3649); // 9: PYTHIA
  sampleXsec.push_back(5.8123); // 12

  std::vector<std::vector<double> > signalProcId;
  for(int imass=0; imass<27;imass++) {
    std::vector<double> massId;
    massId.push_back(9000+mH[imass]); // ggH->2LAndTau2Nu
    massId.push_back(8000+mH[imass]); // qqH->2LAndTau2Nu
    signalProcId.push_back(massId);
  }

  std::vector<int> sampleProcessId; 
  // ids are taken from: https://docs.google.com/spreadsheet/ccc?key=0Ankm0DuoD0h0dHRhV1VNSlV1NEdhNFdyOXh3eFpSMHc&hl=en_US#gid=31
  sampleProcessId.push_back(80); // Wjets
  sampleProcessId.push_back(36); // Z->ll [10-50]
  sampleProcessId.push_back(37); // Z->ll [50-inf]
  sampleProcessId.push_back(19); // tW
  sampleProcessId.push_back(20); // tbarW
  sampleProcessId.push_back(10); // ttbar
  sampleProcessId.push_back(0); // qqWW Mad
  sampleProcessId.push_back(1); // ggWW
  sampleProcessId.push_back(74); // WZ->3LNu
  sampleProcessId.push_back(71); // ZZ->anything
  sampleProcessId.push_back(6); // qqWW PYTHIA
  sampleProcessId.push_back(81); // Vgamma madgraph inclusive

  // signal samples
  float nEvH[27][2];
  for(int imass=0; imass<27; imass++) {
    for(int i=0; i<2; i++) {
      nEvH[imass][i] = 0.0;
    }
  }

  for(int imass=0; imass<27; imass++) {
    std::vector<TChain*> massChain = signalChains[imass];
    for(int i=0; i<2; i++) {

      cout << "\tProcessing signal sample mass # " << mH[imass] << "..." << endl;
      
      Int_t           nCuts;
      Float_t         nSel[29];   //[nCuts]
      
      // List of branches
      TBranch        *b_nCuts;   //!
      TBranch        *b_nSel;   //!
      
      massChain[i]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
      massChain[i]->SetBranchAddress("nSel", nSel, &b_nSel);
    
      Long64_t nentries = massChain[i]->GetEntries();

      Long64_t nbytes = 0, nb = 0;
      // loop over files (>1 if VecBos in batch, splitted in many jobs)
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        nb = massChain[i]->GetEntry(jentry);   nbytes += nb;

        nEvH[imass][i] += nSel[0];
      }
    }
  }

  cout << "debug: start bkg" << endl;

  // backgrounds
  float nEv[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<NSAMPLES; isample++) {

    cout << "processing sample " << isample << endl;
    Int_t           nCuts;
    Float_t         nSel[29];   //[nCuts]
    
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

  cout << "Processed all the chains." << endl; 

  if (sampleXsec.size() != sampleName.size() ) cout << "nasty error! check sizes..." << endl;

  cout << signalSampleName.size() << " " << signalProcId.size() << " " << signalXSec.size() << Higgs_8TevOver7TeV_ratio.size() << endl; 

  std::ofstream weightsFile;
  weightsFile.open("weightTrees.sh");
  weightsFile << "#! /bin/sh\n\n" << std::endl;
  weightsFile << "mkdir -p results/merged_skim" << std::endl;
  weightsFile << "lumiEE=$1" << std::endl;
  weightsFile << "lumiMM=$2" << std::endl;
  weightsFile << "lumiEM=$3" << std::endl;
  
  // now write ee
  weightsFile << "echo \"Adding weights for ee datasets for \" $lumiEE \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<27; imass++) {
    std::vector<double> massXsec = signalXSec[imass];
    std::vector<TString> massSampleName = signalSampleName[imass];
    std::vector<double> massId = signalProcId[imass];
    for(int i=0; i<2; i++) {
      int release = 0;
      float w = weight(nEvH[imass][i], massXsec[i], 1., 1.);
      weightsFile << "addWeights(\"" << massSampleName[i].Data() << "\", " << w << "*$lumiEE, " << massId[i] << " ,1, " << release << ");" << std::endl;
    }
  }
  for(int isample=0; isample<NSAMPLES; isample++) {
    int release = 0;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    weightsFile << "addWeights(\"" << sampleName[isample].Data() << "\", " << w << "*$lumiEE, " << sampleProcessId[isample] << " ,1, " << release << ");" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;
  
  // now write mm
  weightsFile << "echo \"Adding weights for mm datasets for \" $lumiMM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<27; imass++) {
    std::vector<double> massXsec = signalXSec[imass];
    std::vector<TString> massSampleName = signalSampleName[imass];
    std::vector<double> massId = signalProcId[imass];
    for(int i=0; i<2; i++) {
      cout << "Events processed for sample: " << massSampleName[i] << " = " << nEvH[imass][i] << endl;
      int release = 0;
      float w = weight(nEvH[imass][i], massXsec[i], 1., 1.);
      TString massSampleNameMM = massSampleName[i].ReplaceAll("_ee","_mm");
      weightsFile << "addWeights(\"" << massSampleNameMM.Data() << "\", " << w << "*$lumiMM, " << massId[i] << " ,0, " << release << ");" << std::endl;
    }
  }
  for(int isample=0; isample<NSAMPLES; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    int release = 0;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMM = sampleName[isample].ReplaceAll("_ee","_mm");
    weightsFile << "addWeights(\"" << sampleNameMM.Data() << "\", " << w << "*$lumiMM, " << sampleProcessId[isample] << " ,0, " << release << ");" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;


  // now write em
  weightsFile << "echo \"Adding weights for em datasets for \" $lumiEM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<27; imass++) {
    std::vector<double> massXsec = signalXSec[imass];
    std::vector<TString> massSampleName = signalSampleName[imass];
    std::vector<double> massId = signalProcId[imass];
    for(int i=0; i<2; i++) {
      int release = 0;
      float w = weight(nEvH[imass][i], massXsec[i], 1., 1.);
      TString massSampleNameEM = massSampleName[i].ReplaceAll("_ee","_em");
      weightsFile << "addWeights(\"" << massSampleNameEM.Data() << "\", " << w << "*$lumiEM, " << massId[i] << " ,2, " << release << ");" << std::endl;
    }
  }
  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    int release = 0;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEM = sampleName[isample].ReplaceAll("_mm","_em");
    weightsFile << "addWeights(\"" << sampleNameEM.Data() << "\", " << w << "*$lumiEM, " << sampleProcessId[isample] << " ,2, " << release << ");" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;


  // now write me
  weightsFile << "echo \"Adding weights for me datasets for \" $lumiEM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<27; imass++) {
    std::vector<double> massXsec = signalXSec[imass];
    std::vector<TString> massSampleName = signalSampleName[imass];
    std::vector<double> massId = signalProcId[imass];
    for(int i=0; i<2; i++) {
      int release = 0;
      float w = weight(nEvH[imass][i], massXsec[i], 1., 1.);
      TString massSampleNameME = massSampleName[i].ReplaceAll("_ee","_me");
      weightsFile << "addWeights(\"" << massSampleNameME.Data() << "\", " << w << "*$lumiEM, " << massId[i] << " ,3, " << release << ");" << std::endl;
    }
  }
  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    int release = 0;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEM = sampleName[isample].ReplaceAll("_em","_me");
    weightsFile << "addWeights(\"" << sampleNameEM.Data() << "\", " << w << "*$lumiEM, " << sampleProcessId[isample] << " ,3, " << release << ");" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;

  weightsFile << "echo \"done weighting.\"" << std::endl;

}

double weight(double ngen, double xsec, double filtereff, double lumi) {

  if(ngen==0) return 0;
  return xsec * filtereff * lumi / ngen;

}

