#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void countEvents() {

  char nametree[200];
  sprintf(nametree,"PRESELECTION_EVENT_COUNTER");

  cout << "nametree = " << nametree << endl;

  TChain *chains[38];
  for(int isample=0; isample<38; isample++) {
    chains[isample] = new TChain(nametree);
  }

  chains[0]->Add("results/WJetsMADGRAPH/WJets-madgraph/*Counters.root");

  chains[1]->Add("results/ZJetsMADGRAPH/ZJets-madgraph/*Counters.root");

  chains[2]->Add("results/TTbar/TTbarJets-madgraph/*Counters.root");

  chains[3]->Add("results/QCD/QCD_EMEnriched_Pt20to30/*Counters.root");
  chains[4]->Add("results/QCD/QCD_EMEnriched_Pt30to80/*Counters.root");
  chains[5]->Add("results/QCD/QCD_EMEnriched_Pt80to170/*Counters.root");

  chains[6]->Add("results/QCD/QCD_BCtoE_Pt20to30/*Counters.root");
  chains[7]->Add("results/QCD/QCD_BCtoE_Pt30to80/*Counters.root");
  chains[8]->Add("results/QCD/QCD_BCtoE_Pt80to170/*Counters.root");

  chains[9]->Add("results/SingleTop/SingleTop_sChannel-madgraph/*Counters.root");
  chains[10]->Add("results/SingleTop/SingleTop_tChannel-madgraph/*Counters.root");
  chains[11]->Add("results/SingleTop/SingleTop_tWChannel-madgraph/*Counters.root");

  chains[12]->Add("results/PhotonJet/PhotonJet_Pt0to15/*Counters.root");
  chains[13]->Add("results/PhotonJet/PhotonJet_Pt15to20/*Counters.root");
  chains[14]->Add("results/PhotonJet/PhotonJet_Pt20to30/*Counters.root");
  chains[15]->Add("results/PhotonJet/PhotonJet_Pt30to50/*Counters.root");
  chains[16]->Add("results/PhotonJet/PhotonJet_Pt50to80/*Counters.root");
  chains[17]->Add("results/PhotonJet/PhotonJet_Pt80to120/*Counters.root");
  chains[18]->Add("results/PhotonJet/PhotonJet_Pt120to170/*Counters.root");
  chains[19]->Add("results/PhotonJet/PhotonJet_Pt170to300/*Counters.root");
  chains[20]->Add("results/PhotonJet/PhotonJet_Pt300to500/*Counters.root");
  chains[21]->Add("results/PhotonJet/PhotonJet_Pt500toInf/*Counters.root");

  chains[22]->Add("results/WJetsALPGEN/W1Jets_Pt0to100-alpgen/*Counters.root");
  chains[23]->Add("results/WJetsALPGEN/W1Jets_Pt100to300-alpgen/*Counters.root");
  chains[24]->Add("results/WJetsALPGEN/W1Jets_Pt300to800-alpgen/*Counters.root");
  chains[25]->Add("results/WJetsALPGEN/W1Jets_Pt800to1600-alpgen/*Counters.root");

  chains[26]->Add("results/HiggsWW/H120_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[27]->Add("results/HiggsWW/H130_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[28]->Add("results/HiggsWW/H140_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[29]->Add("results/HiggsWW/H150_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[30]->Add("results/HiggsWW/H160_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[31]->Add("results/HiggsWW/H170_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[32]->Add("results/HiggsWW/H180_2W_2lnu_gluonfusion_7TeV/*Counters.root");
  chains[33]->Add("results/HiggsWW/H190_2W_2lnu_gluonfusion_7TeV/*Counters.root");

  chains[34]->Add("results/DiBosons/Wgamma/*Counters.root");
  chains[35]->Add("results/DiBosons/WW_2l_7TeV/*Counters.root");
  chains[36]->Add("results/DiBosons/WZ_3l_7TeV/*Counters.root");
  chains[37]->Add("results/DiBosons/ZZ_2l2nu/*Counters.root");

  cout << "chains added. " << endl;

  std::vector<std::string> sampleName;

  sampleName.push_back("WJetsMADGRAPH");

  sampleName.push_back("ZJetsMADGRAPH");

  sampleName.push_back("TTbar_mcatnlo");

  sampleName.push_back("QCD_QCD_EMEnriched_Pt20to30");
  sampleName.push_back("QCD_QCD_EMEnriched_Pt30to80");
  sampleName.push_back("QCD_QCD_EMEnriched_Pt80to170");

  sampleName.push_back("QCD_QCD_BCtoE_Pt20to30");
  sampleName.push_back("QCD_QCD_BCtoE_Pt30to80");
  sampleName.push_back("QCD_QCD_BCtoE_Pt80to170");

  sampleName.push_back("SingleTop_SingleTop_sChannel_madgraph");
  sampleName.push_back("SingleTop_SingleTop_tChannel_madgraph");
  sampleName.push_back("SingleTop_SingleTop_tWChannel_madgraph");

  sampleName.push_back("PhotonJet_Pt0to15");
  sampleName.push_back("PhotonJet_Pt15to20");
  sampleName.push_back("PhotonJet_Pt20to30");
  sampleName.push_back("PhotonJet_Pt30to50");
  sampleName.push_back("PhotonJet_Pt50to80");
  sampleName.push_back("PhotonJet_Pt80to120");
  sampleName.push_back("PhotonJet_Pt120to170");
  sampleName.push_back("PhotonJet_Pt170to300");
  sampleName.push_back("PhotonJet_Pt300to500");
  sampleName.push_back("PhotonJet_Pt500toInf");

  sampleName.push_back("W1Jets_Pt0to100-alpgen");
  sampleName.push_back("W1Jets_Pt100to300-alpgen");
  sampleName.push_back("W1Jets_Pt300to800-alpgen");
  sampleName.push_back("W1Jets_Pt800to1600-alpgen");

  sampleName.push_back("H120_2W_2lnu");
  sampleName.push_back("H130_2W_2lnu");
  sampleName.push_back("H140_2W_2lnu");
  sampleName.push_back("H150_2W_2lnu");
  sampleName.push_back("H160_2W_2lnu");
  sampleName.push_back("H170_2W_2lnu");
  sampleName.push_back("H180_2W_2lnu");
  sampleName.push_back("H190_2W_2lnu");

  sampleName.push_back("Wgamma");
  sampleName.push_back("WW_2l");
  sampleName.push_back("WZ_3l");
  sampleName.push_back("ZZ_2l2nu");

  float nEv[38];
  for(int isample=0; isample<38; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<38; isample++) {

    cout << "\tProcessing sample # " << isample << "..." << endl;

    Int_t           nCuts;
    Float_t         nSel[20];   //[nCuts]
    
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

  for(int isample=0; isample<38; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
  }
  
}

