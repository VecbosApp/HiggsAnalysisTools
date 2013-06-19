// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <algorithm>
#include <iostream>

// HiggsAnalysisTools includes
#include "LumiReweightingStandAlone.h"
#include "HiggsAnalysisTools/include/HWWKinematics.hh"

#include "addWeightsToTreetH.hh"

using namespace std;
using namespace reweight;

addWeightsToTreetH::addWeightsToTreetH(const char* filename, float baseW, int processId, int finalstate, int release) {
  filename_ = std::string(filename);
  baseW_ = baseW;
  processId_ = processId;
  finalstate_ = finalstate;
  release_ = release;
}

void addWeightsToTreetH::addWeights() {

  // 2011 sub-rns lumi
  Float_t lumiA = 2.1;
  Float_t lumiB = 2.5;

  cout << "Adding weight branch to file " << filename_ << " with weight " << baseW_ << endl;
  if (release_==0) cout << "Offline efficiency computed using 52X samples" << endl;
  if (release_==1) cout << "Offline efficiency computed using 42X samples" << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename_.c_str());
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("latino");
  } else {
    cout << "File " << filename_ << " not existing !" << endl;
    return;
  }

  LumiReWeighting LumiWeights( "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/s7pileup200.root",
                               "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/puRun2012_5100ipb_71.root",
                               "hNPU","pileup");
    
  // reading root files with electrons and muons efficiencies

  // Muons
  TFile fileSFmuons52("/afs/cern.ch/work/e/emanuele/public/effsf/muon_scale_factors_52X.root");
  TH2F *histoSFmuons52 = (TH2F*)fileSFmuons52.Get("muonDATAMCratio")->Clone("effSFmuons52");

  // Electrons
  TFile fileSFEle52("/afs/cern.ch/work/e/emanuele/public/effsf/electron_scale_factors_52X.root");
  TH2F *histoSFele52 = (TH2F*)fileSFEle52.Get("newhwwWP_ratio")->Clone("effSFele52");

  fileOrig->cd();

  if ( treeOrig ) {

    int nentriesOrig = treeOrig->GetEntries();
    
    TFile *fileNew = TFile::Open(filename_.c_str(),"recreate");
    TTree *treeNew = new TTree("latino","tree with 2 lepton selection");

    std::vector<TTree*> trees; 
    trees.push_back(treeNew);
    
    // Add also a branch with jet category (1 for njets=0, -1 for njets=1: useful for the fit)
    // and a branch with float final selection bool (for roofit)
    Int_t           run;
    Int_t           lumi;
    ULong64_t       event;
    Float_t         puweight;
    Float_t         met;
    Float_t         pfMet;
    Float_t         caloMet;
    Float_t         projMet;
    Float_t         deltaPhi;
    Float_t         deltaR;
    Float_t         transvMass, transvMassUp, transvMassDown;
    Float_t         eleInvMass;
    Float_t         maxPtEle;
    Float_t         minPtEle;
    Float_t         detaLeptons;
    Bool_t          finalLeptons;
    Bool_t          jetVeto;
    Bool_t          uncorrJetVeto;
    Bool_t          preDeltaPhi;
    Bool_t          finalSelection;
    Bool_t          step[29];
    Float_t         npu[3];
    Bool_t          hlt;
    Float_t         KFactor;
    Bool_t          promptDecay;

    Int_t           njets;
    Int_t           ncbIDjets;
    Int_t           nuncorrjets;

    Int_t numbtagCSVMmvaIDcentraljets;//MVA ID for jets
    Int_t numbtagCSVLmvaIDcentraljets;
    Int_t numbtagCSVTmvaIDcentraljets;

    Int_t numbtagCSVMmvaIDforwardjets;//MVA ID for jets
    Int_t numbtagCSVLmvaIDforwardjets;
    Int_t numbtagCSVTmvaIDforwardjets;

    Int_t numbtagCSVMcbIDcentraljets;//Cut Based ID for jets
    Int_t numbtagCSVLcbIDcentraljets;
    Int_t numbtagCSVTcbIDcentraljets;
    
    Int_t numbtagCSVMcbIDforwardjets;//Cut Based ID for jets
    Int_t numbtagCSVLcbIDforwardjets;
    Int_t numbtagCSVTcbIDforwardjets;

    Int_t           nVtx;
    Float_t         dxyEVT;
    Float_t         dszEVT;
    Float_t         softtche;
    Float_t         hardbjpb;
    Float_t         bTagSecVertex;
    Float_t         leadingJetBTagTrackCount;
    Float_t         subleadingJetBTagTrackCount;
    Float_t         subleadingJetsMaxBTagTrackCount;
    Float_t         leadingJetBTagBProb;
    Float_t         subleadingJetBTagBProb;
    Float_t         subleadingJetsMaxBTagBProb;
    Float_t         pt[3];
    Float_t         eta[3];
    Float_t         scEnergy[3];
    Float_t         R9[3];
    Float_t         pxChMet;
    Float_t         pyChMet;
    Float_t         pzChMet;
    Float_t         dymva1;
    Float_t         pxLeadJet[3];
    Float_t         pyLeadJet[3];
    Float_t         pzLeadJet[3];
    Float_t         pxSecondJet[3];
    Float_t         pySecondJet[3];
    Float_t         pzSecondJet[3];
    Float_t         jetmva1, jetmva2;
    Int_t           jetid1, jetid2;
    Float_t         pxL1;
    Float_t         pyL1;
    Float_t         pzL1;
    Float_t         pxL2;
    Float_t         pyL2;
    Float_t         pzL2;
    Float_t         pxL3;
    Float_t         pyL3;
    Float_t         pzL3;
    Float_t         eneL1;
    Float_t         eneL2;
    Float_t         eneL3;
    Int_t           typeL1;
    Int_t           typeL2;
    Int_t           typeL3;
    Int_t           nSoftMu;
    Int_t           nSoftMuNoJets;
    Int_t           nsoftjet;
    Int_t           nsoftbjet;
    Int_t           numExtraLep;

    // lepton quantities. Depending on the channel ...
    Int_t           ch[3];
    Float_t         lepiso[3];

    TLorentzVector *sumJetsV4 = 0;
    TLorentzVector *uncorrSumJetsV4 = 0;
    TVector3       *pfmetV = 0;
    float          pfmetX, pfmetY;

    Float_t         signPFMet;
    Float_t         signPFChargedMet;
    Float_t         mtrchargedMet;

    // Fake rate weights (only for W+jets selection)
    Float_t         weightFP;
    Float_t         weightStatFP;
    Float_t         weightFF;
    Float_t         weightStatFF;
    Float_t         weightPP;
    Float_t         weightStatPP;
    //
    Float_t         weightFP15;
    Float_t         weightStatFP15;
    Float_t         weightFF15;
    Float_t         weightStatFF15;
    Float_t         weightPP15;
    Float_t         weightStatPP15;
    //
    Float_t         weightFP50;
    Float_t         weightStatFP50;
    Float_t         weightFF50;
    Float_t         weightStatFF50;
    Float_t         weightPP50;
    Float_t         weightStatPP50;
    //
    Float_t         weightFPQCD;
    Float_t         weightStatFPQCD;
    Float_t         weightFFQCD;
    Float_t         weightStatFFQCD;
    Float_t         weightPPQCD;
    Float_t         weightStatPPQCD;
    //
    Int_t           tight;

    // DY generator level quantities
    Float_t         genmll;
    Float_t         genptll;
    Float_t         genyll;

    // PDF weights
    //     Double_t       cteq66W[45];
    //     Double_t       mstwW[31];
    //     Double_t       nnpdfW[101];
    
    treeOrig->SetBranchAddress("run", &run);
    treeOrig->SetBranchAddress("lumi", &lumi);
    treeOrig->SetBranchAddress("event", &event);
    treeOrig->SetBranchAddress("puweight", &puweight);
    treeOrig->SetBranchAddress("genmll", &genmll);
    treeOrig->SetBranchAddress("genptll", &genptll);
    treeOrig->SetBranchAddress("genyll", &genyll);
    treeOrig->SetBranchAddress("hlt", &hlt);
    treeOrig->SetBranchAddress("met", &met);
    treeOrig->SetBranchAddress("pfMet", &pfMet);
    treeOrig->SetBranchAddress("caloMet", &caloMet);
    treeOrig->SetBranchAddress("projMet", &projMet);
    treeOrig->SetBranchAddress("deltaPhi", &deltaPhi);
    treeOrig->SetBranchAddress("deltaR", &deltaR);
    treeOrig->SetBranchAddress("transvMass", &transvMass);
    treeOrig->SetBranchAddress("transvMassUp", &transvMassUp);
    treeOrig->SetBranchAddress("transvMassDown", &transvMassDown);
    treeOrig->SetBranchAddress("eleInvMass", &eleInvMass);
    treeOrig->SetBranchAddress("maxPtEle", &maxPtEle);
    treeOrig->SetBranchAddress("minPtEle", &minPtEle);
    treeOrig->SetBranchAddress("detaLeptons", &detaLeptons);
    treeOrig->SetBranchAddress("nVtx", &nVtx);
    treeOrig->SetBranchAddress("finalLeptons", &finalLeptons);
    treeOrig->SetBranchAddress("jetVeto", &jetVeto);
    treeOrig->SetBranchAddress("uncorrJetVeto", &uncorrJetVeto);
    treeOrig->SetBranchAddress("preDeltaPhi", &preDeltaPhi);
    treeOrig->SetBranchAddress("finalSelection", &finalSelection);
    treeOrig->SetBranchAddress("step", step);
    treeOrig->SetBranchAddress("npu", npu);
    treeOrig->SetBranchAddress("KFactor", &KFactor);
    treeOrig->SetBranchAddress("promptDecay", &promptDecay);
    treeOrig->SetBranchAddress("njets", &njets);
    treeOrig->SetBranchAddress("ncbIDjets", &ncbIDjets);
    treeOrig->SetBranchAddress("nuncorrjets", &nuncorrjets);

    treeOrig->SetBranchAddress("numbtagCSVMmvaIDcentraljets", &numbtagCSVMmvaIDcentraljets);
    treeOrig->SetBranchAddress("numbtagCSVLmvaIDcentraljets", &numbtagCSVLmvaIDcentraljets);
    treeOrig->SetBranchAddress("numbtagCSVTmvaIDcentraljets", &numbtagCSVTmvaIDcentraljets);

    treeOrig->SetBranchAddress("numbtagCSVMmvaIDforwardjets", &numbtagCSVMmvaIDforwardjets);
    treeOrig->SetBranchAddress("numbtagCSVLmvaIDforwardjets", &numbtagCSVLmvaIDforwardjets);
    treeOrig->SetBranchAddress("numbtagCSVTmvaIDforwardjets", &numbtagCSVTmvaIDforwardjets);

    treeOrig->SetBranchAddress("numbtagCSVMcbIDcentraljets", &numbtagCSVMcbIDcentraljets);
    treeOrig->SetBranchAddress("numbtagCSVLcbIDcentraljets", &numbtagCSVLcbIDcentraljets);
    treeOrig->SetBranchAddress("numbtagCSVTcbIDcentraljets", &numbtagCSVTcbIDcentraljets);

    treeOrig->SetBranchAddress("numbtagCSVMcbIDforwardjets", &numbtagCSVMcbIDforwardjets);
    treeOrig->SetBranchAddress("numbtagCSVLcbIDforwardjets", &numbtagCSVLcbIDforwardjets);
    treeOrig->SetBranchAddress("numbtagCSVTcbIDforwardjets", &numbtagCSVTcbIDforwardjets);

    treeOrig->SetBranchAddress("dxyEVT", &dxyEVT);
    treeOrig->SetBranchAddress("dszEVT", &dszEVT);
    treeOrig->SetBranchAddress("softtche", &softtche);
    treeOrig->SetBranchAddress("hardbjpb", &hardbjpb);
    treeOrig->SetBranchAddress("bTagSecVertex", &bTagSecVertex);
    treeOrig->SetBranchAddress("leadingJetBTagTrackCount", &leadingJetBTagTrackCount);
    treeOrig->SetBranchAddress("subleadingJetBTagTrackCount", &subleadingJetBTagTrackCount);
    treeOrig->SetBranchAddress("subleadingJetsMaxBTagTrackCount", &subleadingJetsMaxBTagTrackCount);
    treeOrig->SetBranchAddress("leadingJetBTagJetBProb", &leadingJetBTagBProb);
    treeOrig->SetBranchAddress("subleadingJetBTagJetBProb", &subleadingJetBTagBProb);
    treeOrig->SetBranchAddress("subleadingJetsMaxBTagJetBProb", &subleadingJetsMaxBTagBProb);
    treeOrig->SetBranchAddress("pt", pt);
    treeOrig->SetBranchAddress("eta", eta);
    treeOrig->SetBranchAddress("scEnergy", scEnergy);
    treeOrig->SetBranchAddress("R9", R9);
    treeOrig->SetBranchAddress("pxChMet", &pxChMet);
    treeOrig->SetBranchAddress("pyChMet", &pyChMet);
    treeOrig->SetBranchAddress("pzChMet", &pzChMet); 
    treeOrig->SetBranchAddress("dymva1", &dymva1);
    treeOrig->SetBranchAddress("pxLeadJet", pxLeadJet);
    treeOrig->SetBranchAddress("pyLeadJet", pyLeadJet);
    treeOrig->SetBranchAddress("pzLeadJet", pzLeadJet);
    treeOrig->SetBranchAddress("pxSecondJet", pxSecondJet);
    treeOrig->SetBranchAddress("pySecondJet", pySecondJet);
    treeOrig->SetBranchAddress("pzSecondJet", pzSecondJet);
    treeOrig->SetBranchAddress("jetmva1", &jetmva1);
    treeOrig->SetBranchAddress("jetmva2", &jetmva2);
    treeOrig->SetBranchAddress("jetid1", &jetid1);
    treeOrig->SetBranchAddress("jetid2", &jetid2);
    treeOrig->SetBranchAddress("pxL1", &pxL1);
    treeOrig->SetBranchAddress("pyL1", &pyL1);
    treeOrig->SetBranchAddress("pzL1", &pzL1);
    treeOrig->SetBranchAddress("pxL2", &pxL2);
    treeOrig->SetBranchAddress("pyL2", &pyL2);
    treeOrig->SetBranchAddress("pzL2", &pzL2);
    treeOrig->SetBranchAddress("pxL3", &pxL3);
    treeOrig->SetBranchAddress("pyL3", &pyL3);
    treeOrig->SetBranchAddress("pzL3", &pzL3);
    treeOrig->SetBranchAddress("eneL1", &eneL1);
    treeOrig->SetBranchAddress("eneL2", &eneL2);
    treeOrig->SetBranchAddress("eneL3", &eneL3);
    treeOrig->SetBranchAddress("typeL1", &typeL1);
    treeOrig->SetBranchAddress("typeL2", &typeL2);
    treeOrig->SetBranchAddress("typeL3", &typeL3);
    treeOrig->SetBranchAddress("nSoftMu", &nSoftMu);
    treeOrig->SetBranchAddress("nSoftMuNoJets", &nSoftMuNoJets);
    treeOrig->SetBranchAddress("nsoftjet", &nsoftjet);
    treeOrig->SetBranchAddress("nsoftbjet", &nsoftbjet);
    treeOrig->SetBranchAddress("numExtraLep", &numExtraLep);
    treeOrig->SetBranchAddress("sumJetsV4", &sumJetsV4);
    treeOrig->SetBranchAddress("uncorrSumJetsV4", &uncorrSumJetsV4);
    treeOrig->SetBranchAddress("pfmetV", &pfmetV);
    treeOrig->SetBranchAddress("ch", ch);
    treeOrig->SetBranchAddress("lepiso", lepiso);
    treeOrig->SetBranchAddress("signPFMet", &signPFMet);
    treeOrig->SetBranchAddress("signPFChargedMet", &signPFChargedMet);
    treeOrig->SetBranchAddress("mtrchargedMet", &mtrchargedMet);
    //     treeOrig->SetBranchAddress("cteq66W", cteq66W);
    //     treeOrig->SetBranchAddress("mstwW",mstwW);
    //     treeOrig->SetBranchAddress("nnpdfW",nnpdfW);
    
    // electron ID (only filled for electrons)
    Float_t pt_1,       pt_2;
    Float_t eta_1,      eta_2;
    Float_t R9_1,       R9_2;
    Float_t scEnergy_1, scEnergy_2;

    // convert the booleans into integers (to insert in RooDataset)
    Int_t         i_jetVeto;
    Int_t         i_uncorrJetVeto;
    Int_t         i_preDeltaPhi;
    Int_t         i_finalSelection;
    Int_t         i_promptDecay;
    Int_t         i_WWSel0j;
    Int_t         i_WWSel1j;
    Int_t         i_hlt;

    // other kinematics
    float deltaPhi_LL_MET;
    float deltaPhi_LLJ1_MET;
    float deltaPhi_LL_JET1;
    float deltaPhi_LL_JET2;
    float deltaPhi_MET_JET1;
    float deltaPhi_MET_JET2;
    float deltaPhi_LL_JJ;
    float dphilmet1, dphilmet2;
    float jetpt1, jeteta1, jetphi1;
    float jetpt2, jeteta2, jetphi2;
    float chmet;
    float pmet, pmet2;
    float L1pt, L1eta, L1phi;
    float L2pt, L2eta, L2phi;
    float L3pt, L3eta, L3phi;
    float dileptonPt; 
    float trileptonPt; 
    float mtw1, mtw2;
    float mr,mrnew,dphillr,ddphillr;
    float jetcat = 0;
    float consecevent = -1;
    float iMet;
    float dphill; // deltaphi in radians
    
    // variables to be converted in float...
    float f_run, f_lumi;
    float f_hlt, f_nVtx, f_njets, f_ncbIDjets, f_nuncorrjets;

    float f_numbtagCSVMmvaIDcentraljets, f_numbtagCSVLmvaIDcentraljets, f_numbtagCSVTmvaIDcentraljets;
    float f_numbtagCSVMmvaIDforwardjets, f_numbtagCSVLmvaIDforwardjets, f_numbtagCSVTmvaIDforwardjets;  
    float f_numbtagCSVMcbIDcentraljets, f_numbtagCSVLcbIDcentraljets, f_numbtagCSVTcbIDcentraljets;
    float f_numbtagCSVMcbIDforwardjets, f_numbtagCSVLcbIDforwardjets, f_numbtagCSVTcbIDforwardjets;
 
    float f_zveto, f_bveto_ip, f_bveto_mu, f_bveto_munj, f_bveto, f_dphiveto, f_typeL1, f_typeL2, f_typeL3,
      f_nSoftMu, f_nSoftMuNoJets, f_numExtraLep, f_finalstate, f_processId, sameflav, f_nsoftbjet, f_nsoftjet;
    float f_ch[3];
    
    // vetoes
    int zveto, bveto_ip, bveto_mu, bveto_munj, bveto;

    // additional (dummy for the moment)
    Float_t effW   = 1.0;   
    Float_t triggW = 1.0;

    // pileup
    float puW;
    
    for(int i=0; i<(int)trees.size();i++) {
      TTree *theTreeNew = trees[i];

      // the selected final state: eee = 0, mmm = 1, eem = 2, mme = 3
      theTreeNew->Branch("channel",  &f_finalstate, "channel/F");
      theTreeNew->Branch("sameflav", &sameflav, "sameflav/F");

      // one integer containing the process identifier (for MC, 0 for data)
      theTreeNew->Branch("dataset", &f_processId, "dataset/F");


      // Copys branches
      theTreeNew->Branch("run", &f_run, "run/F");
      theTreeNew->Branch("lumi", &f_lumi, "lumi/F");
      theTreeNew->Branch("event", &event, "event/l");
      theTreeNew->Branch("genmll", &genmll, "genmll/F");
      theTreeNew->Branch("genptll", &genptll, "genptll/F");
      theTreeNew->Branch("genyll", &genyll, "genyll/F");
      theTreeNew->Branch("puW", &puW, "puW/F");
      theTreeNew->Branch("effW", &effW, "effW/F");
      theTreeNew->Branch("triggW", &triggW, "triggW/F");
      theTreeNew->Branch("pfmet", &pfMet, "pfmet/F");
      theTreeNew->Branch("chmet", &chmet, "chmet/F");  
      theTreeNew->Branch("dymva1", &dymva1, "dymva1/F");
      theTreeNew->Branch("trigger", &f_hlt, "trigger/F");
      theTreeNew->Branch("caloMet", &caloMet, "caloMet/F");
      theTreeNew->Branch("ppfmet", &pmet, "ppfmet/F");
      theTreeNew->Branch("pchmet", &pmet2, "pchmet/F");
      theTreeNew->Branch("mpmet", &projMet, "mpmet/F");
      theTreeNew->Branch("dphill", &dphill, "dphill/F");
      theTreeNew->Branch("drll", &deltaR, "drll/F");
      theTreeNew->Branch("mll", &eleInvMass, "mll/F");
      theTreeNew->Branch("mth", &transvMass, "mth/F");
      theTreeNew->Branch("mthUp", &transvMassUp, "mthUp/F");
      theTreeNew->Branch("mthDown", &transvMassDown, "mthDown/F");
      theTreeNew->Branch("mtw1", &mtw1, "mtw1/F");
      theTreeNew->Branch("mtw2", &mtw2, "mtw2/F");
      theTreeNew->Branch("pt1" , &L1pt, "pt1/F");
      theTreeNew->Branch("pt2",  &L2pt, "pt2/F");
      theTreeNew->Branch("pt3",  &L3pt, "pt3/F");
      theTreeNew->Branch("eta1", &L1eta, "eta1/F");
      theTreeNew->Branch("eta2", &L2eta, "eta2/F");
      theTreeNew->Branch("eta3", &L3eta, "eta3/F");
      theTreeNew->Branch("phi1", &L1phi, "phi1/F");
      theTreeNew->Branch("phi2", &L2phi, "phi2/F");
      theTreeNew->Branch("phi3", &L3phi, "phi3/F");
      theTreeNew->Branch("detaLeptons", &detaLeptons, "detaLeptons/F");
      theTreeNew->Branch("nvtx", &f_nVtx, "nvtx/F");
      theTreeNew->Branch("jetVeto", &i_jetVeto, "jetVeto/I");
      theTreeNew->Branch("uncorrJetVeto", &i_uncorrJetVeto, "uncorrJetVeto/I");
      theTreeNew->Branch("preDeltaPhi", &i_preDeltaPhi, "preDeltaPhi/I");
      theTreeNew->Branch("finalSelection", &i_finalSelection, "finalSelection/I");
      theTreeNew->Branch("WWSel0j", &i_WWSel0j, "WWSel0j/I");
      theTreeNew->Branch("WWSel1j", &i_WWSel1j, "WWSel1j/I");
      theTreeNew->Branch("kfW", &KFactor, "kfW/F");
      theTreeNew->Branch("npu", npu, "npu[3]/F");
      theTreeNew->Branch("promptDecay", &i_promptDecay, "promptDecay/I");
      theTreeNew->Branch("njet", &f_njets, "njet/F");
      theTreeNew->Branch("ncbIDjet", &f_ncbIDjets, "ncbIDjet/F");
      
      theTreeNew->Branch("numbtagCSVMmvaIDcentraljets", &f_numbtagCSVMmvaIDcentraljets, "numbtagCSVMmvaIDcentraljets/F");
      theTreeNew->Branch("numbtagCSVLmvaIDcentraljets", &f_numbtagCSVLmvaIDcentraljets, "numbtagCSVLmvaIDcentraljets/F");
      theTreeNew->Branch("numbtagCSVTmvaIDcentraljets", &f_numbtagCSVTmvaIDcentraljets, "numbtagCSVTmvaIDcentraljets/F");

      theTreeNew->Branch("numbtagCSVMmvaIDforwardjets", &f_numbtagCSVMmvaIDforwardjets, "numbtagCSVMmvaIDforwardjets/F");
      theTreeNew->Branch("numbtagCSVLmvaIDforwardjets", &f_numbtagCSVLmvaIDforwardjets, "numbtagCSVLmvaIDforwardjets/F");
      theTreeNew->Branch("numbtagCSVTmvaIDforwardjets", &f_numbtagCSVTmvaIDforwardjets, "numbtagCSVTmvaIDforwardjets/F");

      theTreeNew->Branch("numbtagCSVMcbIDcentraljets", &f_numbtagCSVMcbIDcentraljets, "numbtagCSVMcbIDcentraljets/F");
      theTreeNew->Branch("numbtagCSVLcbIDcentraljets", &f_numbtagCSVLcbIDcentraljets, "numbtagCSVLcbIDcentraljets/F");
      theTreeNew->Branch("numbtagCSVTcbIDcentraljets", &f_numbtagCSVTcbIDcentraljets, "numbtagCSVTcbIDcentraljets/F");

      theTreeNew->Branch("numbtagCSVMcbIDforwardjets", &f_numbtagCSVMcbIDforwardjets, "numbtagCSVMcbIDforwardjets/F");
      theTreeNew->Branch("numbtagCSVLcbIDforwardjets", &f_numbtagCSVLcbIDforwardjets, "numbtagCSVLcbIDforwardjets/F");
      theTreeNew->Branch("numbtagCSVTcbIDforwardjets", &f_numbtagCSVTcbIDforwardjets, "numbtagCSVTcbIDforwardjets/F");

      theTreeNew->Branch("nuncorrjets", &f_nuncorrjets, "nuncorrjets/F");
      theTreeNew->Branch("dxyEVT", &dxyEVT, "dxyEVT/F");
      theTreeNew->Branch("dszEVT", &dszEVT, "dszEVT/F");
      theTreeNew->Branch("softtche", &softtche, "softtche/F");
      theTreeNew->Branch("hardbjpb", &hardbjpb, "hardbjpb/F");
      theTreeNew->Branch("bTagSecVertex", &bTagSecVertex, "bTagSecVertex/F");
      theTreeNew->Branch("leadingJetBTagTrackCount", &leadingJetBTagTrackCount, "leadingJetBTagTrackCount/F");
      theTreeNew->Branch("subleadingJetBTagTrackCount", &subleadingJetBTagTrackCount, "subleadingJetBTagTrackCount/F");
      theTreeNew->Branch("subleadingJetsMaxBTagTrackCount", &subleadingJetsMaxBTagTrackCount, "subleadingJetsMaxBTagTrackCount/F");
      theTreeNew->Branch("leadingJetBTagBProb", &leadingJetBTagBProb, "leadingJetBTagBProb/F");
      theTreeNew->Branch("subleadingJetBTagBProb", &subleadingJetBTagBProb, "subleadingJetBTagBProb/F");
      theTreeNew->Branch("subleadingJetsMaxBTagBProb", &subleadingJetsMaxBTagBProb, "subleadingJetsMaxBTagBProb/F");
      theTreeNew->Branch("step", step, "step[29]/O");
      theTreeNew->Branch("zveto", &f_zveto, "zveto/F");
      theTreeNew->Branch("bveto_ip", &f_bveto_ip, "bveto_ip/F");
      theTreeNew->Branch("bveto_mu", &f_bveto_mu, "bveto_mu/F");
      theTreeNew->Branch("bveto_munj", &f_bveto_munj, "bveto_munj/F");
      theTreeNew->Branch("bveto", &f_bveto, "bveto/F");
      theTreeNew->Branch("dphiveto", &f_dphiveto, "dphiveto/F");
      theTreeNew->Branch("gammaMRStar", &mr, "gammaMRStar/F");
      theTreeNew->Branch("mr", &mrnew, "mr/F");
      theTreeNew->Branch("dphillr", &dphillr, "dphillr/F");
      theTreeNew->Branch("ddphillr", &ddphillr, "ddphillr/F");
      theTreeNew->Branch("pfmetsign", &signPFMet, "pfmetsign/F");
      theTreeNew->Branch("chmetsign", &signPFChargedMet, "chmetsign/F");
      theTreeNew->Branch("pt1_eleid", &pt_1, "pt1_eleid/F");
      theTreeNew->Branch("eta1_eleid", &eta_1, "eta1_eleid/F");
      theTreeNew->Branch("R91_eleid", &R9_1, "R91_eleid/F");
      theTreeNew->Branch("scEnergy1_eleid", &scEnergy_1, "scEnergy1_eleid/F");
      theTreeNew->Branch("pt2_eleid", &pt_2, "pt2_eleid/F");
      theTreeNew->Branch("eta2_eleid", &eta_2, "eta2_eleid/F");
      theTreeNew->Branch("R92_eleid", &R9_2, "R92_eleid/F");
      theTreeNew->Branch("scEnergy2_eleid", &scEnergy_2, "scEnergy2_eleid/F");
      theTreeNew->Branch("pxChMet", &pxChMet, "pxChMet/F");
      theTreeNew->Branch("pyChMet", &pyChMet, "pyChMet/F");
      theTreeNew->Branch("pzChMet", &pzChMet, "pzChMet/F");
      theTreeNew->Branch("pxLeadJet", pxLeadJet, "pxLeadJet[3]/F");
      theTreeNew->Branch("pyLeadJet", pyLeadJet, "pyLeadJet[3]/F");
      theTreeNew->Branch("pzLeadJet", pzLeadJet, "pzLeadJet[3]/F");
      theTreeNew->Branch("pxSecondJet", pxSecondJet, "pxSecondJet[3]/F");
      theTreeNew->Branch("pySecondJet", pySecondJet, "pySecondJet[3]/F");
      theTreeNew->Branch("pzSecondJet", pzSecondJet, "pzSecondJet[3]/F");
      theTreeNew->Branch("jetmva1", &jetmva1, "jetmva1/F");
      theTreeNew->Branch("jetmva2", &jetmva2, "jetmva2/F");
      theTreeNew->Branch("jetid1", &jetid1, "jetid1/I");
      theTreeNew->Branch("jetid2", &jetid2, "jetid2/I");
      theTreeNew->Branch("pxL1", &pxL1, "pxL1/F");
      theTreeNew->Branch("pyL1", &pyL1, "pyL1/F");
      theTreeNew->Branch("pzL1", &pzL1, "pzL1/F");
      theTreeNew->Branch("eneL1",  &eneL1,  "eneL1/F");
      theTreeNew->Branch("typeL1", &f_typeL1, "typeL1/F");
      theTreeNew->Branch("pxL2", &pxL2, "pxL2/F");
      theTreeNew->Branch("pyL2", &pyL2, "pyL2/F");
      theTreeNew->Branch("pzL2", &pzL2, "pzL2/F");
      theTreeNew->Branch("eneL2",  &eneL2,  "eneL2/F");
      theTreeNew->Branch("typeL2", &f_typeL2, "typeL2/F");
      theTreeNew->Branch("dphilmet1", &dphilmet1, "dphilmet1/F");
      theTreeNew->Branch("dphilmet2", &dphilmet2, "dphilmet2/F");
      theTreeNew->Branch("dphillmet", &deltaPhi_LL_MET, "dphillmet/F");
      theTreeNew->Branch("deltaPhi_LLJ1_MET", &deltaPhi_LLJ1_MET, "deltaPhi_LLJ1_MET/F");
      theTreeNew->Branch("dphilljet", &deltaPhi_LL_JET1, "dphilljet/F");
      theTreeNew->Branch("dphilljet2", &deltaPhi_LL_JET2, "dphilljet2/F");
      theTreeNew->Branch("deltaPhi_MET_JET1", &deltaPhi_MET_JET1, "deltaPhi_MET_JET1/F");
      theTreeNew->Branch("deltaPhi_MET_JET2", &deltaPhi_MET_JET2, "deltaPhi_MET_JET2/F");
      theTreeNew->Branch("deltaPhi_LL_JJ", &deltaPhi_LL_JJ, "deltaPhi_LL_JJ/F");
      theTreeNew->Branch("ptll" , &dileptonPt, "ptll/F");
      theTreeNew->Branch("ptlll", &trileptonPt,"ptlll/F");
      theTreeNew->Branch("jetpt1", &jetpt1, "jetpt1/F");
      theTreeNew->Branch("jeteta1", &jeteta1, "jeteta1/F");
      theTreeNew->Branch("jetphi1", &jetphi1, "jetphi1/F");
      theTreeNew->Branch("jetpt2", &jetpt2, "jetpt2/F");
      theTreeNew->Branch("jeteta2", &jeteta2, "jeteta2/F");
      theTreeNew->Branch("jetphi2", &jetphi2, "jetphi2/F");
      theTreeNew->Branch("nSoftMu", &f_nSoftMu, "nSoftMu/I");
      theTreeNew->Branch("nSoftMuNoJets", &f_nSoftMuNoJets, "nSoftMuNoJets/F");
      theTreeNew->Branch("nsoftjet", &f_nsoftjet, "nsoftjet/F");
      theTreeNew->Branch("nsoftbjet", &f_nsoftbjet, "nsoftbjet/F");
      theTreeNew->Branch("nextra", &f_numExtraLep, "nextra/F");
      theTreeNew->Branch("sumJetsV4", "TLorentzVector", &sumJetsV4);
      theTreeNew->Branch("uncorrSumJetsV4", "TLorentzVector", &uncorrSumJetsV4);
      theTreeNew->Branch("pfmetV", "TVector3", &pfmetV);
      theTreeNew->Branch("pfmetX", &pfmetX, "pfmetX/F");
      theTreeNew->Branch("pfmetY", &pfmetY, "pfmetY/F");
      theTreeNew->Branch("ch1", &(f_ch[0]), "ch1/F");
      theTreeNew->Branch("ch2", &(f_ch[1]), "ch2/F");
      theTreeNew->Branch("ch3", &(f_ch[2]), "ch3/F");
      theTreeNew->Branch("iso1", &(lepiso[0]), "iso1/F");
      theTreeNew->Branch("iso2", &(lepiso[1]), "iso2/F");
      theTreeNew->Branch("iso3", &(lepiso[2]), "iso3/F");
      theTreeNew->Branch("jetcat", &jetcat,  "jetcat/F");

      theTreeNew->Branch("consecevent", &consecevent, "consecevent/F");
      theTreeNew->Branch("baseW", &baseW_,  "baseW/F");
      
      //       theTreeNew->Branch("cteq66W", cteq66W, "cteq66W[45]/D");
      //       theTreeNew->Branch("mstwW",   mstwW,   "mstwW[31]/D");
      //       theTreeNew->Branch("nnpdfW",  nnpdfW,  "nnpdfW[101]/D");

      theTreeNew->Branch("iMet", &iMet, "iMet/F");

    }

    int j =0;

    // event container. It happens that some events are duplicate, and are consecutive in the same run
    //     std::vector<ULong64_t> eventsInRun;
    //     int lastrun=0;

    for(int i=0; i<nentriesOrig; i++) {
      if (i%10000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
      treeOrig->GetEntry(i);
      
      /// do this check only in data
      //if(processId>=100 && processId<1000) {
      //  if(run!=lastrun) {
      //    eventsInRun.clear();
      //    lastrun=run;
      //  }
      //  vector<ULong64_t>::iterator it;
      //  it = find(eventsInRun.begin(), eventsInRun.end(), event);
      //  if(it!=eventsInRun.end()) {
      //    cout << "Found dup. Run = " << run << "   event = " << event << ". Skipping it." << endl;
      //    continue;
      //  }
      //  eventsInRun.push_back(event);
      //}

      if     (njets==0) jetcat = 1;
      else if(njets==1) jetcat = -1;
      else jetcat = -2;

      pt_1       = pt[0];
      eta_1      = eta[0];
      scEnergy_1 = scEnergy[0];
      R9_1       = R9[0];
      pt_2       = pt[1];
      eta_2      = eta[1];
      scEnergy_2 = scEnergy[1];
      R9_2       = R9[1];

      // consider only events with 0 or 1 jet
      // and fit variables within fit range
      //     if (njets<=1 && 
      //         met>=0 && met<=200 &&
      //         deltaPhi>=0 && deltaPhi<=180 &&
      //         maxPtEle>=20 && maxPtEle<=200 &&
      //         eleInvMass>=12 && eleInvMass<=150 &&
      //         bTagImpPar>=-1001 && bTagImpPar<=2 &&
      //         finalSelection) {
      
      i_jetVeto = (jetVeto) ? 1 : 0;
      i_uncorrJetVeto = (uncorrJetVeto) ? 1 : 0;
      i_preDeltaPhi = (preDeltaPhi) ? 1 : 0;
      i_finalSelection = (finalSelection) ? 1 : 0;
      i_promptDecay = (promptDecay) ? 1 : 0;
      i_hlt = (hlt) ? 1 : 0;

      pfmetX = pfmetV->Px();
      pfmetY = pfmetV->Py();

      zveto = bveto_ip = bveto_mu = bveto_munj = bveto = -999;

      TVector3 TV_L1( pxL1, pyL1, pzL1 );
      TVector3 TV_L2( pxL2, pyL2, pzL2 );
      TVector3 TV_L3( pxL3, pyL3, pzL3 );
      TVector3 TV_L1p2   = TV_L1 + TV_L2;
      TVector3 TV_L1p2p3 = TV_L1 + TV_L2 + TV_L3;
      TVector3 TV_chmet( pxChMet, pyChMet, pzChMet );
      TVector3 TV_jet1( pxLeadJet[0], pyLeadJet[0], pzLeadJet[0] );
      TVector3 TV_jet2( pxSecondJet[0], pySecondJet[0], pzSecondJet[0] );
      deltaPhi_LL_MET   = (180./3.14) * TV_chmet.DeltaPhi(TV_L1p2);
      float l1eta = TV_L1.Eta();
      float l2eta = TV_L2.Eta();
      float l3eta = TV_L3.Eta();
      float l1pt  = TV_L1.Pt();
      float l2pt  = TV_L2.Pt();
      float l3pt  = TV_L3.Pt();

      zveto = (fabs(eleInvMass-91.1876)>15) ? 1 : 0;
      bveto_ip = 1;
      if(softtche>=2.1) bveto_ip = 0;
      // if(hardbjpb>=1.05) bveto_ip = 0; // revert to 2011 definition
      // and fix the temporary bug in the definition of maxTCHE for 1 jet bin
      if(njets==1 && leadingJetBTagTrackCount>2.1) {
        bveto_ip=0;
        for(int s=13; s<29; s++) step[s]=0;
      }
      bveto_mu = (nSoftMu==0) ? 1 : 0;
      bveto = (bveto_ip && bveto_mu) ? 1 : 0;
      bveto_munj = (nSoftMuNoJets==0) ? 1 : 0;

      dphilmet1 = fabs(pfmetV->DeltaPhi(TV_L1));
      dphilmet2 = fabs(pfmetV->DeltaPhi(TV_L2));

      chmet = TV_chmet.Pt();    
      pmet = GetProjectedMet(*pfmetV,TV_L1,TV_L2);
      pmet2 = GetProjectedMet(TV_chmet,TV_L1,TV_L2);

      mtw1 = calcMT(*pfmetV,TV_L1);
      mtw2 = calcMT(*pfmetV,TV_L2);
        
      iMet = projMet * cos(TV_chmet.Angle(*pfmetV));

      dphill = deltaPhi * TMath::Pi() / 180.;
      //! razor-like variables
      TLorentzVector FV_L1(TV_L1,eneL1);
      TLorentzVector FV_L2(TV_L2,eneL2);
      HWWKinematics kine(FV_L1,FV_L2,*pfmetV);
      mr = 2*kine.CalcMR();
      mrnew = 2*kine.CalcMRNEW();
      dphillr = fabs(kine.CalcDeltaPhiRFRAME());
      ddphillr = fabs(kine.CalcDoubleDphiRFRAME());

      // PU weights
      puW = LumiWeights.weight(npu[1]);

      //  offline efficiency scale factors
      //  for 3 leptons for 52X only
      Float_t eff1=1.; 
      Float_t eff2=1.;
      Float_t eff3=1.;
      Float_t effA1, effA2, effB1, effB2;
      effA1 = effA2 = effB1 = effB2 = 1.;
      if (processId_>0) { // MC => apply scale factors
        if (finalstate_ == 0) {   // eee
	  eff1 = getOfflineEff(l1pt, l1eta, histoSFele52);    
	  eff2 = getOfflineEff(l2pt, l2eta, histoSFele52);    
	  eff3 = getOfflineEff(l3pt, l3eta, histoSFele52);    
        } else if (finalstate_ == 1) { // mmm
	  eff1 = getOfflineEff(l1pt, l1eta, histoSFmuons52);
	  eff2 = getOfflineEff(l2pt, l2eta, histoSFmuons52);
	  eff3 = getOfflineEff(l3pt, l3eta, histoSFmuons52);
        } else if (finalstate_ == 2) { // eem
	  eff1 = getOfflineEff(l1pt, l1eta, histoSFele52);
	  eff2 = getOfflineEff(l2pt, l2eta, histoSFele52);
	  eff3 = getOfflineEff(l3pt, l3eta, histoSFmuons52);
        } else if (finalstate_ == 3) { // mme
	  eff1 = getOfflineEff(l1pt, l1eta, histoSFmuons52);
	  eff2 = getOfflineEff(l2pt, l2eta, histoSFele52);
	  eff3 = getOfflineEff(l3pt, l3eta, histoSFele52);
        } 
        effW = eff1*eff2*eff3;
      } else { // data
        effW = 1.;
      }

      if(sqrt(pow(pxLeadJet[0],2)+pow(pyLeadJet[0],2))>15) {
        jetpt1 = TV_jet1.Pt();
        jeteta1 = TV_jet1.Eta();
        jetphi1 = TV_jet1.Phi();
        TVector3 TV_L12pJ1 = TV_L1p2 + TV_jet1;
        deltaPhi_LLJ1_MET = (180./TMath::Pi()) * fabs(TV_chmet.DeltaPhi(TV_L12pJ1));   
        deltaPhi_MET_JET1 = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_chmet));
        deltaPhi_LL_JET1  = (180./TMath::Pi()) * fabs(TV_jet1.DeltaPhi(TV_L1p2));
      } else {
        deltaPhi_LLJ1_MET = -1.;
        deltaPhi_MET_JET1 = -1.;
        deltaPhi_LL_JET1 = -1.;
      }
      
      if(sqrt(pow(pxSecondJet[0],2)+pow(pySecondJet[0],2))>15) {
        jetpt2 = TV_jet2.Pt();
        jeteta2 = TV_jet2.Eta();
        jetphi2 = TV_jet2.Phi();
        deltaPhi_MET_JET2 = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_chmet));
        deltaPhi_LL_JET2  = (180./TMath::Pi()) * fabs(TV_jet2.DeltaPhi(TV_L1p2));
      } else {
        deltaPhi_MET_JET2 = -1.;
        deltaPhi_LL_JET2 = -1.;
      }

      if(sqrt(pow(pxLeadJet[0],2)+pow(pyLeadJet[0],2))>15 && sqrt(pow(pxSecondJet[0],2)+pow(pySecondJet[0],2))>15) {
        TVector3 TV_J1p2  = TV_jet1 + TV_jet2;
        deltaPhi_LL_JJ    = (180./3.14) * TV_L1p2.DeltaPhi(TV_J1p2);
      } else {
        deltaPhi_LL_JJ = -1.;
      }

      dileptonPt   = TV_L1p2.Pt();
      trileptonPt  = TV_L1p2p3.Pt();

      L1pt  = TV_L1.Pt();
      L2pt  = TV_L2.Pt();
      L3pt  = TV_L3.Pt();

      L1eta = TV_L1.Eta();
      L2eta = TV_L2.Eta();
      L3eta = TV_L3.Eta();

      L1phi = TV_L1.Phi();
      L2phi = TV_L2.Phi();
      L3phi = TV_L3.Phi();

      consecevent = (float)j;
      
      sameflav = (finalstate_<2) ? 1. : 0;

      i_WWSel0j = (step[14] && (dymva1>0.6 || !sameflav) && njets==0) ? 1 : 0;
      i_WWSel1j = (step[14] && (dymva1>0.3 || !sameflav) && njets==1) ? 1 : 0;
      
      // change the format of the integers -> float
      f_run = (float)run;
      f_lumi = (float)lumi;
      f_hlt = (float)i_hlt;
      f_nVtx = (float)nVtx;
      f_njets = (float)njets;
      f_ncbIDjets = (float)ncbIDjets;

      f_numbtagCSVMmvaIDcentraljets = (float)numbtagCSVMmvaIDcentraljets;
      f_numbtagCSVLmvaIDcentraljets = (float)numbtagCSVLmvaIDcentraljets;
      f_numbtagCSVTmvaIDcentraljets = (float)numbtagCSVTmvaIDcentraljets;

      f_numbtagCSVMmvaIDforwardjets = (float)numbtagCSVMmvaIDforwardjets;
      f_numbtagCSVLmvaIDforwardjets = (float)numbtagCSVLmvaIDforwardjets;
      f_numbtagCSVTmvaIDforwardjets = (float)numbtagCSVTmvaIDforwardjets;

      f_numbtagCSVMcbIDcentraljets = (float)numbtagCSVMcbIDcentraljets;
      f_numbtagCSVLcbIDcentraljets = (float)numbtagCSVLcbIDcentraljets;
      f_numbtagCSVTcbIDcentraljets = (float)numbtagCSVTcbIDcentraljets;

      f_numbtagCSVMcbIDforwardjets = (float)numbtagCSVMcbIDforwardjets;
      f_numbtagCSVLcbIDforwardjets = (float)numbtagCSVLcbIDforwardjets;
      f_numbtagCSVTcbIDforwardjets = (float)numbtagCSVTcbIDforwardjets;

      f_nuncorrjets = (float)nuncorrjets;
      f_zveto = (float)zveto;
      f_bveto_ip = (float)bveto_ip;
      f_bveto_mu = (float)bveto_mu;
      f_bveto_munj = (float)bveto_munj;
      f_bveto = (float)bveto;
      f_dphiveto = (float) ((jetpt1>15 && deltaPhi_LL_JET1<165) || jetpt1<=15);
      f_typeL1 = (float)typeL1;
      f_typeL2 = (float)typeL2;
      f_typeL3 = (float)typeL3;
      f_nSoftMuNoJets = (float)nSoftMuNoJets;
      f_nSoftMu = (float)nSoftMu;
      f_numExtraLep = (float)numExtraLep;
      f_nsoftjet = (float)nsoftjet;
      f_nsoftbjet = (float)nsoftbjet;
      f_ch[0] = (float)ch[0];
      f_ch[1] = (float)ch[1];
      f_ch[2] = (float)ch[2];
      f_finalstate = (float)finalstate_;
      f_processId = (float)processId_;

      bool fitSel = pfMet>20. && eleInvMass>12 && projMet>20. && (njets==0 || njets==1 || (f_dphiveto || !sameflav )  ) && bveto_mu==1 && numExtraLep==0 && bveto_ip==1 && ( !sameflav || (pfMet>45.0 ) );

      if(processId_ < 100 || processId_ >= 1000 ) { // MC
        treeNew->Fill();
        
      } else { // data: apply the trigger 
        if(hlt) {
          treeNew->Fill();
        }
      }
      j++;
      //    }
    }
  
    fileNew->cd();
    treeNew->Write();
    fileNew->Close();

    fileOrig->cd();
    fileOrig->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename_ << endl;
    return;
  }
}

float addWeightsToTreetH::GetProjectedMet(TVector3 met, TVector3 p1, TVector3 p2) {

  float projMET = 0.0;
  float deltaPhi1 = fabs(p1.DeltaPhi(met));
  float deltaPhi2 = fabs(p2.DeltaPhi(met));
  float deltaphi = TMath::Min(deltaPhi1,deltaPhi2);
  if(deltaphi<TMath::Pi()/2.) projMET = met.Mag() * sin(deltaphi);
  else projMET = met.Mag();

  return projMET;
}

float addWeightsToTreetH::calcMT(TVector3 met, TVector3 lepton) {
  return sqrt( 2.*(lepton.Pt())*(met.Pt())*( 1 - cos(met.DeltaPhi(lepton))) );
}

float addWeightsToTreetH::getOfflineEff(float pT, float eta, TH2F *myH) {

  float theEff=-1.;

  int   xBins = myH->GetXaxis()->GetNbins();
  float xMin  = myH->GetXaxis()->GetBinLowEdge(1);
  float xMax  = myH->GetXaxis()->GetBinUpEdge(xBins);
  int   yBins = myH->GetYaxis()->GetNbins();
  float yMin  = myH->GetYaxis()->GetBinLowEdge(1);
  float yMax  = myH->GetYaxis()->GetBinUpEdge(yBins);
  //  int theBin = myH->FindBin(pT, fabs(eta));  // 42X Clara's maps are pt(x);eta(y)
  int theBin = myH->FindBin(fabs(eta), pT);
  //  if (pT>xMin && pT<xMax && fabs(eta)>yMin && fabs(eta)<yMax) {
  if (pT>yMin && pT<yMax && fabs(eta)>xMin && fabs(eta)<xMax) { 
    theEff = myH->GetBinContent(theBin);
  } else {
    theEff = 1.;
    //    cout << "pT = " << pT << ", eta = " << eta << ": pT or eta out of histo bounds. Put SF = 1" << endl;
  }

  // cout << "pT = " << pT << ", eta = " << eta << ", eff = " << theEff << endl;

  return theEff;
}


