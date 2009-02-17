#include <string>

#include <TTree.h>

#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/ElectronBestCandidateSelector.hh"
#include "HiggsAnalysisTools/include/BestLeptonSelectorWjets.hh"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

BestLeptonSelectorWjets::BestLeptonSelectorWjets(TTree *tree)
  : HiggsBase(tree) {

}

BestLeptonSelectorWjets::~BestLeptonSelectorWjets() { }

void BestLeptonSelectorWjets::Loop() {

  if(fChain == 0) return;

  float nTot = 0;
  float nByPt = 0;
  float nBySCenergy = 0;
  float nByTrackerIsolation = 0;
  float nByHcalIsolation = 0;
  float nByElectronIdLH = 0;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // find electron from W+jets mc truth
    int idxGen = indexMcTruthElectronWToENu();
    
    // trigger
    Utils anaUtils;
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

    if(!passedHLT) continue;

    _acceptanceElectrons.clear();
    // now apply acceptance / pt cuts and fill the acceptance electrons
    for(int iele=0; iele<nEle; iele++) {
      if(fabs(etaEle[iele])<2.5 && etEle[iele]>10) {
        _acceptanceElectrons.push_back(iele);
      }
    }

    // select the best electron among the accepted with different criteria
    _bestByPt = -1; _bestBySCenergy = -1; _bestByTrackerIsolation = -1; _bestByHcalIsolation = -1; _bestByElectronIdLH = -1;
    getBestElectronFunny(_acceptanceElectrons);

    if(_acceptanceElectrons.size() >= 2) {


      // look for at least one reco ele matching mc ele (for normalization)
      for(unsigned int i=0; i<_acceptanceElectrons.size(); i++) {
        int iele = _acceptanceElectrons[i];
        if(deltaR_MCmatch(idxGen,iele)<0.3) {
          nTot++;
          break;
        }
      }

      if(deltaR_MCmatch(idxGen,_bestByPt)<0.3) nByPt++;
      if(deltaR_MCmatch(idxGen,_bestBySCenergy)<0.3) nBySCenergy++;
      if(deltaR_MCmatch(idxGen,_bestByTrackerIsolation)<0.3) nByTrackerIsolation++;
      if(deltaR_MCmatch(idxGen,_bestByHcalIsolation)<0.3) nByHcalIsolation++;
      if(deltaR_MCmatch(idxGen,_bestByElectronIdLH)<0.3) nByElectronIdLH++;
    }

  } // end loop on events

  cout << "=== display results ===" << endl;
  cout << " nTot = " << nTot << endl;
  cout << " nByPt eff = " << nByPt << endl;
  cout << " nBySCenergy eff = " << nBySCenergy << endl;
  cout << " nByTrackerIsolation eff = " << nByTrackerIsolation << endl;
  cout << " nByHcalIsolation eff = " << nByHcalIsolation << endl;
  cout << " nByElectronIdLH eff = " << nByElectronIdLH << endl;
  cout << "========================" << endl;

}

int BestLeptonSelectorWjets::indexMcTruthElectronWToENu() {

  int idx=-1;
  for(int imc=0; imc<50; imc++) {
    //    if ( fabs(idMc[imc])==11 && fabs(idMc[mothMc[imc]])==24 ) {
    if ( fabs(idMc[imc])==11 &&  statusMc[imc]==3 ) {
      idx=imc;
      break;
    }
  }

  return idx;

}

float BestLeptonSelectorWjets::deltaR_MCmatch(int iMc, int iEle) {

  TVector3 pMcParticle;
  pMcParticle.SetMagThetaPhi(pMc[iMc], thetaMc[iMc], phiMc[iMc]);

  TVector3 pRecoElectron(pxEle[iEle],pyEle[iEle],pzEle[iEle]);

  float dr = pMcParticle.DeltaR(pRecoElectron);

  return dr;

}

void BestLeptonSelectorWjets::getBestElectronFunny(std::vector<int> goodElectrons) {

  std::vector<ElectronQualityData> electronQual;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    ElectronQualityData quality;
    quality.index = eleIndex;
    quality.ptAtVtx = etEle[eleIndex];
    quality.SCenergy = ecalEle[eleIndex];
    quality.trackerSumPt = eleSumPt04Ele[eleIndex];
    quality.hcalSumEt = eleSumHadEt04Ele[eleIndex];
    quality.electronIdLH = eleLikelihoodEle[eleIndex];
    electronQual.push_back(quality);
  }
  ElectronBestCandidateSelector selector(electronQual);

  _bestByPt = selector.bestByPt().first;
  _bestBySCenergy = selector.bestBySCenergy().first;
  _bestByTrackerIsolation = selector.bestByTrackerIsolation().first;
  _bestByHcalIsolation = selector.bestByEcalIsolation().first;
  _bestByElectronIdLH = selector.bestByElectronIdLH().first;

}
