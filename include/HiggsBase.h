/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  3 10:19:15 2010 by ROOT version 5.22/00d
// from TTree ntp1/ntp1
// found on file: /tmp/crovelli/default_MC_9_1.root
//////////////////////////////////////////////////////////

#ifndef HiggsBase_h
#define HiggsBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <string>

class HiggsBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nl1Technical;
   Int_t           l1Technical[3];   //[nl1Technical]
   Int_t           nl1Global;
   Int_t           l1Global[5];   //[nl1Global]
   Float_t         rhoFastjet;
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Int_t           bunchCrossing;
   Int_t           orbitNumber;
   Int_t           nPU;
   Int_t           nMc;
   Float_t         pMc[101];   //[nMc]
   Float_t         thetaMc[101];   //[nMc]
   Float_t         etaMc[101];   //[nMc]
   Float_t         phiMc[101];   //[nMc]
   Float_t         energyMc[101];   //[nMc]
   Int_t           idMc[101];   //[nMc]
   Int_t           mothMc[101];   //[nMc]
   Int_t           statusMc[101];   //[nMc]
   Int_t           nTrg;
   Int_t           firedTrg[100];   //[nTrg]
   Int_t           nHLT;
   Int_t           indexHLT[1000];   //[nHLT]
   std::vector<std::string>  *nameHLT;
   Int_t           nEle;
   Int_t           chargeEle[100];   //[nEle]
   Float_t         energyEle[100];   //[nEle]
   Float_t         thetaEle[100];   //[nEle]
   Float_t         etaEle[100];   //[nEle]
   Float_t         phiEle[100];   //[nEle]
   Float_t         pxEle[100];   //[nEle]
   Float_t         pyEle[100];   //[nEle]
   Float_t         pzEle[100];   //[nEle]
   Float_t         vertexXEle[100];   //[nEle]
   Float_t         vertexYEle[100];   //[nEle]
   Float_t         vertexZEle[100];   //[nEle]
   Int_t           fiducialFlagsEle[100];   //[nEle]
   Int_t           recoFlagsEle[100];   //[nEle]
   Int_t           energyCorrectionsEle[100];   //[nEle]
   Int_t           superClusterIndexEle[100];   //[nEle]
   Int_t           PFsuperClusterIndexEle[100];   //[nEle]
   Int_t           trackIndexEle[100];   //[nEle]
   Int_t           gsfTrackIndexEle[100];   //[nEle]
   Float_t         convDistEle[100];   //[nEle]
   Float_t         convDcotEle[100];   //[nEle]
   Float_t         convRadiusEle[100];   //[nEle]
   Int_t           convTrackIndexEle[100];   //[nEle]
   Int_t           scPixChargeEle[100];   //[nEle]
   Bool_t          hasMatchedConversionEle[100];   //[nEle]
   Int_t           classificationEle[100];   //[nEle]
   Int_t           standardClassificationEle[100];   //[nEle]
   Float_t         fbremEle[100];   //[nEle]
   Int_t           nbremsEle[100];   //[nEle]
   Float_t         hOverEEle[100];   //[nEle]
   Float_t         eSuperClusterOverPEle[100];   //[nEle]
   Float_t         eSeedOverPoutEle[100];   //[nEle]
   Float_t         deltaEtaAtVtxEle[100];   //[nEle]
   Float_t         deltaPhiAtVtxEle[100];   //[nEle]
   Float_t         deltaEtaAtCaloEle[100];   //[nEle]
   Float_t         deltaPhiAtCaloEle[100];   //[nEle]
   Float_t         dr03TkSumPtEle[100];   //[nEle]
   Float_t         dr03EcalRecHitSumEtEle[100];   //[nEle]
   Float_t         dr03HcalTowerSumEtEle[100];   //[nEle]
   Float_t         dr04TkSumPtEle[100];   //[nEle]
   Float_t         dr04EcalRecHitSumEtEle[100];   //[nEle]
   Float_t         dr04HcalTowerSumEtEle[100];   //[nEle]
   Float_t         scBasedEcalSum03Ele[100];   //[nEle]
   Float_t         scBasedEcalSum04Ele[100];   //[nEle]
   Float_t         dr03HcalTowerSumEtFullConeEle[5];   //[nEle]
   Float_t         dr04HcalTowerSumEtFullConeEle[5];   //[nEle]
   Float_t         eleIdLikelihoodEle[100];   //[nEle]
   Float_t         pflowMVAEle[100];   //[nEle]
   Float_t         pfChargedIsoEle[100];   //[nEle]
   Float_t         pfNeutralIsoEle[100];   //[nEle]
   Float_t         pfPhotonIsoEle[100];   //[nEle]
   Float_t         pfGenericChargedIsoEle[100];   //[nEle]
   Float_t         pfGenericNeutralIsoEle[100];   //[nEle]
   Float_t         pfGenericPhotonIsoEle[100];   //[nEle]
   Float_t         pfGenericNoOverChargedIsoEle[100];   //[nEle]
   Float_t         pfGenericNoOverNeutralIsoEle[100];   //[nEle]
   Float_t         pfGenericNoOverPhotonIsoEle[100];   //[nEle]
   Float_t         pfCombinedIsoEle[100];   //[nEle]
   Int_t           nPho;
   Int_t           chargePho[100];   //[nPho]
   Float_t         energyPho[100];   //[nPho]
   Float_t         thetaPho[100];   //[nPho]
   Float_t         etaPho[100];   //[nPho]
   Float_t         phiPho[100];   //[nPho]
   Float_t         pxPho[100];   //[nPho]
   Float_t         pyPho[100];   //[nPho]
   Float_t         pzPho[100];   //[nPho]
   Float_t         vertexXPho[100];   //[nPho]
   Float_t         vertexYPho[100];   //[nPho]
   Float_t         vertexZPho[100];   //[nPho]
   Int_t           fiducialFlagsPho[100];   //[nPho]
   Int_t           recoFlagsPho[100];   //[nPho]
   Int_t           superClusterIndexPho[100];   //[nPho]
   Int_t           PFsuperClusterIndexPho[100];   //[nPho]
   Float_t         hOverEPho[100];   //[nPho]
   Float_t         dr03TkSumPtPho[100];   //[nPho]
   Float_t         dr03HollowTkSumPtPho[100];   //[nPho]
   Float_t         dr03EcalRecHitSumEtPho[100];   //[nPho]
   Float_t         dr03HcalTowerSumEtPho[100];   //[nPho]
   Float_t         dr04TkSumPtPho[100];   //[nPho]
   Float_t         dr04HollowTkSumPtPho[100];   //[nPho]
   Float_t         dr04EcalRecHitSumEtPho[100];   //[nPho]
   Float_t         dr04HcalTowerSumEtPho[100];   //[nPho]
   Float_t         chargedHadronIsoPho[100];   //[nPho]
   Float_t         neutralHadronIsoPho[100];   //[nPho]
   Float_t         photonIsoPho[100];   //[nPho]
   Int_t           hasPixelSeedPho[100];   //[nPho]
   Bool_t          hasMatchedConversionPho[100];   //[nPho]
   Int_t           nSC;
   Int_t           nBCSC[100];   //[nSC]
   Int_t           nCrystalsSC[100];   //[nSC]
   Float_t         rawEnergySC[100];   //[nSC]
   Float_t         energySC[100];   //[nSC]
   Float_t         etaSC[100];   //[nSC]
   Float_t         thetaSC[100];   //[nSC]
   Float_t         phiSC[100];   //[nSC]
   Float_t         phiWidthSC[100];   //[nSC]
   Float_t         etaWidthSC[100];   //[nSC]
   Float_t         e3x3SC[100];   //[nSC]
   Float_t         e5x5SC[100];   //[nSC]
   Float_t         eMaxSC[100];   //[nSC]
   Float_t         e2x2SC[100];   //[nSC]
   Float_t         e2ndSC[100];   //[nSC]
   Float_t         e1x5SC[100];   //[nSC]
   Float_t         e2x5MaxSC[100];   //[nSC]
   Float_t         e4SwissCrossSC[100];   //[nSC]
   Float_t         covIEtaIEtaSC[100];   //[nSC]
   Float_t         covIEtaIPhiSC[100];   //[nSC]
   Float_t         covIPhiIPhiSC[100];   //[nSC]
   Float_t         hOverESC[100];   //[nSC]
   Int_t           recoFlagSC[100];   //[nSC]
   Float_t         timeSC[100];   //[nSC]
   Float_t         chi2SC[100];   //[nSC]
   Float_t         seedEnergySC[100];   //[nSC]
   Int_t           nPFSC;
   Int_t           nBCPFSC[100];   //[nPFSC]
   Int_t           nCrystalsPFSC[100];   //[nPFSC]
   Float_t         rawEnergyPFSC[100];   //[nPFSC]
   Float_t         energyPFSC[100];   //[nPFSC]
   Float_t         etaPFSC[100];   //[nPFSC]
   Float_t         thetaPFSC[100];   //[nPFSC]
   Float_t         phiPFSC[100];   //[nPFSC]
   Float_t         phiWidthPFSC[100];   //[nPFSC]
   Float_t         etaWidthPFSC[100];   //[nPFSC]
   Float_t         e3x3PFSC[100];   //[nPFSC]
   Float_t         e5x5PFSC[100];   //[nPFSC]
   Float_t         eMaxPFSC[100];   //[nPFSC]
   Float_t         e2x2PFSC[100];   //[nPFSC]
   Float_t         e2ndPFSC[100];   //[nPFSC]
   Float_t         e1x5PFSC[100];   //[nPFSC]
   Float_t         e2x5MaxPFSC[100];   //[nPFSC]
   Float_t         e4SwissCrossPFSC[100];   //[nPFSC]
   Float_t         covIEtaIEtaPFSC[100];   //[nPFSC]
   Float_t         covIEtaIPhiPFSC[100];   //[nPFSC]
   Float_t         covIPhiIPhiPFSC[100];   //[nPFSC]
   Float_t         hOverEPFSC[100];   //[nPFSC]
   Int_t           recoFlagPFSC[100];   //[nPFSC]
   Float_t         timePFSC[100];   //[nPFSC]
   Float_t         chi2PFSC[100];   //[nPFSC]
   Float_t         seedEnergyPFSC[100];   //[nPFSC]
   Int_t           nTrack;
   Float_t         pxTrack[2000];   //[nTrack]
   Float_t         pyTrack[2000];   //[nTrack]
   Float_t         pzTrack[2000];   //[nTrack]
   Int_t           vtxIndexTrack[2000];   //[nTrack]
   Float_t         vtxWeightTrack[2000];   //[nTrack]
   Float_t         chargeTrack[2000];   //[nTrack]
   Float_t         ptErrorTrack[2000];   //[nTrack]
   Float_t         trackValidHitsTrack[2000];   //[nTrack]
   Float_t         trackLostHitsTrack[2000];   //[nTrack]
   Float_t         trackNormalizedChi2Track[2000];   //[nTrack]
   Int_t           qualityMaskTrack[2000];   //[nTrack]
   Float_t         impactPar3DTrack[2000];   //[nTrack]
   Float_t         impactPar3DErrorTrack[2000];   //[nTrack]
   Float_t         transvImpactParTrack[2000];   //[nTrack]
   Float_t         transvImpactParErrorTrack[2000];   //[nTrack]
   Float_t         impactPar3DBiasedTrack[2000];   //[nTrack]
   Float_t         impactPar3DBiasedErrorTrack[2000];   //[nTrack]
   Float_t         transvImpactParBiasedTrack[2000];   //[nTrack]
   Float_t         transvImpactParBiasedErrorTrack[2000];   //[nTrack]
   Float_t         trackVxTrack[2000];   //[nTrack]
   Float_t         trackVyTrack[2000];   //[nTrack]
   Float_t         trackVzTrack[2000];   //[nTrack]
   Float_t         recHitsSizeTrack[2000];   //[nTrack]
   Int_t           pixelHitsTrack[2000];   //[nTrack]
   Int_t           expInnerLayersTrack[2000];   //[nTrack]
   Int_t           numberOfValidPixelBarrelHitsTrack[2000];   //[nTrack]
   Int_t           numberOfValidPixelEndcapHitsTrack[2000];   //[nTrack]
   Int_t           numberOfValidStripTIBHitsTrack[2000];   //[nTrack]
   Int_t           numberOfValidStripTIDHitsTrack[2000];   //[nTrack]
   Int_t           numberOfValidStripTOBHitsTrack[2000];   //[nTrack]
   Int_t           numberOfValidStripTECHitsTrack[2000];   //[nTrack]
   Int_t           nGsfTrack;
   Float_t         pxGsfTrack[500];   //[nGsfTrack]
   Float_t         pyGsfTrack[500];   //[nGsfTrack]
   Float_t         pzGsfTrack[500];   //[nGsfTrack]
   Int_t           vtxIndexGsfTrack[500];   //[nGsfTrack]
   Float_t         vtxWeightGsfTrack[500];   //[nGsfTrack]
   Float_t         chargeGsfTrack[500];   //[nGsfTrack]
   Float_t         ptErrorGsfTrack[500];   //[nGsfTrack]
   Float_t         trackValidHitsGsfTrack[500];   //[nGsfTrack]
   Float_t         trackLostHitsGsfTrack[500];   //[nGsfTrack]
   Float_t         trackNormalizedChi2GsfTrack[500];   //[nGsfTrack]
   Int_t           qualityMaskGsfTrack[500];   //[nGsfTrack]
   Float_t         impactPar3DGsfTrack[500];   //[nGsfTrack]
   Float_t         impactPar3DErrorGsfTrack[500];   //[nGsfTrack]
   Float_t         transvImpactParGsfTrack[500];   //[nGsfTrack]
   Float_t         transvImpactParErrorGsfTrack[500];   //[nGsfTrack]
   Float_t         impactPar3DBiasedGsfTrack[500];   //[nGsfTrack]
   Float_t         impactPar3DBiasedErrorGsfTrack[500];   //[nGsfTrack]
   Float_t         transvImpactParBiasedGsfTrack[500];   //[nGsfTrack]
   Float_t         transvImpactParBiasedErrorGsfTrack[500];   //[nGsfTrack]
   Float_t         trackVxGsfTrack[500];   //[nGsfTrack]
   Float_t         trackVyGsfTrack[500];   //[nGsfTrack]
   Float_t         trackVzGsfTrack[500];   //[nGsfTrack]
   Float_t         recHitsSizeGsfTrack[500];   //[nGsfTrack]
   Int_t           pixelHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           expInnerLayersGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidPixelBarrelHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidPixelEndcapHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidStripTIBHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidStripTIDHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidStripTOBHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           numberOfValidStripTECHitsGsfTrack[500];   //[nGsfTrack]
   Int_t           chargeModeGsfTrack[500];   //[nGsfTrack]
   Float_t         pxModeGsfTrack[500];   //[nGsfTrack]
   Float_t         pyModeGsfTrack[500];   //[nGsfTrack]
   Float_t         pzModeGsfTrack[500];   //[nGsfTrack]
   Int_t           nGlobalMuonTrack;
   Float_t         pxGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         pyGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         pzGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           vtxIndexGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         vtxWeightGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         chargeGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         ptErrorGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackValidHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackLostHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackNormalizedChi2GlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           qualityMaskGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         impactPar3DGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         impactPar3DErrorGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         transvImpactParGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         transvImpactParErrorGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         impactPar3DBiasedGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         impactPar3DBiasedErrorGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         transvImpactParBiasedGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         transvImpactParBiasedErrorGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackVxGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackVyGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         trackVzGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Float_t         recHitsSizeGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           pixelHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           expInnerLayersGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIBHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIDHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTOBHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTECHitsGlobalMuonTrack[100];   //[nGlobalMuonTrack]
   Int_t           nSTAMuonTrack;
   Float_t         pxSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         pySTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         pzSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         chargeSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         ptErrorSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackValidHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackLostHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackNormalizedChi2STAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           qualityMaskSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackVxSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackVySTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         trackVzSTAMuonTrack[100];   //[nSTAMuonTrack]
   Float_t         recHitsSizeSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           pixelHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           expInnerLayersSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIBHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIDHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTOBHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTECHitsSTAMuonTrack[100];   //[nSTAMuonTrack]
   Int_t           nPV;
   Float_t         PVxPV[100];   //[nPV]
   Float_t         PVyPV[100];   //[nPV]
   Float_t         PVzPV[100];   //[nPV]
   Float_t         PVErrxPV[100];   //[nPV]
   Float_t         PVErryPV[100];   //[nPV]
   Float_t         PVErrzPV[100];   //[nPV]
   Float_t         SumPtPV[100];   //[nPV]
   Float_t         ndofPV[100];   //[nPV]
   Float_t         chi2PV[100];   //[nPV]
   Float_t         pxChMetPV[100];   //[nPV]
   Float_t         pyChMetPV[100];   //[nPV]
   Float_t         pzChMetPV[100];   //[nPV]
   Int_t           isFakePV[100];  // [nPV]
   Int_t           nMuon;
   Int_t           chargeMuon[100];   //[nMuon]
   Float_t         energyMuon[100];   //[nMuon]
   Float_t         thetaMuon[100];   //[nMuon]
   Float_t         etaMuon[100];   //[nMuon]
   Float_t         phiMuon[100];   //[nMuon]
   Float_t         pxMuon[100];   //[nMuon]
   Float_t         pyMuon[100];   //[nMuon]
   Float_t         pzMuon[100];   //[nMuon]
   Float_t         vertexXMuon[100];   //[nMuon]
   Float_t         vertexYMuon[100];   //[nMuon]
   Float_t         vertexZMuon[100];   //[nMuon]
   Int_t           trackIndexMuon[100];   //[nMuon]
   Int_t           standAloneTrackIndexMuon[100];   //[nMuon]
   Int_t           combinedTrackIndexMuon[100];   //[nMuon]
   Int_t           muonIdMuon[100];   //[nMuon]
   Int_t           typeMuon[100];   //[nMuon]
   Int_t           numberOfMatchesMuon[100];   //[nMuon]
   Float_t         sumPt03Muon[100];   //[nMuon]
   Float_t         emEt03Muon[100];   //[nMuon]
   Float_t         hadEt03Muon[100];   //[nMuon]
   Float_t         hoEt03Muon[100];   //[nMuon]
   Float_t         nTrk03Muon[100];   //[nMuon]
   Float_t         nJets03Muon[100];   //[nMuon]
   Float_t         sumPt05Muon[100];   //[nMuon]
   Float_t         emEt05Muon[100];   //[nMuon]
   Float_t         hadEt05Muon[100];   //[nMuon]
   Float_t         hoEt05Muon[100];   //[nMuon]
   Float_t         nTrk05Muon[100];   //[nMuon]
   Float_t         nJets05Muon[100];   //[nMuon]
   Float_t         pfChargedIsoMuon[100];   //[nMuon]
   Float_t         pfNeutralIsoMuon[100];   //[nMuon]
   Float_t         pfPhotonIsoMuon[100];   //[nMuon]
   Float_t         pfGenericChargedIsoMuon[100];   //[nMuon]
   Float_t         pfGenericNeutralIsoMuon[100];   //[nMuon]
   Float_t         pfGenericPhotonIsoMuon[100];   //[nMuon]
   Float_t         pfGenericNoOverChargedIsoMuon[100];   //[nMuon]
   Float_t         pfGenericNoOverNeutralIsoMuon[100];   //[nMuon]
   Float_t         pfGenericNoOverPhotonIsoMuon[100];   //[nMuon]
   Float_t         pfCombinedIsoMuon[100];   //[nMuon]
   Float_t         EcalExpDepoMuon[100];   //[nMuon]
   Float_t         HcalExpDepoMuon[100];   //[nMuon]
   Float_t         HoExpDepoMuon[100];   //[nMuon]
   Float_t         emS9Muon[100];   //[nMuon]
   Float_t         hadS9Muon[100];   //[nMuon]
   Float_t         hoS9Muon[100];   //[nMuon]
   Float_t         CaloCompMuon[100];   //[nMuon]
   Int_t           nReducedPFCand;
   Int_t           chargeReducedPFCand[100];   //[nReducedPFCand]
   Float_t         energyReducedPFCand[100];   //[nReducedPFCand]
   Float_t         thetaReducedPFCand[100];   //[nReducedPFCand]
   Float_t         etaReducedPFCand[100];   //[nReducedPFCand]
   Float_t         phiReducedPFCand[100];   //[nReducedPFCand]
   Float_t         pxReducedPFCand[100];   //[nReducedPFCand]
   Float_t         pyReducedPFCand[100];   //[nReducedPFCand]
   Float_t         pzReducedPFCand[100];   //[nReducedPFCand]
   Float_t         vertexXReducedPFCand[100];   //[nReducedPFCand]
   Float_t         vertexYReducedPFCand[100];   //[nReducedPFCand]
   Float_t         vertexZReducedPFCand[100];   //[nReducedPFCand]
   Int_t           nMet;
   Int_t           chargeMet[1];   //[nMet]
   Float_t         energyMet[1];   //[nMet]
   Float_t         thetaMet[1];   //[nMet]
   Float_t         etaMet[1];   //[nMet]
   Float_t         phiMet[1];   //[nMet]
   Float_t         pxMet[1];   //[nMet]
   Float_t         pyMet[1];   //[nMet]
   Float_t         pzMet[1];   //[nMet]
   Float_t         vertexXMet[1];   //[nMet]
   Float_t         vertexYMet[1];   //[nMet]
   Float_t         vertexZMet[1];   //[nMet]
   Int_t           nTCMet;
   Int_t           chargeTCMet[1];   //[nTCMet]
   Float_t         energyTCMet[1];   //[nTCMet]
   Float_t         thetaTCMet[1];   //[nTCMet]
   Float_t         etaTCMet[1];   //[nTCMet]
   Float_t         phiTCMet[1];   //[nTCMet]
   Float_t         pxTCMet[1];   //[nTCMet]
   Float_t         pyTCMet[1];   //[nTCMet]
   Float_t         pzTCMet[1];   //[nTCMet]
   Float_t         vertexXTCMet[1];   //[nTCMet]
   Float_t         vertexYTCMet[1];   //[nTCMet]
   Float_t         vertexZTCMet[1];   //[nTCMet]
   Int_t           nPFMet;
   Int_t           chargePFMet[1];   //[nPFMet]
   Float_t         energyPFMet[1];   //[nPFMet]
   Float_t         thetaPFMet[1];   //[nPFMet]
   Float_t         etaPFMet[1];   //[nPFMet]
   Float_t         phiPFMet[1];   //[nPFMet]
   Float_t         pxPFMet[1];   //[nPFMet]
   Float_t         pyPFMet[1];   //[nPFMet]
   Float_t         pzPFMet[1];   //[nPFMet]
   Float_t         vertexXPFMet[1];   //[nPFMet]
   Float_t         vertexYPFMet[1];   //[nPFMet]
   Float_t         vertexZPFMet[1];   //[nPFMet]
   Float_t         sumEtPFMet[1];   //[nPFMet]
   Float_t         mEtSigPFMet[1];   //[nPFMet]
   Float_t         significancePFMet[1];   //[nPFMet]
   Int_t           nPFChMet;
   Int_t           chargePFChMet[1];   //[nPFChMet]
   Float_t         energyPFChMet[1];   //[nPFChMet]
   Float_t         thetaPFChMet[1];   //[nPFChMet]
   Float_t         etaPFChMet[1];   //[nPFChMet]
   Float_t         phiPFChMet[1];   //[nPFChMet]
   Float_t         pxPFChMet[1];   //[nPFChMet]
   Float_t         pyPFChMet[1];   //[nPFChMet]
   Float_t         pzPFChMet[1];   //[nPFChMet]
   Float_t         vertexXPFChMet[1];   //[nPFChMet]
   Float_t         vertexYPFChMet[1];   //[nPFChMet]
   Float_t         vertexZPFChMet[1];   //[nPFChMet]
   Float_t         sumEtPFChMet[1];   //[nPFChMet]
   Float_t         mEtSigPFChMet[1];   //[nPFChMet]
   Float_t         significancePFChMet[1];   //[nPFChMet]
   Int_t           nGenMet;
   Int_t           chargeGenMet[1];   //[nGenMet]
   Float_t         energyGenMet[1];   //[nGenMet]
   Float_t         thetaGenMet[1];   //[nGenMet]
   Float_t         etaGenMet[1];   //[nGenMet]
   Float_t         phiGenMet[1];   //[nGenMet]
   Float_t         pxGenMet[1];   //[nGenMet]
   Float_t         pyGenMet[1];   //[nGenMet]
   Float_t         pzGenMet[1];   //[nGenMet]
   Float_t         vertexXGenMet[1];   //[nGenMet]
   Float_t         vertexYGenMet[1];   //[nGenMet]
   Float_t         vertexZGenMet[1];   //[nGenMet]
   Int_t           nAK5Jet;
   Int_t           chargeAK5Jet[200];   //[nAK5Jet]
   Float_t         energyAK5Jet[200];   //[nAK5Jet]
   Float_t         thetaAK5Jet[200];   //[nAK5Jet]
   Float_t         etaAK5Jet[200];   //[nAK5Jet]
   Float_t         phiAK5Jet[200];   //[nAK5Jet]
   Float_t         pxAK5Jet[200];   //[nAK5Jet]
   Float_t         pyAK5Jet[200];   //[nAK5Jet]
   Float_t         pzAK5Jet[200];   //[nAK5Jet]
   Float_t         vertexXAK5Jet[200];   //[nAK5Jet]
   Float_t         vertexYAK5Jet[200];   //[nAK5Jet]
   Float_t         vertexZAK5Jet[200];   //[nAK5Jet]
   Float_t         emFracAK5Jet[200];   //[nAK5Jet]
   Float_t         hadFracAK5Jet[200];   //[nAK5Jet]
   Int_t           IdAK5Jet[200];   //[nAK5Jet]
   Int_t           nHitAK5Jet[200];   //[nAK5Jet]
   Int_t           nHit90AK5Jet[200];   //[nAK5Jet]
   Float_t         fHPDAK5Jet[200];   //[nAK5Jet]
   Float_t         covEtaEtaAK5Jet[200];   //[nAK5Jet]
   Float_t         covPhiPhiAK5Jet[200];   //[nAK5Jet]
   Float_t         fLSAK5Jet[200];   //[nAK5Jet]
   Float_t         fOOTAK5Jet[200];   //[nAK5Jet]
   Float_t         combinedSecondaryVertexBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         combinedSecondaryVertexMVABJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         jetBProbabilityBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         jetProbabilityBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         simpleSecondaryVertexHighEffBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         simpleSecondaryVertexHighPurBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         softMuonBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         softMuonByIP3dBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         softMuonByPtBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         softElectronByIP3dBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         softElectronByPtBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         trackCountingHighPurBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         trackCountingHighEffBJetTagsAK5Jet[200];   //[nAK5Jet]
   Float_t         uncorrEnergyAK5Jet[200];   //[nAK5Jet]
   Int_t           nAK5PFPUcorrJet;
   Int_t           chargeAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         energyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         thetaAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         etaAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         phiAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         pxAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         pyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         pzAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         vertexXAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         vertexYAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         vertexZAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         chargedHadronEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         neutralHadronEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         photonEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         electronEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         muonEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         HFHadronEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         HFEMEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           chargedHadronMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           neutralHadronMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           photonMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           electronMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           muonMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           HFHadronMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           HFEMMultiplicityAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         chargedEmEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         neutralEmEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         jetBProbabilityBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         jetProbabilityBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softMuonBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softMuonByIP3dBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softMuonByPtBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softElectronBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softElectronByIP3dBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         softElectronByPtBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         trackCountingHighPurBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         trackCountingHighEffBJetTagsAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         uncorrEnergyAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         rmsCandAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Float_t         ptDAK5PFPUcorrJet[200];   //[nAK5PFPUcorrJet]
   Int_t           nAK5GenJet;
   Int_t           chargeAK5GenJet[200];   //[nAK5GenJet]
   Float_t         energyAK5GenJet[200];   //[nAK5GenJet]
   Float_t         thetaAK5GenJet[200];   //[nAK5GenJet]
   Float_t         etaAK5GenJet[200];   //[nAK5GenJet]
   Float_t         phiAK5GenJet[200];   //[nAK5GenJet]
   Float_t         pxAK5GenJet[200];   //[nAK5GenJet]
   Float_t         pyAK5GenJet[200];   //[nAK5GenJet]
   Float_t         pzAK5GenJet[200];   //[nAK5GenJet]
   Float_t         vertexXAK5GenJet[200];   //[nAK5GenJet]
   Float_t         vertexYAK5GenJet[200];   //[nAK5GenJet]
   Float_t         vertexZAK5GenJet[200];   //[nAK5GenJet]
   Double_t        genPtHat;
   Double_t        genProcessId;
   Double_t        genWeight;
   Double_t        genAlphaQCD;
   Double_t        genAlphaQED;
   Int_t           nTriggerObsPassing;
   Int_t           sizePassing[1000];
   Int_t           indexPassingPerPath[1000];
   Int_t           indexPassing[1000];
   Int_t           nTriggerObs;
   Float_t         triggerObsPt[4000]; 
   Float_t         triggerObsEta[4000]; 
   Float_t         triggerObsPhi[4000]; 
   Float_t         triggerObsMass[4000]; 
  
   // List of branches
   TBranch        *b_nl1Technical;   //!
   TBranch        *b_l1Technical;   //!
   TBranch        *b_nl1Global;   //!
   TBranch        *b_l1Global;   //!
   TBranch        *b_rhoFastjet;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_energyMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_statusMc;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_firedTrg;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_nameHLT;   //!
   TBranch        *b_indexHLT;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_thetaEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_pxEle;   //!
   TBranch        *b_pyEle;   //!
   TBranch        *b_pzEle;   //!
   TBranch        *b_vertexXEle;   //!
   TBranch        *b_vertexYEle;   //!
   TBranch        *b_vertexZEle;   //!
   TBranch        *b_fiducialFlagsEle;   //!
   TBranch        *b_recoFlagsEle;   //!
   TBranch        *b_energyCorrectionsEle;   //!
   TBranch        *b_superClusterIndexEle;   //!
   TBranch        *b_PFsuperClusterIndexEle;   //!
   TBranch        *b_trackIndexEle;   //!
   TBranch        *b_gsfTrackIndexEle;   //!
   TBranch        *b_convDistEle;   //!
   TBranch        *b_convDcotEle;   //!
   TBranch        *b_convRadiusEle;   //!
   TBranch        *b_convTrackIndexEle;   //!
   TBranch        *b_scPixChargeEle;   //!
   TBranch        *b_hasMatchedConversionEle;   //!
   TBranch        *b_classificationEle;   //!
   TBranch        *b_standardClassificationEle;   //!
   TBranch        *b_fbremEle;   //!
   TBranch        *b_nbremsEle;   //!
   TBranch        *b_hOverEEle;   //!
   TBranch        *b_eSuperClusterOverPEle;   //!
   TBranch        *b_eSeedOverPoutEle;   //!
   TBranch        *b_deltaEtaAtVtxEle;   //!
   TBranch        *b_deltaPhiAtVtxEle;   //!
   TBranch        *b_deltaEtaAtCaloEle;   //!
   TBranch        *b_deltaPhiAtCaloEle;   //!
   TBranch        *b_dr03TkSumPtEle;   //!
   TBranch        *b_dr03EcalRecHitSumEtEle;   //!
   TBranch        *b_dr03HcalTowerSumEtEle;   //!
   TBranch        *b_dr04TkSumPtEle;   //!
   TBranch        *b_dr04EcalRecHitSumEtEle;   //!
   TBranch        *b_dr04HcalTowerSumEtEle;   //!
   TBranch        *b_scBasedEcalSum03Ele;   //!
   TBranch        *b_scBasedEcalSum04Ele;   //!
   TBranch        *b_dr03HcalTowerSumEtFullConeEle;   //!
   TBranch        *b_dr04HcalTowerSumEtFullConeEle;   //!
   TBranch        *b_eleIdLikelihoodEle;   //!
   TBranch        *b_pflowMVAEle;   //!
   TBranch        *b_pfChargedIsoEle;   //!
   TBranch        *b_pfNeutralIsoEle;   //!
   TBranch        *b_pfPhotonIsoEle;   //!
   TBranch        *b_pfGenericChargedIsoEle;   //!
   TBranch        *b_pfGenericNeutralIsoEle;   //!
   TBranch        *b_pfGenericPhotonIsoEle;   //!
   TBranch        *b_pfGenericNoOverChargedIsoEle;   //!
   TBranch        *b_pfGenericNoOverNeutralIsoEle;   //!
   TBranch        *b_pfGenericNoOverPhotonIsoEle;   //!
   TBranch        *b_pfCombinedIsoEle;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_chargePho;   //!
   TBranch        *b_energyPho;   //!
   TBranch        *b_thetaPho;   //!
   TBranch        *b_etaPho;   //!
   TBranch        *b_phiPho;   //!
   TBranch        *b_pxPho;   //!
   TBranch        *b_pyPho;   //!
   TBranch        *b_pzPho;   //!
   TBranch        *b_vertexXPho;   //!
   TBranch        *b_vertexYPho;   //!
   TBranch        *b_vertexZPho;   //!
   TBranch        *b_fiducialFlagsPho;   //!
   TBranch        *b_recoFlagsPho;   //!
   TBranch        *b_superClusterIndexPho;   //!
   TBranch        *b_PFsuperClusterIndexPho;   //!
   TBranch        *b_hOverEPho;   //!
   TBranch        *b_dr03TkSumPtPho;   //!
   TBranch        *b_dr03HollowTkSumPtPho;   //!
   TBranch        *b_dr03EcalRecHitSumEtPho;   //!
   TBranch        *b_dr03HcalTowerSumEtPho;   //!
   TBranch        *b_dr04TkSumPtPho;   //!
   TBranch        *b_dr04HollowTkSumPtPho;   //!
   TBranch        *b_dr04EcalRecHitSumEtPho;   //!
   TBranch        *b_dr04HcalTowerSumEtPho;   //!
   TBranch        *b_chargedHadronIsoPho;   //!
   TBranch        *b_neutralHadronIsoPho;   //!
   TBranch        *b_photonIsoPho;   //!
   TBranch        *b_hasPixelSeedPho;   //!
   TBranch        *b_hasMatchedConversionPho;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_nBCSC;   //!
   TBranch        *b_nCrystalsSC;   //!
   TBranch        *b_rawEnergySC;   //!
   TBranch        *b_energySC;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_thetaSC;   //!
   TBranch        *b_phiSC;   //!
   TBranch        *b_phiWidthSC;   //!
   TBranch        *b_etaWidthSC;   //!
   TBranch        *b_e3x3SC;   //!
   TBranch        *b_e5x5SC;   //!
   TBranch        *b_eMaxSC;   //!
   TBranch        *b_e2x2SC;   //!
   TBranch        *b_e2ndSC;   //!
   TBranch        *b_e1x5SC;   //!
   TBranch        *b_e2x5MaxSC;   //!
   TBranch        *b_e4SwissCrossSC;   //!
   TBranch        *b_covIEtaIEtaSC;   //!
   TBranch        *b_covIEtaIPhiSC;   //!
   TBranch        *b_covIPhiIPhiSC;   //!
   TBranch        *b_hOverESC;   //!
   TBranch        *b_recoFlagSC;   //!
   TBranch        *b_timeSC;   //!
   TBranch        *b_chi2SC;   //!
   TBranch        *b_seedEnergySC;   //!
   TBranch        *b_nPFSC;   //!
   TBranch        *b_nBCPFSC;   //!
   TBranch        *b_nCrystalsPFSC;   //!
   TBranch        *b_rawEnergyPFSC;   //!
   TBranch        *b_energyPFSC;   //!
   TBranch        *b_etaPFSC;   //!
   TBranch        *b_thetaPFSC;   //!
   TBranch        *b_phiPFSC;   //!
   TBranch        *b_phiWidthPFSC;   //!
   TBranch        *b_etaWidthPFSC;   //!
   TBranch        *b_e3x3PFSC;   //!
   TBranch        *b_e5x5PFSC;   //!
   TBranch        *b_eMaxPFSC;   //!
   TBranch        *b_e2x2PFSC;   //!
   TBranch        *b_e2ndPFSC;   //!
   TBranch        *b_e1x5PFSC;   //!
   TBranch        *b_e2x5MaxPFSC;   //!
   TBranch        *b_e4SwissCrossPFSC;   //!
   TBranch        *b_covIEtaIEtaPFSC;   //!
   TBranch        *b_covIEtaIPhiPFSC;   //!
   TBranch        *b_covIPhiIPhiPFSC;   //!
   TBranch        *b_hOverEPFSC;   //!
   TBranch        *b_recoFlagPFSC;   //!
   TBranch        *b_timePFSC;   //!
   TBranch        *b_chi2PFSC;   //!
   TBranch        *b_seedEnergyPFSC;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_pxTrack;   //!
   TBranch        *b_pyTrack;   //!
   TBranch        *b_pzTrack;   //!
   TBranch        *b_vtxIndexTrack;   //!
   TBranch        *b_vtxWeightTrack;   //!
   TBranch        *b_chargeTrack;   //!
   TBranch        *b_ptErrorTrack;   //!
   TBranch        *b_trackValidHitsTrack;   //!
   TBranch        *b_trackLostHitsTrack;   //!
   TBranch        *b_trackNormalizedChi2Track;   //!
   TBranch        *b_qualityMaskTrack;   //!
   TBranch        *b_impactPar3DTrack;   //!
   TBranch        *b_impactPar3DErrorTrack;   //!
   TBranch        *b_transvImpactParTrack;   //!
   TBranch        *b_transvImpactParErrorTrack;   //!
   TBranch        *b_impactPar3DBiasedTrack;   //!
   TBranch        *b_impactPar3DBiasedErrorTrack;   //!
   TBranch        *b_transvImpactParBiasedTrack;   //!
   TBranch        *b_transvImpactParBiasedErrorTrack;   //!
   TBranch        *b_trackVxTrack;   //!
   TBranch        *b_trackVyTrack;   //!
   TBranch        *b_trackVzTrack;   //!
   TBranch        *b_pixelHitsTrack;   //!
   TBranch        *b_expInnerLayersTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsTrack;   //!
   TBranch        *b_nGsfTrack;   //!
   TBranch        *b_pxGsfTrack;   //!
   TBranch        *b_pyGsfTrack;   //!
   TBranch        *b_pzGsfTrack;   //!
   TBranch        *b_vtxIndexGsfTrack;   //!
   TBranch        *b_vtxWeightGsfTrack;   //!
   TBranch        *b_chargeGsfTrack;   //!
   TBranch        *b_ptErrorGsfTrack;   //!
   TBranch        *b_trackValidHitsGsfTrack;   //!
   TBranch        *b_trackLostHitsGsfTrack;   //!
   TBranch        *b_trackNormalizedChi2GsfTrack;   //!
   TBranch        *b_qualityMaskGsfTrack;   //!
   TBranch        *b_impactPar3DGsfTrack;   //!
   TBranch        *b_impactPar3DErrorGsfTrack;   //!
   TBranch        *b_transvImpactParGsfTrack;   //!
   TBranch        *b_transvImpactParErrorGsfTrack;   //!
   TBranch        *b_impactPar3DBiasedGsfTrack;   //!
   TBranch        *b_impactPar3DBiasedErrorGsfTrack;   //!
   TBranch        *b_transvImpactParBiasedGsfTrack;   //!
   TBranch        *b_transvImpactParBiasedErrorGsfTrack;   //!
   TBranch        *b_trackVxGsfTrack;   //!
   TBranch        *b_trackVyGsfTrack;   //!
   TBranch        *b_trackVzGsfTrack;   //!
   TBranch        *b_pixelHitsGsfTrack;   //!
   TBranch        *b_expInnerLayersGsfTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGsfTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGsfTrack;   //!
   TBranch        *b_chargeModeGsfTrack;   //!
   TBranch        *b_pxModeGsfTrack;   //!
   TBranch        *b_pyModeGsfTrack;   //!
   TBranch        *b_pzModeGsfTrack;   //!
   TBranch        *b_nGlobalMuonTrack;   //!
   TBranch        *b_pxGlobalMuonTrack;   //!
   TBranch        *b_pyGlobalMuonTrack;   //!
   TBranch        *b_pzGlobalMuonTrack;   //!
   TBranch        *b_vtxIndexGlobalMuonTrack;   //!
   TBranch        *b_vtxWeightGlobalMuonTrack;   //!
   TBranch        *b_chargeGlobalMuonTrack;   //!
   TBranch        *b_ptErrorGlobalMuonTrack;   //!
   TBranch        *b_trackValidHitsGlobalMuonTrack;   //!
   TBranch        *b_trackLostHitsGlobalMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2GlobalMuonTrack;   //!
   TBranch        *b_qualityMaskGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DErrorGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParErrorGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DBiasedGlobalMuonTrack;   //!
   TBranch        *b_impactPar3DBiasedErrorGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParBiasedGlobalMuonTrack;   //!
   TBranch        *b_transvImpactParBiasedErrorGlobalMuonTrack;   //!
   TBranch        *b_trackVxGlobalMuonTrack;   //!
   TBranch        *b_trackVyGlobalMuonTrack;   //!
   TBranch        *b_trackVzGlobalMuonTrack;   //!
   TBranch        *b_pixelHitsGlobalMuonTrack;   //!
   TBranch        *b_expInnerLayersGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGlobalMuonTrack;   //!
   TBranch        *b_nSTAMuonTrack;   //!
   TBranch        *b_pxSTAMuonTrack;   //!
   TBranch        *b_pySTAMuonTrack;   //!
   TBranch        *b_pzSTAMuonTrack;   //!
   TBranch        *b_chargeSTAMuonTrack;   //!
   TBranch        *b_ptErrorSTAMuonTrack;   //!
   TBranch        *b_trackValidHitsSTAMuonTrack;   //!
   TBranch        *b_trackLostHitsSTAMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2STAMuonTrack;   //!
   TBranch        *b_qualityMaskSTAMuonTrack;   //!
   TBranch        *b_trackVxSTAMuonTrack;   //!
   TBranch        *b_trackVySTAMuonTrack;   //!
   TBranch        *b_trackVzSTAMuonTrack;   //!
   TBranch        *b_pixelHitsSTAMuonTrack;   //!
   TBranch        *b_expInnerLayersSTAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsSTAMuonTrack;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVxPV;   //!
   TBranch        *b_PVyPV;   //!
   TBranch        *b_PVzPV;   //!
   TBranch        *b_PVErrxPV;   //!
   TBranch        *b_PVErryPV;   //!
   TBranch        *b_PVErrzPV;   //!
   TBranch        *b_SumPtPV;   //!
   TBranch        *b_ndofPV;   //!
   TBranch        *b_chi2PV;   //!
   TBranch        *b_pxChMetPV;   //!
   TBranch        *b_pyChMetPV;   //!
   TBranch        *b_pzChMetPV;   //!
   TBranch        *b_isFakePV;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_thetaMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_pxMuon;   //!
   TBranch        *b_pyMuon;   //!
   TBranch        *b_pzMuon;   //!
   TBranch        *b_vertexXMuon;   //!
   TBranch        *b_vertexYMuon;   //!
   TBranch        *b_vertexZMuon;   //!
   TBranch        *b_trackIndexMuon;   //!
   TBranch        *b_standAloneTrackIndexMuon;   //!
   TBranch        *b_combinedTrackIndexMuon;   //!
   TBranch        *b_typeMuon;   //!
   TBranch        *b_numberOfMatchesMuon;   //!
   TBranch        *b_muonIdMuon;   //!
   TBranch        *b_sumPt03Muon;   //!
   TBranch        *b_emEt03Muon;   //!
   TBranch        *b_hadEt03Muon;   //!
   TBranch        *b_hoEt03Muon;   //!
   TBranch        *b_nTrk03Muon;   //!
   TBranch        *b_nJets03Muon;   //!
   TBranch        *b_sumPt05Muon;   //!
   TBranch        *b_emEt05Muon;   //!
   TBranch        *b_hadEt05Muon;   //!
   TBranch        *b_hoEt05Muon;   //!
   TBranch        *b_nTrk05Muon;   //!
   TBranch        *b_nJets05Muon;   //!
   TBranch        *b_pfChargedIsoMuon;   //!
   TBranch        *b_pfNeutralIsoMuon;   //!
   TBranch        *b_pfPhotonIsoMuon;   //!
   TBranch        *b_pfGenericChargedIsoMuon;   //!
   TBranch        *b_pfGenericNeutralIsoMuon;   //!
   TBranch        *b_pfGenericPhotonIsoMuon;   //!
   TBranch        *b_pfGenericNoOverChargedIsoMuon;   //!
   TBranch        *b_pfGenericNoOverNeutralIsoMuon;   //!
   TBranch        *b_pfGenericNoOverPhotonIsoMuon;   //!
   TBranch        *b_pfCombinedIsoMuon;   //!
   TBranch        *b_EcalExpDepoMuon;   //!
   TBranch        *b_HcalExpDepoMuon;   //!
   TBranch        *b_HoExpDepoMuon;   //!
   TBranch        *b_emS9Muon;   //!
   TBranch        *b_hadS9Muon;   //!
   TBranch        *b_hoS9Muon;   //!
   TBranch        *b_CaloCompMuon;   //!
   TBranch        *b_nReducedPFCand;   //!
   TBranch        *b_chargeReducedPFCand;   //!
   TBranch        *b_energyReducedPFCand;   //!
   TBranch        *b_thetaReducedPFCand;   //!
   TBranch        *b_etaReducedPFCand;   //!
   TBranch        *b_phiReducedPFCand;   //!
   TBranch        *b_pxReducedPFCand;   //!
   TBranch        *b_pyReducedPFCand;   //!
   TBranch        *b_pzReducedPFCand;   //!
   TBranch        *b_vertexXReducedPFCand;   //!
   TBranch        *b_vertexYReducedPFCand;   //!
   TBranch        *b_vertexZReducedPFCand;   //!
   TBranch        *b_nMet;   //!
   TBranch        *b_chargeMet;   //!
   TBranch        *b_energyMet;   //!
   TBranch        *b_thetaMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_vertexXMet;   //!
   TBranch        *b_vertexYMet;   //!
   TBranch        *b_vertexZMet;   //!
   TBranch        *b_nTCMet;   //!
   TBranch        *b_chargeTCMet;   //!
   TBranch        *b_energyTCMet;   //!
   TBranch        *b_thetaTCMet;   //!
   TBranch        *b_etaTCMet;   //!
   TBranch        *b_phiTCMet;   //!
   TBranch        *b_pxTCMet;   //!
   TBranch        *b_pyTCMet;   //!
   TBranch        *b_pzTCMet;   //!
   TBranch        *b_vertexXTCMet;   //!
   TBranch        *b_vertexYTCMet;   //!
   TBranch        *b_vertexZTCMet;   //!
   TBranch        *b_nPFMet;   //!
   TBranch        *b_chargePFMet;   //!
   TBranch        *b_energyPFMet;   //!
   TBranch        *b_thetaPFMet;   //!
   TBranch        *b_etaPFMet;   //!
   TBranch        *b_phiPFMet;   //!
   TBranch        *b_pxPFMet;   //!
   TBranch        *b_pyPFMet;   //!
   TBranch        *b_pzPFMet;   //!
   TBranch        *b_vertexXPFMet;   //!
   TBranch        *b_vertexYPFMet;   //!
   TBranch        *b_vertexZPFMet;   //!
   TBranch        *b_sumEtPFMet;   //!
   TBranch        *b_mEtSigPFMet;   //!
   TBranch        *b_significancePFMet;   //!
   TBranch        *b_nPFChMet;   //!
   TBranch        *b_chargePFChMet;   //!
   TBranch        *b_energyPFChMet;   //!
   TBranch        *b_thetaPFChMet;   //!
   TBranch        *b_etaPFChMet;   //!
   TBranch        *b_phiPFChMet;   //!
   TBranch        *b_pxPFChMet;   //!
   TBranch        *b_pyPFChMet;   //!
   TBranch        *b_pzPFChMet;   //!
   TBranch        *b_vertexXPFChMet;   //!
   TBranch        *b_vertexYPFChMet;   //!
   TBranch        *b_vertexZPFChMet;   //!
   TBranch        *b_sumEtPFChMet;   //!
   TBranch        *b_mEtSigPFChMet;   //!
   TBranch        *b_significancePFChMet;   //!
   TBranch        *b_nGenMet;   //!
   TBranch        *b_chargeGenMet;   //!
   TBranch        *b_energyGenMet;   //!
   TBranch        *b_thetaGenMet;   //!
   TBranch        *b_etaGenMet;   //!
   TBranch        *b_phiGenMet;   //!
   TBranch        *b_pxGenMet;   //!
   TBranch        *b_pyGenMet;   //!
   TBranch        *b_pzGenMet;   //!
   TBranch        *b_vertexXGenMet;   //!
   TBranch        *b_vertexYGenMet;   //!
   TBranch        *b_vertexZGenMet;   //!
   TBranch        *b_nAK5Jet;   //!
   TBranch        *b_chargeAK5Jet;   //!
   TBranch        *b_energyAK5Jet;   //!
   TBranch        *b_thetaAK5Jet;   //!
   TBranch        *b_etaAK5Jet;   //!
   TBranch        *b_phiAK5Jet;   //!
   TBranch        *b_pxAK5Jet;   //!
   TBranch        *b_pyAK5Jet;   //!
   TBranch        *b_pzAK5Jet;   //!
   TBranch        *b_vertexXAK5Jet;   //!
   TBranch        *b_vertexYAK5Jet;   //!
   TBranch        *b_vertexZAK5Jet;   //!
   TBranch        *b_emFracAK5Jet;   //!
   TBranch        *b_hadFracAK5Jet;   //!
   TBranch        *b_IdAK5Jet;   //!
   TBranch        *b_nHitAK5Jet;   //!
   TBranch        *b_nHit90AK5Jet;   //!
   TBranch        *b_fHPDAK5Jet;   //!
   TBranch        *b_covEtaEtaAK5Jet;   //!
   TBranch        *b_covPhiPhiAK5Jet;   //!
   TBranch        *b_fLSAK5Jet;   //!
   TBranch        *b_fOOTAK5Jet;   //!
   TBranch        *b_combinedSecondaryVertexBJetTagsAK5Jet;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5Jet;   //!
   TBranch        *b_jetBProbabilityBJetTagsAK5Jet;   //!
   TBranch        *b_jetProbabilityBJetTagsAK5Jet;   //!
   TBranch        *b_simpleSecondaryVertexHighEffBJetTagsAK5Jet;   //!
   TBranch        *b_simpleSecondaryVertexHighPurBJetTagsAK5Jet;   //!
   TBranch        *b_softMuonBJetTagsAK5Jet;   //!
   TBranch        *b_softMuonByIP3dBJetTagsAK5Jet;   //!
   TBranch        *b_softMuonByPtBJetTagsAK5Jet;   //!
   TBranch        *b_softElectronByIP3dBJetTagsAK5Jet;   //!
   TBranch        *b_softElectronByPtBJetTagsAK5Jet;   //!
   TBranch        *b_trackCountingHighPurBJetTagsAK5Jet;   //!
   TBranch        *b_trackCountingHighEffBJetTagsAK5Jet;   //!
   TBranch        *b_uncorrEnergyAK5Jet;   //!
   TBranch        *b_nAK5PFPUcorrJet;   //!
   TBranch        *b_chargeAK5PFPUcorrJet;   //!
   TBranch        *b_energyAK5PFPUcorrJet;   //!
   TBranch        *b_thetaAK5PFPUcorrJet;   //!
   TBranch        *b_etaAK5PFPUcorrJet;   //!
   TBranch        *b_phiAK5PFPUcorrJet;   //!
   TBranch        *b_pxAK5PFPUcorrJet;   //!
   TBranch        *b_pyAK5PFPUcorrJet;   //!
   TBranch        *b_pzAK5PFPUcorrJet;   //!
   TBranch        *b_vertexXAK5PFPUcorrJet;   //!
   TBranch        *b_vertexYAK5PFPUcorrJet;   //!
   TBranch        *b_vertexZAK5PFPUcorrJet;   //!
   TBranch        *b_chargedHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_neutralHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_photonEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_electronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_muonEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_HFHadronEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_HFEMEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_chargedHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_neutralHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_photonMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_electronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_muonMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_HFHadronMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_HFEMMultiplicityAK5PFPUcorrJet;   //!
   TBranch        *b_chargedEmEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_neutralEmEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_combinedSecondaryVertexBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_jetBProbabilityBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_jetProbabilityBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonByIP3dBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softMuonByPtBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronByIP3dBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_softElectronByPtBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_trackCountingHighPurBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_trackCountingHighEffBJetTagsAK5PFPUcorrJet;   //!
   TBranch        *b_uncorrEnergyAK5PFPUcorrJet;   //!
   TBranch        *b_rmsCandAK5PFPUcorrJet;   //!
   TBranch        *b_ptDAK5PFPUcorrJet;   //!
   TBranch        *b_nAK5GenJet;   //!
   TBranch        *b_chargeAK5GenJet;   //!
   TBranch        *b_energyAK5GenJet;   //!
   TBranch        *b_thetaAK5GenJet;   //!
   TBranch        *b_etaAK5GenJet;   //!
   TBranch        *b_phiAK5GenJet;   //!
   TBranch        *b_pxAK5GenJet;   //!
   TBranch        *b_pyAK5GenJet;   //!
   TBranch        *b_pzAK5GenJet;   //!
   TBranch        *b_vertexXAK5GenJet;   //!
   TBranch        *b_vertexYAK5GenJet;   //!
   TBranch        *b_vertexZAK5GenJet;   //!
   TBranch        *b_genPtHat;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genAlphaQCD;   //!
   TBranch        *b_genAlphaQED;   //!
   TBranch        *b_nTriggerObsPassing;
   TBranch        *b_sizePassing;
   TBranch        *b_indexPassingPerPath;
   TBranch        *b_indexPassing;
   TBranch        *b_nTriggerObs;
   TBranch        *b_triggerObsPt; 
   TBranch        *b_triggerObsEta; 
   TBranch        *b_triggerObsPhi; 
   TBranch        *b_triggerObsMass; 

   HiggsBase(TTree *tree=0);
   virtual ~HiggsBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HiggsBase_cxx
HiggsBase::HiggsBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/crovelli/default_MC_9_1.root");
      if (!f) {
         f = new TFile("/tmp/crovelli/default_MC_9_1.root");
      }
      tree = (TTree*)gDirectory->Get("ntp1");

   }
   Init(tree);
}

HiggsBase::~HiggsBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HiggsBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HiggsBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HiggsBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   nameHLT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nl1Technical", &nl1Technical, &b_nl1Technical);
   fChain->SetBranchAddress("l1Technical", l1Technical, &b_l1Technical);
   fChain->SetBranchAddress("nl1Global", &nl1Global, &b_nl1Global);
   fChain->SetBranchAddress("l1Global", l1Global, &b_l1Global);
   fChain->SetBranchAddress("rhoFastjet", &rhoFastjet, &b_rhoFastjet);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("nameHLT", &nameHLT, &b_nameHLT);
   fChain->SetBranchAddress("indexHLT", indexHLT, &b_indexHLT);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("fiducialFlagsEle", fiducialFlagsEle, &b_fiducialFlagsEle);
   fChain->SetBranchAddress("recoFlagsEle", recoFlagsEle, &b_recoFlagsEle);
   fChain->SetBranchAddress("energyCorrectionsEle", energyCorrectionsEle, &b_energyCorrectionsEle);
   fChain->SetBranchAddress("superClusterIndexEle", superClusterIndexEle, &b_superClusterIndexEle);
   fChain->SetBranchAddress("PFsuperClusterIndexEle", PFsuperClusterIndexEle, &b_PFsuperClusterIndexEle);
   fChain->SetBranchAddress("trackIndexEle", trackIndexEle, &b_trackIndexEle);
   fChain->SetBranchAddress("gsfTrackIndexEle", gsfTrackIndexEle, &b_gsfTrackIndexEle);
   fChain->SetBranchAddress("convDistEle", convDistEle, &b_convDistEle);
   fChain->SetBranchAddress("convDcotEle", convDcotEle, &b_convDcotEle);
   fChain->SetBranchAddress("convRadiusEle", convRadiusEle, &b_convRadiusEle);
   fChain->SetBranchAddress("convTrackIndexEle", convTrackIndexEle, &b_convTrackIndexEle);
   fChain->SetBranchAddress("scPixChargeEle", scPixChargeEle, &b_scPixChargeEle);
   fChain->SetBranchAddress("hasMatchedConversionEle", hasMatchedConversionEle, &b_hasMatchedConversionEle);
   fChain->SetBranchAddress("classificationEle", classificationEle, &b_classificationEle);
   fChain->SetBranchAddress("standardClassificationEle", standardClassificationEle, &b_standardClassificationEle);
   fChain->SetBranchAddress("fbremEle", fbremEle, &b_fbremEle);
   fChain->SetBranchAddress("nbremsEle", nbremsEle, &b_nbremsEle);
   fChain->SetBranchAddress("hOverEEle", hOverEEle, &b_hOverEEle);
   fChain->SetBranchAddress("eSuperClusterOverPEle", eSuperClusterOverPEle, &b_eSuperClusterOverPEle);
   fChain->SetBranchAddress("eSeedOverPoutEle", eSeedOverPoutEle, &b_eSeedOverPoutEle);
   fChain->SetBranchAddress("deltaEtaAtVtxEle", deltaEtaAtVtxEle, &b_deltaEtaAtVtxEle);
   fChain->SetBranchAddress("deltaPhiAtVtxEle", deltaPhiAtVtxEle, &b_deltaPhiAtVtxEle);
   fChain->SetBranchAddress("deltaEtaAtCaloEle", deltaEtaAtCaloEle, &b_deltaEtaAtCaloEle);
   fChain->SetBranchAddress("deltaPhiAtCaloEle", deltaPhiAtCaloEle, &b_deltaPhiAtCaloEle);
   fChain->SetBranchAddress("dr03TkSumPtEle", dr03TkSumPtEle, &b_dr03TkSumPtEle);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtEle", dr03EcalRecHitSumEtEle, &b_dr03EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr03HcalTowerSumEtEle", dr03HcalTowerSumEtEle, &b_dr03HcalTowerSumEtEle);
   fChain->SetBranchAddress("dr04TkSumPtEle", dr04TkSumPtEle, &b_dr04TkSumPtEle);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtEle", dr04EcalRecHitSumEtEle, &b_dr04EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr04HcalTowerSumEtEle", dr04HcalTowerSumEtEle, &b_dr04HcalTowerSumEtEle);
   fChain->SetBranchAddress("scBasedEcalSum03Ele", scBasedEcalSum03Ele, &b_scBasedEcalSum03Ele);
   fChain->SetBranchAddress("scBasedEcalSum04Ele", scBasedEcalSum04Ele, &b_scBasedEcalSum04Ele);
   fChain->SetBranchAddress("dr03HcalTowerSumEtFullConeEle", dr03HcalTowerSumEtFullConeEle, &b_dr03HcalTowerSumEtFullConeEle);
   fChain->SetBranchAddress("dr04HcalTowerSumEtFullConeEle", dr04HcalTowerSumEtFullConeEle, &b_dr04HcalTowerSumEtFullConeEle);
   fChain->SetBranchAddress("eleIdLikelihoodEle", eleIdLikelihoodEle, &b_eleIdLikelihoodEle);
   fChain->SetBranchAddress("pflowMVAEle", pflowMVAEle, &b_pflowMVAEle);
   fChain->SetBranchAddress("pfChargedIsoEle", pfChargedIsoEle, &b_pfChargedIsoEle);
   fChain->SetBranchAddress("pfNeutralIsoEle", pfNeutralIsoEle, &b_pfNeutralIsoEle);
   fChain->SetBranchAddress("pfPhotonIsoEle", pfPhotonIsoEle, &b_pfPhotonIsoEle);
   fChain->SetBranchAddress("pfGenericChargedIsoEle", pfGenericChargedIsoEle, &b_pfGenericChargedIsoEle);
   fChain->SetBranchAddress("pfGenericNeutralIsoEle", pfGenericNeutralIsoEle, &b_pfGenericNeutralIsoEle);
   fChain->SetBranchAddress("pfGenericPhotonIsoEle", pfGenericPhotonIsoEle, &b_pfGenericPhotonIsoEle);
   fChain->SetBranchAddress("pfGenericNoOverChargedIsoEle", pfGenericNoOverChargedIsoEle, &b_pfGenericNoOverChargedIsoEle);
   fChain->SetBranchAddress("pfGenericNoOverNeutralIsoEle", pfGenericNoOverNeutralIsoEle, &b_pfGenericNoOverNeutralIsoEle);
   fChain->SetBranchAddress("pfGenericNoOverPhotonIsoEle", pfGenericNoOverPhotonIsoEle, &b_pfGenericNoOverPhotonIsoEle);
   fChain->SetBranchAddress("pfCombinedIsoEle", pfCombinedIsoEle, &b_pfCombinedIsoEle);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("chargePho", chargePho, &b_chargePho);
   fChain->SetBranchAddress("energyPho", energyPho, &b_energyPho);
   fChain->SetBranchAddress("thetaPho", thetaPho, &b_thetaPho);
   fChain->SetBranchAddress("etaPho", etaPho, &b_etaPho);
   fChain->SetBranchAddress("phiPho", phiPho, &b_phiPho);
   fChain->SetBranchAddress("pxPho", pxPho, &b_pxPho);
   fChain->SetBranchAddress("pyPho", pyPho, &b_pyPho);
   fChain->SetBranchAddress("pzPho", pzPho, &b_pzPho);
   fChain->SetBranchAddress("vertexXPho", vertexXPho, &b_vertexXPho);
   fChain->SetBranchAddress("vertexYPho", vertexYPho, &b_vertexYPho);
   fChain->SetBranchAddress("vertexZPho", vertexZPho, &b_vertexZPho);
   fChain->SetBranchAddress("fiducialFlagsPho", fiducialFlagsPho, &b_fiducialFlagsPho);
   fChain->SetBranchAddress("recoFlagsPho", recoFlagsPho, &b_recoFlagsPho);
   fChain->SetBranchAddress("superClusterIndexPho", superClusterIndexPho, &b_superClusterIndexPho);
   fChain->SetBranchAddress("PFsuperClusterIndexPho", PFsuperClusterIndexPho, &b_PFsuperClusterIndexPho);
   fChain->SetBranchAddress("hOverEPho", hOverEPho, &b_hOverEPho);
   fChain->SetBranchAddress("dr03TkSumPtPho", dr03TkSumPtPho, &b_dr03TkSumPtPho);
   fChain->SetBranchAddress("dr03HollowTkSumPtPho", dr03HollowTkSumPtPho, &b_dr03HollowTkSumPtPho);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtPho", dr03EcalRecHitSumEtPho, &b_dr03EcalRecHitSumEtPho);
   fChain->SetBranchAddress("dr03HcalTowerSumEtPho", dr03HcalTowerSumEtPho, &b_dr03HcalTowerSumEtPho);
   fChain->SetBranchAddress("dr04TkSumPtPho", dr04TkSumPtPho, &b_dr04TkSumPtPho);
   fChain->SetBranchAddress("dr04HollowTkSumPtPho", dr04HollowTkSumPtPho, &b_dr04HollowTkSumPtPho);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtPho", dr04EcalRecHitSumEtPho, &b_dr04EcalRecHitSumEtPho);
   fChain->SetBranchAddress("dr04HcalTowerSumEtPho", dr04HcalTowerSumEtPho, &b_dr04HcalTowerSumEtPho);
   fChain->SetBranchAddress("chargedHadronIsoPho", chargedHadronIsoPho, &b_chargedHadronIsoPho);
   fChain->SetBranchAddress("neutralHadronIsoPho", neutralHadronIsoPho, &b_neutralHadronIsoPho);
   fChain->SetBranchAddress("photonIsoPho", photonIsoPho, &b_photonIsoPho);
   fChain->SetBranchAddress("hasPixelSeedPho", hasPixelSeedPho, &b_hasPixelSeedPho);
   fChain->SetBranchAddress("hasMatchedConversionPho", hasMatchedConversionPho, &b_hasMatchedConversionPho);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("nBCSC", nBCSC, &b_nBCSC);
   fChain->SetBranchAddress("nCrystalsSC", nCrystalsSC, &b_nCrystalsSC);
   fChain->SetBranchAddress("rawEnergySC", rawEnergySC, &b_rawEnergySC);
   fChain->SetBranchAddress("energySC", energySC, &b_energySC);
   fChain->SetBranchAddress("etaSC", etaSC, &b_etaSC);
   fChain->SetBranchAddress("thetaSC", thetaSC, &b_thetaSC);
   fChain->SetBranchAddress("phiSC", phiSC, &b_phiSC);
   fChain->SetBranchAddress("phiWidthSC", phiWidthSC, &b_phiWidthSC);
   fChain->SetBranchAddress("etaWidthSC", etaWidthSC, &b_etaWidthSC);
   fChain->SetBranchAddress("e3x3SC", e3x3SC, &b_e3x3SC);
   fChain->SetBranchAddress("e5x5SC", e5x5SC, &b_e5x5SC);
   fChain->SetBranchAddress("eMaxSC", eMaxSC, &b_eMaxSC);
   fChain->SetBranchAddress("e2x2SC", e2x2SC, &b_e2x2SC);
   fChain->SetBranchAddress("e2ndSC", e2ndSC, &b_e2ndSC);
   fChain->SetBranchAddress("e1x5SC", e1x5SC, &b_e1x5SC);
   fChain->SetBranchAddress("e2x5MaxSC", e2x5MaxSC, &b_e2x5MaxSC);
   fChain->SetBranchAddress("e4SwissCrossSC", e4SwissCrossSC, &b_e4SwissCrossSC);
   fChain->SetBranchAddress("covIEtaIEtaSC", covIEtaIEtaSC, &b_covIEtaIEtaSC);
   fChain->SetBranchAddress("covIEtaIPhiSC", covIEtaIPhiSC, &b_covIEtaIPhiSC);
   fChain->SetBranchAddress("covIPhiIPhiSC", covIPhiIPhiSC, &b_covIPhiIPhiSC);
   fChain->SetBranchAddress("hOverESC", hOverESC, &b_hOverESC);
   fChain->SetBranchAddress("recoFlagSC", recoFlagSC, &b_recoFlagSC);
   fChain->SetBranchAddress("timeSC", timeSC, &b_timeSC);
   fChain->SetBranchAddress("chi2SC", chi2SC, &b_chi2SC);
   fChain->SetBranchAddress("seedEnergySC", seedEnergySC, &b_seedEnergySC);
   fChain->SetBranchAddress("nPFSC", &nPFSC, &b_nPFSC);
   fChain->SetBranchAddress("nBCPFSC", nBCPFSC, &b_nBCPFSC);
   fChain->SetBranchAddress("nCrystalsPFSC", nCrystalsPFSC, &b_nCrystalsPFSC);
   fChain->SetBranchAddress("rawEnergyPFSC", rawEnergyPFSC, &b_rawEnergyPFSC);
   fChain->SetBranchAddress("energyPFSC", energyPFSC, &b_energyPFSC);
   fChain->SetBranchAddress("etaPFSC", etaPFSC, &b_etaPFSC);
   fChain->SetBranchAddress("thetaPFSC", thetaPFSC, &b_thetaPFSC);
   fChain->SetBranchAddress("phiPFSC", phiPFSC, &b_phiPFSC);
   fChain->SetBranchAddress("phiWidthPFSC", phiWidthPFSC, &b_phiWidthPFSC);
   fChain->SetBranchAddress("etaWidthPFSC", etaWidthPFSC, &b_etaWidthPFSC);
   fChain->SetBranchAddress("e3x3PFSC", e3x3PFSC, &b_e3x3PFSC);
   fChain->SetBranchAddress("e5x5PFSC", e5x5PFSC, &b_e5x5PFSC);
   fChain->SetBranchAddress("eMaxPFSC", eMaxPFSC, &b_eMaxPFSC);
   fChain->SetBranchAddress("e2x2PFSC", e2x2PFSC, &b_e2x2PFSC);
   fChain->SetBranchAddress("e2ndPFSC", e2ndPFSC, &b_e2ndPFSC);
   fChain->SetBranchAddress("e1x5PFSC", e1x5PFSC, &b_e1x5PFSC);
   fChain->SetBranchAddress("e2x5MaxPFSC", e2x5MaxPFSC, &b_e2x5MaxPFSC);
   fChain->SetBranchAddress("e4SwissCrossPFSC", e4SwissCrossPFSC, &b_e4SwissCrossPFSC);
   fChain->SetBranchAddress("covIEtaIEtaPFSC", covIEtaIEtaPFSC, &b_covIEtaIEtaPFSC);
   fChain->SetBranchAddress("covIEtaIPhiPFSC", covIEtaIPhiPFSC, &b_covIEtaIPhiPFSC);
   fChain->SetBranchAddress("covIPhiIPhiPFSC", covIPhiIPhiPFSC, &b_covIPhiIPhiPFSC);
   fChain->SetBranchAddress("hOverEPFSC", hOverEPFSC, &b_hOverEPFSC);
   fChain->SetBranchAddress("recoFlagPFSC", recoFlagPFSC, &b_recoFlagPFSC);
   fChain->SetBranchAddress("timePFSC", timePFSC, &b_timePFSC);
   fChain->SetBranchAddress("chi2PFSC", chi2PFSC, &b_chi2PFSC);
   fChain->SetBranchAddress("seedEnergyPFSC", seedEnergyPFSC, &b_seedEnergyPFSC);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("pxTrack", pxTrack, &b_pxTrack);
   fChain->SetBranchAddress("pyTrack", pyTrack, &b_pyTrack);
   fChain->SetBranchAddress("pzTrack", pzTrack, &b_pzTrack);
   fChain->SetBranchAddress("vtxIndexTrack", vtxIndexTrack, &b_vtxIndexTrack);
   fChain->SetBranchAddress("vtxWeightTrack", vtxWeightTrack, &b_vtxWeightTrack);
   fChain->SetBranchAddress("chargeTrack", chargeTrack, &b_chargeTrack);
   fChain->SetBranchAddress("ptErrorTrack", ptErrorTrack, &b_ptErrorTrack);
   fChain->SetBranchAddress("trackValidHitsTrack", trackValidHitsTrack, &b_trackValidHitsTrack);
   fChain->SetBranchAddress("trackLostHitsTrack", trackLostHitsTrack, &b_trackLostHitsTrack);
   fChain->SetBranchAddress("trackNormalizedChi2Track", trackNormalizedChi2Track, &b_trackNormalizedChi2Track);
   fChain->SetBranchAddress("qualityMaskTrack", qualityMaskTrack, &b_qualityMaskTrack);
   fChain->SetBranchAddress("impactPar3DTrack", impactPar3DTrack, &b_impactPar3DTrack);
   fChain->SetBranchAddress("impactPar3DErrorTrack", impactPar3DErrorTrack, &b_impactPar3DErrorTrack);
   fChain->SetBranchAddress("transvImpactParTrack", transvImpactParTrack, &b_transvImpactParTrack);
   fChain->SetBranchAddress("transvImpactParErrorTrack", transvImpactParErrorTrack, &b_transvImpactParErrorTrack);
   fChain->SetBranchAddress("impactPar3DBiasedTrack", impactPar3DBiasedTrack, &b_impactPar3DBiasedTrack);
   fChain->SetBranchAddress("impactPar3DBiasedErrorTrack", impactPar3DBiasedErrorTrack, &b_impactPar3DBiasedErrorTrack);
   fChain->SetBranchAddress("transvImpactParBiasedTrack", transvImpactParBiasedTrack, &b_transvImpactParBiasedTrack);
   fChain->SetBranchAddress("transvImpactParBiasedErrorTrack", transvImpactParBiasedErrorTrack, &b_transvImpactParBiasedErrorTrack);
   fChain->SetBranchAddress("trackVxTrack", trackVxTrack, &b_trackVxTrack);
   fChain->SetBranchAddress("trackVyTrack", trackVyTrack, &b_trackVyTrack);
   fChain->SetBranchAddress("trackVzTrack", trackVzTrack, &b_trackVzTrack);
   fChain->SetBranchAddress("pixelHitsTrack", pixelHitsTrack, &b_pixelHitsTrack);
   fChain->SetBranchAddress("expInnerLayersTrack", expInnerLayersTrack, &b_expInnerLayersTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsTrack", numberOfValidPixelBarrelHitsTrack, &b_numberOfValidPixelBarrelHitsTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsTrack", numberOfValidPixelEndcapHitsTrack, &b_numberOfValidPixelEndcapHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsTrack", numberOfValidStripTIBHitsTrack, &b_numberOfValidStripTIBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsTrack", numberOfValidStripTIDHitsTrack, &b_numberOfValidStripTIDHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsTrack", numberOfValidStripTOBHitsTrack, &b_numberOfValidStripTOBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsTrack", numberOfValidStripTECHitsTrack, &b_numberOfValidStripTECHitsTrack);
   fChain->SetBranchAddress("nGsfTrack", &nGsfTrack, &b_nGsfTrack);
   fChain->SetBranchAddress("pxGsfTrack", pxGsfTrack, &b_pxGsfTrack);
   fChain->SetBranchAddress("pyGsfTrack", pyGsfTrack, &b_pyGsfTrack);
   fChain->SetBranchAddress("pzGsfTrack", pzGsfTrack, &b_pzGsfTrack);
   fChain->SetBranchAddress("vtxIndexGsfTrack", vtxIndexGsfTrack, &b_vtxIndexGsfTrack);
   fChain->SetBranchAddress("vtxWeightGsfTrack", vtxWeightGsfTrack, &b_vtxWeightGsfTrack);
   fChain->SetBranchAddress("chargeGsfTrack", chargeGsfTrack, &b_chargeGsfTrack);
   fChain->SetBranchAddress("ptErrorGsfTrack", ptErrorGsfTrack, &b_ptErrorGsfTrack);
   fChain->SetBranchAddress("trackValidHitsGsfTrack", trackValidHitsGsfTrack, &b_trackValidHitsGsfTrack);
   fChain->SetBranchAddress("trackLostHitsGsfTrack", trackLostHitsGsfTrack, &b_trackLostHitsGsfTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GsfTrack", trackNormalizedChi2GsfTrack, &b_trackNormalizedChi2GsfTrack);
   fChain->SetBranchAddress("qualityMaskGsfTrack", qualityMaskGsfTrack, &b_qualityMaskGsfTrack);
   fChain->SetBranchAddress("impactPar3DGsfTrack", impactPar3DGsfTrack, &b_impactPar3DGsfTrack);
   fChain->SetBranchAddress("impactPar3DErrorGsfTrack", impactPar3DErrorGsfTrack, &b_impactPar3DErrorGsfTrack);
   fChain->SetBranchAddress("transvImpactParGsfTrack", transvImpactParGsfTrack, &b_transvImpactParGsfTrack);
   fChain->SetBranchAddress("transvImpactParErrorGsfTrack", transvImpactParErrorGsfTrack, &b_transvImpactParErrorGsfTrack);
   fChain->SetBranchAddress("impactPar3DBiasedGsfTrack", impactPar3DBiasedGsfTrack, &b_impactPar3DBiasedGsfTrack);
   fChain->SetBranchAddress("impactPar3DBiasedErrorGsfTrack", impactPar3DBiasedErrorGsfTrack, &b_impactPar3DBiasedErrorGsfTrack);
   fChain->SetBranchAddress("transvImpactParBiasedGsfTrack", transvImpactParBiasedGsfTrack, &b_transvImpactParBiasedGsfTrack);
   fChain->SetBranchAddress("transvImpactParBiasedErrorGsfTrack", transvImpactParBiasedErrorGsfTrack, &b_transvImpactParBiasedErrorGsfTrack);
   fChain->SetBranchAddress("trackVxGsfTrack", trackVxGsfTrack, &b_trackVxGsfTrack);
   fChain->SetBranchAddress("trackVyGsfTrack", trackVyGsfTrack, &b_trackVyGsfTrack);
   fChain->SetBranchAddress("trackVzGsfTrack", trackVzGsfTrack, &b_trackVzGsfTrack);
   fChain->SetBranchAddress("pixelHitsGsfTrack", pixelHitsGsfTrack, &b_pixelHitsGsfTrack);
   fChain->SetBranchAddress("expInnerLayersGsfTrack", expInnerLayersGsfTrack, &b_expInnerLayersGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGsfTrack", numberOfValidPixelBarrelHitsGsfTrack, &b_numberOfValidPixelBarrelHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGsfTrack", numberOfValidPixelEndcapHitsGsfTrack, &b_numberOfValidPixelEndcapHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGsfTrack", numberOfValidStripTIBHitsGsfTrack, &b_numberOfValidStripTIBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGsfTrack", numberOfValidStripTIDHitsGsfTrack, &b_numberOfValidStripTIDHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGsfTrack", numberOfValidStripTOBHitsGsfTrack, &b_numberOfValidStripTOBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGsfTrack", numberOfValidStripTECHitsGsfTrack, &b_numberOfValidStripTECHitsGsfTrack);
   fChain->SetBranchAddress("chargeModeGsfTrack", chargeModeGsfTrack, &b_chargeModeGsfTrack);
   fChain->SetBranchAddress("pxModeGsfTrack", pxModeGsfTrack, &b_pxModeGsfTrack);
   fChain->SetBranchAddress("pyModeGsfTrack", pyModeGsfTrack, &b_pyModeGsfTrack);
   fChain->SetBranchAddress("pzModeGsfTrack", pzModeGsfTrack, &b_pzModeGsfTrack);
   fChain->SetBranchAddress("nGlobalMuonTrack", &nGlobalMuonTrack, &b_nGlobalMuonTrack);
   fChain->SetBranchAddress("pxGlobalMuonTrack", pxGlobalMuonTrack, &b_pxGlobalMuonTrack);
   fChain->SetBranchAddress("pyGlobalMuonTrack", pyGlobalMuonTrack, &b_pyGlobalMuonTrack);
   fChain->SetBranchAddress("pzGlobalMuonTrack", pzGlobalMuonTrack, &b_pzGlobalMuonTrack);
   fChain->SetBranchAddress("vtxIndexGlobalMuonTrack", vtxIndexGlobalMuonTrack, &b_vtxIndexGlobalMuonTrack);
   fChain->SetBranchAddress("vtxWeightGlobalMuonTrack", vtxWeightGlobalMuonTrack, &b_vtxWeightGlobalMuonTrack);
   fChain->SetBranchAddress("chargeGlobalMuonTrack", chargeGlobalMuonTrack, &b_chargeGlobalMuonTrack);
   fChain->SetBranchAddress("ptErrorGlobalMuonTrack", ptErrorGlobalMuonTrack, &b_ptErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackValidHitsGlobalMuonTrack", trackValidHitsGlobalMuonTrack, &b_trackValidHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackLostHitsGlobalMuonTrack", trackLostHitsGlobalMuonTrack, &b_trackLostHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GlobalMuonTrack", trackNormalizedChi2GlobalMuonTrack, &b_trackNormalizedChi2GlobalMuonTrack);
   fChain->SetBranchAddress("qualityMaskGlobalMuonTrack", qualityMaskGlobalMuonTrack, &b_qualityMaskGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DGlobalMuonTrack", impactPar3DGlobalMuonTrack, &b_impactPar3DGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DErrorGlobalMuonTrack", impactPar3DErrorGlobalMuonTrack, &b_impactPar3DErrorGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParGlobalMuonTrack", transvImpactParGlobalMuonTrack, &b_transvImpactParGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParErrorGlobalMuonTrack", transvImpactParErrorGlobalMuonTrack, &b_transvImpactParErrorGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DBiasedGlobalMuonTrack", impactPar3DBiasedGlobalMuonTrack, &b_impactPar3DBiasedGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DBiasedErrorGlobalMuonTrack", impactPar3DBiasedErrorGlobalMuonTrack, &b_impactPar3DBiasedErrorGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParBiasedGlobalMuonTrack", transvImpactParBiasedGlobalMuonTrack, &b_transvImpactParBiasedGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParBiasedErrorGlobalMuonTrack", transvImpactParBiasedErrorGlobalMuonTrack, &b_transvImpactParBiasedErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackVxGlobalMuonTrack", trackVxGlobalMuonTrack, &b_trackVxGlobalMuonTrack);
   fChain->SetBranchAddress("trackVyGlobalMuonTrack", trackVyGlobalMuonTrack, &b_trackVyGlobalMuonTrack);
   fChain->SetBranchAddress("trackVzGlobalMuonTrack", trackVzGlobalMuonTrack, &b_trackVzGlobalMuonTrack);
   fChain->SetBranchAddress("pixelHitsGlobalMuonTrack", pixelHitsGlobalMuonTrack, &b_pixelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("expInnerLayersGlobalMuonTrack", expInnerLayersGlobalMuonTrack, &b_expInnerLayersGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGlobalMuonTrack", numberOfValidPixelBarrelHitsGlobalMuonTrack, &b_numberOfValidPixelBarrelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGlobalMuonTrack", numberOfValidPixelEndcapHitsGlobalMuonTrack, &b_numberOfValidPixelEndcapHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGlobalMuonTrack", numberOfValidStripTIBHitsGlobalMuonTrack, &b_numberOfValidStripTIBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGlobalMuonTrack", numberOfValidStripTIDHitsGlobalMuonTrack, &b_numberOfValidStripTIDHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGlobalMuonTrack", numberOfValidStripTOBHitsGlobalMuonTrack, &b_numberOfValidStripTOBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGlobalMuonTrack", numberOfValidStripTECHitsGlobalMuonTrack, &b_numberOfValidStripTECHitsGlobalMuonTrack);
   fChain->SetBranchAddress("nSTAMuonTrack", &nSTAMuonTrack, &b_nSTAMuonTrack);
   fChain->SetBranchAddress("pxSTAMuonTrack", pxSTAMuonTrack, &b_pxSTAMuonTrack);
   fChain->SetBranchAddress("pySTAMuonTrack", pySTAMuonTrack, &b_pySTAMuonTrack);
   fChain->SetBranchAddress("pzSTAMuonTrack", pzSTAMuonTrack, &b_pzSTAMuonTrack);
   fChain->SetBranchAddress("chargeSTAMuonTrack", chargeSTAMuonTrack, &b_chargeSTAMuonTrack);
   fChain->SetBranchAddress("ptErrorSTAMuonTrack", ptErrorSTAMuonTrack, &b_ptErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackValidHitsSTAMuonTrack", trackValidHitsSTAMuonTrack, &b_trackValidHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackLostHitsSTAMuonTrack", trackLostHitsSTAMuonTrack, &b_trackLostHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2STAMuonTrack", trackNormalizedChi2STAMuonTrack, &b_trackNormalizedChi2STAMuonTrack);
   fChain->SetBranchAddress("qualityMaskSTAMuonTrack", qualityMaskSTAMuonTrack, &b_qualityMaskSTAMuonTrack);
   fChain->SetBranchAddress("trackVxSTAMuonTrack", trackVxSTAMuonTrack, &b_trackVxSTAMuonTrack);
   fChain->SetBranchAddress("trackVySTAMuonTrack", trackVySTAMuonTrack, &b_trackVySTAMuonTrack);
   fChain->SetBranchAddress("trackVzSTAMuonTrack", trackVzSTAMuonTrack, &b_trackVzSTAMuonTrack);
   fChain->SetBranchAddress("pixelHitsSTAMuonTrack", pixelHitsSTAMuonTrack, &b_pixelHitsSTAMuonTrack);
   fChain->SetBranchAddress("expInnerLayersSTAMuonTrack", expInnerLayersSTAMuonTrack, &b_expInnerLayersSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsSTAMuonTrack", numberOfValidPixelBarrelHitsSTAMuonTrack, &b_numberOfValidPixelBarrelHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsSTAMuonTrack", numberOfValidPixelEndcapHitsSTAMuonTrack, &b_numberOfValidPixelEndcapHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsSTAMuonTrack", numberOfValidStripTIBHitsSTAMuonTrack, &b_numberOfValidStripTIBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsSTAMuonTrack", numberOfValidStripTIDHitsSTAMuonTrack, &b_numberOfValidStripTIDHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsSTAMuonTrack", numberOfValidStripTOBHitsSTAMuonTrack, &b_numberOfValidStripTOBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsSTAMuonTrack", numberOfValidStripTECHitsSTAMuonTrack, &b_numberOfValidStripTECHitsSTAMuonTrack);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVxPV", PVxPV, &b_PVxPV);
   fChain->SetBranchAddress("PVyPV", PVyPV, &b_PVyPV);
   fChain->SetBranchAddress("PVzPV", PVzPV, &b_PVzPV);
   fChain->SetBranchAddress("PVErrxPV", PVErrxPV, &b_PVErrxPV);
   fChain->SetBranchAddress("PVErryPV", PVErryPV, &b_PVErryPV);
   fChain->SetBranchAddress("PVErrzPV", PVErrzPV, &b_PVErrzPV);
   fChain->SetBranchAddress("SumPtPV", SumPtPV, &b_SumPtPV);
   fChain->SetBranchAddress("ndofPV", ndofPV, &b_ndofPV);
   fChain->SetBranchAddress("chi2PV", chi2PV, &b_chi2PV);
   fChain->SetBranchAddress("pxChMetPV", pxChMetPV, &b_pxChMetPV);
   fChain->SetBranchAddress("pyChMetPV", pyChMetPV, &b_pyChMetPV);
   fChain->SetBranchAddress("pzChMetPV", pzChMetPV, &b_pzChMetPV);
   fChain->SetBranchAddress("isFakePV", isFakePV, &b_isFakePV);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("chargeMuon", chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("energyMuon", energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("thetaMuon", thetaMuon, &b_thetaMuon);
   fChain->SetBranchAddress("etaMuon", etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("pxMuon", pxMuon, &b_pxMuon);
   fChain->SetBranchAddress("pyMuon", pyMuon, &b_pyMuon);
   fChain->SetBranchAddress("pzMuon", pzMuon, &b_pzMuon);
   fChain->SetBranchAddress("vertexXMuon", vertexXMuon, &b_vertexXMuon);
   fChain->SetBranchAddress("vertexYMuon", vertexYMuon, &b_vertexYMuon);
   fChain->SetBranchAddress("vertexZMuon", vertexZMuon, &b_vertexZMuon);
   fChain->SetBranchAddress("trackIndexMuon", trackIndexMuon, &b_trackIndexMuon);
   fChain->SetBranchAddress("standAloneTrackIndexMuon", standAloneTrackIndexMuon, &b_standAloneTrackIndexMuon);
   fChain->SetBranchAddress("combinedTrackIndexMuon", combinedTrackIndexMuon, &b_combinedTrackIndexMuon);
   fChain->SetBranchAddress("muonIdMuon", muonIdMuon, &b_muonIdMuon);
   fChain->SetBranchAddress("typeMuon", typeMuon, &b_typeMuon);
   fChain->SetBranchAddress("numberOfMatchesMuon", numberOfMatchesMuon, &b_numberOfMatchesMuon);
   fChain->SetBranchAddress("sumPt03Muon", sumPt03Muon, &b_sumPt03Muon);
   fChain->SetBranchAddress("emEt03Muon", emEt03Muon, &b_emEt03Muon);
   fChain->SetBranchAddress("hadEt03Muon", hadEt03Muon, &b_hadEt03Muon);
   fChain->SetBranchAddress("hoEt03Muon", hoEt03Muon, &b_hoEt03Muon);
   fChain->SetBranchAddress("nTrk03Muon", nTrk03Muon, &b_nTrk03Muon);
   fChain->SetBranchAddress("nJets03Muon", nJets03Muon, &b_nJets03Muon);
   fChain->SetBranchAddress("sumPt05Muon", sumPt05Muon, &b_sumPt05Muon);
   fChain->SetBranchAddress("emEt05Muon", emEt05Muon, &b_emEt05Muon);
   fChain->SetBranchAddress("hadEt05Muon", hadEt05Muon, &b_hadEt05Muon);
   fChain->SetBranchAddress("hoEt05Muon", hoEt05Muon, &b_hoEt05Muon);
   fChain->SetBranchAddress("nTrk05Muon", nTrk05Muon, &b_nTrk05Muon);
   fChain->SetBranchAddress("nJets05Muon", nJets05Muon, &b_nJets05Muon);
   fChain->SetBranchAddress("pfChargedIsoMuon", pfChargedIsoMuon, &b_pfChargedIsoMuon);
   fChain->SetBranchAddress("pfNeutralIsoMuon", pfNeutralIsoMuon, &b_pfNeutralIsoMuon);
   fChain->SetBranchAddress("pfPhotonIsoMuon", pfPhotonIsoMuon, &b_pfPhotonIsoMuon);
   fChain->SetBranchAddress("pfGenericChargedIsoMuon", pfGenericChargedIsoMuon, &b_pfGenericChargedIsoMuon);
   fChain->SetBranchAddress("pfGenericNeutralIsoMuon", pfGenericNeutralIsoMuon, &b_pfGenericNeutralIsoMuon);
   fChain->SetBranchAddress("pfGenericPhotonIsoMuon", pfGenericPhotonIsoMuon, &b_pfGenericPhotonIsoMuon);
   fChain->SetBranchAddress("pfGenericNoOverChargedIsoMuon", pfGenericNoOverChargedIsoMuon, &b_pfGenericNoOverChargedIsoMuon);
   fChain->SetBranchAddress("pfGenericNoOverNeutralIsoMuon", pfGenericNoOverNeutralIsoMuon, &b_pfGenericNoOverNeutralIsoMuon);
   fChain->SetBranchAddress("pfGenericNoOverPhotonIsoMuon", pfGenericNoOverPhotonIsoMuon, &b_pfGenericNoOverPhotonIsoMuon);
   fChain->SetBranchAddress("pfCombinedIsoMuon", pfCombinedIsoMuon, &b_pfCombinedIsoMuon);
   fChain->SetBranchAddress("EcalExpDepoMuon", EcalExpDepoMuon, &b_EcalExpDepoMuon);
   fChain->SetBranchAddress("HcalExpDepoMuon", HcalExpDepoMuon, &b_HcalExpDepoMuon);
   fChain->SetBranchAddress("HoExpDepoMuon", HoExpDepoMuon, &b_HoExpDepoMuon);
   fChain->SetBranchAddress("emS9Muon", emS9Muon, &b_emS9Muon);
   fChain->SetBranchAddress("hadS9Muon", hadS9Muon, &b_hadS9Muon);
   fChain->SetBranchAddress("hoS9Muon", hoS9Muon, &b_hoS9Muon);
   fChain->SetBranchAddress("CaloCompMuon", CaloCompMuon, &b_CaloCompMuon);
   fChain->SetBranchAddress("nReducedPFCand", &nReducedPFCand, &b_nReducedPFCand);
   fChain->SetBranchAddress("chargeReducedPFCand", chargeReducedPFCand, &b_chargeReducedPFCand);
   fChain->SetBranchAddress("energyReducedPFCand", energyReducedPFCand, &b_energyReducedPFCand);
   fChain->SetBranchAddress("thetaReducedPFCand", thetaReducedPFCand, &b_thetaReducedPFCand);
   fChain->SetBranchAddress("etaReducedPFCand", etaReducedPFCand, &b_etaReducedPFCand);
   fChain->SetBranchAddress("phiReducedPFCand", phiReducedPFCand, &b_phiReducedPFCand);
   fChain->SetBranchAddress("pxReducedPFCand", pxReducedPFCand, &b_pxReducedPFCand);
   fChain->SetBranchAddress("pyReducedPFCand", pyReducedPFCand, &b_pyReducedPFCand);
   fChain->SetBranchAddress("pzReducedPFCand", pzReducedPFCand, &b_pzReducedPFCand);
   fChain->SetBranchAddress("vertexXReducedPFCand", vertexXReducedPFCand, &b_vertexXReducedPFCand);
   fChain->SetBranchAddress("vertexYReducedPFCand", vertexYReducedPFCand, &b_vertexYReducedPFCand);
   fChain->SetBranchAddress("vertexZReducedPFCand", vertexZReducedPFCand, &b_vertexZReducedPFCand);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("nTCMet", &nTCMet, &b_nTCMet);
   fChain->SetBranchAddress("chargeTCMet", chargeTCMet, &b_chargeTCMet);
   fChain->SetBranchAddress("energyTCMet", energyTCMet, &b_energyTCMet);
   fChain->SetBranchAddress("thetaTCMet", thetaTCMet, &b_thetaTCMet);
   fChain->SetBranchAddress("etaTCMet", etaTCMet, &b_etaTCMet);
   fChain->SetBranchAddress("phiTCMet", phiTCMet, &b_phiTCMet);
   fChain->SetBranchAddress("pxTCMet", pxTCMet, &b_pxTCMet);
   fChain->SetBranchAddress("pyTCMet", pyTCMet, &b_pyTCMet);
   fChain->SetBranchAddress("pzTCMet", pzTCMet, &b_pzTCMet);
   fChain->SetBranchAddress("vertexXTCMet", vertexXTCMet, &b_vertexXTCMet);
   fChain->SetBranchAddress("vertexYTCMet", vertexYTCMet, &b_vertexYTCMet);
   fChain->SetBranchAddress("vertexZTCMet", vertexZTCMet, &b_vertexZTCMet);
   fChain->SetBranchAddress("nPFMet", &nPFMet, &b_nPFMet);
   fChain->SetBranchAddress("chargePFMet", chargePFMet, &b_chargePFMet);
   fChain->SetBranchAddress("energyPFMet", energyPFMet, &b_energyPFMet);
   fChain->SetBranchAddress("thetaPFMet", thetaPFMet, &b_thetaPFMet);
   fChain->SetBranchAddress("etaPFMet", etaPFMet, &b_etaPFMet);
   fChain->SetBranchAddress("phiPFMet", phiPFMet, &b_phiPFMet);
   fChain->SetBranchAddress("pxPFMet", pxPFMet, &b_pxPFMet);
   fChain->SetBranchAddress("pyPFMet", pyPFMet, &b_pyPFMet);
   fChain->SetBranchAddress("pzPFMet", pzPFMet, &b_pzPFMet);
   fChain->SetBranchAddress("vertexXPFMet", vertexXPFMet, &b_vertexXPFMet);
   fChain->SetBranchAddress("vertexYPFMet", vertexYPFMet, &b_vertexYPFMet);
   fChain->SetBranchAddress("vertexZPFMet", vertexZPFMet, &b_vertexZPFMet);
   fChain->SetBranchAddress("sumEtPFMet", sumEtPFMet, &b_sumEtPFMet);
   fChain->SetBranchAddress("mEtSigPFMet", mEtSigPFMet, &b_mEtSigPFMet);
   fChain->SetBranchAddress("significancePFMet", significancePFMet, &b_significancePFMet);
   fChain->SetBranchAddress("nPFChMet", &nPFChMet, &b_nPFChMet);
   fChain->SetBranchAddress("chargePFChMet", chargePFChMet, &b_chargePFChMet);
   fChain->SetBranchAddress("energyPFChMet", energyPFChMet, &b_energyPFChMet);
   fChain->SetBranchAddress("thetaPFChMet", thetaPFChMet, &b_thetaPFChMet);
   fChain->SetBranchAddress("etaPFChMet", etaPFChMet, &b_etaPFChMet);
   fChain->SetBranchAddress("phiPFChMet", phiPFChMet, &b_phiPFChMet);
   fChain->SetBranchAddress("pxPFChMet", pxPFChMet, &b_pxPFChMet);
   fChain->SetBranchAddress("pyPFChMet", pyPFChMet, &b_pyPFChMet);
   fChain->SetBranchAddress("pzPFChMet", pzPFChMet, &b_pzPFChMet);
   fChain->SetBranchAddress("vertexXPFChMet", vertexXPFChMet, &b_vertexXPFChMet);
   fChain->SetBranchAddress("vertexYPFChMet", vertexYPFChMet, &b_vertexYPFChMet);
   fChain->SetBranchAddress("vertexZPFChMet", vertexZPFChMet, &b_vertexZPFChMet);
   fChain->SetBranchAddress("sumEtPFChMet", sumEtPFChMet, &b_sumEtPFChMet);
   fChain->SetBranchAddress("mEtSigPFChMet", mEtSigPFChMet, &b_mEtSigPFChMet);
   fChain->SetBranchAddress("significancePFChMet", significancePFChMet, &b_significancePFChMet);
   fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
   fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
   fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
   fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
   fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
   fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
   fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
   fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
   fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
   fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
   fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
   fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   fChain->SetBranchAddress("nAK5Jet", &nAK5Jet, &b_nAK5Jet);
   fChain->SetBranchAddress("chargeAK5Jet", chargeAK5Jet, &b_chargeAK5Jet);
   fChain->SetBranchAddress("energyAK5Jet", energyAK5Jet, &b_energyAK5Jet);
   fChain->SetBranchAddress("thetaAK5Jet", thetaAK5Jet, &b_thetaAK5Jet);
   fChain->SetBranchAddress("etaAK5Jet", etaAK5Jet, &b_etaAK5Jet);
   fChain->SetBranchAddress("phiAK5Jet", phiAK5Jet, &b_phiAK5Jet);
   fChain->SetBranchAddress("pxAK5Jet", pxAK5Jet, &b_pxAK5Jet);
   fChain->SetBranchAddress("pyAK5Jet", pyAK5Jet, &b_pyAK5Jet);
   fChain->SetBranchAddress("pzAK5Jet", pzAK5Jet, &b_pzAK5Jet);
   fChain->SetBranchAddress("vertexXAK5Jet", vertexXAK5Jet, &b_vertexXAK5Jet);
   fChain->SetBranchAddress("vertexYAK5Jet", vertexYAK5Jet, &b_vertexYAK5Jet);
   fChain->SetBranchAddress("vertexZAK5Jet", vertexZAK5Jet, &b_vertexZAK5Jet);
   fChain->SetBranchAddress("emFracAK5Jet", emFracAK5Jet, &b_emFracAK5Jet);
   fChain->SetBranchAddress("hadFracAK5Jet", hadFracAK5Jet, &b_hadFracAK5Jet);
   fChain->SetBranchAddress("IdAK5Jet", IdAK5Jet, &b_IdAK5Jet);
   fChain->SetBranchAddress("nHitAK5Jet", nHitAK5Jet, &b_nHitAK5Jet);
   fChain->SetBranchAddress("nHit90AK5Jet", nHit90AK5Jet, &b_nHit90AK5Jet);
   fChain->SetBranchAddress("fHPDAK5Jet", fHPDAK5Jet, &b_fHPDAK5Jet);
   fChain->SetBranchAddress("covEtaEtaAK5Jet", covEtaEtaAK5Jet, &b_covEtaEtaAK5Jet);
   fChain->SetBranchAddress("covPhiPhiAK5Jet", covPhiPhiAK5Jet, &b_covPhiPhiAK5Jet);
   fChain->SetBranchAddress("fLSAK5Jet", fLSAK5Jet, &b_fLSAK5Jet);
   fChain->SetBranchAddress("fOOTAK5Jet", fOOTAK5Jet, &b_fOOTAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5Jet", combinedSecondaryVertexBJetTagsAK5Jet, &b_combinedSecondaryVertexBJetTagsAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5Jet", combinedSecondaryVertexMVABJetTagsAK5Jet, &b_combinedSecondaryVertexMVABJetTagsAK5Jet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5Jet", jetBProbabilityBJetTagsAK5Jet, &b_jetBProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5Jet", jetProbabilityBJetTagsAK5Jet, &b_jetProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5Jet", simpleSecondaryVertexHighEffBJetTagsAK5Jet, &b_simpleSecondaryVertexHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5Jet", simpleSecondaryVertexHighPurBJetTagsAK5Jet, &b_simpleSecondaryVertexHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5Jet", softMuonBJetTagsAK5Jet, &b_softMuonBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonByIP3dBJetTagsAK5Jet", softMuonByIP3dBJetTagsAK5Jet, &b_softMuonByIP3dBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonByPtBJetTagsAK5Jet", softMuonByPtBJetTagsAK5Jet, &b_softMuonByPtBJetTagsAK5Jet);
   fChain->SetBranchAddress("softElectronByIP3dBJetTagsAK5Jet", softElectronByIP3dBJetTagsAK5Jet, &b_softElectronByIP3dBJetTagsAK5Jet);
   fChain->SetBranchAddress("softElectronByPtBJetTagsAK5Jet", softElectronByPtBJetTagsAK5Jet, &b_softElectronByPtBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5Jet", trackCountingHighPurBJetTagsAK5Jet, &b_trackCountingHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5Jet", trackCountingHighEffBJetTagsAK5Jet, &b_trackCountingHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5Jet", combinedSecondaryVertexBJetTagsAK5Jet, &b_combinedSecondaryVertexBJetTagsAK5Jet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5Jet", combinedSecondaryVertexMVABJetTagsAK5Jet, &b_combinedSecondaryVertexMVABJetTagsAK5Jet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5Jet", jetBProbabilityBJetTagsAK5Jet, &b_jetBProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5Jet", jetProbabilityBJetTagsAK5Jet, &b_jetProbabilityBJetTagsAK5Jet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5Jet", simpleSecondaryVertexHighEffBJetTagsAK5Jet, &b_simpleSecondaryVertexHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5Jet", simpleSecondaryVertexHighPurBJetTagsAK5Jet, &b_simpleSecondaryVertexHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5Jet", softMuonBJetTagsAK5Jet, &b_softMuonBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5Jet", trackCountingHighPurBJetTagsAK5Jet, &b_trackCountingHighPurBJetTagsAK5Jet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5Jet", trackCountingHighEffBJetTagsAK5Jet, &b_trackCountingHighEffBJetTagsAK5Jet);
   fChain->SetBranchAddress("uncorrEnergyAK5Jet", uncorrEnergyAK5Jet, &b_uncorrEnergyAK5Jet);
   fChain->SetBranchAddress("nAK5PFPUcorrJet", &nAK5PFPUcorrJet, &b_nAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargeAK5PFPUcorrJet", chargeAK5PFPUcorrJet, &b_chargeAK5PFPUcorrJet);
   fChain->SetBranchAddress("energyAK5PFPUcorrJet", energyAK5PFPUcorrJet, &b_energyAK5PFPUcorrJet);
   fChain->SetBranchAddress("thetaAK5PFPUcorrJet", thetaAK5PFPUcorrJet, &b_thetaAK5PFPUcorrJet);
   fChain->SetBranchAddress("etaAK5PFPUcorrJet", etaAK5PFPUcorrJet, &b_etaAK5PFPUcorrJet);
   fChain->SetBranchAddress("phiAK5PFPUcorrJet", phiAK5PFPUcorrJet, &b_phiAK5PFPUcorrJet);
   fChain->SetBranchAddress("pxAK5PFPUcorrJet", pxAK5PFPUcorrJet, &b_pxAK5PFPUcorrJet);
   fChain->SetBranchAddress("pyAK5PFPUcorrJet", pyAK5PFPUcorrJet, &b_pyAK5PFPUcorrJet);
   fChain->SetBranchAddress("pzAK5PFPUcorrJet", pzAK5PFPUcorrJet, &b_pzAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexXAK5PFPUcorrJet", vertexXAK5PFPUcorrJet, &b_vertexXAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexYAK5PFPUcorrJet", vertexYAK5PFPUcorrJet, &b_vertexYAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexZAK5PFPUcorrJet", vertexZAK5PFPUcorrJet, &b_vertexZAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFPUcorrJet", chargedHadronEnergyAK5PFPUcorrJet, &b_chargedHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFPUcorrJet", neutralHadronEnergyAK5PFPUcorrJet, &b_neutralHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonEnergyAK5PFPUcorrJet", photonEnergyAK5PFPUcorrJet, &b_photonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronEnergyAK5PFPUcorrJet", electronEnergyAK5PFPUcorrJet, &b_electronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonEnergyAK5PFPUcorrJet", muonEnergyAK5PFPUcorrJet, &b_muonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronEnergyAK5PFPUcorrJet", HFHadronEnergyAK5PFPUcorrJet, &b_HFHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMEnergyAK5PFPUcorrJet", HFEMEnergyAK5PFPUcorrJet, &b_HFEMEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronMultiplicityAK5PFPUcorrJet", chargedHadronMultiplicityAK5PFPUcorrJet, &b_chargedHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronMultiplicityAK5PFPUcorrJet", neutralHadronMultiplicityAK5PFPUcorrJet, &b_neutralHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonMultiplicityAK5PFPUcorrJet", photonMultiplicityAK5PFPUcorrJet, &b_photonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronMultiplicityAK5PFPUcorrJet", electronMultiplicityAK5PFPUcorrJet, &b_electronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFPUcorrJet", muonMultiplicityAK5PFPUcorrJet, &b_muonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronMultiplicityAK5PFPUcorrJet", HFHadronMultiplicityAK5PFPUcorrJet, &b_HFHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMMultiplicityAK5PFPUcorrJet", HFEMMultiplicityAK5PFPUcorrJet, &b_HFEMMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFPUcorrJet", chargedEmEnergyAK5PFPUcorrJet, &b_chargedEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFPUcorrJet", neutralEmEnergyAK5PFPUcorrJet, &b_neutralEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5PFPUcorrJet", combinedSecondaryVertexBJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet", combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5PFPUcorrJet", jetBProbabilityBJetTagsAK5PFPUcorrJet, &b_jetBProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5PFPUcorrJet", jetProbabilityBJetTagsAK5PFPUcorrJet, &b_jetProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5PFPUcorrJet", softMuonBJetTagsAK5PFPUcorrJet, &b_softMuonBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByIP3dBJetTagsAK5PFPUcorrJet", softMuonByIP3dBJetTagsAK5PFPUcorrJet, &b_softMuonByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByPtBJetTagsAK5PFPUcorrJet", softMuonByPtBJetTagsAK5PFPUcorrJet, &b_softMuonByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronBJetTagsAK5PFPUcorrJet", softElectronBJetTagsAK5PFPUcorrJet, &b_softElectronBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByIP3dBJetTagsAK5PFPUcorrJet", softElectronByIP3dBJetTagsAK5PFPUcorrJet, &b_softElectronByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByPtBJetTagsAK5PFPUcorrJet", softElectronByPtBJetTagsAK5PFPUcorrJet, &b_softElectronByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5PFPUcorrJet", trackCountingHighPurBJetTagsAK5PFPUcorrJet, &b_trackCountingHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5PFPUcorrJet", trackCountingHighEffBJetTagsAK5PFPUcorrJet, &b_trackCountingHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("uncorrEnergyAK5PFPUcorrJet", uncorrEnergyAK5PFPUcorrJet, &b_uncorrEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("rmsCandAK5PFPUcorrJet", rmsCandAK5PFPUcorrJet, &b_rmsCandAK5PFPUcorrJet);
   fChain->SetBranchAddress("ptDAK5PFPUcorrJet", ptDAK5PFPUcorrJet, &b_ptDAK5PFPUcorrJet);
   fChain->SetBranchAddress("nAK5GenJet", &nAK5GenJet, &b_nAK5GenJet);
   fChain->SetBranchAddress("chargeAK5GenJet", chargeAK5GenJet, &b_chargeAK5GenJet);
   fChain->SetBranchAddress("energyAK5GenJet", energyAK5GenJet, &b_energyAK5GenJet);
   fChain->SetBranchAddress("thetaAK5GenJet", thetaAK5GenJet, &b_thetaAK5GenJet);
   fChain->SetBranchAddress("etaAK5GenJet", etaAK5GenJet, &b_etaAK5GenJet);
   fChain->SetBranchAddress("phiAK5GenJet", phiAK5GenJet, &b_phiAK5GenJet);
   fChain->SetBranchAddress("pxAK5GenJet", pxAK5GenJet, &b_pxAK5GenJet);
   fChain->SetBranchAddress("pyAK5GenJet", pyAK5GenJet, &b_pyAK5GenJet);
   fChain->SetBranchAddress("pzAK5GenJet", pzAK5GenJet, &b_pzAK5GenJet);
   fChain->SetBranchAddress("vertexXAK5GenJet", vertexXAK5GenJet, &b_vertexXAK5GenJet);
   fChain->SetBranchAddress("vertexYAK5GenJet", vertexYAK5GenJet, &b_vertexYAK5GenJet);
   fChain->SetBranchAddress("vertexZAK5GenJet", vertexZAK5GenJet, &b_vertexZAK5GenJet);
   fChain->SetBranchAddress("genPtHat", &genPtHat, &b_genPtHat);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
   fChain->SetBranchAddress("nTriggerObsPassing",&nTriggerObsPassing,&b_nTriggerObsPassing);
   fChain->SetBranchAddress("sizePassing",sizePassing,&b_sizePassing);
   fChain->SetBranchAddress("indexPassingPerPath",indexPassingPerPath,&b_indexPassingPerPath);
   fChain->SetBranchAddress("indexPassing",indexPassing,&b_indexPassing);
   fChain->SetBranchAddress("nTriggerObs",&nTriggerObs,&b_nTriggerObs);
   fChain->SetBranchAddress("triggerObsPt",triggerObsPt,&b_triggerObsPt);
   fChain->SetBranchAddress("triggerObsEta",triggerObsEta,&b_triggerObsEta);
   fChain->SetBranchAddress("triggerObsPhi",triggerObsPhi,&b_triggerObsPhi);
   fChain->SetBranchAddress("triggerObsMass",triggerObsMass,&b_triggerObsMass);

   Notify();
}

Bool_t HiggsBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HiggsBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HiggsBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HiggsBase_cxx
