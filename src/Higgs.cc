#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <fstream>
#include <sstream>

#include "TString.h"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "HiggsAnalysisTools/include/Higgs.hh"

using namespace bits;

Higgs::Higgs(TTree *tree) : HiggsBase(tree)
{
  jsonFile = "";
  lastFile = "";

  //initialize the JES objects for calo and PF
  jecUnc_calo = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/GR_R_42_V19_AK5Calo_Uncertainty.txt");
  jecUnc_PF = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/GR_R_42_V19_AK5PF_Uncertainty.txt");

}

Higgs::~Higgs()
{
  // By this time, the destructor of HiggsBase has not yet been called.
  // This means that the tree has not yet been deleted.
  // So, we do nothing here.
}

void Higgs::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool Higgs::hasPassedHLT() {
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

void Higgs::setJsonGoodRunList(const std::string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void Higgs::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool Higgs::isGoodRunLS()
{
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //  std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //      std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
	return true;
    }
  return false;
}

bool Higgs::reloadTriggerMask(bool newVersion)
{
  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {   
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
            if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
              {
                triggerMask.push_back( indexHLT[i] ) ;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {

        std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
        TTree *treeCond;
        std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
        std::vector<std::string>  *nameHLT_;
        std::vector<unsigned int> *indexHLT_;

        //To get the pointers for the vectors
        nameHLT_=0;
        indexHLT_=0;

        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);

        std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++) 
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
                  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (int i=0;i<m_requiredTriggers.size();++i)
          std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
}

float Higgs::mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met) {
  float pTll = (pl1.Vect() + pl2.Vect()).Pt();
  float mll = (pl1 + pl2).M();
  float El = sqrt(pTll*pTll + mll*mll);
  float pTnu = met.Pt();
  float Enu = sqrt(pTnu*pTnu + mll*mll);
  float Ex = (pl1+pl2).X() + met.X();
  float Ey = (pl1+pl2).Y() + met.Y();
  float mnu = mll;

  return sqrt(mll*mll + mnu*mnu + 2*(El*Enu-Ex*Ex-Ey*Ey));
}

float Higgs::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  Utils anaUtils;
  bool inEB=anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB);
  measurements.pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
  if(inEB && fabs(etaEle[eleIndex])<1.0) measurements.subdet = 0;
  else if (inEB && fabs(etaEle[eleIndex])>=1.0) measurements.subdet = 1;
  else measurements.subdet = 2;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 pIn(pxGsfTrack[gsftrack],pyGsfTrack[gsftrack],pzGsfTrack[gsftrack]);
  measurements.OneOverEMinusOneOverP = 1./(eSuperClusterOverPEle[eleIndex]*pIn.Mag()) - 1./pIn.Mag();
  return lh.resultLog(measurements);
}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float Higgs::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float Higgs::SigmaiPiP(int electron) {
  float spp;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    spp = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      spp = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      spp = 999.;
    }
  }
  return spp;
}

bool Higgs::isPFJetID(float eta, float neutralHadFrac, float neutralEmFraction, int nConstituents, float chargedHadFraction, 
                      float chargedMultiplicity, float chargedEmFraction, int WP) {
  switch(WP) {
  case none:
    return true;
    break;
  case loose:
    if(neutralHadFrac>=0.99 || neutralEmFraction>=0.99 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case medium:
    if(neutralHadFrac>=0.95 || neutralEmFraction>=0.95 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case tight:
    if(neutralHadFrac>=0.90 || neutralEmFraction>=0.90 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  default:
    std::cout << "Jet::isPFJetID(nt WP). Requested wrong Working point. Available are loose, medium, tight." << std::endl;
    return false;
  }
  return true;
}

/// *****************************************
/// From C. Rogan, M. Pierini, M. Spiropulu 
/// *****************************************

//
// this is the 'new' MRstar
//
double Higgs::CalcMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  return temp;
}



//
// this is the 'new' MRstar, times 'gamma_{R*}' - I would recommend making 'R*' with this as 
// the denominator and 'M_{T}^{R}' as the numerator (the next function in this file)
//
double Higgs::CalcGammaMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar
  temp *= mygamma;

  return temp;
}

//
// This is 'M_{T}^{R}', the guy that should be used in the numerator of 'R' or 'R*'
//
double Higgs::CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){

  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;

  temp = sqrt(temp);

  return temp;
}

std::vector<int> Higgs::sortElectronsByPt(std::vector<int> electrons) {
  int tmp;
  int max;
  for(int i=0;i<(int)electrons.size();i++) {
      max = i;
      float maxelePt = GetPt(pxEle[i],pyEle[i]);
      for(int x=i; x<(int)electrons.size(); x++) {
        float xelePt = GetPt(pxEle[x],pyEle[x]);
        if(xelePt > maxelePt) {
          max = x;
        }
      }
      tmp = electrons[i];
      electrons[i] = electrons[max];
      electrons[max] = tmp;
  }
  return electrons;
}

std::vector<int> Higgs::sortMuonsByPt(std::vector<int> muons) {
  int tmp;
  int max;
  for(int i=0;i<(int)muons.size();i++) {
      max = i;
      float maxmuonPt = GetPt(pxMuon[i],pyMuon[i]);
      for(int x=i; x<(int)muons.size(); x++) {
        float xmuonPt = GetPt(pxMuon[x],pyMuon[x]);
        if(xmuonPt > maxmuonPt) {
          max = x;
        }
      }
      tmp = muons[i];
      muons[i] = muons[max];
      muons[max] = tmp;
  }
  return muons;
}

TLorentzVector Higgs::GetJESCorrected(TLorentzVector p4jet, const char *ScaleDirection) {

  float mass = p4jet.M();
  float ptUnscaled = p4jet.Pt();
  
  // estimate the uncertainty
  jecUnc_PF->setJetEta(p4jet.Eta());
  jecUnc_PF->setJetPt(ptUnscaled);

  int scaleEnergy = 0;
  if(TString(ScaleDirection).Contains("Up")) scaleEnergy = 1.0;
  if(TString(ScaleDirection).Contains("Down")) scaleEnergy = -1.0;

  // apply the uncertainty
  float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
  float p = pt/fabs(sin(p4jet.Theta()));
  float energy = sqrt(p*p+mass*mass);

  TLorentzVector p4Scaled;
  p4Scaled.SetPtEtaPhiE(pt,p4jet.Eta(),p4jet.Phi(),energy);

  return p4Scaled;
}

TVector3 Higgs::pfChargedMet(TVector3 lep1, TVector3 lep2) {

  float chMetP3x = pxPFChMet[0];
  float chMetP3y = pyPFChMet[0];

  // charged PF MET has been computed with all the PF cands (inverted -p)
  // first remove the contribution in dR = 0.1 to avoid double counting
  for(int i=0; i<nReducedPFCand; i++) {
    TVector3 pfCandP3(pxReducedPFCand[i],pyReducedPFCand[i],pzReducedPFCand[i]);
    if(pfCandP3.DeltaR(lep1)<=0.1 || pfCandP3.DeltaR(lep2)<=0.1) {
      chMetP3x += pxReducedPFCand[i];
      chMetP3y += pyReducedPFCand[i];
    }
  }

  // then add back the RECO leptons
  chMetP3x -= (lep1.Px() + lep2.Px());
  chMetP3y -= (lep1.Py() + lep2.Py());
  
  return TVector3(chMetP3x,chMetP3y,0.0);

}

std::string Higgs::getHLTPathForRun(int runN, std::string fullname) {
  TString fullName = TString(fullname.c_str());
  TObjArray* selectionTokens = fullName.Tokenize(":");
  if (selectionTokens->GetEntries()!=2) {
    std::cout << "Wrong trigger strings " << selectionTokens->GetEntries() << std::endl;
    return std::string("NOPATH");
  }
  TString RunRange =((TObjString*)(*selectionTokens)[0])->GetString();
  TString HLTPathName =((TObjString*)(*selectionTokens)[1])->GetString();
  
  TObjArray* runs = RunRange.Tokenize("-");
  if (runs->GetEntries()!=2) {
    std::cout << "Wrong trigger run range strings " << runs->GetEntries() << std::endl;
    return std::string("NOPATH");    
  }
  
  const char *minStr = (((TObjString*)(*runs)[0])->GetString()).Data();
  const char *maxStr = (((TObjString*)(*runs)[1])->GetString()).Data();

  int min = atoi(minStr);
  int max = atoi(maxStr);

  if(runN>=min && runN<=max) return std::string(HLTPathName.Data());
  else return std::string("NOPATH");
}
