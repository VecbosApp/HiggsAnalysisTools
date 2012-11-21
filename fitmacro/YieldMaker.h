#ifndef YIELDMAKER_H
#define YIELDMAKER_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include "HWWKinematics.hh"

enum channels { of0j, of1j, sf0j, sf1j };

float getFitChannel(float channel, float njet) {
  if((int)channel>=2 && (int)njet==0) return of0j;
  if((int)channel>=2 && (int)njet==1) return of1j;
  if((int)channel<2 && (int)njet==0) return sf0j;
  if((int)channel<2 && (int)njet==1) return sf1j;
  return -1;
}

class YieldMaker {

 protected:
  RooRealVar rrchannel  ;	
  RooRealVar rrmr     ;
  RooRealVar rrD     ;
  RooRealVar rrmt    ;
  RooRealVar rrweight   ;
  RooRealVar rrprocess  ;
  RooArgSet argset      ;
  RooDataSet dataset    ;
  
 public :        

 YieldMaker():
  rrchannel  (RooRealVar("channel",  "channel",   0., 10.)), 
    rrmr     (RooRealVar("mr",       "mr",      0., 10000000.)),
    rrD      (RooRealVar("dphillr",  "dphillr",      0., TMath::Pi())),
    rrmt     (RooRealVar("mt",       "mt",           0., 10000000.)),
    rrweight (RooRealVar("weight",   "weight",    -100000., 10000000.)),
    rrprocess(RooRealVar("dataset",  "dataset",   0., 10000000.)),
    argset(RooArgSet(rrchannel, rrmr, rrD, rrmt, rrweight, rrprocess, "argset")),
    dataset(RooDataSet("dataset", "dataset", argset))
      {}

   float getYield(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, float procmin=-100000, float procmax=100000) {
    
    float yield = 0.0;
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && mt>mtmin && mt<mtmax && proc>=procmin && proc<=procmax) yield += weight;
    }

    return yield;

  }

   float getCount(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, float procmin=-100000, float procmax=100000) {
    
    float yield = 0.0;
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && mt>mtmin && mt<mtmax && proc>=procmin && proc<=procmax) yield += 1.0;
    }

    return yield;

  }

   void get1DHist(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, float procmin, float procmax, TH1* hist) {
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && mt>mtmin && mt<mtmax && proc>=procmin && proc<=procmax) hist->Fill(mr,weight);
    }

  }


   void get2DHist(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, float procmin, float procmax, TH2* hist) {
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && mt>mtmin && mt<mtmax && proc>=procmin && proc<=procmax) hist->Fill(mr,D,weight);
    }

  }


   RooDataSet getFitDataSet(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, float procmin=-100000, float procmax=100000) {
    RooRealVar mr("mr",   "mr",            100, 50, 1000,      "GeV/c^{2}");
    RooRealVar D("dphillr",   "dphillr",   100, 0,  TMath::Pi());
    RooRealVar w("weight", "weight", 0.,  -10.,  10000.);
    RooArgSet  aset(mr, D, w, "aset");
    RooDataSet dset("dataset","", aset,RooFit::WeightVar("weight"));
    

    for (int i = 0; i < dataset.numEntries(); i++) {
      float m      = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && m>mrmin && m<mrmax && dphi>dphimin && dphi<dphimax && mt>mtmin && mt<mtmax && proc>=procmin && proc<=procmax) {
	aset.setRealValue("mr", m);
	aset.setRealValue("D",  dphi);
        aset.setRealValue("mt", mt);
	aset.setRealValue("weight", weight);
	dset.add(aset);
      }
    }

    return dset;

  }

   void getDataSet1D(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, RooDataSet& dset, RooRealVar& m, RooRealVar& w) {

    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float ch     = dataset.get(i)->getRealValue("channel");
      float weight = dataset.get(i)->getRealValue("weight");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && dphi>dphimin && dphi<dphimax && mt>mtmin && mt<mtmax) {
	m.setVal(mr);
        w.setVal(weight);
        RooArgSet aset(m, "argset_obs");
	dset.add(aset,weight);
      }
    }
  }

   void getDataSet2D(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float mtmin, float mtmax, RooDataSet& dset, RooRealVar& m, RooRealVar& D, RooRealVar& w) {

    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float mt     = dataset.get(i)->getRealValue("mt");
      float ch     = dataset.get(i)->getRealValue("channel");
      float weight = dataset.get(i)->getRealValue("weight");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && dphi>dphimin && dphi<dphimax && mt>mtmin && mt<mtmax) {
	m.setVal(mr);
	D.setVal(dphi);
        w.setVal(weight);
	RooArgSet aset(m, D, "argset_obs");
	dset.add(aset,weight);
      }
    }
  }


  void fill(std::string filepath) {

    cout << "\tYieldMaker. Filling dataset from: " << filepath << endl;

    TFile *file = TFile::Open(filepath.c_str());
    TTree *tree = (TTree*)file->Get("latino");

    float eta1       = 0.0;
    float eta2       = 0.0;
    float phi1       = 0.0;
    float phi2       = 0.0;
    float pt1        = 0.0;
    float pt2        = 0.0;
    float pfmet      = 0.0;
    float pfmetphi   = 0.0;
    float mth        = 0.0;
    float boostedMR  = 0.0;
    float baseweight = 0.0;
    float effweight  = 0.0;
    float puweight   = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("eta1",    &eta1);
    tree->SetBranchAddress("eta2",    &eta2);
    tree->SetBranchAddress("phi1",    &phi1);
    tree->SetBranchAddress("phi2",    &phi2);
    tree->SetBranchAddress("pt1",     &pt1);
    tree->SetBranchAddress("pt2",     &pt2);
    tree->SetBranchAddress("pfmet",   &pfmet);
    tree->SetBranchAddress("pfmetphi",&pfmetphi);
    tree->SetBranchAddress("mth",     &mth);
    tree->SetBranchAddress("boostedMR",&boostedMR);
    tree->SetBranchAddress("baseW",   &baseweight);
    tree->SetBranchAddress("puW",     &puweight);
    tree->SetBranchAddress("effW",    &effweight);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);
    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      float mr         = 0.0;
      float dphillr    = 0.0;

      int l1id,l2id;
      if(ch==0) { l1id=13; l2id=13; }
      if(ch==1) { l1id=11; l2id=11; }
      if(ch==2) { l1id=11; l2id=13; }
      if(ch==3) { l1id=13; l2id=11; }

      float elmass=0.51E-03;
      float mumass=105.7E-03;
      float l1m = (l1id==11) ? elmass : mumass;
      float l2m = (l2id==11) ? elmass : mumass;
      TLorentzVector l1,l2;
      TVector3 met;
      l1.SetPtEtaPhiM(pt1,eta1,phi1,l1m);
      l2.SetPtEtaPhiM(pt2,eta2,phi2,l2m);
      met.SetPtEtaPhi(pfmet,0.0,pfmetphi);

      HWWKinematics kine(l1,l2,met);
      mr = 2*kine.CalcMRNEW(); // met is the smeared one...
      // mr = 2*boostedMR;
      dphillr = fabs(kine.CalcDeltaPhiRFRAME());

      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphillr);
      argset.setRealValue("mt", mth);
      float channel = getFitChannel(ch,njet);
      argset.setRealValue("channel", channel);
      argset.setRealValue("dataset", proc);
      float weight = baseweight*effweight*puweight;
      argset.setRealValue("weight", weight);
      if(channel<0) continue;
      dataset.add(argset);
    }

  }

};



class WJetsYieldMaker : public YieldMaker {

 public :        

  WJetsYieldMaker():YieldMaker(){}

  void fill(std::string filepath) {

    cout << "\tWJetsYieldMaker. Filling dataset from: " << filepath << endl;

    TFile *file = TFile::Open(filepath.c_str());
    TTree *tree = (TTree*)file->Get("latino");

    float eta1       = 0.0;
    float eta2       = 0.0;
    float phi1       = 0.0;
    float phi2       = 0.0;
    float pt1        = 0.0;
    float pt2        = 0.0;
    float pfmet      = 0.0;
    float pfmetphi   = 0.0;
    float mth        = 0.0;
    float boostedMR  = 0.0;
    float fakeweight = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("eta1",    &eta1);
    tree->SetBranchAddress("eta2",    &eta2);
    tree->SetBranchAddress("phi1",    &phi1);
    tree->SetBranchAddress("phi2",    &phi2);
    tree->SetBranchAddress("pt1",     &pt1);
    tree->SetBranchAddress("pt2",     &pt2);
    tree->SetBranchAddress("pfmet",   &pfmet);
    tree->SetBranchAddress("pfmetphi",&pfmetphi);
    tree->SetBranchAddress("mth",     &mth);
    tree->SetBranchAddress("boostedMR",&boostedMR);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);
    tree->SetBranchAddress("fake2W",  &fakeweight);

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      float mr         = 0.0;
      float dphillr    = 0.0;

      int l1id,l2id;
      if(ch==0) { l1id=13; l2id=13; }
      if(ch==1) { l1id=11; l2id=11; }
      if(ch==2) { l1id=11; l2id=13; }
      if(ch==3) { l1id=13; l2id=11; }

      float elmass=0.51E-03;
      float mumass=105.7E-03;
      float l1m = (l1id==11) ? elmass : mumass;
      float l2m = (l2id==11) ? elmass : mumass;
      TLorentzVector l1,l2;
      TVector3 met;
      l1.SetPtEtaPhiM(pt1,eta1,phi1,l1m);
      l2.SetPtEtaPhiM(pt2,eta2,phi2,l2m);
      met.SetPtEtaPhi(pfmet,0.0,pfmetphi);

      HWWKinematics kine(l1,l2,met);
      mr = 2*kine.CalcMRNEW(); // met is the smeared one... 
      // mr = 2*boostedMR;
      dphillr = fabs(kine.CalcDeltaPhiRFRAME());

      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphillr);
      argset.setRealValue("mt", mth);
      float channel = getFitChannel(ch,njet);
      argset.setRealValue("channel", channel);
      argset.setRealValue("dataset", proc);
      float weight = fakeweight;
      argset.setRealValue("weight", weight);
      if(channel<0) continue;
      dataset.add(argset);
    }

  }


};


class DataYieldMaker : public YieldMaker {

 public :        

  DataYieldMaker():YieldMaker(){}

  void fill(std::string filepath) {

    cout << "\tDataYieldMaker. Filling dataset from: " << filepath << endl;

    TFile *file = TFile::Open(filepath.c_str());
    TTree *tree = (TTree*)file->Get("latino");

    float eta1       = 0.0;
    float eta2       = 0.0;
    float phi1       = 0.0;
    float phi2       = 0.0;
    float pt1        = 0.0;
    float pt2        = 0.0;
    float pfmet      = 0.0;
    float pfmetphi   = 0.0;
    float mth        = 0.0;
    float boostedMR  = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("eta1",    &eta1);
    tree->SetBranchAddress("eta2",    &eta2);
    tree->SetBranchAddress("phi1",    &phi1);
    tree->SetBranchAddress("phi2",    &phi2);
    tree->SetBranchAddress("pt1",     &pt1);
    tree->SetBranchAddress("pt2",     &pt2);
    tree->SetBranchAddress("pfmet",   &pfmet);
    tree->SetBranchAddress("pfmetphi",&pfmetphi);
    tree->SetBranchAddress("mth",     &mth);
    tree->SetBranchAddress("boostedMR",&boostedMR);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);

      float mr         = 0.0;
      float dphillr    = 0.0;

      int l1id,l2id;
      if(ch==0) { l1id=13; l2id=13; }
      if(ch==1) { l1id=11; l2id=11; }
      if(ch==2) { l1id=11; l2id=13; }
      if(ch==3) { l1id=13; l2id=11; }

      float elmass=0.51E-03;
      float mumass=105.7E-03;
      float l1m = (l1id==11) ? elmass : mumass;
      float l2m = (l2id==11) ? elmass : mumass;
      TLorentzVector l1,l2;
      TVector3 met;
      l1.SetPtEtaPhiM(pt1,eta1,phi1,l1m);
      l2.SetPtEtaPhiM(pt2,eta2,phi2,l2m);
      met.SetPtEtaPhi(pfmet,0.0,pfmetphi);

      HWWKinematics kine(l1,l2,met);
      mr = 2*kine.CalcMRNEW(); // met is the smeared one...
      // mr = 2*boostedMR;
      dphillr = fabs(kine.CalcDeltaPhiRFRAME());

      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphillr);
      argset.setRealValue("mt", mth);
      float channel = getFitChannel(ch,njet);
      argset.setRealValue("channel", channel);
      argset.setRealValue("dataset", proc);
      argset.setRealValue("weight",  1.0);
      if(channel<0) continue;
      dataset.add(argset);
    }

  }


};

#endif

