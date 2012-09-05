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
  RooRealVar rrweight   ;
  RooRealVar rrprocess  ;
  RooArgSet argset      ;
  RooDataSet dataset    ;
  
 public :        

 YieldMaker():
  rrchannel  (RooRealVar("channel",  "channel",   0., 10.)), 
    rrmr     (RooRealVar("mr",       "mr",      0., 10000000.)),
    rrD      (RooRealVar("dphillr",  "dphillr",      0., TMath::Pi())),
    rrweight (RooRealVar("weight",   "weight",    -100000., 10000000.)),
    rrprocess(RooRealVar("dataset",  "dataset",   0., 10000000.)),
    argset(RooArgSet(rrchannel, rrmr, rrD, rrweight, rrprocess, "argset")),
    dataset(RooDataSet("dataset", "dataset", argset))
      {}

  float getYield(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float procmin=-100000, float procmax=100000) {
    
    float yield = 0.0;
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && proc>=procmin && proc<=procmax) yield += weight;
    }

    return yield;

  }

  float getCount(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float procmin=-100000, float procmax=100000) {
    
    float yield = 0.0;
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && proc>=procmin && proc<=procmax) yield += 1.0;
    }

    return yield;

  }

  void get1DHist(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float procmin, float procmax, TH1* hist) {
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && proc>=procmin && proc<=procmax) hist->Fill(mr,weight);
    }

  }


  void get2DHist(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float procmin, float procmax, TH2* hist) {
    
    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float D      = dataset.get(i)->getRealValue("dphillr");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && D>dphimin && D<dphimax && proc>=procmin && proc<=procmax) hist->Fill(mr,D,weight);
    }

  }


  RooDataSet getFitDataSet(int channel, float mrmin, float mrmax, float dphimin, float dphimax, float procmin=-100000, float procmax=100000) {
    RooRealVar mr("mr",   "mr",            100, 50, 1000,      "GeV/c^{2}");
    RooRealVar D("dphillr",   "dphillr",   100, 0,  TMath::Pi());
    RooRealVar w("weight", "weight", 0.,  -10.,  10000.);
    RooArgSet  aset(mr, D, w, "aset");
    RooDataSet dset("dataset","", aset);
    

    for (int i = 0; i < dataset.numEntries(); i++) {
      float m      = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float weight = dataset.get(i)->getRealValue("weight");
      float ch     = dataset.get(i)->getRealValue("channel");
      float proc   = dataset.get(i)->getRealValue("dataset");

      if (channel == (int)ch && m>mrmin && m<mrmax && dphi>dphimin && dphi<dphimax && proc>=procmin && proc<=procmax) {
	aset.setRealValue("mr", m);
	aset.setRealValue("D",  dphi);
	aset.setRealValue("weight", weight);
	dset.add(aset);
      }
    }

    return dset;

  }

  void fill(std::string filepath) {

    cout << "\tYieldMaker. Filling dataset from: " << filepath << endl;

    TFile *file = TFile::Open(filepath.c_str());
    TTree *tree = (TTree*)file->Get("latinoFitSkim");

    float mr         = 0.0;
    float dphi       = 0.0;
    float baseweight = 0.0;
    float effweight  = 0.0;
    float puweight   = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("mr",      &mr);
    tree->SetBranchAddress("dphillr", &dphi);
    tree->SetBranchAddress("baseW",   &baseweight);
    tree->SetBranchAddress("puW",     &puweight);
    tree->SetBranchAddress("effW",    &effweight);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);
    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphi);
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
    TTree *tree = (TTree*)file->Get("latinoFitSkim");

    float mr         = 0.0;
    float dphi       = 0.0;
    float fakeweight = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("mr",      &mr);
    tree->SetBranchAddress("dphillr", &dphi);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);
    tree->SetBranchAddress("fake2W",  &fakeweight);

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphi);
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
    TTree *tree = (TTree*)file->Get("latinoFitSkim");

    float mr         = 0.0;
    float dphi       = 0.0;
    float ch         = 0.0;
    float proc       = 0.0;
    float njet       = 0.0;

    tree->SetBranchAddress("mr",      &mr);
    tree->SetBranchAddress("dphillr", &dphi);
    tree->SetBranchAddress("channel", &ch);
    tree->SetBranchAddress("dataset", &proc);
    tree->SetBranchAddress("njet",    &njet);

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      argset.setRealValue("mr",      mr);
      argset.setRealValue("dphillr", dphi);
      float channel = getFitChannel(ch,njet);
      argset.setRealValue("channel", channel);
      argset.setRealValue("dataset", proc);
      argset.setRealValue("weight",  1.0);
      if(channel<0) continue;
      dataset.add(argset);
    }

  }

  void getDataSet1D(int channel, float mrmin, float mrmax, float dphimin, float dphimax, RooDataSet& dset, RooRealVar& m) {

    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float ch     = dataset.get(i)->getRealValue("channel");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && dphi>dphimin && dphi<dphimax) {
	m.setVal(mr);
	RooArgSet aset(m, "argset_obs");
	dset.add(aset);
      }
    }
  }

  void getDataSet2D(int channel, float mrmin, float mrmax, float dphimin, float dphimax, RooDataSet& dset, RooRealVar& m, RooRealVar& D) {

    for (int i = 0; i < dataset.numEntries(); i++) {
      float mr     = dataset.get(i)->getRealValue("mr");
      float dphi   = dataset.get(i)->getRealValue("dphillr");
      float ch     = dataset.get(i)->getRealValue("channel");

      if (channel == (int)ch && mr>mrmin && mr<mrmax && dphi>dphimin && dphi<dphimax) {
	m.setVal(mr);
	D.setVal(dphi);
	RooArgSet aset(m, D, "argset_obs");
	dset.add(aset);
      }
    }
  }


};

#endif

