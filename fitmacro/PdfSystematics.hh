#include <string>
#include <sstream>

#include <RooRealVar.h>
#include <RooArgList.h>

using namespace std;

class WWSystematics {
public:
  WWSystematics(string process) { _process=process; }
  ~WWSystematics() {}

  string getFormulaSyst() {
    stringstream fss;
    fss << "@0*(1+@1+@2+@3+@4+@5+@6+@7+@8+@9+@10+@11+@12";
    if(_process.find("qqww")!=string::npos) fss << "+@13+@14+@15";
    fss << ")";
    return fss.str();
  }

  RooArgList getParSystematics(string par, string chstr, string tevstr, float centralval) {
    RooRealVar *WW_par = new RooRealVar((_process+"_"+chstr+tevstr+"_"+par).c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_MC = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_MC").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaleup = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaleup").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaledn = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaledn").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_resmet = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_resmet").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_rese = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_rese").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_resmu = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_resmu").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaleupe = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaleupe").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaledne = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaledne").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaleupmu = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaleupmu").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scalednmu = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scalednmu").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scaleupj = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scaleupj").c_str(), "", 0., -10., 10.);
    RooRealVar *WW_par_err_scalednj = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_scalednj").c_str(), "", 0., -10., 10.);

    WW_par->setVal(centralval);
    WW_par->setConstant(kTRUE);

    RooArgList WW_par_err(*WW_par,
                          *WW_par_err_resmet,    *WW_par_err_rese, *WW_par_err_resmu, 
                          *WW_par_err_scaleupe,  *WW_par_err_scaledne,
                          *WW_par_err_scaleupmu, *WW_par_err_scalednmu);
    WW_par_err.add(*WW_par_err_scaleupj); 
    WW_par_err.add(*WW_par_err_scalednj);

    if(_process.find("qqww")!=string::npos) {
      WW_par_err.add(*WW_par_err_MC);
      WW_par_err.add(*WW_par_err_scaleup);
      WW_par_err.add(*WW_par_err_scaledn);
    }
    return WW_par_err;
  }

private:
  string _process;
};


class WJetsSystematics {
public:
  WJetsSystematics(string process) { _process=process; }
  ~WJetsSystematics() {}

  string getFormulaSyst() {
    stringstream fss;
    fss << "@0*(1+@1+@2)";
    return fss.str();
  }

  RooArgList getParSystematics(string par, string chstr, string tevstr, float centralval) {
    RooRealVar *WJets_par = new RooRealVar((_process+"_"+chstr+tevstr+"_"+par).c_str(), "", 0., -10., 10.);
    RooRealVar *WJets_par_err_fakerateup = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_fakerateup").c_str(), "", 0., -10., 10.);
    RooRealVar *WJets_par_err_fakeratedn = new RooRealVar((_process+"_"+chstr+tevstr+par+"_err_fakeratedn").c_str(), "", 0., -10., 10.);

    WJets_par->setVal(centralval);
    WJets_par->setConstant(kTRUE);

    RooArgList WJets_par_err(*WJets_par,
                             *WJets_par_err_fakerateup, *WJets_par_err_fakeratedn);
    return WJets_par_err;
  }

private:
  string _process;
};

