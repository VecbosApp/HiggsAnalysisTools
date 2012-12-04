#include <RooArgList.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooHistFunc.h>
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"

#include "YieldMaker.h"
#include "CardTemplate.h"
#include "findAndReplace.h"
#include "yields.h"
#include "XSecProvider.h"
#include "FitSelection.hh"
#include "ConfigParser.hh"
#include "PdfSystematics.hh"

using namespace RooFit;
using namespace std;

class HiggsMassPointInfo {

public:
  HiggsMassPointInfo() {}
  ~HiggsMassPointInfo() {}

  float lumi;
  bool  do1D;
  bool  doFFT;
  std::string treeFolder;
  std::string hww2DShapesfilename;

  DataYieldMaker  ymaker_data;
  YieldMaker      ymaker_hi;
  YieldMaker      ymaker_qqww;
  YieldMaker      ymaker_ggww;
  YieldMaker      ymaker_top;
  YieldMaker      ymaker_dysf;
  YieldMaker      ymaker_dyof;
  YieldMaker      ymaker_others;
  WJetsYieldMaker ymaker_wj;
  XSecProvider xsecProvider;
  FitSelection sel;
  
 std::string getSignalCBMeanString(int ch) {
    stringstream fss;
    fss << "( ";  
    if (!doFFT) fss << "@0 + ";
    if (ch == of0j) fss << "26.86 - 0.24065*@0";
    if (ch == of1j) fss << "29.16 - 0.24410*@0";
    if (ch == sf0j) fss << "49.28 - 0.36044*@0";
    if (ch == sf1j) fss << "29.18 - 0.23906*@0";
    fss << " ) + @0*@1";
    return fss.str();
  }

  std::string getSignalCBSigmaString(int ch) {
    stringstream fss;
    fss << "( ";  
    if (ch == of0j) fss << "-1.878 + 0.1887*@0";
    if (ch == of1j) fss << " 3.154 + 0.1657*@0";
    if (ch == sf0j) fss << "-4.480 + 0.1816*@0";
    if (ch == sf1j) fss << "-0.771 + 0.1839*@0";
    fss << " ) * (1+@1)";
    return fss.str();
  }

  std::string getSignalCBAlphaString(int ch) {
    stringstream fss;
    if (ch == of0j) fss << "9.572";
    if (ch == of1j) fss << "8.718";
    if (ch == sf0j) fss << "8.304";
    if (ch == sf1j) fss << "4.050";
    return fss.str();
  }

  std::string getSignalCBNString(int ch) {
    stringstream fss;
    if (ch == of0j) fss << "16.9 - 0.0637*@0";
    if (ch == of1j) fss << "15.8 - 0.00064*@0";
    if (ch == sf0j) fss << "-13.2 + 0.1470*@0";
    if (ch == sf1j) fss << "17.7 - 0.0595*@0";
    return fss.str();
  }

  void createCard(float mass, float mrMin, float mrMax, int ch) {

    ConfigParser cp("config/mrfit_of_hcp.txt");

    std::string tevstr = "_8TeV"; 

    std::string chstr;
    if (ch == of0j) chstr = "of_0j";
    if (ch == of1j) chstr = "of_1j";
    if (ch == sf0j) chstr = "sf_0j";
    if (ch == sf1j) chstr = "sf_1j";

    stringstream mass_str_ss;
    mass_str_ss << mass;
    std::string mass_str = mass_str_ss.str();
        
    std::cout << "Creating datacard for " << mass_str << " GeV mass point and channel " << chstr << " ... " << std::endl;
    
    std::string card_name   = do1D ? (std::string("card_1D_m")+mass_str+"_8TeV_") : (std::string("card_2D_m")+mass_str+"_8TeV_");
    card_name += chstr;
    std::string workspace = card_name+"_workspace.root";

    ScaleFactors sf(ch);

    float yield_data   = ymaker_data   .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax);
    float yield_qqww   = ymaker_qqww   .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax) * sf.getWW() * lumi;
    float yield_ggww   = ymaker_ggww   .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax) * sf.getWW() * lumi;
    float yield_top    = ymaker_top    .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax) * sf.getTop() * lumi;
    float yield_dy     = 0.0; 
    if(ch==of0j || ch==of1j) yield_dy = ymaker_dyof .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax);
    else                     yield_dy = ymaker_dysf .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax) * sf.getDY() * lumi;
    float yield_others = ymaker_others .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax) * lumi;
    float yield_wj     = ymaker_wj     .getYield(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax);

    std::string card   = createCardTemplate(ch, do1D, workspace.c_str());

    std::string binname;
    if (ch == of0j) binname = "of_0j";
    if (ch == of1j) binname = "of_1j";
    if (ch == sf0j) binname = "sf_0j";
    if (ch == sf1j) binname = "sf_1j";

    int jet = (ch == of0j || ch == sf0j) ? 0 : 1;
    card = findAndReplace(card, "GGH_PDF"         , xsecProvider.get8TeVggHPDFDown(mass),            xsecProvider.get8TeVggHPDFUp(mass));
    card = findAndReplace(card, "VBF_PDF"         , xsecProvider.get8TeVVBFPDFDown(mass),            xsecProvider.get8TeVVBFPDFUp(mass));
    card = findAndReplace(card, "GGH_QCD"         , xsecProvider.get8TeVggHExclQCDDown(mass),        xsecProvider.get8TeVggHExclQCDUp(mass));
    card = findAndReplace(card, "GGH1_QCD"        , xsecProvider.get8TeVggH1inExclQCDDown(mass,jet), xsecProvider.get8TeVggH1inExclQCDUp(mass,jet));
    card = findAndReplace(card, "GGH2_QCD"        , xsecProvider.get8TeVggH2inExclQCDDown(mass,jet), xsecProvider.get8TeVggH2inExclQCDUp(mass,jet));
    card = findAndReplace(card, "VBF_QCD"         , xsecProvider.get8TeVVBFQCDDown(mass),            xsecProvider.get8TeVVBFQCDUp(mass));
    card = findAndReplace(card, "SIG_GGH_YIELD"   , 1);
    card = findAndReplace(card, "SIG_VBF_YIELD"   , 1);
    card = findAndReplace(card, "BKG_qqWW_YIELD"  , yield_qqww);
    card = findAndReplace(card, "BKG_ggWW_YIELD"  , yield_ggww);
    card = findAndReplace(card, "BKG_TOP_YIELD"   , yield_top);
    card = findAndReplace(card, "BKG_DY_YIELD"    , yield_dy);
    card = findAndReplace(card, "BKG_OTHERS_YIELD", yield_others);
    card = findAndReplace(card, "BKG_WJETS_YIELD" , yield_wj);
    card = findAndReplace(card, "BIN"             , binname);
    card = findAndReplace(card, "OBS"             , yield_data);
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_MC"),            cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","qqWWMCatNLONom","me","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_MC"),           cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","qqWWMCatNLONom","si","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-qcd"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMCatNLONom","qqWWMCatNLOUp","me","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-qcd"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMCatNLONom","qqWWMCatNLOUp","si","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-qcd"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMCatNLONom","qqWWMCatNLODown","me","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-qcd"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMCatNLONom","qqWWMCatNLODown","si","nominals","nominals"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-met"),       cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-met"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-met"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-met"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-j"));

    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_res-met"),       cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","res-met"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_res-met"),      cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","res-met"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","res-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","res-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","res-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","res-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaleup-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaleup-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaleup-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaleup-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaledn-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaledn-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaledn-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaledn-e"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaleup-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaleup-m"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaleup-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaleup-m"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaledn-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaledn-m"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaledn-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaledn-m"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaleup-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaleup-j"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaleup-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaleup-j"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_mean_err_scaledn-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","me","scaledn-j"));
    card = findAndReplace(card, ("BKG_GGWW_"+chstr+tevstr+"_sigma_err_scaledn-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"ggWW","si","scaledn-j"));

    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-met"),       cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-met"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-met"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-met"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_res-e"),         cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_res-e"),        cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","res-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-e"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-m"),    cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-m"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaleup-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaleup-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaleup-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaleup-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_mean_err_scaledn-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","me","scaledn-j"));
    card = findAndReplace(card, ("BKG_QQWW_"+chstr+tevstr+"_sigma_err_scaledn-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"qqWWMadgraph","si","scaledn-j"));

    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_mean_err_res-met"),       cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","res-met"));
    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_sigma_err_res-met"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","res-met"));
    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_mean_err_scaleup-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaleup-j"));
    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_sigma_err_scaleup-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaleup-j"));
    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_mean_err_scaledn-j"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaledn-j"));
    card = findAndReplace(card, ("BKG_TOP_"+chstr+tevstr+"_sigma_err_scaledn-j"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaledn-j"));

    card = findAndReplace(card, ("BKG_WJ_"+chstr+tevstr+"_mean_err_fakerateup"),       cp.getRelUncertainty(getStringFitChannel(ch),"wjets","me","fakerateup"));
    card = findAndReplace(card, ("BKG_WJ_"+chstr+tevstr+"_sigma_err_fakerateup"),      cp.getRelUncertainty(getStringFitChannel(ch),"wjets","si","fakerateup"));
    card = findAndReplace(card, ("BKG_WJ_"+chstr+tevstr+"_mean_err_fakeratedn"),       cp.getRelUncertainty(getStringFitChannel(ch),"wjets","me","fakeratedn"));
    card = findAndReplace(card, ("BKG_WJ_"+chstr+tevstr+"_sigma_err_fakeratedn"),      cp.getRelUncertainty(getStringFitChannel(ch),"wjets","si","fakeratedn"));

    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_mean_err_res-met"),       cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","res-met"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_sigma_err_res-met"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","res-met"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_mean_err_scaleup-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaleup-e"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_sigma_err_scaleup-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaleup-e"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_mean_err_scaledn-e"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaledn-e"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_sigma_err_scaledn-e"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaledn-e"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_mean_err_scaleup-m"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaleup-m"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_sigma_err_scaleup-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaleup-m"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_mean_err_scaledn-m"),      cp.getRelUncertainty(getStringFitChannel(ch),"Top","me","scaledn-m"));
    card = findAndReplace(card, ("BKG_OTHERS_"+chstr+tevstr+"_sigma_err_scaledn-m"),     cp.getRelUncertainty(getStringFitChannel(ch),"Top","si","scaledn-m"));


    ofstream file;
    file.open ((card_name +".txt").c_str());
    file << card;
    file.close();
    
    RooWorkspace w("w", "");
    
    RooRealVar CMS_ww2l_dphi ("CMS_ww2l_dphi" , "#Delta #Phi" ,   0,       TMath::Pi(), "");
    RooRealVar CMS_ww2l_mr_1D("CMS_ww2l_mr_1D", "M_{R}",          mrMin,   mrMax,       "GeV/c^{2}");
    RooRealVar weight("weight","weight",1.0,-1000.,1000.);
    if (doFFT) CMS_ww2l_mr_1D.setBins(100000, "fft");
        
    if (do1D) {
      RooArgSet argset_obs(CMS_ww2l_mr_1D, "argset_obs");
      RooDataSet data_obs("data_obs", "data_obs", argset_obs);
      
      ymaker_data.getDataSet1D(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax, data_obs, CMS_ww2l_mr_1D, weight);
    
      w.import(data_obs);
    }

    else {
      RooArgSet argset_obs(CMS_ww2l_mr_1D, CMS_ww2l_dphi, "argset_obs");
      RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
      ymaker_data.getDataSet2D(ch, mrMin, mrMax, sel.dphimin, sel.dphimax, sel.mtmin, sel.mtmax, data_obs, CMS_ww2l_mr_1D, CMS_ww2l_dphi, weight);

      w.import(data_obs);
    }
    

    ///////////////////// Define parameters //////////////////////////////////

    float qqWWme = 0.;
    float qqWWsi = 0.;
    float ggWWme = 0.;
    float ggWWsi = 0.;
    float Topme = 0.;
    float Topsi = 0.;
    float DYofme1 = 0.;
    float DYofsi1 = 0.;
    float DYofme2 = 0.;
    float DYofsi2 = 0.;
    float DYoffrac = 0.;
    float DYsfa0  = 0.;
    float DYsfa1  = 0.;
    float DYsfa2  = 0.;
    float DYsfa3  = 0.;
    float DYsfa4  = 0.;
    float DYsfa5  = 0.;
    float DYsfa6  = 0.;
    float DYsfa7  = 0.;
    float DYsfa8  = 0.;
    float DYsfa9  = 0.;
    float DYsfa10 = 0.;
    float DYsfa11 = 0.;
    float DYsfa12 = 0.;
    float DYsfa13 = 0.;
    float WJme    = 0.;
    float WJsi    = 0.;
    float Otofme  = 0.;
    float Otofsi  = 0.;
    float Otsfa0  = 0.;
    float Otsfa1  = 0.;
    float Otsfa2  = 0.;
    float Otsfa3  = 0.;
    float Otsfa4  = 0.;
    float Otsfa5  = 0.;
    float Otsfa6  = 0.;
    float Otsfa7  = 0.;
    float Otsfa8  = 0.;
    float Otsfa9  = 0.;
    float Otsfa10 = 0.;
    float Otsfa11 = 0.;
    float Otsfa12 = 0.;
    float Otsfa13 = 0.;

    if(ch == of0j || ch == of1j) {
      qqWWme = cp.getLandau(getStringFitChannel(ch),"qqWWMadgraph","nominals")[0];
      qqWWsi = cp.getLandau(getStringFitChannel(ch),"qqWWMadgraph","nominals")[1];
      
      ggWWme = cp.getLandau(getStringFitChannel(ch),"ggWW","nominals")[0];
      ggWWsi = cp.getLandau(getStringFitChannel(ch),"ggWW","nominals")[1];
      
      Topme = cp.getLandau(getStringFitChannel(ch),"Top","nominals")[0];
      Topsi = cp.getLandau(getStringFitChannel(ch),"Top","nominals")[1];
      
      DYofme1 = cp.getDoubleGaussian("of0j","embeddedtt","nominals")[0];
      DYofsi1 = cp.getDoubleGaussian("of0j","embeddedtt","nominals")[1];
      DYofme2 = cp.getDoubleGaussian("of0j","embeddedtt","nominals")[2];
      DYofsi2 = cp.getDoubleGaussian("of0j","embeddedtt","nominals")[3];
      DYoffrac = cp.getDoubleGaussian("of0j","embeddedtt","nominals")[4];
      
      WJme = cp.getLandau(getStringFitChannel(ch),"wjets","nominals")[0];
      WJsi = cp.getLandau(getStringFitChannel(ch),"wjets","nominals")[1]; 
      
      Otofme = cp.getLandau(getStringFitChannel(ch),"Ot","nominals")[0];
      Otofsi = cp.getLandau(getStringFitChannel(ch),"Ot","nominals")[1];
    }

    else if(ch == sf0j) {
      qqWWme = 166.184;
      qqWWsi = 31.4472;

      ggWWme = 165.034;
      ggWWsi = 26.5768;

      Topme = 198.501;
      Topsi = 43.3945;

      DYsfa0  = 84.7481;
      DYsfa1  = 13.0831;
      DYsfa2  = 56.1839;
      DYsfa3  = 0.0362868;
      DYsfa4  = 200.527;
      DYsfa5  = 17.7071;
      DYsfa6  = 13.9596;
      DYsfa7  = 0.0518413;
      DYsfa8  = 5.53392;
      DYsfa9  = 0.91199;
      DYsfa10  = 107.116;
      DYsfa11  = -8.36242;
      DYsfa12  = 4064.28;
      DYsfa13  = 0.00377819;

      WJme = 126.356;
      WJsi = 20.1992;

      Otsfa0  = 87.1441;
      Otsfa1  = 8.19755;
      Otsfa2  = 58.9242;
      Otsfa3  = 0.0224312;
      Otsfa4  = 207.433;
      Otsfa5  = 19.2934;
      Otsfa6  = 11.915;
      Otsfa7  = 0.166626;
      Otsfa8  = 65.7923;
      Otsfa9  = 0.0639726;
      Otsfa10  = 122.923;
      Otsfa11  = 10.95;
      Otsfa12  = 932.008;
      Otsfa13  = 1.43236e-08;
    }

    else if(ch == sf1j) {
      qqWWme = 173.467;
      qqWWsi = 37.7352;
      
      ggWWme = 168.458;
      ggWWsi = 30.4483;

      Topme = 185.053;
      Topsi = 40.3746;

      DYsfa0  = 79.0588;
      DYsfa1  = 30.0091;
      DYsfa2  = 84.2549;
      DYsfa3  = 0.0155325;
      DYsfa4  = 197.92;
      DYsfa5  = 18.9961;
      DYsfa6  = 24.912;
      DYsfa7  = 0.098307;
      DYsfa8  = 11.9076;
      DYsfa9  = 0.199545;
      DYsfa10  = 21.3897;
      DYsfa11  = -18.2451;
      DYsfa12  = 754.799;
      DYsfa13  = 0.0958732;

      WJme = 135.239;
      WJsi = 26.0509;

      Otsfa0  = 192.072;
      Otsfa1  = 10.1597;
      Otsfa2  = 56.8007;
      Otsfa3  = 0.701599;
      Otsfa4  = 228.752;
      Otsfa5  = 87.2529;
      Otsfa6  = 58.8586;
      Otsfa7  = 0.18201;
      Otsfa8  = 5.15327;
      Otsfa9  = 0.318942;
      Otsfa10  = 86.3412;
      Otsfa11  = 19.9938;
      Otsfa12  = 213.419;
      Otsfa13  = 0.0378662;
    }


    stringstream lumiss;
    lumiss << lumi;
    std::string lumistr = lumiss.str(); 

    WWSystematics qqwwsyst("bkg_qqww");
    RooArgList qqww_mean_al  = qqwwsyst.getParSystematics("mean",chstr,tevstr,qqWWme);
    RooArgList qqww_sigma_al = qqwwsyst.getParSystematics("sigma",chstr,tevstr,qqWWsi); 
    RooFormulaVar qqww_mean (("bkg_qqww_"+chstr+tevstr+"_mean" ).c_str(), qqwwsyst.getFormulaSyst().c_str(), qqww_mean_al);
    RooFormulaVar qqww_sigma(("bkg_qqww_"+chstr+tevstr+"_sigma").c_str(), qqwwsyst.getFormulaSyst().c_str(), qqww_sigma_al);

    WWSystematics ggwwsyst("bkg_ggww");
    RooArgList ggww_mean_al  = ggwwsyst.getParSystematics("mean",chstr,tevstr,ggWWme);
    RooArgList ggww_sigma_al = ggwwsyst.getParSystematics("sigma",chstr,tevstr,ggWWsi);
    RooFormulaVar ggww_mean (("bkg_ggww_"+chstr+tevstr+"_mean" ).c_str(), ggwwsyst.getFormulaSyst().c_str(), ggww_mean_al);
    RooFormulaVar ggww_sigma(("bkg_ggww_"+chstr+tevstr+"_sigma").c_str(), ggwwsyst.getFormulaSyst().c_str(), ggww_sigma_al);

    WWSystematics topsyst("bkg_top");
    RooArgList top_mean_al  = topsyst.getParSystematics("mean",chstr,tevstr,Topme);
    RooArgList top_sigma_al = topsyst.getParSystematics("sigma",chstr,tevstr,Topsi);
    RooFormulaVar top_mean (("bkg_top_"+chstr+tevstr+"_mean" ).c_str(), topsyst.getFormulaSyst().c_str(), top_mean_al);
    RooFormulaVar top_sigma(("bkg_top_"+chstr+tevstr+"_sigma").c_str(), topsyst.getFormulaSyst().c_str(), top_sigma_al);

    RooRealVar dyof_mean1 (("bkg_dy_"+chstr+tevstr+"_mean1" ).c_str(), "", DYofme1);
    RooRealVar dyof_sigma1(("bkg_dy_"+chstr+tevstr+"_sigma1").c_str(), "", DYofsi1);
    RooRealVar dyof_mean2 (("bkg_dy_"+chstr+tevstr+"_mean2" ).c_str(), "", DYofme2);
    RooRealVar dyof_sigma2(("bkg_dy_"+chstr+tevstr+"_sigma2").c_str(), "", DYofsi2);
    RooRealVar dyof_frac  (("bkg_dy_"+chstr+tevstr+"_frac").c_str(), "", DYoffrac);

    WJetsSystematics wjsyst("bkg_wj");
    RooArgList wj_mean_al  = wjsyst.getParSystematics("mean",chstr,tevstr,WJme);
    RooArgList wj_sigma_al = wjsyst.getParSystematics("sigma",chstr,tevstr,WJsi);    
    RooFormulaVar wj_mean (("bkg_wj_"+chstr+tevstr+"_mean" ).c_str(), wjsyst.getFormulaSyst().c_str(), wj_mean_al);
    RooFormulaVar wj_sigma(("bkg_wj_"+chstr+tevstr+"_sigma").c_str(), wjsyst.getFormulaSyst().c_str(), wj_sigma_al);

    WWSystematics otherssyst("bkg_others");
    RooArgList others_mean_al  = otherssyst.getParSystematics("mean",chstr,tevstr,Otofme);
    RooArgList others_sigma_al = otherssyst.getParSystematics("sigma",chstr,tevstr,Otofsi);    
    RooFormulaVar othersof_mean (("bkg_others_"+chstr+tevstr+"_mean" ).c_str(), otherssyst.getFormulaSyst().c_str(), others_mean_al);
    RooFormulaVar othersof_sigma(("bkg_others_"+chstr+tevstr+"_sigma").c_str(), otherssyst.getFormulaSyst().c_str(), others_sigma_al);

    RooRealVar dysf_a0 (("bkg_dy_"+chstr+tevstr+"_a0" ).c_str(), "", DYsfa0);
    RooRealVar dysf_a1 (("bkg_dy_"+chstr+tevstr+"_a1" ).c_str(), "", DYsfa1);
    RooRealVar dysf_a2 (("bkg_dy_"+chstr+tevstr+"_a2" ).c_str(), "", DYsfa2);
    RooRealVar dysf_a3 (("bkg_dy_"+chstr+tevstr+"_a3" ).c_str(), "", DYsfa3);
    RooRealVar dysf_a4 (("bkg_dy_"+chstr+tevstr+"_a4" ).c_str(), "", DYsfa4);
    RooRealVar dysf_a5 (("bkg_dy_"+chstr+tevstr+"_a5" ).c_str(), "", DYsfa5);
    RooRealVar dysf_a6 (("bkg_dy_"+chstr+tevstr+"_a6" ).c_str(), "", DYsfa6);
    RooRealVar dysf_a7 (("bkg_dy_"+chstr+tevstr+"_a7" ).c_str(), "", DYsfa7);
    RooRealVar dysf_a8 (("bkg_dy_"+chstr+tevstr+"_a8" ).c_str(), "", DYsfa8);
    RooRealVar dysf_a9 (("bkg_dy_"+chstr+tevstr+"_a9" ).c_str(), "", DYsfa9);
    RooRealVar dysf_a10 (("bkg_dy_"+chstr+tevstr+"_a10" ).c_str(), "", DYsfa10);
    RooRealVar dysf_a11 (("bkg_dy_"+chstr+tevstr+"_a11" ).c_str(), "", DYsfa11);
    RooRealVar dysf_a12 (("bkg_dy_"+chstr+tevstr+"_a12" ).c_str(), "", DYsfa12);
    RooRealVar dysf_a13 (("bkg_dy_"+chstr+tevstr+"_a13" ).c_str(), "", DYsfa13);
	
    RooRealVar otherssf_a0 (("bkg_others_"+chstr+tevstr+"_a0" ).c_str(), "", Otsfa0);
    RooRealVar otherssf_a1 (("bkg_others_"+chstr+tevstr+"_a1" ).c_str(), "", Otsfa1);
    RooRealVar otherssf_a2 (("bkg_others_"+chstr+tevstr+"_a2" ).c_str(), "", Otsfa2);
    RooRealVar otherssf_a3 (("bkg_others_"+chstr+tevstr+"_a3" ).c_str(), "", Otsfa3);
    RooRealVar otherssf_a4 (("bkg_others_"+chstr+tevstr+"_a4" ).c_str(), "", Otsfa4);
    RooRealVar otherssf_a5 (("bkg_others_"+chstr+tevstr+"_a5" ).c_str(), "", Otsfa5);
    RooRealVar otherssf_a6 (("bkg_others_"+chstr+tevstr+"_a6" ).c_str(), "", Otsfa6);
    RooRealVar otherssf_a7 (("bkg_others_"+chstr+tevstr+"_a7" ).c_str(), "", Otsfa7);
    RooRealVar otherssf_a8 (("bkg_others_"+chstr+tevstr+"_a8" ).c_str(), "", Otsfa8);
    RooRealVar otherssf_a9 (("bkg_others_"+chstr+tevstr+"_a9" ).c_str(), "", Otsfa9);
    RooRealVar otherssf_a10 (("bkg_others_"+chstr+tevstr+"_a10" ).c_str(), "", Otsfa10);
    RooRealVar otherssf_a11 (("bkg_others_"+chstr+tevstr+"_a11" ).c_str(), "", Otsfa11);
    RooRealVar otherssf_a12 (("bkg_others_"+chstr+tevstr+"_a12" ).c_str(), "", Otsfa12);
    RooRealVar otherssf_a13 (("bkg_others_"+chstr+tevstr+"_a13" ).c_str(), "", Otsfa13);
	
    RooRealVar sig_mean_err    (("sig_"+chstr+"_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
    RooRealVar sig_sigma_err   (("sig_"+chstr+"_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);
    
    RooRealVar masshiggs       ("MH", "", mass);

    InterpolatedYields ggHy(ch,ggH);
    InterpolatedYields qqHy(ch,vbfH);

    RooRealVar ggh_gamma_BW    (("sig_ggh_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
    RooFormulaVar ggh_mean_CB  (("sig_ggh_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (ch).c_str() , RooArgList(masshiggs, sig_mean_err));
    RooFormulaVar ggh_sigma_CB (("sig_ggh_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
    RooFormulaVar ggh_alpha    (("sig_ggh_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(ch).c_str() , RooArgList(masshiggs));
    RooFormulaVar ggh_n        (("sig_ggh_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(ch).c_str()     , RooArgList(masshiggs));
    RooFormulaVar ggh_norm     ("ggH_norm"                                    , (ggHy.getInterpolatedYieldString(lumi)).c_str()          , RooArgList(masshiggs));
    
    RooRealVar vbf_gamma_BW    (("sig_vbf_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
    RooFormulaVar vbf_mean_CB  (("sig_vbf_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (ch).c_str() , RooArgList(masshiggs, sig_mean_err));
    RooFormulaVar vbf_sigma_CB (("sig_vbf_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
    RooFormulaVar vbf_alpha    (("sig_vbf_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(ch).c_str() , RooArgList(masshiggs));
    RooFormulaVar vbf_n        (("sig_vbf_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(ch).c_str()     , RooArgList(masshiggs));
    RooFormulaVar vbf_norm     ("qqH_norm"                                    , (qqHy.getInterpolatedYieldString(lumi)).c_str()          , RooArgList(masshiggs));
      
    /////////// Set parameters to constant //////////////////

    dyof_mean1  .setConstant(kTRUE);
    dyof_sigma1 .setConstant(kTRUE);
    dyof_mean2  .setConstant(kTRUE);
    dyof_sigma2 .setConstant(kTRUE);
    dyof_frac   .setConstant(kTRUE);

    dysf_a0  .setConstant(kTRUE);
    dysf_a1  .setConstant(kTRUE);
    dysf_a2  .setConstant(kTRUE);
    dysf_a3  .setConstant(kTRUE);
    dysf_a4  .setConstant(kTRUE);
    dysf_a5  .setConstant(kTRUE);
    dysf_a6  .setConstant(kTRUE);
    dysf_a7  .setConstant(kTRUE);
    dysf_a8  .setConstant(kTRUE);
    dysf_a9  .setConstant(kTRUE);
    dysf_a10 .setConstant(kTRUE);
    dysf_a11 .setConstant(kTRUE);
    dysf_a12 .setConstant(kTRUE);
    dysf_a13 .setConstant(kTRUE);
    
    otherssf_a0  .setConstant(kTRUE);
    otherssf_a1  .setConstant(kTRUE);
    otherssf_a2  .setConstant(kTRUE);
    otherssf_a3  .setConstant(kTRUE);
    otherssf_a4  .setConstant(kTRUE);
    otherssf_a5  .setConstant(kTRUE);
    otherssf_a6  .setConstant(kTRUE);
    otherssf_a7  .setConstant(kTRUE);
    otherssf_a8  .setConstant(kTRUE);
    otherssf_a9  .setConstant(kTRUE);
    otherssf_a10 .setConstant(kTRUE);
    otherssf_a11 .setConstant(kTRUE);
    otherssf_a12 .setConstant(kTRUE);
    otherssf_a13 .setConstant(kTRUE);

    masshiggs     .setConstant(kTRUE);
    sig_mean_err  .setConstant(kTRUE);
    sig_sigma_err .setConstant(kTRUE);
    ggh_gamma_BW  .setConstant(kTRUE);
    vbf_gamma_BW  .setConstant(kTRUE);
    
    ////////////////// Define the PDFs /////////////////////////////////

    const char* bkg_qqww_pdf_name   = do1D ? "bkg_qqww"   : "bkg_qqww_1D" ;
    const char* bkg_ggww_pdf_name   = do1D ? "bkg_ggww"   : "bkg_ggww_1D" ;
    const char* bkg_top_pdf_name    = do1D ? "bkg_top"    : "bkg_top_1D" ;
    const char* bkg_dy_pdf_name     = do1D ? "bkg_dy"     : "bkg_dy_1D" ;
    const char* bkg_wj_pdf_name     = do1D ? "bkg_wj"     : "bkg_wj_1D" ;
    const char* bkg_others_pdf_name = do1D ? "bkg_others" : "bkg_others_1D" ;
    const char* sig_ggH_pdf_name    = do1D ? "ggH"        : "ggH_1D"  ;
    const char* sig_qqH_pdf_name    = do1D ? "qqH"        : "qqH_1D"  ;

    RooLandau *bkg_qqww_pdf = new RooLandau(bkg_qqww_pdf_name,"",CMS_ww2l_mr_1D,qqww_mean,qqww_sigma);

    RooLandau *bkg_ggww_pdf = new RooLandau(bkg_ggww_pdf_name,"",CMS_ww2l_mr_1D,ggww_mean,ggww_sigma);

    RooLandau *bkg_top_pdf = new RooLandau(bkg_top_pdf_name,"",CMS_ww2l_mr_1D,top_mean,top_sigma);

    RooLandau *bkg_wj_pdf = new RooLandau(bkg_wj_pdf_name,"",CMS_ww2l_mr_1D,wj_mean,wj_sigma);

    RooAbsPdf *bkg_dy_pdf, *bkg_others_pdf;
    bkg_dy_pdf=bkg_others_pdf=0;
    if(ch == of0j || ch == of1j) {
      RooGaussian *dygauss1 = new RooGaussian("dygauss1","dygauss1",CMS_ww2l_mr_1D,dyof_mean1,dyof_sigma1);
      RooGaussian *dygauss2 = new RooGaussian("dygauss2","dygauss2",CMS_ww2l_mr_1D,dyof_mean2,dyof_sigma2);
      bkg_dy_pdf = new RooAddPdf(bkg_dy_pdf_name,"",*dygauss1,*dygauss2,dyof_frac);
      bkg_others_pdf = new RooLandau(bkg_others_pdf_name,"",CMS_ww2l_mr_1D,othersof_mean,othersof_sigma);
    }
    if(ch == sf0j || ch == sf1j) {
      bkg_dy_pdf = new RooqqZZPdf_v2(bkg_dy_pdf_name,"",CMS_ww2l_mr_1D,
				     dysf_a0,
				     dysf_a1,
				     dysf_a2,
				     dysf_a3,
				     dysf_a4,
				     dysf_a5,
				     dysf_a6,
				     dysf_a7,
				     dysf_a8,
				     dysf_a9,
				     dysf_a10,
				     dysf_a11,
				     dysf_a12,
				     dysf_a13);

      bkg_others_pdf = new RooqqZZPdf_v2(bkg_others_pdf_name,"",CMS_ww2l_mr_1D,
					 otherssf_a0,
					 otherssf_a1,
					 otherssf_a2,
					 otherssf_a3,
					 otherssf_a4,
					 otherssf_a5,
					 otherssf_a6,
					 otherssf_a7,
					 otherssf_a8,
					 otherssf_a9,
					 otherssf_a10,
					 otherssf_a11,
					 otherssf_a12,
					 otherssf_a13);
    }

    RooCBShape      signalCB_ggH  (doFFT ? "signalCB_ggH" : sig_ggH_pdf_name, "", CMS_ww2l_mr_1D, ggh_mean_CB,ggh_sigma_CB,ggh_alpha,ggh_n);
    RooBreitWigner  signalBW_ggH   ("signalBW_ggH", "", CMS_ww2l_mr_1D, masshiggs,ggh_gamma_BW);
    RooFFTConvPdf*  sig_ggH_pdf = new RooFFTConvPdf(sig_ggH_pdf_name, "", CMS_ww2l_mr_1D, signalBW_ggH,signalCB_ggH,2);
    sig_ggH_pdf->setBufferFraction(0.2);

    RooCBShape      signalCB_qqH   (doFFT? "signalCB_qqH" : sig_qqH_pdf_name, "", CMS_ww2l_mr_1D, vbf_mean_CB,vbf_sigma_CB,vbf_alpha,vbf_n);
    RooBreitWigner  signalBW_qqH   ("signalBW_qqH", "", CMS_ww2l_mr_1D, masshiggs,vbf_gamma_BW);
    RooFFTConvPdf*  sig_qqH_pdf = new RooFFTConvPdf(sig_qqH_pdf_name, "", CMS_ww2l_mr_1D, signalBW_qqH,signalCB_qqH,2);
    sig_qqH_pdf->setBufferFraction(0.2);

    w.import(ggh_norm);
    w.import(vbf_norm);

    if (do1D) { 
      w.import(*bkg_qqww_pdf);
      w.import(*bkg_ggww_pdf);
      w.import(*bkg_top_pdf);
      w.import(*bkg_dy_pdf);
      w.import(*bkg_wj_pdf);
      w.import(*bkg_others_pdf);
      if(doFFT) {
        w.import(*sig_ggH_pdf);
        w.import(*sig_qqH_pdf);
      } else {
        w.import(signalCB_ggH);
        w.import(signalCB_qqH);
      }
    } 

    else {
      TFile *shapes2Dfile = TFile::Open(hww2DShapesfilename.c_str());
      RooArgList v2dList(CMS_ww2l_mr_1D, CMS_ww2l_dphi);
      RooArgSet  v2dSet (CMS_ww2l_mr_1D, CMS_ww2l_dphi);
      
      TH2F* dphishape_qqww = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_qqww_"+chstr).c_str()));
      TH2F* dphishape_ggww = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_ggww_"+chstr).c_str()));
      TH2F* dphishape_top = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_top_"+chstr).c_str()));
      TH2F* dphishape_dy = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_dy_"+chstr).c_str()));
      TH2F* dphishape_wj = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_wj_"+chstr).c_str()));
      TH2F* dphishape_others = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_others_"+chstr).c_str()));
      TH2F* dphishape_sig = (TH2F*)(shapes2Dfile->Get(("hist2D_sig_"+chstr).c_str()));

      // fake for the moment
      // TH2F* dphishape_ww_up = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_ww_"+chstr).c_str()));
      // TH2F* dphishape_top = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_top_"+chstr).c_str()));
      // TH2F* dphishape_dy = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_dy_"+chstr).c_str()));
      // TH2F* dphishape_wj = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_wj_"+chstr).c_str()));
      // TH2F* dphishape_others = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_others_"+chstr).c_str()));
      // TH2F* dphishape_sig = (TH2F*)(shapes2Dfile->Get(("hist2D_bkg_sig_"+chstr).c_str()));

      RooDataHist rhist_qqww   (("rhist_qqww_" +chstr+tevstr).c_str(), "", v2dList, dphishape_qqww);
      RooDataHist rhist_ggww   (("rhist_ggww_" +chstr+tevstr).c_str(), "", v2dList, dphishape_ggww);
      RooDataHist rhist_top    (("rhist_top_" +chstr+tevstr).c_str(), "", v2dList, dphishape_top);
      RooDataHist rhist_dy     (("rhist_dy_" +chstr+tevstr).c_str(), "", v2dList, dphishape_dy);
      RooDataHist rhist_wj     (("rhist_wj_" +chstr+tevstr).c_str(), "", v2dList, dphishape_wj);
      RooDataHist rhist_others (("rhist_others_" +chstr+tevstr).c_str(), "", v2dList, dphishape_others);
      RooDataHist rhist_ggH    (("rhist_ggH_" +chstr+tevstr).c_str(), "", v2dList, dphishape_sig);
      RooDataHist rhist_qqH    (("rhist_qqH_" +chstr+tevstr).c_str(), "", v2dList, dphishape_sig);
	    
      RooHistPdf rpdf_qqww   (("bkg_qqww_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_qqww);
      RooHistPdf rpdf_ggww   (("bkg_ggww_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_ggww);
      RooHistPdf rpdf_top    (("bkg_top_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_top);
      RooHistPdf rpdf_dy     (("bkg_dy_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_dy);
      RooHistPdf rpdf_wj     (("bkg_wj_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_wj);
      RooHistPdf rpdf_others (("bkg_others_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_wj);
      RooHistPdf rpdf_ggH    (("bkg_ggH_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_ggH);
      RooHistPdf rpdf_qqH    (("bkg_qqH_dphi2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_qqH);
	    
      // will be used for syst
      // RooRealVar CMS_ww2l_bkg("CMS_ww2l_bkgDPHI" ,"" ,0,-10,10); 
	
      FastVerticalInterpHistPdf2D plpdf_qqww   (("bkg_qqww_FVIHP_" +chstr+tevstr).c_str(),   "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_qqww)   ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_ggww   (("bkg_ggww_FVIHP_" +chstr+tevstr).c_str(),   "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_ggww)   ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_top    (("bkg_top_FVIHP_" +chstr+tevstr).c_str(),    "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_top)    ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_dy     (("bkg_dy_FVIHP_" +chstr+tevstr).c_str(),     "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_dy)     ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_wj     (("bkg_wj_FVIHP_" +chstr+tevstr).c_str(),     "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_wj)     ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_others (("bkg_others_FVIHP_" +chstr+tevstr).c_str(), "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_others) ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_ggH    (("sig_ggH_FVIHP_"  +chstr+tevstr).c_str(),   "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_ggH)    ,RooArgList()                ,1.0,1);
      FastVerticalInterpHistPdf2D plpdf_qqH    (("sig_qqH_FVIHP_"  +chstr+tevstr).c_str(),   "",CMS_ww2l_mr_1D,CMS_ww2l_dphi,true,RooArgList(rpdf_qqH)    ,RooArgList()                ,1.0,1);

      RooProdPdf bkg_qqww_pdf_2D   ("bkg_qqww"   , "", *bkg_qqww_pdf   ,Conditional(plpdf_qqww   , RooArgSet(CMS_ww2l_dphi))); 
      RooProdPdf bkg_ggww_pdf_2D   ("bkg_ggww"   , "", *bkg_ggww_pdf   ,Conditional(plpdf_ggww   , RooArgSet(CMS_ww2l_dphi))); 
      RooProdPdf bkg_top_pdf_2D    ("bkg_top"    , "", *bkg_top_pdf    ,Conditional(plpdf_top    , RooArgSet(CMS_ww2l_dphi))); 
      RooProdPdf bkg_dy_pdf_2D     ("bkg_dy"     , "", *bkg_dy_pdf     ,Conditional(plpdf_dy     , RooArgSet(CMS_ww2l_dphi))); 
      RooProdPdf bkg_wj_pdf_2D     ("bkg_wj"     , "", *bkg_wj_pdf     ,Conditional(plpdf_wj     , RooArgSet(CMS_ww2l_dphi))); 
      RooProdPdf bkg_others_pdf_2D ("bkg_others" , "", *bkg_others_pdf ,Conditional(plpdf_others , RooArgSet(CMS_ww2l_dphi))); 
	
      if(doFFT) {
        RooProdPdf sig_ggH_pdf_2D  ("ggH", "", *sig_ggH_pdf   ,Conditional(plpdf_ggH  , RooArgSet(CMS_ww2l_dphi))); 
        RooProdPdf sig_qqH_pdf_2D  ("qqH", "", *sig_qqH_pdf   ,Conditional(plpdf_qqH  , RooArgSet(CMS_ww2l_dphi))); 
        
        w.import(bkg_qqww_pdf_2D); 
        w.import(bkg_ggww_pdf_2D); 
        w.import(bkg_top_pdf_2D); 
        w.import(bkg_dy_pdf_2D); 
        w.import(bkg_wj_pdf_2D); 
        w.import(bkg_others_pdf_2D); 
        w.import(sig_ggH_pdf_2D); 
        w.import(sig_qqH_pdf_2D); 
      } else {
        RooProdPdf sig_ggH_pdf_2D  ("ggH", "", signalCB_ggH   ,Conditional(plpdf_ggH  , RooArgSet(CMS_ww2l_dphi))); 
        RooProdPdf sig_qqH_pdf_2D  ("qqH", "", signalCB_qqH   ,Conditional(plpdf_qqH  , RooArgSet(CMS_ww2l_dphi))); 
        
        w.import(bkg_qqww_pdf_2D); 
        w.import(bkg_ggww_pdf_2D); 
        w.import(bkg_top_pdf_2D); 
        w.import(bkg_dy_pdf_2D); 
        w.import(bkg_wj_pdf_2D); 
        w.import(bkg_others_pdf_2D); 
        w.import(sig_ggH_pdf_2D); 
        w.import(sig_qqH_pdf_2D); 
      }

    }

    w.writeToFile(workspace.c_str());

  }

};

void createWorkspace() {

  HiggsMassPointInfo hmpi8;
  hmpi8.lumi = 5.26;
  hmpi8.do1D = true;
  hmpi8.doFFT = false;
  hmpi8.treeFolder = "latinos_tree_skim_of/";
  hmpi8.hww2DShapesfilename = "hww2DShapes.root";
  hmpi8.xsecProvider.initXsec();
  hmpi8.xsecProvider.initQCDScale();
  hmpi8.xsecProvider.initPDF();
  hmpi8.xsecProvider.initJetBinFracs();

  /*
  hmpi8.ymaker_data   .fill(hmpi8.treeFolder+"/data/latino_RunA_892pbinv.root");
  hmpi8.ymaker_data   .fill(hmpi8.treeFolder+"/data/latino_RunB_4404pbinv.root");
  hmpi8.ymaker_data   .fill(hmpi8.treeFolder+"/data/latino_RunC_6807pbinv.root");
  hmpi8.ymaker_qqww   .fill(hmpi8.treeFolder+"/nominals/latino_000_WWJets2LMad.root");
  hmpi8.ymaker_ggww   .fill(hmpi8.treeFolder+"/nominals/latino_001_GluGluToWWTo4L.root");
  hmpi8.ymaker_top    .fill(hmpi8.treeFolder+"/nominals/latino_019_TTTo2L2Nu2B.root");
  hmpi8.ymaker_top    .fill(hmpi8.treeFolder+"/nominals/latino_011_TtWFullDR.root");
  hmpi8.ymaker_top    .fill(hmpi8.treeFolder+"/nominals/latino_012_TbartWFullDR.root");
  hmpi8.ymaker_dysf   .fill(hmpi8.treeFolder+"/nominals/latino_037_DY50toLLMad.root");
  hmpi8.ymaker_dysf   .fill(hmpi8.treeFolder+"/nominals/latino_036_DY10toLLMad.root");
  hmpi8.ymaker_dyof   .fill(hmpi8.treeFolder+"/nominals/latino_RunABC_DYtt_8fb.root");
  hmpi8.ymaker_wj     .fill(hmpi8.treeFolder+"/wj/latino_RunABC_LooseLoose_skimww.root");
  hmpi8.ymaker_others .fill(hmpi8.treeFolder+"/nominals/latino_074_WZJetsMad.root");
  hmpi8.ymaker_others .fill(hmpi8.treeFolder+"/nominals/latino_075_ZZJetsMad.root");
  */

  // for (float i = 114.; i <= 180.; i += 1.) {
  for (float i = 125.; i <= 125.; i += 1.) {  
    for(int j = 0; j < 4; ++j) hmpi8.createCard(i, 50, 500, j);
  }

  hmpi8.do1D = false;
  // for (float i = 114.; i <= 180.; i += 1.) {
  for (float i = 125.; i <= 125.; i += 1.) {  
    for(int j = 0; j < 4; ++j) hmpi8.createCard(i, 50, 500, j);
  }

}
