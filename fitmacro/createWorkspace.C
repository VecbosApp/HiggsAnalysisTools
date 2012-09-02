#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>

#include "YieldMaker.h"
#include "CardTemplate.h"
#include "findAndReplace.h"

using namespace RooFit;
using namespace std;

enum channels { of0j, of1j, sf0j, sf1j };

struct HiggsMassPointInfo {

    float lumi;
    float dphiMin;
    float dphiMax;
    bool  do1D;
    std::string treeFolder;
    std::string melafilename;
    
    std::map<std::string, std::vector<float> > interpolationmap;

    DataYieldMaker  ymaker_dt;
    YieldMaker      ymaker_hi;
    YieldMaker      ymaker_ww;
    YieldMaker      ymaker_to;
    YieldMaker      ymaker_dy;
    YieldMaker      ymaker_ot;
    WJetsYieldMaker ymaker_wj;

 std::string getSignalCBMeanString(int ch) {
   stringstream fss;
   fss << "( ";  
   if (ch == of0j) "26.86 - 0.24065*@0";
   if (ch == of1j) "29.16 - 0.24410*@0";
   if (ch == sf0j) "49.28 - 0.36044*@0";
   if (ch == sf1j) "29.18 - 0.23906*@0";
   fss << " ) + @0*@1";
   return fss.str();
}

 std::string getSignalCBSigmaString(int ch) {
   stringstream fss;
   fss << "( ";  
   if (ch == of0j) "-1.878 + 0.1887*@0";
   if (ch == of1j) " 3.154 + 0.1657*@0";
   if (ch == sf0j) "-4.480 + 0.1816*@0";
   if (ch == sf1j) "-0.771 + 0.1839*@0";
   fss << " ) + @0*@1";
   return fss.str();
}

 std::string getSignalCBAlphaString(int ch) {
   stringstream fss;
   fss << "( ";  
   if (ch == of0j) "9.572";
   if (ch == of1j) "8.718";
   if (ch == sf0j) "8.304";
   if (ch == sf1j) "4.050";
   return fss.str();
}

 std::string getSignalCBNString(int ch) {
   stringstream fss;
   fss << "( ";  
   if (ch == of0j) "16.9 - 0.0637*@0";
   if (ch == of1j) "15.8 - 0.00064*@0";
   if (ch == sf0j) "-13.2 + 0.1470*@0";
   if (ch == sf1j) "17.7 - 0.0595*@0";
   return fss.str();
}

  void createCard(float mass, float mrMin, float mrMax, int ch) {

    std::string chstr;
    if (ch == of0j) chstr = "of_0j";
    if (ch == of1j) chstr = "of_1j";
    if (ch == sf0j) chstr = "sf_0j";
    if (ch == sf1j) chstr = "sf_1j";

    stringstream mass_str_ss;
    mass_str_ss << mass;
    std::string mass_str = mass_str_ss.str();
        
    std::cout << "Creating datacard for " << mass_str << " GeV mass point and channel " << chstr << " ... " << std::endl;
    
    std::string card_name   = do1D ? (std::string("card_1D_m")+mass_str+"_8TeV_") : (std::string("card_2D_m")+mass_str+"_8TeV");
    card_name += chstr;
    std::string workspace = card_name+"_workspace.root";

    float yield_dt = ymaker_dt.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_hi = ymaker_hi.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_ww = ymaker_ww.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_to = ymaker_to.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_dy = ymaker_dy.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_ot = ymaker_ot.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);
    float yield_wj = ymaker_wj.getYield(ch, mrMin, mrMax, dphiMin, dphiMax);

    std::string card   = createCardTemplate(ch, do1D, workspace.c_str());

    std::string binname;
    if (ch == of0j) binname = "of_0j";
    if (ch == of1j) binname = "of_1j";
    if (ch == sf0j) binname = "sf_0j";
    if (ch == sf1j) binname = "sf_1j";

    card = findAndReplace(card, "SIG_GGH_YIELD"   , hi);
    card = findAndReplace(card, "SIG_VBF_YIELD"   , 0.); // to be changed
    card = findAndReplace(card, "BKG_WW_YIELD"    , yield_ww);
    card = findAndReplace(card, "BKG_TOP_YIELD"   , yield_top);
    card = findAndReplace(card, "BKG_DY_YIELD"    , yield_dy);
    card = findAndReplace(card, "BKG_OTHERS_YIELD", yield_ot);
    card = findAndReplace(card, "BKG_WJETS_YIELD" , yield_wj);
    card = findAndReplace(card, "BIN"            , binname);
    card = findAndReplace(card, "OBS"            , yield_dt);

    ofstream file;
    file.open ((card_name +".txt").c_str());
    file << card;
    file.close();
    
    RooWorkspace w("w", "");
    
    RooRealVar CMS_ww2l_dphi ("CMS_ww2l_dphi" , "#Delta #Phi" ,   0,       TMath::Pi(), "");
    RooRealVar CMS_ww2l_mr_1D("CMS_ww2l_mr_1D", "M_{R}",          mrMin,   mrMax,       "GeV/c^{2}");
    if (doFFT) CMS_zz4l_mass_1D.setBins(100000, "fft");
        
    if (do1D) {
      RooArgSet argset_obs(CMS_ww2l_mr_1D, "argset_obs");
      RooDataSet data_obs("data_obs", "data_obs", argset_obs);
      
      ymaker_data.getDataSet1D(ch, mrMin, mrMax, dphiMin, dphiMax, data_obs, CMS_ww2l_mr_1D);
    
      w.import(data_obs);
    }

    else {
      RooArgSet argset_obs(CMS_ww2l_mr_1D, CMS_ww2l_dphi, "argset_obs");
      RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
      ymaker_data.getDataSet2D(ch, z1min, z2min, mrMin, mrMax, dphiMin, dphiMax, data_obs, CMS_ww2l_mr_1D, CMS_ww2l_dphi);

      w.import(data_obs);
    }
    

    ///////////////////// Define parameters //////////////////////////////////
    
    float WWme = 0.;
    float WWsi = 0.;
    float Topme = 0.;
    float Topsi = 0.;
    float DYofme = 0.;
    float DYofsi = 0.;
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
    float WJofme = 0.;
    float WJofsi = 0.;
    float WJsfa0  = 0.;
    float WJsfa1  = 0.;
    float WJsfa2  = 0.;
    float WJsfa3  = 0.;
    float WJsfa4  = 0.;
    float WJsfa5  = 0.;
    float WJsfa6  = 0.;
    float WJsfa7  = 0.;
    float WJsfa8  = 0.;
    float WJsfa9  = 0.;
    float WJsfa10 = 0.;
    float WJsfa11 = 0.;
    float WJsfa12 = 0.;
    float WJsfa13 = 0.;
    float Otofme = 0.;
    float Otofsi = 0.;
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

    if(ch == of0j) {
      WWme = 158.109;
      WWsi = 33.6625;
      
      Topme = 193.609;
      Topsi = 46.473;
      
      DYofme = 130.131;
      DYofsi = 33.3977;

      WJofme = 125.842;
      WJofsi = 23.9167;

      Otofme = 100.089;
      Otofsi = 19.2616;
    }

    else if(ch == of1j) {
      WWme = 170.113;
      WWsi = 38.5888;
      
      Topme = 183.549;
      Topsi = 42.8287;

      DYofme = 100.007;
      DYofsi = 53.6901;

      WJofme = 143.556;
      WJofsi = 29.4511;
      
      Otofme = 123.646;
      Otofsi = 27.2714;
    }

    else if(ch == sf0j) {
      WWme = 166.081;
      WWsi = 31.1662;

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

      WJsfa0  = 98.8649;
      WJsfa1  = 17.5007;
      WJsfa2  = 72.1404;
      WJsfa3  = 0.44289;
      WJsfa4  = 209.894;
      WJsfa5  = 52.3488;
      WJsfa6  = 7.55082;
      WJsfa7  = 0.241892;
      WJsfa8  = 0.0841174;
      WJsfa9  = 0.490682;
      WJsfa10  = 195.594;
      WJsfa11  = -0.68302;
      WJsfa12  = 9927.89;
      WJsfa13  = 0.00424511;

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
      WWme = 173.194;
      WWsi = 37.3955;

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

      WJsfa0  = 104.474;
      WJsfa1  = 24.6371;
      WJsfa2  = 82.8418;
      WJsfa3  = 0.0475532;
      WJsfa4  = 184.969;
      WJsfa5  = 10.9556;
      WJsfa6  = 25.4052;
      WJsfa7  = 0.0580383;
      WJsfa8  = 21.3376;
      WJsfa9  = 0.0290728;
      WJsfa10  = 48.1927;
      WJsfa11  = -1.69871;
      WJsfa12  = 9935.12;
      WJsfa13  = 0.326497;

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


    std::string tevstr = "_8TeV"; 
    stringstream lumiss;
    lumiss << lumi;
    std::string lumistr = lumiss.str(); 

    RooRealVar ww_mean (("bkg_ww_"+chstr+tevstr+"_mean" ).c_str(), "", WWme);
    RooRealVar ww_sigma(("bkg_ww_"+chstr+tevstr+"_sigma").c_str(), "", WWsi);

    RooRealVar top_mean (("bkg_top_"+chstr+tevstr+"_mean" ).c_str(), "", Topme);
    RooRealVar top_sigma(("bkg_top_"+chstr+tevstr+"_sigma").c_str(), "", Topsi);

    RooRealVar dyof_mean (("bkg_dy_"+chstr+tevstr+"_mean" ).c_str(), "", DYofme);
    RooRealVar dyof_sigma(("bkg_dy_"+chstr+tevstr+"_sigma").c_str(), "", DYofsi);

    RooRealVar wjof_mean (("bkg_wj_"+chstr+tevstr+"_mean" ).c_str(), "", WJofme);
    RooRealVar wjof_sigma(("bkg_wj_"+chstr+tevstr+"_sigma").c_str(), "", WJofsi);

    RooRealVar others_mean (("bkg_others_"+chstr+tevstr+"_mean" ).c_str(), "", Otofme);
    RooRealVar others_sigma(("bkg_others_"+chstr+tevstr+"_sigma").c_str(), "", Otofsi);

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
	
    RooRealVar wjsf_a0 (("bkg_wj_"+chstr+tevstr+"_a0" ).c_str(), "", WJsfa0);
    RooRealVar wjsf_a1 (("bkg_wj_"+chstr+tevstr+"_a1" ).c_str(), "", WJsfa1);
    RooRealVar wjsf_a2 (("bkg_wj_"+chstr+tevstr+"_a2" ).c_str(), "", WJsfa2);
    RooRealVar wjsf_a3 (("bkg_wj_"+chstr+tevstr+"_a3" ).c_str(), "", WJsfa3);
    RooRealVar wjsf_a4 (("bkg_wj_"+chstr+tevstr+"_a4" ).c_str(), "", WJsfa4);
    RooRealVar wjsf_a5 (("bkg_wj_"+chstr+tevstr+"_a5" ).c_str(), "", WJsfa5);
    RooRealVar wjsf_a6 (("bkg_wj_"+chstr+tevstr+"_a6" ).c_str(), "", WJsfa6);
    RooRealVar wjsf_a7 (("bkg_wj_"+chstr+tevstr+"_a7" ).c_str(), "", WJsfa7);
    RooRealVar wjsf_a8 (("bkg_wj_"+chstr+tevstr+"_a8" ).c_str(), "", WJsfa8);
    RooRealVar wjsf_a9 (("bkg_wj_"+chstr+tevstr+"_a9" ).c_str(), "", WJsfa9);
    RooRealVar wjsf_a10 (("bkg_wj_"+chstr+tevstr+"_a10" ).c_str(), "", WJsfa10);
    RooRealVar wjsf_a11 (("bkg_wj_"+chstr+tevstr+"_a11" ).c_str(), "", WJsfa11);
    RooRealVar wjsf_a12 (("bkg_wj_"+chstr+tevstr+"_a12" ).c_str(), "", WJsfa12);
    RooRealVar wjsf_a13 (("bkg_wj_"+chstr+tevstr+"_a13" ).c_str(), "", WJsfa13);
	
    
