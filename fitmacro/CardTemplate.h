#ifndef CARDTEMPLATE_H
#define CARDTEMPLATE_H

#include "YieldMaker.h"

std::string createCardTemplate(int channel, bool do1D, std::string workspacefilename) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += std::string("shapes * * ") + workspacefilename + " w:$PROCESS\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN           BIN            BIN           BIN           BIN           BIN           BIN            BIN\n";
        card += "process ggH           qqH            bkg_qqww      bkg_ggww      bkg_top       bkg_wj        bkg_others     bkg_dy\n";
        card += "process -5            -4             1             2             3             4             5              6\n";
        card += "rate    SIG_GGH_YIELD SIG_VBF_YIELD  BKG_qqWW_YIELD  BKG_ggWW_YIELD  BKG_TOP_YIELD  BKG_WJETS_YIELD  BKG_OTHERS_YIELD  BKG_DY_YIELD\n";
        card += "------------\n";
        //                                            ggH      qqH     qqWW     ggWW     top      wj       ot       dy
        card += "lumi_8TeV                 lnN        1.05     1.05    -        -        -        -        1.05     -\n";
        card += "pdf_gg                    lnN        GGH_PDF  -       -        1.040    -        -        -        -\n";
        card += "pdf_qqbar                 lnN        -        VBF_PDF 1.040    -        -        -        1.040    -\n";
        card += "pdf_hww2l2nu_accept       lnN        1.02     1.02    -        -        -        -        -        -\n";
      if(channel == of0j || channel == sf0j) {
        card += "UEPS                      lnN        0.943    -       -        -        -        -        -        -\n";  
      }
      else if(channel == of1j || channel == sf1j) {
        card += "UEPS                      lnN        1.084    -       -        -        -        -        -        -\n";  
      }
        card += "QCDscale_ggH              lnN        GGH_QCD  -       -        -        -        -        -        -\n";
        card += "QCDscale_ggH1in           lnN        GGH1_QCD  -       -        -        -        -        -        -\n";
        card += "QCDscale_ggH2in           lnN        GGH2_QCD  -       -       -        -        -        -        -\n";
        card += "QCDscale_qqH              lnN        -        VBF_PDF -        -        -        -        -        -\n";
        card += "QCDscale_WW               lnN        -        -       1.00     -        -        -        -        -\n";
        card += "QCDscale_WW1in            lnN        -        -       1.00     -        -        -        -        -\n";
        card += "QCDscale_WW2in            lnN        -        -       1.00     -        -        -        -        -\n";
        card += "QCDscale_VV               lnN        -        -       -        -        -        -        1.04     -\n";
        card += "QCDscale_ggWW             lnN        -        -       -        1.30     -        -        -        -\n";
        card += "QCDscale_ggH_ACCEPT       lnN        1.02     -       -        -        -        -        -        -\n";
        card += "QCDscale_qqH_ACCEPT       lnN        -        1.02    -        -        -        -        -        -\n";
        card += "BRhiggs_WW2l              lnN        1.02     1.02    -        -        -        -        -        -\n";
        card += "FakeRate                  lnN        -        -       -        -        -        1.360    -        -\n";
        card += "CMS_eff_m                 lnN        1.015    1.015   -        -        -        -        1.015    -\n";
        card += "CMS_eff_e                 lnN        1.01     1.01    -        -        -        -        1.010    -\n";
        //                                            ggH      qqH     qqWW     ggWW     top      wj       ot       dy
    if(channel == of0j || channel == sf0j) {
        card += "CMS_hww_0j_ttbar          lnN        -         -       -       -        1.210    -        -        -\n";
        card += "CMS_hww_0j_WW             lnN        -         -       1.080   1.080    -        -        -        -\n";
    } 
    else if(channel == of1j || channel == sf1j) {
        card += "CMS_hww_1j_ttbar          lnN        -         -       -       -        1.061    -        -        -\n";
        card += "CMS_hww_1j_WW             lnN        -         -       1.130   1.130    -        -        -        -\n";
    } 
    if (channel == of0j) {
        card += "sig_of_0j_mean_err_7TeV   param      0        0.005                     \n";
        card += "sig_of_0j_sigma_err_7TeV  param      0        0.3                       \n";
        card += "CMS_hwwof_0j_Z            lnN        -         -       -        -       -         -       -        1.10\n";
    }
    else if (channel == of1j) {
        card += "sig_of_1j_mean_err_7TeV      param      0        0.004                     \n";
        card += "sig_of_1j_sigma_err_7TeV     param      0        0.3                       \n";
        card += "CMS_hwwof_1j_Z            lnN        -         -       -        -       -         -       -        1.10\n";
    }
    else if (channel == sf0j) {
        card += "sig_sf_0j_mean_err_7TeV   param      0        0.005                      \n";
        card += "sig_sf_0j_sigma_err_7TeV  param      0        0.3                        \n";
        card += "CMS_hwwsf_0j_Z            lnN        -         -       -        -       -         -       -        1.45\n";
    }
    else if (channel == sf1j) {
        card += "sig_sf_1j_mean_err_7TeV   param      0        0.005                      \n";
        card += "sig_sf_1j_sigma_err_7TeV  param      0        0.3                        \n";
        card += "CMS_hwwsf_1j_Z            lnN        -         -       -        -       -         -       -        1.185\n";
    }
    if (!do1D) {
      // TODO: shape syst on deltaphi
      //        card += "CMS_ww2l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

#endif
