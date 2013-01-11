#ifndef CARDTEMPLATE_H
#define CARDTEMPLATE_H

#include "YieldMaker.h"

std::string createCardTemplate(int channel, bool do1D, std::string workspacefilename) {

    std::string chstr;
    if (channel == of0j) chstr = "of_0j";
    if (channel == of1j) chstr = "of_1j";
    if (channel == sf0j) chstr = "sf_0j";
    if (channel == sf1j) chstr = "sf_1j";

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
        card += "bin     BIN         BIN          BIN         BIN         BIN         BIN         BIN         BIN          BIN          BIN\n";
        card += "process ggH           vbfH           wzttH         bkg_qqww      bkg_ggww      bkg_top       bkg_wj        bkg_others     bkg_dy         bkg_wgstar\n";
        card += "process -1            0              -2            1             2             3             4             5              6              7\n";
        card += "rate    SIG_GGH_YIELD SIG_VBF_YIELD  SIG_WZTT_YIELD  BKG_qqWW_YIELD  BKG_ggWW_YIELD  BKG_TOP_YIELD  BKG_WJETS_YIELD  BKG_OTHERS_YIELD  BKG_DY_YIELD  BKG_WGSTAR_YIELD\n";
        card += "------------\n";
        //                                            ggH      qqH     wzttH    qqWW     ggWW     top      wj       ot       dy       wgstar
        card += "lumi_8TeV                 lnN        1.044    1.044   1.044    -        -        -        -        1.044    -        1.044\n";
        card += "pdf_gg                    lnN        GGH_PDF  -       -        -        1.040    -        -        -        -        -\n";
        card += "pdf_qqbar                 lnN        -        VBF_PDF VBF_PDF  1.000    -        -        -        1.040    -        1.040\n";
        card += "pdf_hww2l2nu_accept       lnN        1.02     1.02    -        -        -        -        -        -        -        -\n";
      if(channel == of0j || channel == sf0j) {
        card += "UEPS                      lnN        0.943    -       -        -        -        -        -        -        -        -\n";  
      }
      else if(channel == of1j || channel == sf1j) {
        card += "UEPS                      lnN        1.084    -       -        -        -        -        -        -        -        -\n";  
      }
        card += "QCDscale_ggH              lnN        GGH_QCD  -       -        -        -        -        -        -        -        -\n";
        card += "QCDscale_ggH1in           lnN        GGH1_QCD  -      -        -        -        -        -        -        -        -\n";
        card += "QCDscale_ggH2in           lnN        GGH2_QCD  -      -        -        -        -        -        -        -        -\n";
        card += "QCDscale_qqH              lnN        -        VBF_PDF -        -        -        -        -        -        -        -\n";
        card += "QCDscale_WW               lnN        -        -       -        1.00     -        -        -        -        -        -\n";
        card += "QCDscale_WW1in            lnN        -        -       -        1.00     -        -        -        -        -        -\n";
        card += "QCDscale_WW2in            lnN        -        -       -        1.00     -        -        -        -        -        -\n";
        card += "QCDscale_VV               lnN        -        -       -        -        -        -        -        1.030    -        -\n";
        card += "QCDscale_V                lnN        -        -       -        -        -        -        -        -        -        1.300\n";
        card += "QCDscale_ggWW             lnN        -        -       -        -        1.30     -        -        -        -        -\n";
        card += "QCDscale_ggH_ACCEPT       lnN        1.02     -       -        -        -        -        -        -        -        -\n";
        card += "QCDscale_qqH_ACCEPT       lnN        -        1.02    -        -        -        -        -        -        -        -\n";
        card += "BRhiggs_WW2l              lnN        1.02     1.02    1.02     -        -        -        -        -        -        -\n";
        card += "FakeRate                  lnN        -        -       -        -        -        -        1.360    -        -        -\n";
        card += "CMS_eff_m                 lnN        1.030    1.030   1.030    -        -        -        -        1.030    1.030    1.030\n";
        card += "CMS_eff_e                 lnN        1.040    1.040   1.040    -        -        -        -        1.040    1.040    1.040\n";
        //                                            ggH      qqH     wzttH    qqWW     ggWW     top      wj       ot       dy       wgstar
    if(channel == of0j || channel == sf0j) {
        card += "CMS_hww_0j_ttbar          lnN        -         -       -       -        -        1.210    -        -        -        -\n";
        card += "CMS_hww_0j_WW             lnN        -         -       -       1.076    1.076    -        -        -        -        -\n";
    } 
    else if(channel == of1j || channel == sf1j) {
        card += "CMS_hww_1j_ttbar          lnN        -         -       -       -       -        1.032    -        -        -        -\n";
        card += "CMS_hww_1j_WW             lnN        -         -       -       1.150   1.150    -        -        -        -        -\n";
    } 
    if (channel == of0j) {
        card += "CMS_hwwof_0j_Z            lnN        -         -       -       -        -       -         -       -        1.151    -\n";
    }
    else if (channel == of1j) {
        card += "CMS_hwwof_1j_Z            lnN        -         -       -        -       -         -       -       -        1.151    -\n";
    }
    else if (channel == sf0j) {
        card += "CMS_hwwsf_0j_Z            lnN        -         -       -        -       -         -       -       -        1.45     -\n";
    }
    else if (channel == sf1j) {
        card += "CMS_hwwsf_1j_Z            lnN        -         -       -        -       -         -       -       -        1.185    -\n";
    }
    if (!do1D) {
      // TODO: shape syst on deltaphi
      //        card += "CMS_ww2l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    // shape systematics here. Implemented only for OF 0/1 jet
    if (channel == of0j || channel == of1j) {
      card += "sig_"+chstr+"_8TeV_mean_err_res_met        param    0    SIG_"+chstr+"_8TeV_mean_err_res_met      [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_res_met       param    0    SIG_"+chstr+"_8TeV_sigma_err_res_met     [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_res_met       param    0    SIG_"+chstr+"_8TeV_alpha_err_res_met     [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_res_met           param    0    SIG_"+chstr+"_8TeV_n_err_res_met         [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_res_e          param    0    SIG_"+chstr+"_8TeV_mean_err_res_e        [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_res_e         param    0    SIG_"+chstr+"_8TeV_sigma_err_res_e       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_res_e         param    0    SIG_"+chstr+"_8TeV_alpha_err_res_e       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_res_e             param    0    SIG_"+chstr+"_8TeV_n_err_res_e           [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaleup_e      param    0    SIG_"+chstr+"_8TeV_mean_err_scaleup_e    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaleup_e     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaleup_e   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaleup_e     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaleup_e   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaleup_e         param    0    SIG_"+chstr+"_8TeV_n_err_scaleup_e       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaledn_e      param    0    SIG_"+chstr+"_8TeV_mean_err_scaledn_e    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaledn_e     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaledn_e   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaledn_e     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaledn_e   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaledn_e         param    0    SIG_"+chstr+"_8TeV_n_err_scaledn_e       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaleup_m      param    0    SIG_"+chstr+"_8TeV_mean_err_scaleup_m    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaleup_m     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaleup_m   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaleup_m     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaleup_m   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaleup_m         param    0    SIG_"+chstr+"_8TeV_n_err_scaleup_m       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaledn_m      param    0    SIG_"+chstr+"_8TeV_mean_err_scaledn_m    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaledn_m     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaledn_m   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaledn_m     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaledn_m   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaledn_m         param    0    SIG_"+chstr+"_8TeV_n_err_scaledn_m       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaleup_j      param    0    SIG_"+chstr+"_8TeV_mean_err_scaleup_j    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaleup_j     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaleup_j   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaleup_j     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaleup_j   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaleup_j         param    0    SIG_"+chstr+"_8TeV_n_err_scaleup_j       [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_mean_err_scaledn_j      param    0    SIG_"+chstr+"_8TeV_mean_err_scaledn_j    [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_sigma_err_scaledn_j     param    0    SIG_"+chstr+"_8TeV_sigma_err_scaledn_j   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_alpha_err_scaledn_j     param    0    SIG_"+chstr+"_8TeV_alpha_err_scaledn_j   [-3,3]    \n";
      card += "sig_"+chstr+"_8TeV_n_err_scaledn_j         param    0    SIG_"+chstr+"_8TeV_n_err_scaledn_j       [-3,3]    \n";

      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_MC           param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_MC           [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_MC          param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_MC          [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaleup_qcd  param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaleup_qcd  [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaleup_qcd param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaleup_qcd [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaledn_qcd  param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaledn_qcd  [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaledn_qcd param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaledn_qcd [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_res_met      param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_res_met      [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_res_met     param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_res_met     [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_res_e        param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_res_e        [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_res_e       param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_res_e       [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaleup_e    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaleup_e    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaleup_e   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaleup_e   [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaledn_e    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaledn_e    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaledn_e   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaledn_e   [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaleup_m    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaleup_m    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaleup_m   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaleup_m   [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaledn_m    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaledn_m    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaledn_m   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaledn_m   [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaleup_j    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaleup_j    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaleup_j   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaleup_j   [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_mean_err_scaledn_j    param    0    BKG_QQWW_"+chstr+"_8TeV_mean_err_scaledn_j    [-3,3]    \n";
      card += "bkg_qqww_"+chstr+"_8TeV_sigma_err_scaledn_j   param    0    BKG_QQWW_"+chstr+"_8TeV_sigma_err_scaledn_j   [-3,3]    \n";

      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_res_met      param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_res_met     [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_res_met     param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_res_met    [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_res_e        param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_res_e       [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_res_e       param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_res_e      [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaleup_e    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaleup_e   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaleup_e   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaleup_e  [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaledn_e    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaledn_e   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaledn_e   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaledn_e  [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaleup_m    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaleup_m   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaleup_m   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaleup_m  [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaledn_m    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaledn_m   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaledn_m   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaledn_m  [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaleup_j    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaleup_j   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaleup_j   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaleup_j  [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_mean_err_scaledn_j    param    0    BKG_GGWW_"+chstr+"_8TeV_mean_err_scaledn_j   [-3,3]    \n";
      card += "bkg_ggww_"+chstr+"_8TeV_sigma_err_scaledn_j   param    0    BKG_GGWW_"+chstr+"_8TeV_sigma_err_scaledn_j  [-3,3]    \n";

      card += "bkg_top_"+chstr+"_8TeV_mean_err_res_met      param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_res_met     [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_res_met     param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_res_met    [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_res_e        param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_res_e       [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_res_e       param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_res_e      [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaleup_e    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaleup_e   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaleup_e   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaleup_e  [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaledn_e    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaledn_e   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaledn_e   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaledn_e  [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaleup_m    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaleup_m   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaleup_m   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaleup_m  [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaledn_m    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaledn_m   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaledn_m   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaledn_m  [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaleup_j    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaleup_j   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaleup_j   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaleup_j  [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_mean_err_scaledn_j    param    0    BKG_TOP_"+chstr+"_8TeV_mean_err_scaledn_j   [-3,3]    \n";
      card += "bkg_top_"+chstr+"_8TeV_sigma_err_scaledn_j   param    0    BKG_TOP_"+chstr+"_8TeV_sigma_err_scaledn_j  [-3,3]    \n";

      card += "bkg_wj_"+chstr+"_8TeV_mean_err_fakerateup     param    0    BKG_WJ_"+chstr+"_8TeV_mean_err_fakerateup    [-3,3]    \n";
      card += "bkg_wj_"+chstr+"_8TeV_sigma_err_fakerateup    param    0    BKG_WJ_"+chstr+"_8TeV_sigma_err_fakerateup   [-3,3]    \n";
      card += "bkg_wj_"+chstr+"_8TeV_mean_err_fakeratedn     param    0    BKG_WJ_"+chstr+"_8TeV_mean_err_fakeratedn    [-3,3]    \n";
      card += "bkg_wj_"+chstr+"_8TeV_sigma_err_fakeratedn    param    0    BKG_WJ_"+chstr+"_8TeV_sigma_err_fakeratedn   [-3,3]    \n";

      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_res_met      param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_res_met     [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_res_met     param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_res_met    [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_res_e        param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_res_e       [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_res_e       param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_res_e      [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaleup_e    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaleup_e   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaleup_e   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaleup_e  [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaledn_e    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaledn_e   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaledn_e   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaledn_e  [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaleup_m    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaleup_m   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaleup_m   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaleup_m  [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaledn_m    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaledn_m   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaledn_m   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaledn_m  [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaleup_j    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaleup_j   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaleup_j   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaleup_j  [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_mean_err_scaledn_j    param    0    BKG_WGSTAR_"+chstr+"_8TeV_mean_err_scaledn_j   [-3,3]    \n";
      card += "bkg_wgstar_"+chstr+"_8TeV_sigma_err_scaledn_j   param    0    BKG_WGSTAR_"+chstr+"_8TeV_sigma_err_scaledn_j  [-3,3]    \n";

      card += "bkg_others_"+chstr+"_8TeV_mean_err_res_met      param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_res_met     [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_res_met     param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_res_met    [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_res_e        param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_res_e       [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_res_e       param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_res_e      [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaleup_e    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaleup_e   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaleup_e   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaleup_e  [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaledn_e    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaledn_e   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaledn_e   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaledn_e  [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaleup_m    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaleup_m   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaleup_m   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaleup_m  [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaledn_m    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaledn_m   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaledn_m   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaledn_m  [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaleup_j    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaleup_j   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaleup_j   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaleup_j  [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_mean_err_scaledn_j    param    0    BKG_OTHERS_"+chstr+"_8TeV_mean_err_scaledn_j   [-3,3]    \n";
      card += "bkg_others_"+chstr+"_8TeV_sigma_err_scaledn_j   param    0    BKG_OTHERS_"+chstr+"_8TeV_sigma_err_scaledn_j  [-3,3]    \n";

    }
    return card;
}

#endif
