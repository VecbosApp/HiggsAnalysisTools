#ifndef CARDTEMPLATE_H
#define CARDTEMPLATE_H

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
        card += "bin     BIN           BIN            BIN           BIN            BIN            BIN             BIN\n";
        card += "process ggH           qqH            bkg_ww        bkg_top        bkg_wjets      bkg_others      bkg_dy\n";
        card += "process -5            -4             1             2              3              4               5\n";
        card += "rate    SIG_GGH_YIELD SIG_VBF_YIELD  BKG_WW_YIELD  BKG_TOP_YIELD  BKG_WJETS_YIELD  BKG_OTHERS_YIELD  BKG_DY_YIELD\n";
        card += "------------\n";
        card += "lumi_8TeV                 lnN        1.05    1.05    1.05     1.05     -            1.05          -\n";
        card += "pdf_gg                    lnN        GGH_PDF -       -        TOP_PDF  WJETS_PDF    OTHERS_PDF    DY_PDF\n";
        card += "pdf_qqbar                 lnN         -      VBF_PDF WW_PDF   -        -            -             -\n";
        card += "pdf_hww2l2nu_accept       lnN        1.02    1.02    -        -        -            -             -\n";
        card += "QCDscale_ggH              lnN        GGH_QCD -       -        -        -            -             -\n";
        card += "QCDscale_qqH              lnN        -       VBF_PDF -        -        -            -             -\n";
	//        card += "QCDscale_ggVV             lnN TODO
        card += "QCDscale_VV               lnN        -       -       WW_QCD   -        -            -             -\n";
        card += "BRhiggs_WW2l              lnN        1.02    1.02    -        -        -            -             -\n";
        card += "CMS_Wjets                 lnN        -       -       -        -        -            -             0.36\n";
        card += "CMS_eff_m                 lnN        1.015   1.015   1.015    1.015    -            1.015         -\n";
        card += "CMS_eff_e                 lnN        1.01    1.01    1.01     1.01     -            1.01          -\n";
    if (channel == 0) {
        card += "sig_of_0j_mean_err_7TeV   param      0        0.005                     \n";
        card += "sig_of_0j_sigma_err_7TeV  param      0        0.3                       \n";
    }
    else if (channel == 1) {
        card += "sig_of_1j_mean_err_7TeV      param      0        0.004                     \n";
        card += "sig_of_1j_sigma_err_7TeV     param      0        0.3                       \n";
    }
    else if (channel == 2) {
        card += "sig_sf_0j_mean_err_7TeV   param      0        0.005                      \n";
        card += "sig_sf_0j_sigma_err_7TeV  param      0        0.3                        \n";
    }
    else if (channel == 3) {
        card += "sig_sf_1j_mean_err_7TeV   param      0        0.005                      \n";
        card += "sig_sf_1j_sigma_err_7TeV  param      0        0.3                        \n";
    }
    if (!do1D) {
        card += "CMS_ww2l_bkgMELA          param      0       1       [-3,3]             \n";
    }
    return card;
}

#endif
