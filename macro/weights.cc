#include <iostream>
#include <string>

using namespace std;

double weight(string sample, double ngen, double xsec, double filtereff, double lumi = 10, double prescale = 1);

int main(int argc, char* argv[]) {
  
  double wantLumi = 1000.0; //pb-1  
  double w;

  std::cout << "sample" << "\t" << "xsec" << "\t" << "weight" << std:: endl;
  
  w = weight("H120", 103700, 0.247143, 1.0, wantLumi);
  w = weight("H130", 108200, 0.452859, 1.0, wantLumi);
  w = weight("H140", 109550, 0.649260, 1.0, wantLumi);
  w = weight("H150", 109550, 0.787871, 1.0, wantLumi);
  w = weight("H155", 97850, 0.842093, 1.0, wantLumi);
  w = weight("H160", 109550, 0.897043, 1.0, wantLumi);
  w = weight("H165", 101600, 0.867591, 1.0, wantLumi);
  w = weight("H170", 108150, 0.808914, 1.0, wantLumi);
  w = weight("H175", 117450, 0.751628, 1.0, wantLumi);
  w = weight("H180", 105450, 0.685617, 1.0, wantLumi);
  w = weight("H190", 103550, 0.503611, 1.0, wantLumi);

  // MCFM NLO
  w = weight("Wgamma", 107050, 23.2*1.8, 1.0, wantLumi);
  w = weight("WW_2l", 221100, 42.9*0.324*0.324, 1.0, wantLumi);
  w = weight("WZ_3l", 109000, 18.3*(3*0.108)*(3*0.0337), 1.0, wantLumi);
  w = weight("ZZ_2l2nu", 109074, 5.9*(2*0.107*0.20), 1.0, wantLumi);

  // MCFM NLO
  w = weight("Wjets", 9.7689e+06, 3*9679.9, 1.0, wantLumi);
  w = weight("Zjets",  1.03492e+06,  3*1606.6, 1.0, wantLumi);
  w = weight("TTjets", 1.2834e+06,   165, 1.0, wantLumi);

  w = weight("SingleTop_sChannel", 312055, 4.6*(3*0.1080), 1.0, wantLumi);
  w = weight("SingleTop_tChannel", 478593, 63*(3*0.1080), 1.0, wantLumi);
  w = weight("SingleTop_tWChannel", 16437, 10.6, 1.0, wantLumi);

  w = weight("QCD_EMenriched_Pt20to30", 2.27098e+07, 235.5E+06, 0.0073, wantLumi);
  w = weight("QCD_EMenriched_Pt30to80", 3.11427e+07, 59.3E+06, 0.059, wantLumi);
  w = weight("QCD_EMenriched_Pt80to170", 2.17402e+06, 0.906E+06, 0.148, wantLumi);
  w = weight("QCD_BCtoE_Pt20to30", 1.94146e+06, 235.5E+06, 0.00046, wantLumi);
  w = weight("QCD_BCtoE_Pt30to80", 0, 59.3E+06, 0.00234, wantLumi);
  w = weight("QCD_BCtoE_Pt80to170", 8674, 0.906E+06, 0.0104, wantLumi);
  w = weight("PhotonJet_Pt0to15", 115265, 84.46E+06, 1.0, wantLumi);
  w = weight("PhotonJet_Pt15to20", 108560, 114700, 1.0, wantLumi);
  w = weight("PhotonJet_Pt20to30", 0, 57180, 1.0, wantLumi);
  w = weight("PhotonJet_Pt30to50", 110000, 16520, 1.0, wantLumi);
  w = weight("PhotonJet_Pt50to80", 109730, 2723, 1.0, wantLumi);
  w = weight("PhotonJet_Pt80to120", 10827, 446.2, 1.0, wantLumi);
  w = weight("PhotonJet_Pt120to170", 122281, 84.43, 1.0, wantLumi);
  w = weight("PhotonJet_Pt170to300", 125128, 22.55, 1.0, wantLumi);
  w = weight("PhotonJet_Pt300to500", 7606, 1.545, 1.0, wantLumi);
  w = weight("PhotonJet_Pt500toInf", 6895, 0.0923, 1.0, wantLumi);

  w = weight("InclusiveMu", 6570971, 0.5091E+09, 0.0002881, wantLumi, 4); // run on SD prescaled

  //  w = weight("Wmunu", 75065, 11840, 0.691, wantLumi);
  //  w = weight("Zmumu", 11690, 11840, 0.691, wantLumi);
  

  return 0;

}

double weight(string sample, double ngen, double xsec, double filtereff, 
	   double lumi, double prescale) {

  double W = xsec * filtereff * lumi * prescale / ngen;

  std::cout << sample << "\t" << xsec << "\t" << W << std:: endl;

  return W;

}
