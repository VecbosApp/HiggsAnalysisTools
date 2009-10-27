#include <iostream>
#include <string>

using namespace std;

int weight(string sample, double ngen, double xsec, double filtereff, double lumi = 10, double prescale = 1);

int main(int argc, char* argv[]) {
  
  double wantLumi = 100.0; //pb-1  
  int w;

  std::cout << "sample" << "\t" << "xsec" << "\t" << "prescale" << "\t" << "weight" << std:: endl;
  
  w = weight("H120", 109550, 0.429291, 1.0, wantLumi);
  w = weight("H130", 108500, 0.798372, 1.0, wantLumi);
  w = weight("H140", 109700, 1.160314, 1.0, wantLumi);
  w = weight("H150", 109550, 1.427122, 1.0, wantLumi);
  w = weight("H155", 106550, 1.535030, 1.0, wantLumi);
  w = weight("H160", 108450, 1.645316, 1.0, wantLumi);
  w = weight("H165",  98550, 1.601116, 1.0, wantLumi);
  w = weight("H170", 108000, 1.501921, 1.0, wantLumi);
  w = weight("H175",  98400, 1.403838, 1.0, wantLumi);
  w = weight("H180", 108000, 1.288077, 1.0, wantLumi);
  w = weight("H190", 108000, 0.956866, 1.0, wantLumi);
  w = weight("H200", 105650, 0.811859, 1.0, wantLumi);
  w = weight("H210", 108000, 0.720998, 1.0, wantLumi);
  w = weight("H220", 107900, 0.652178, 1.0, wantLumi);
  w = weight("H230", 105350, 0.596385, 1.0, wantLumi);
  w = weight("H240", 105500, 0.549726, 1.0, wantLumi);
  w = weight("H250", 107300, 0.510142, 1.0, wantLumi);
  w = weight("H275", 113150, 0.434763, 1.0, wantLumi);
  w = weight("H300", 100500, 0.385833, 1.0, wantLumi);
  w = weight("H350",  94500, 0.382384, 1.0, wantLumi);
  w = weight("H400", 101250, 0.294458, 1.0, wantLumi);
  w = weight("H450", 103400, 0.194491, 1.0, wantLumi);
  w = weight("H500", 105650, 0.126372, 1.0, wantLumi);
  w = weight("H550",  98750, 0.082537, 1.0, wantLumi);
  w = weight("H600", 100400, 0.054517, 1.0, wantLumi);

  w = weight("Wgamma", 103720, 11960.0, 1.0, wantLumi);
  w = weight("WW", 8600326, 44.8, 1.0, wantLumi);
  w = weight("WZ", 5115685, 17.4, 1.0, wantLumi);
  w = weight("ZZ", 3056640, 7.1, 1.0, wantLumi);

  // crosse sections LO at 10 TeV (even if the sample is at 7 TeV...)
  w = weight("Wjets", 1204434, 46050, 1.0, wantLumi);
  w = weight("Zjets",  181376,  7164, 1.0, wantLumi);
  w = weight("TTjets", 221331,   415, 1.0, wantLumi);

  w = weight("TTbar", 528940 , 375.0, 1.0, wantLumi);
  //  w = weight("Wenu", 2157227, 11840, 0.738, wantLumi);
  //  w = weight("Zee", 2675110, 1944, 1.0, wantLumi);
  w = weight("QCD_EMenriched_Pt20to30", 33.853E+06, 400E+06, 0.008, wantLumi);
  w = weight("QCD_EMenriched_Pt30to80", 42.284E+06, 100E+06, 0.047, wantLumi);
  w = weight("QCD_EMenriched_Pt80to170", 5.375E+06, 1.9E+06, 0.15, wantLumi);
  w = weight("QCD_BCtoE_Pt20to30", 2.468E+06, 400E+06, 0.00048, wantLumi);
  w = weight("QCD_BCtoE_Pt30to80", 2.041E+06, 100E+06, 0.0024, wantLumi);
  w = weight("QCD_BCtoE_Pt80to170", 1.042E+06, 1.9E+06, 0.012, wantLumi);
  w = weight("PhotonJet_Pt0to15", 109595, 100.5E+06, 1.0, wantLumi);
  w = weight("PhotonJet_Pt15to20", 107480, 168000, 1.0, wantLumi);
  w = weight("PhotonJet_Pt20to30", 108084, 84870, 1.0, wantLumi);
  w = weight("PhotonJet_Pt30to50", 124040, 26320, 1.0, wantLumi);
  w = weight("PhotonJet_Pt50to80", 112160, 4589, 1.0, wantLumi);
  w = weight("PhotonJet_Pt80to120", 110000, 786.4, 1.0, wantLumi);
  w = weight("PhotonJet_Pt120to170", 110000, 164.8, 1.0, wantLumi);
  w = weight("PhotonJet_Pt170to300", 106143, 45.96, 1.0, wantLumi);
  w = weight("PhotonJet_Pt300to500", 108290, 3.708, 1.0, wantLumi);
  w = weight("PhotonJet_Pt500toInf", 104735, 0.3285, 1.0, wantLumi);

  w = weight("InclusiveMu", 6570971, 0.5091E+09, 0.0002881, wantLumi, 4); // run on SD prescaled

  //  w = weight("Wmunu", 75065, 11840, 0.691, wantLumi);
  //  w = weight("Zmumu", 11690, 11840, 0.691, wantLumi);
  

  return 0;

}

int weight(string sample, double ngen, double xsec, double filtereff, 
	   double lumi, double prescale) {

  double W = xsec * filtereff * lumi * prescale / ngen;

  std::cout << sample << "\t" << xsec << "\t" << prescale << "\t" << W << std:: endl;

  return W;

}
