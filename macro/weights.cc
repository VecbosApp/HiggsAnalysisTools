#include <iostream>
#include <string>

using namespace std;

double weight(string sample, double ngen, double xsec, double filtereff, double lumi = 10, double prescale = 1);

int main(int argc, char* argv[]) {
  
  double wantLumi = 1.0; //pb-1  
  double w;

  std::cout << "sample" << "\t" << "xsec" << "\t" << "weight" << std:: endl;
  
  w = weight("H120", 9269, 0.247143, 1.0, wantLumi);
  w = weight("H130", 109995, 0.452859, 1.0, wantLumi);
  w = weight("H140", 109550, 0.649260, 1.0, wantLumi);
  w = weight("H150", 109550, 0.787871, 1.0, wantLumi);
  w = weight("H155", 97850, 0.842093, 1.0, wantLumi);
  w = weight("H160", 147903, 0.897043 * 4./9., 1.0, wantLumi);
  w = weight("H165", 101600, 0.867591, 1.0, wantLumi);
  w = weight("H170", 108150, 0.808914, 1.0, wantLumi);
  w = weight("H175", 117450, 0.751628, 1.0, wantLumi);
  w = weight("H180", 105450, 0.685617, 1.0, wantLumi);
  w = weight("H190", 103550, 0.503611, 1.0, wantLumi);

  // MCFM NLO
  //  w = weight("Wgamma", 107050, 23.2*1.8, 1.0, wantLumi);
  w = weight("WW", 110000, 4.50347, 1.0, wantLumi);
  w = weight("WZ", 110000, 0.599442, 1.0, wantLumi);
  w = weight("ZZ", 3.11337e+06, 5.9, 1.0, wantLumi);

  // MCFM NLO
  w = weight("Wenu", 2.79853e+06, 31314./3. * 0.742, 1.0, wantLumi);
  w = weight("Wmunu", 5.33094e+06, 31314./3. * 0.742, 1.0, wantLumi);
  w = weight("Wtaunu", 3.91925e+06, 31314./3. * 0.742, 1.0, wantLumi);

  w = weight("DY(ee) > 20 GeV",  2.04899e+06,  4998/3., 1.0, wantLumi);
  w = weight("DY(mm) > 20 GeV",  1.89893e+06,  4998/3., 1.0, wantLumi);
  w = weight("DY(tautau) > 20 GeV", 1.84822e+06, 4998./3., 1.0, wantLumi);

  w = weight("TTjets", 2.01075e+06,   165, 1.0, wantLumi);

  w = weight("SingleTop_sChannel", 535570, 4.21 * (0.1080*3), 1.0, wantLumi);
  w = weight("SingleTop_tChannel", 926480, 64.6*(3*0.1080), 1.0, wantLumi);
  w = weight("SingleTop_tWChannel", 494961, 10.6 * (0.1080*3), 1.0, wantLumi);

  return 0;

}

double weight(string sample, double ngen, double xsec, double filtereff, 
	   double lumi, double prescale) {

  double W = xsec * filtereff * lumi * prescale / ngen;

  std::cout << sample << "\t" << xsec << "\t" << W << std:: endl;

  return W;

}
