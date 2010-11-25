#include <iostream>
#include <string>

using namespace std;

double weight(string sample, double ngen, double xsec, double filtereff, double lumi = 10, double prescale = 1);

int main(int argc, char* argv[]) {
  
  double wantLumi = 1.0; //pb-1  
  double w;

  std::cout << "sample" << "\t" << "xsec" << "\t" << "weight" << std:: endl;
  
  w = weight("H120", 9269, 0.247143, 1.0, wantLumi);
  w = weight("H130", 100000, 0.452859, 1.0, wantLumi);
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
  w = weight("ZZ", 2.11337e+06, 5.9, 1.0, wantLumi);

  // MCFM NLO
  w = weight("W(e+mu+tau)nu", 2.14853e+06 + 5.33094e+06 + 1.133e+06, 31314./3. * 0.742, 1.0, wantLumi);

  w = weight("DY 10-50 GeV",  188145,  1950, 1.0, wantLumi); // not good: is >50 xsec - 10-50 xsec
  w = weight("DY >50 GeV",  2027,  3048, 1.0, wantLumi); // mll > 50 GeV

  w = weight("TTjets", 2.01075e+06,   165, 1.0, wantLumi);

  w = weight("SingleTop_sChannel", 394967, 4.21 * (0.1080*3), 1.0, wantLumi);
  w = weight("SingleTop_tChannel", 484060, 64.6*(3*0.1080), 1.0, wantLumi);
  w = weight("SingleTop_tWChannel", 494961, 10.6 * (0.1080*3), 1.0, wantLumi);

  return 0;

}

double weight(string sample, double ngen, double xsec, double filtereff, 
	   double lumi, double prescale) {

  double W = xsec * filtereff * lumi * prescale / ngen;

  std::cout << sample << "\t" << xsec << "\t" << W << std:: endl;

  return W;

}
