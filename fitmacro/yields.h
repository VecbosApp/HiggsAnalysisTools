#ifndef YIELDS_H
#define YIELDS_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std;
enum production { ggH, vbfH, wzttH };

class InterpolatedYields {

 protected:
  int channel_;
  int production_;
  bool do7TeV_;

 public:
 
  InterpolatedYields(int channel, int production, bool do7TeV):
    channel_(channel), 
    production_(production),
    do7TeV_(do7TeV) {}
  
  virtual ~InterpolatedYields() {}

  std::string getInterpolatedYieldString(float lumi) {
      stringstream fss;
      fss << "( ";  
      if(do7TeV_) {
        if (channel_ == of0j) {
          if(production_ == ggH) fss << "(1143.04) + (-26.5719*@0) + (0.200632*@0*@0) + (-0.000484228*@0*@0*@0)";
          if(production_ == vbfH) fss << "(13.7678) + (-0.310861*@0) + (0.00228496*@0*@0) + (-5.38769e-06*@0*@0*@0)";
          if(production_ == wzttH) fss << "(-5.71273) + (0.0416904*@0) + (0.000319699*@0*@0) + (-1.97016e-06*@0*@0*@0)";
        }
        if (channel_ == of1j) {
          if(production_ == ggH) fss << "(567.858) + (-12.8277*@0) + (0.0943411*@0*@0) + (-0.000222659*@0*@0*@0)";
          if(production_ == vbfH) fss << "(54.4908) + (-1.22711*@0) + (0.00898224*@0*@0) + (-2.10377e-05*@0*@0*@0)";
          if(production_ == wzttH) fss << "(19.2559) + (-0.515524*@0) + (0.0043563*@0*@0) + (-1.13483e-05*@0*@0*@0)";
        }
        if (channel_ == sf0j) {
          if(production_ == ggH) fss << "(1245.84) + (-28.2856*@0) + (0.209094*@0*@0) + (-0.000496139*@0*@0*@0)";
          if(production_ == vbfH) fss << "(12.4096) + (-0.276765*@0) + (0.00200813*@0*@0) + (-4.67134e-06*@0*@0*@0)";
          if(production_ == wzttH) fss << "(-10.1661) + (0.143923*@0) + (-0.000445489*@0*@0) + (-1.20312e-07*@0*@0*@0)";
        }
        if (channel_ == sf1j) {
          if(production_ == ggH) fss << " (602.641) + (-13.3449*@0) + (0.0963668*@0*@0) + (-0.000224109*@0*@0*@0)";
          if(production_ == vbfH) fss << "(60.4911) + (-1.32972*@0) + (0.00952514*@0*@0) + (-2.19434e-05*@0*@0*@0)";
          if(production_ == wzttH) fss << "(17.6056) + (-0.486441*@0) + (0.00418762*@0*@0) + (-1.10529e-05*@0*@0*@0)";
        }
        fss << " )*" << lumi;        
      } else {
        if (channel_ == of0j) {
          if(production_ == ggH) fss << "(1217.28) + (-28.4603*@0) + (0.215522*@0*@0) + (-0.000519669*@0*@0*@0)";
          if(production_ == vbfH) fss << "(20.1416) + (-0.458206*@0) + (0.00337207*@0*@0) + (-7.88492e-06*@0*@0*@0)";
          if(production_ == wzttH) fss << "(-3.71819) + (0.00350664*@0) + (0.000561685*@0*@0) + (-2.47807e-06*@0*@0*@0)";
        }
        if (channel_ == of1j) {
          if(production_ == ggH) fss << "(571.628) + (-13.1265*@0) + (0.0977026*@0*@0) + (-0.000231797*@0*@0*@0)";
          if(production_ == vbfH) fss << "(69.228) + (-1.57485*@0) + (0.0115823*@0*@0) + (-2.70207e-05*@0*@0*@0)";
          if(production_ == wzttH) fss << "(16.5684) + (-0.462684*@0) + (0.00401544*@0*@0) + (-1.06248e-05*@0*@0*@0)";
        }
        if (channel_ == sf0j) {
          if(production_ == ggH) fss << "(1183.86) + (-27.3423*@0) + (0.204699*@0*@0) + (-0.000488713*@0*@0*@0)";
          if(production_ == vbfH) fss << "(17.1506) + (-0.386847*@0) + (0.00281906*@0*@0) + (-6.51418e-06*@0*@0*@0)";
          if(production_ == wzttH) fss << "(-4.74726) + (0.0341373*@0) + (0.000289565*@0*@0) + (-1.7456e-06*@0*@0*@0)";
        }
        if (channel_ == sf1j) {
          if(production_ == ggH) fss << "(461.347) + (-10.5413*@0) + (0.0779025*@0*@0) + (-0.000182962*@0*@0*@0)";
          if(production_ == vbfH) fss << "(53.9046) + (-1.21768*@0) + (0.00886474*@0*@0) + (-2.03634e-05*@0*@0*@0)";
          if(production_ == wzttH) fss << "(29.2128) + (-0.726066*@0) + (0.00581923*@0*@0) + (-1.47149e-05*@0*@0*@0)";
        }
        fss << " )*" << lumi;
      }
      return fss.str();
  }

};

class ScaleFactors {

 protected:
  int channel_;
  bool do7TeV_;

 public:
  ScaleFactors(int channel, bool do7TeV):
    channel_(channel),
    do7TeV_(do7TeV) {}

    virtual ~ScaleFactors() {}

    // from AN-13-022 v6
    float getWW() {
      if(do7TeV_) {
        if(channel_==of0j || channel_==sf0j) return 1.20; 
        if(channel_==of1j || channel_==sf1j) return 1.10;        
      } else {
        if(channel_==of0j || channel_==sf0j) return 1.05; 
        if(channel_==of1j || channel_==sf1j) return 0.90;
      }
      return 1.0;
    }
    float getTop() { 
      if(do7TeV_) {
        if(channel_==of0j || channel_==sf0j) return 1.31; 
        if(channel_==of1j || channel_==sf1j) return 1.08;
      } else {
        if(channel_==of0j || channel_==sf0j) return 0.98; 
        if(channel_==of1j || channel_==sf1j) return 1.08;
      }
      return 1.0;
    }
    float getDY() {
      if(channel_==of0j || channel_==of1j) return 1.0; // Dy->tt emb data sample used
      if(channel_==sf0j) return 2.2; // re-calculated
      if(channel_==sf1j) return 1.6; // re-calculated 
      return 1.0;
    }    

};

#endif

