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

 public:
 
 InterpolatedYields(int channel, int production):
  channel_(channel), 
    production_(production) {}
  
  virtual ~InterpolatedYields() {}

  std::string getInterpolatedYieldString(float lumi) {
      stringstream fss;
      fss << "( ";  
      if (channel_ == of0j) {
	if(production_ == ggH) fss << "21821.6  - 503.955*@0 + 3.771120*@0*@0 - 0.008995620*@0*@0*@0";
	if(production_ == vbfH) fss << "342.89  - 7.72813*@0 + 0.056292*@0*@0 - 0.000130115*@0*@0*@0";
        if(production_ == wzttH) fss << "-21.08 - 0.83254*@0 + 0.015693*@0*@0 - 5.59119e-05*@0*@0*@0";
      }
      if (channel_ == of1j) {
	if(production_ == ggH) fss << "9729.7  - 221.505*@0 + 1.633110*@0*@0 - 0.003834270*@0*@0*@0";
        if(production_ == vbfH) fss << "1191.9 - 26.8418*@0 + 0.195277*@0*@0 - 0.000450233*@0*@0*@0";
        if(production_ == wzttH) fss << "345.8 - 9.15234*@0 + 0.076578*@0*@0 - 0.000197964*@0*@0*@0";
      }
      if (channel_ == sf0j) {
	if(production_ == ggH) fss << "   22086.7 - 508.521*@0 + 3.79508*@0*@0 - 0.00903367*@0*@0*@0";
	if(production_ == vbfH) fss << "   315.24 - 7.09028*@0 + 0.05150*@0*@0 - 0.00011855*@0*@0*@0";
        if(production_ == wzttH) fss << "-51.4374 - 0.15662*@0 + 0.01095*@0*@0 - 4.5266e-05*@0*@0*@0";
      }
      if (channel_ == sf1j) {
	if(production_ == ggH) fss << "8444.36 - 192.373*@0 + 1.4169100*@0*@0 - 0.003315*@0*@0*@0";
	if(production_ == vbfH) fss << "968.52 - 21.7376*@0 + 0.1570710*@0*@0 - 0.000358*@0*@0*@0";
        if(production_ == wzttH) fss << "572.1 - 14.1389*@0 + 0.1127810*@0*@0 - 0.000284*@0*@0*@0";
      }
      fss << " )*" << lumi;
      return fss.str();
  }

};

class ScaleFactors {

 protected:
  int channel_;

 public:
  ScaleFactors(int channel):
    channel_(channel) {}

    virtual ~ScaleFactors() {}

    float getWW() {
      if(channel_==of0j || channel_==sf0j) return 1.05; 
      if(channel_==of1j || channel_==sf1j) return 0.90;
      return 1.0;
    }
    float getTop() { 
      // from AN-12-378 v5
      if(channel_==of0j || channel_==sf0j) return 0.98; 
      if(channel_==of1j || channel_==sf1j) return 1.08;
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

