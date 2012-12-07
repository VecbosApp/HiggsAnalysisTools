#ifndef YIELDS_H
#define YIELDS_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std;
enum production { ggH, vbfH };

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
	if(production_ == ggH) fss << "1071.72 - 24.4336*@0 + 0.180539*@0*@0 - 0.000425689*@0*@0*@0";
	if(production_ == vbfH) fss << "14.81  - 0.33084*@0 + 0.002387*@0*@0 - 5.45931e-06*@0*@0*@0";
      }
      if (channel_ == of1j) {
	if(production_ == ggH) fss << "466.12 - 10.5387*@0 + 0.077163*@0*@0 - 0.000179912*@0*@0*@0";
	if(production_ == vbfH) fss << "52.18 - 1.17388*@0 + 0.008527*@0*@0 - 1.96117e-05*@0*@0*@0";
      }
      if (channel_ == sf0j) {
	if(production_ == ggH) fss << "1153.17 - 26.4545*@0 + 0.196776*@0*@0 - 0.000467208*@0*@0*@0";
	if(production_ == vbfH) fss << "32.95  - 0.7147*@0  + 0.005067*@0*@0 - 1.16101e-05*@0*@0*@0";
      }
      if (channel_ == sf1j) {
	if(production_ == ggH) fss << "517.843 - 11.7562*@0 + 0.0865359*@0*@0 - 0.000203278*@0*@0*@0";
	if(production_ == vbfH) fss << "103.60 - 2.25182*@0 + 0.0159676*@0*@0 - 3.64287e-05*@0*@0*@0";
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
      if(channel_==of0j || channel_==sf0j) return 1.14; 
      if(channel_==of1j || channel_==sf1j) return 0.90;
      return 1.0;
    }
    float getTop() { 
      // from AN-12-378 v5
      if(channel_==of0j || channel_==sf0j) return 0.95; 
      if(channel_==of1j || channel_==sf1j) return 1.09;
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

