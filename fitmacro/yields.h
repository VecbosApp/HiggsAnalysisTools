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
          if(production_ == ggH) fss << " 1178.32   - 26.89600*@0 + 0.1997410*@0*@0 - 0.000475800*@0*@0*@0";
          if(production_ == vbfH) fss << "13.6107   - 0.304446*@0 + 0.0022179*@0*@0 - 5.19069e-06*@0*@0*@0";
          if(production_ == wzttH) fss << "0.147542 - 0.071180*@0 + 0.0010193*@0*@0 - 3.38809e-06*@0*@0*@0";
        }
        if (channel_ == of1j) {
          if(production_ == ggH) fss << "564.704  - 12.5945*@0  + 0.0915010*@0*@0 - 0.000213701*@0*@0*@0";
          if(production_ == vbfH) fss << "53.5115 - 1.19516*@0  + 0.0086761*@0*@0 - 2.01715e-05*@0*@0*@0";
          if(production_ == wzttH) fss << "23.549 - 0.587955*@0 + 0.0047250*@0*@0 - 1.19170e-05*@0*@0*@0";
        }
        if (channel_ == sf0j) {
          if(production_ == ggH) fss << "1251.46  - 28.29400*@0 + 0.208324000*@0*@0 - 0.000492681*@0*@0*@0";
          if(production_ == vbfH) fss << "12.737  - 0.282539*@0 + 0.002040230*@0*@0 - 4.72842e-06*@0*@0*@0";
          if(production_ == wzttH) fss << "-7.588 + 0.091194*@0 - 9.30171e-05*@0*@0 - 8.96112e-07*@0*@0*@0";
        }
        if (channel_ == sf1j) {
          if(production_ == ggH) fss << "608.34  - 13.4254*@0 + 0.096626*@0*@0 - 0.000224077*@0*@0*@0";
          if(production_ == vbfH) fss << "59.322 - 1.30037*@0 + 0.009286*@0*@0 - 2.13219e-05*@0*@0*@0";
          if(production_ == wzttH) fss << "21.06 - 0.55640*@0 + 0.004651*@0*@0 - 1.20667e-05*@0*@0*@0";
        }
        // they were calculated for 1 fb-1 @7TeV
        fss << " )*" << lumi;        
      } else {
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
        // they were calculated for 19.47 fb-1 @8TeV
        fss << " )/19.47*" << lumi;
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

