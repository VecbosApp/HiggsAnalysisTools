#include <TString.h>
#include <map>

TString higgsCuts(int mH, bool out) {
  
  std::map<int,TString> cuts;
  cuts.insert(std::make_pair(115,TString("maxPtEle>20 && minPtEle>10 && deltaPhi<115 && transvMass>70 && transvMass<110")));
  cuts.insert(std::make_pair(120,TString("maxPtEle>20 && minPtEle>10 && deltaPhi<115 && transvMass>70 && transvMass<120")));
  cuts.insert(std::make_pair(130,TString("maxPtEle>25 && minPtEle>10 && deltaPhi<90 && transvMass>75 && transvMass<125")));
  cuts.insert(std::make_pair(140,TString("maxPtEle>25 && minPtEle>15 && deltaPhi<90 && transvMass>80 && transvMass<130")));
  cuts.insert(std::make_pair(150,TString("maxPtEle>27 && minPtEle>25 && deltaPhi<90 && transvMass>80 && transvMass<150")));
  cuts.insert(std::make_pair(160,TString("maxPtEle>30 && minPtEle>25 && deltaPhi<60 && transvMass>90 && transvMass<160")));
  cuts.insert(std::make_pair(170,TString("maxPtEle>34 && minPtEle>25 && deltaPhi<60 && transvMass>110 && transvMass<170")));
  cuts.insert(std::make_pair(180,TString("maxPtEle>36 && minPtEle>25 && deltaPhi<70 && transvMass>120 && transvMass<180")));
  cuts.insert(std::make_pair(190,TString("maxPtEle>38 && minPtEle>25 && deltaPhi<90 && transvMass>120 && transvMass<190")));
  cuts.insert(std::make_pair(200,TString("maxPtEle>40 && minPtEle>25 && deltaPhi<100 && transvMass>120 && transvMass<200")));
  cuts.insert(std::make_pair(250,TString("maxPtEle>55 && minPtEle>25 && deltaPhi<140 && transvMass>120 && transvMass<250")));
  cuts.insert(std::make_pair(300,TString("maxPtEle>70 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<300")));
  cuts.insert(std::make_pair(350,TString("maxPtEle>80 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<350")));
  cuts.insert(std::make_pair(400,TString("maxPtEle>90 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<400")));
  cuts.insert(std::make_pair(450,TString("maxPtEle>110 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<450")));
  cuts.insert(std::make_pair(500,TString("maxPtEle>120 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<500")));
  cuts.insert(std::make_pair(550,TString("maxPtEle>130 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<550")));
  cuts.insert(std::make_pair(600,TString("maxPtEle>140 && minPtEle>25 && deltaPhi<175 && transvMass>120 && transvMass<600")));

  std::map<int,TString> cutsOut;
  cutsOut.insert(std::make_pair(115,cuts[115]+TString(" && eleInvMass<40")));
  cutsOut.insert(std::make_pair(120,cuts[120]+TString(" && eleInvMass<40")));
  cutsOut.insert(std::make_pair(130,cuts[130]+TString(" && eleInvMass<45")));
  cutsOut.insert(std::make_pair(140,cuts[140]+TString(" && eleInvMass<45")));
  cutsOut.insert(std::make_pair(150,cuts[150]+TString(" && eleInvMass<50")));
  cutsOut.insert(std::make_pair(160,cuts[160]+TString(" && eleInvMass<50")));
  cutsOut.insert(std::make_pair(170,cuts[170]+TString(" && eleInvMass<50")));
  cutsOut.insert(std::make_pair(180,cuts[180]+TString(" && eleInvMass<60")));
  cutsOut.insert(std::make_pair(190,cuts[190]+TString(" && eleInvMass<80")));
  cutsOut.insert(std::make_pair(200,cuts[200]+TString(" && eleInvMass<90")));
  cutsOut.insert(std::make_pair(250,cuts[250]+TString(" && eleInvMass<150")));
  cutsOut.insert(std::make_pair(300,cuts[300]+TString(" && eleInvMass<200")));
  cutsOut.insert(std::make_pair(350,cuts[350]+TString(" && eleInvMass<250")));
  cutsOut.insert(std::make_pair(400,cuts[400]+TString(" && eleInvMass<300")));
  cutsOut.insert(std::make_pair(450,cuts[450]+TString(" && eleInvMass<350")));
  cutsOut.insert(std::make_pair(500,cuts[500]+TString(" && eleInvMass<400")));
  cutsOut.insert(std::make_pair(550,cuts[550]+TString(" && eleInvMass<450")));
  cutsOut.insert(std::make_pair(600,cuts[600]+TString(" && eleInvMass<500")));

  if(out) return cutsOut[mH];
  else return cuts[mH];

}
