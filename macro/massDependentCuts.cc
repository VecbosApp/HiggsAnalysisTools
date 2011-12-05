#include <TString.h>
#include <map>

TString higgsCuts(int mH, bool out) {
  
  std::map<int,TString> cuts;
  cuts.insert(std::make_pair(110,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115 &&  mth>80  &&  mth<110")));
  cuts.insert(std::make_pair(115,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115 &&  mth>80  &&  mth<110")));
  cuts.insert(std::make_pair(120,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115 &&  mth>80  &&  mth<120")));
  cuts.insert(std::make_pair(130,TString("ptll>45 && pt1>25  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<90  &&  mth>80  &&  mth<125")));
  cuts.insert(std::make_pair(130,TString("ptll>45 && pt1>25  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80  &&  mth<125")));
  cuts.insert(std::make_pair(140,TString("ptll>45 && pt1>25  &&  pt2>15  &&  dphill*180./TMath::Pi()<90  &&  mth>80  &&  mth<130")));
  cuts.insert(std::make_pair(150,TString("ptll>45 && pt1>27  &&  pt2>25  &&  dphill*180./TMath::Pi()<90  &&  mth>80  &&  mth<150")));
  cuts.insert(std::make_pair(160,TString("ptll>45 && pt1>30  &&  pt2>25  &&  dphill*180./TMath::Pi()<60  &&  mth>90  &&  mth<160")));
  cuts.insert(std::make_pair(170,TString("ptll>45 && pt1>34  &&  pt2>25  &&  dphill*180./TMath::Pi()<60  &&  mth>110 &&  mth<170")));
  cuts.insert(std::make_pair(180,TString("ptll>45 && pt1>36  &&  pt2>25  &&  dphill*180./TMath::Pi()<70  &&  mth>120 &&  mth<180")));
  cuts.insert(std::make_pair(190,TString("ptll>45 && pt1>38  &&  pt2>25  &&  dphill*180./TMath::Pi()<90  &&  mth>120 &&  mth<190")));
  cuts.insert(std::make_pair(200,TString("ptll>45 && pt1>40  &&  pt2>25  &&  dphill*180./TMath::Pi()<100 &&  mth>120 &&  mth<200")));
  cuts.insert(std::make_pair(250,TString("ptll>45 && pt1>55  &&  pt2>25  &&  dphill*180./TMath::Pi()<140 &&  mth>120 &&  mth<250")));
  cuts.insert(std::make_pair(300,TString("ptll>45 && pt1>70  &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<300")));
  cuts.insert(std::make_pair(350,TString("ptll>45 && pt1>80  &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<350")));
  cuts.insert(std::make_pair(400,TString("ptll>45 && pt1>90  &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<400")));
  cuts.insert(std::make_pair(450,TString("ptll>45 && pt1>110 &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<450")));
  cuts.insert(std::make_pair(500,TString("ptll>45 && pt1>120 &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<500")));
  cuts.insert(std::make_pair(550,TString("ptll>45 && pt1>130 &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<550")));
  cuts.insert(std::make_pair(600,TString("ptll>45 && pt1>140 &&  pt2>25  &&  dphill*180./TMath::Pi()<175 &&  mth>120 &&  mth<600")));

  std::map<int,TString> cutsOut;
  cutsOut.insert(std::make_pair(110,cuts[110]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(115,cuts[115]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(120,cuts[120]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(130,cuts[130]+TString(" && mll<45")));
  cutsOut.insert(std::make_pair(140,cuts[140]+TString(" && mll<45")));
  cutsOut.insert(std::make_pair(150,cuts[150]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(160,cuts[160]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(170,cuts[170]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(180,cuts[180]+TString(" && mll<60")));
  cutsOut.insert(std::make_pair(190,cuts[190]+TString(" && mll<80")));
  cutsOut.insert(std::make_pair(200,cuts[200]+TString(" && mll<90")));
  cutsOut.insert(std::make_pair(250,cuts[250]+TString(" && mll<150")));
  cutsOut.insert(std::make_pair(300,cuts[300]+TString(" && mll<200")));
  cutsOut.insert(std::make_pair(350,cuts[350]+TString(" && mll<250")));
  cutsOut.insert(std::make_pair(400,cuts[400]+TString(" && mll<300")));
  cutsOut.insert(std::make_pair(450,cuts[450]+TString(" && mll<350")));
  cutsOut.insert(std::make_pair(500,cuts[500]+TString(" && mll<400")));
  cutsOut.insert(std::make_pair(550,cuts[550]+TString(" && mll<450")));
  cutsOut.insert(std::make_pair(600,cuts[600]+TString(" && mll<500")));

  if(out) return cutsOut[mH];
  else return cuts[mH];

}

TString higgsCutsNoMT(int mH, bool out) {
  
  std::map<int,TString> cuts;
  cuts.insert(std::make_pair(110,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115")));
  cuts.insert(std::make_pair(115,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115")));
  cuts.insert(std::make_pair(120,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<115")));
  cuts.insert(std::make_pair(130,TString("ptll>45 && pt1>25  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  dphill*180./TMath::Pi()<90")));
  cuts.insert(std::make_pair(140,TString("ptll>45 && pt1>25  &&  pt2>15  &&  dphill*180./TMath::Pi()<90 ")));
  cuts.insert(std::make_pair(150,TString("ptll>45 && pt1>27  &&  pt2>25  &&  dphill*180./TMath::Pi()<90 ")));
  cuts.insert(std::make_pair(160,TString("ptll>45 && pt1>30  &&  pt2>25  &&  dphill*180./TMath::Pi()<60 ")));
  cuts.insert(std::make_pair(170,TString("ptll>45 && pt1>34  &&  pt2>25  &&  dphill*180./TMath::Pi()<60 ")));
  cuts.insert(std::make_pair(180,TString("ptll>45 && pt1>36  &&  pt2>25  &&  dphill*180./TMath::Pi()<70 ")));
  cuts.insert(std::make_pair(190,TString("ptll>45 && pt1>38  &&  pt2>25  &&  dphill*180./TMath::Pi()<90 ")));
  cuts.insert(std::make_pair(200,TString("ptll>45 && pt1>40  &&  pt2>25  &&  dphill*180./TMath::Pi()<100")));
  cuts.insert(std::make_pair(250,TString("ptll>45 && pt1>55  &&  pt2>25  &&  dphill*180./TMath::Pi()<140")));
  cuts.insert(std::make_pair(300,TString("ptll>45 && pt1>70  &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(350,TString("ptll>45 && pt1>80  &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(400,TString("ptll>45 && pt1>90  &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(450,TString("ptll>45 && pt1>110 &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(500,TString("ptll>45 && pt1>120 &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(550,TString("ptll>45 && pt1>130 &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));
  cuts.insert(std::make_pair(600,TString("ptll>45 && pt1>140 &&  pt2>25  &&  dphill*180./TMath::Pi()<175")));

  std::map<int,TString> cutsOut;
  cutsOut.insert(std::make_pair(110,cuts[110]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(115,cuts[115]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(120,cuts[120]+TString(" && mll<40")));
  cutsOut.insert(std::make_pair(130,cuts[130]+TString(" && mll<45")));
  cutsOut.insert(std::make_pair(140,cuts[140]+TString(" && mll<45")));
  cutsOut.insert(std::make_pair(150,cuts[150]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(160,cuts[160]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(170,cuts[170]+TString(" && mll<50")));
  cutsOut.insert(std::make_pair(180,cuts[180]+TString(" && mll<60")));
  cutsOut.insert(std::make_pair(190,cuts[190]+TString(" && mll<80")));
  cutsOut.insert(std::make_pair(200,cuts[200]+TString(" && mll<90")));
  cutsOut.insert(std::make_pair(250,cuts[250]+TString(" && mll<150")));
  cutsOut.insert(std::make_pair(300,cuts[300]+TString(" && mll<200")));
  cutsOut.insert(std::make_pair(350,cuts[350]+TString(" && mll<250")));
  cutsOut.insert(std::make_pair(400,cuts[400]+TString(" && mll<300")));
  cutsOut.insert(std::make_pair(450,cuts[450]+TString(" && mll<350")));
  cutsOut.insert(std::make_pair(500,cuts[500]+TString(" && mll<400")));
  cutsOut.insert(std::make_pair(550,cuts[550]+TString(" && mll<450")));
  cutsOut.insert(std::make_pair(600,cuts[600]+TString(" && mll<500")));

  if(out) return cutsOut[mH];
  else return cuts[mH];

}

TString higgsCutsBDT(int mH, bool out) {
  
  std::map<int,TString> cuts;
  cuts.insert(std::make_pair(110,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<110")));
  cuts.insert(std::make_pair(115,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<115")));
  cuts.insert(std::make_pair(120,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<120")));
  cuts.insert(std::make_pair(130,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<130")));
  cuts.insert(std::make_pair(140,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<140")));
  cuts.insert(std::make_pair(150,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<150")));
  cuts.insert(std::make_pair(160,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<160")));
  cuts.insert(std::make_pair(170,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<170")));
  cuts.insert(std::make_pair(180,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<180")));
  cuts.insert(std::make_pair(190,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<190")));
  cuts.insert(std::make_pair(200,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<200")));
  cuts.insert(std::make_pair(250,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<250")));
  cuts.insert(std::make_pair(300,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<300")));
  cuts.insert(std::make_pair(350,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<350")));
  cuts.insert(std::make_pair(400,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<400")));
  cuts.insert(std::make_pair(450,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<450")));
  cuts.insert(std::make_pair(500,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<500")));
  cuts.insert(std::make_pair(550,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<550")));
  cuts.insert(std::make_pair(600,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav)) &&  mth>80 && mth<600")));

  std::map<int,TString> cutsOut;
  cutsOut.insert(std::make_pair(110,cuts[110]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(115,cuts[115]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(120,cuts[120]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(130,cuts[130]+TString(" && mll<80")));
  cutsOut.insert(std::make_pair(140,cuts[140]+TString(" && mll<90")));
  cutsOut.insert(std::make_pair(150,cuts[150]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(160,cuts[160]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(170,cuts[170]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(180,cuts[180]+TString(" && mll<110")));
  cutsOut.insert(std::make_pair(190,cuts[190]+TString(" && mll<120")));
  cutsOut.insert(std::make_pair(200,cuts[200]+TString(" && mll<130")));
  cutsOut.insert(std::make_pair(250,cuts[250]+TString(" && mll<250")));
  cutsOut.insert(std::make_pair(300,cuts[300]+TString(" && mll<300")));
  cutsOut.insert(std::make_pair(350,cuts[350]+TString(" && mll<350")));
  cutsOut.insert(std::make_pair(400,cuts[400]+TString(" && mll<400")));
  cutsOut.insert(std::make_pair(450,cuts[450]+TString(" && mll<450")));
  cutsOut.insert(std::make_pair(500,cuts[500]+TString(" && mll<500")));
  cutsOut.insert(std::make_pair(550,cuts[550]+TString(" && mll<550")));
  cutsOut.insert(std::make_pair(600,cuts[600]+TString(" && mll<600")));

  if(out) return cutsOut[mH];
  else return cuts[mH];

}

TString higgsCutsBDTNoMT(int mH, bool out) {
  
  std::map<int,TString> cuts;
  cuts.insert(std::make_pair(110,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(115,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(120,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(130,TString("ptll>45 && pt1>25  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(140,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(150,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(160,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(170,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(180,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(190,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(200,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(250,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(300,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(350,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(400,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(450,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(500,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(550,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));
  cuts.insert(std::make_pair(600,TString("ptll>45 && pt1>20  &&  ((pt2>10 && !sameflav) || (pt2>15 && sameflav))")));

  std::map<int,TString> cutsOut;
  cutsOut.insert(std::make_pair(110,cuts[110]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(115,cuts[115]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(120,cuts[120]+TString(" && mll<70")));
  cutsOut.insert(std::make_pair(130,cuts[130]+TString(" && mll<80")));
  cutsOut.insert(std::make_pair(140,cuts[140]+TString(" && mll<90")));
  cutsOut.insert(std::make_pair(150,cuts[150]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(160,cuts[160]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(170,cuts[170]+TString(" && mll<100")));
  cutsOut.insert(std::make_pair(180,cuts[180]+TString(" && mll<110")));
  cutsOut.insert(std::make_pair(190,cuts[190]+TString(" && mll<120")));
  cutsOut.insert(std::make_pair(200,cuts[200]+TString(" && mll<130")));
  cutsOut.insert(std::make_pair(250,cuts[250]+TString(" && mll<250")));
  cutsOut.insert(std::make_pair(300,cuts[300]+TString(" && mll<300")));
  cutsOut.insert(std::make_pair(350,cuts[350]+TString(" && mll<350")));
  cutsOut.insert(std::make_pair(400,cuts[400]+TString(" && mll<400")));
  cutsOut.insert(std::make_pair(450,cuts[450]+TString(" && mll<450")));
  cutsOut.insert(std::make_pair(500,cuts[500]+TString(" && mll<500")));
  cutsOut.insert(std::make_pair(550,cuts[550]+TString(" && mll<550")));
  cutsOut.insert(std::make_pair(600,cuts[600]+TString(" && mll<600")));

  if(out) return cutsOut[mH];
  else return cuts[mH];

}

