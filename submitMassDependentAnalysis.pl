#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitMassDependentAnalysis.pl -p <prefix>";}

@masses = (120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600);

for($i=0; $i<($#masses+1); $i++) {
    $mass = $masses[$i];
    print "-------------------------->\n";
    print "SUBMITTING MASS SELECTION: mH = $mass ...\n\n";
    open(MASSFILE,">config/higgs/higgsMass.txt");
    print MASSFILE "HiggsMass\t$mass\n";

    $fullprefix=$prefix."/"."OptimMH$mass";

    print "submitting signals...\n";
    $higgsList = "GluGluToHToWWTo2L2Nu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 $higgsList 15 HiggsApp 8nh $fullprefix 1");
    
    $higgsList = "GluGluToHToWWToLNuTauNu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 $higgsList 15 HiggsApp 8nh $fullprefix 1");

    $higgsList = "VBF_HToWWTo2L2Nu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 $higgsList 15 HiggsApp 8nh $fullprefix 1");

    $higgsList = "VBF_HToWWToLNuTauNu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 $higgsList 15 HiggsApp 8nh $fullprefix 1");

    print "done with signals.\n";

    print  "submitting top...\n";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 TToBLNu_TuneZ2_s-channel_7TeV-madgraph 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 TToBLNu_TuneZ2_t-channel_7TeV-madgraph 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 TToBLNu_TuneZ2_tW-channel_7TeV-madgraph 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 TTJets_TuneZ2_7TeV-madgraph-tauola 15 HiggsApp 8nh $fullprefix 1");
    print  "done with top.\n";

    sleep 500;

    print  "submitting V+jets...\n";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 WJetsToLNu_TuneZ2_7TeV-madgraph-tauola 15  HiggsApp 8nh $fullprefix 1");

    sleep 500;

    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToEE_M-20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToMuMu_M-20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToTauTau_M-20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToEE_M-10To20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    print  "done with V+jets.\n";
    
    sleep 500;

    print  "submitting dibosons...\n";
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 WWTo2L2Nu_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 GluGluToWWTo4L_TuneZ2_7TeV-gg2ww-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 WZTo3LNu_TuneZ2_7TeV-pythia6 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 PhotonVJets_7TeV-madgraph 15 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py Spring11_V2 ZZtoAnything_TuneZ2_7TeV-pythia6-tauola 15 HiggsApp 8nh $fullprefix 1");
    print  "done with dibosons.\n";
    
    sleep 500;

    print "\nDONE WITH MASS $mass GeV\n";
    print "<--------------------------\n";
    
}
