#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitMassDependentAnalysis.pl -p <prefix>";}

@masses = (120,130,140,150,160,170,200,300,400,500,600);

for($i=0; $i<($#masses+1); $i++) {
    $mass = $masses[$i];
    print "-------------------------->\n";
    print "SUBMITTING MASS SELECTION: mH = $mass ...\n\n";
    open(MASSFILE,">config/higgs/higgsMass.txt");
    print MASSFILE "HiggsMass\t$mass\n";

    $fullprefix=$prefix."/"."OptimMH$mass";

    print "submitting signals...\n";
    $higgsList = "GluGluToHToWWTo2L2Nu_M-".$mass;
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 1 HiggsApp 8nh $fullprefix 1");

    $higgsList = "GluGluToHToWWToLNuTauNu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 10 HiggsApp 8nh $fullprefix 1");

    $higgsList = "GluGluToHToWWToTauNuQQ_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 10 HiggsApp 8nh $fullprefix 1");

    $higgsList = "VBF_HToWWTo2L2Nu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 10 HiggsApp 8nh $fullprefix 1");

    $higgsList = "VBF_HToWWToLNuTauNu_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 10 HiggsApp 8nh $fullprefix 1");

    $higgsList =  "VBF_HToWWToTauNuQQ_M-".$mass."_7TeV-powheg-pythia6";
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 10 HiggsApp 8nh $fullprefix 1");

    print "done with signals.\n";

    print  "submitting top...\n";
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_s-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_t-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_tW-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py TTbar TTJets_TuneD6T 10 HiggsApp 8nh $fullprefix 1");
    print  "done with top.\n";
    
    print  "submitting V+jets...\n";
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToENu_TuneZ2 10 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToMuNu_TuneZ2 10 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToTauNu_TuneZ2_7TeV-pythia6-tauola 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py ZPYTHIA DYToEE_M-20_CT10_TuneZ2_PU 10 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py ZPYTHIA DYToMuMu_M-20_CT10_TuneZ2_PU 10 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py ZPYTHIA DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola 10 HiggsApp 8nh $fullprefix 1");
    print  "done with V+jets.\n";
    
    print  "submitting dibosons...\n";
    system("python cmst3_submit_manyfilesperjob.py DiBosons WWTo2L2Nu_TuneZ2 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py DiBosons GluGluToWWTo4L_TuneZ2_PU 10 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py DiBosons WZTo3LNu_TuneZ2 5 HiggsApp 8nh $fullprefix 1");
#    system("python cmst3_submit_manyfilesperjob.py DiBosons Wgamma 5 HiggsApp 8nh $fullprefix 1");
    system("python cmst3_submit_manyfilesperjob.py DiBosons ZZtoAnything_TuneZ2 5 HiggsApp 8nh $fullprefix 1");
    print  "done with dibosons.\n";
    
    print "\nDONE WITH MASS $mass GeV\n";
    print "<--------------------------\n";

}
