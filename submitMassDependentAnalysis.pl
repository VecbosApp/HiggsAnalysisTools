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
    system("python cmst3_submit_manyfilesperjob.py HiggsWW $higgsList 1 HiggsApp 8nh $fullprefix");
    print "done with signals.\n";

    print  "submitting top...\n";
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_s-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_t-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py SingleTop TToBLNu_TuneZ2_tW-channel_7TeV-madgraph 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py TTbar TTJets_TuneD6T 10 HiggsApp 8nh $fullprefix");
    print  "done with top.\n";
    
    print  "submitting V+jets...\n";
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToENu_TuneZ2 10 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToMuNu_TuneZ2 10 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py WPYTHIA WToTauNu_TuneZ2_7TeV-pythia6-tauola 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py ZPYTHIA DYJetsToLL_TuneD6T_M-10To50 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py ZPYTHIA DYJetsToLL_TuneD6T_M-50 5 HiggsApp 8nh $fullprefix");
    print  "done with V+jets.\n";
    
#    print  "submitting W+1jet alpgen...\n";
#    system("python cmst3_submit_manyfilesperjob.py WJetsALPGEN W1Jets_Pt0to100-alpgen 5 HiggsApp 8nh $fullprefix");
#    system("python cmst3_submit_manyfilesperjob.py WJetsALPGEN W1Jets_Pt100to300-alpgen 5 HiggsApp 8nh $fullprefix");
#    system("python cmst3_submit_manyfilesperjob.py WJetsALPGEN W1Jets_Pt300to800-alpgen 5 HiggsApp 8nh $fullprefix");
#    system("python cmst3_submit_manyfilesperjob.py WJetsALPGEN W1Jets_Pt800to1600-alpgen 5 HiggsApp 8nh $fullprefix");
#    print  "done with W+1jet alpgen.\n";
    
    print  "submitting dibosons...\n";
    system("python cmst3_submit_manyfilesperjob.py DiBosons WWTo2L2Nu_TuneZ2 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py DiBosons WZTo3LNu_TuneZ2 5 HiggsApp 8nh $fullprefix");
#    system("python cmst3_submit_manyfilesperjob.py DiBosons Wgamma 5 HiggsApp 8nh $fullprefix");
    system("python cmst3_submit_manyfilesperjob.py DiBosons ZZtoAnything_TuneZ2 5 HiggsApp 8nh $fullprefix");
    print  "done with dibosons.\n";
    
    print "\nDONE WITH MASS $mass GeV\n";
    print "<--------------------------";

}
