#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitMassDependentAnalysis.pl -p <prefix>";}

@masses = (120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600);

for($i=0; $i<($#masses+1); $i++) {
    $mass = $masses[$i];
    print "-------------------------->\n";
    print "SUBMITTING MASS SELECTION: mH = $mass ...\n\n";
    open(MASSFILE,">config/higgs/higgsMass.txt");
    print MASSFILE "HiggsMass\t$mass\n";

    $fullprefix=$prefix."/"."OptimMH$mass";

    print  "submitting double ele...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeV DoubleElectron 5 HiggsApp 8nh $fullprefix 0");
    print  "done with double ele.\n";
    
    print  "submitting double mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeV DoubleMu 5 HiggsApp 8nh $fullprefix 0");
    print  "done with double mu.\n";

    print  "submitting mu eg...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeV MuEG 5 HiggsApp 8nh $fullprefix 0");
    print  "done with mu eg.\n";
    
    print "\nDONE WITH MASS $mass GeV\n";
    print "<--------------------------";

}
