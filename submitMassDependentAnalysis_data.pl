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

    print  "submitting E/gamma...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeV dataset_eg_Sep3rdReReco 30 HiggsApp 8nh $fullprefix 0");
    system("python cmst3_submit_manyfilesperjob.py Data7TeV PDElectron_11pbTo40pb 10 HiggsApp 8nh $fullprefix 0");
    print  "done with E/gamma.\n";
    
    print  "submitting mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeV dataset_mu 20 HiggsApp 8nh $fullprefix 0");
    system("python cmst3_submit_manyfilesperjob.py Data7TeV PDMu_11pbTo40pb 20 HiggsApp 8nh $fullprefix 0");
    print  "done with mu.\n";
    
    print "\nDONE WITH MASS $mass GeV\n";
    print "<--------------------------";

}
