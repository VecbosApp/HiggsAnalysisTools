#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('p:');
if($opt_p) {$prefix = $opt_p;}
else { die "usage: ./submitMassDependentAnalysis.pl -p <prefix>";}

#@masses = (120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600);
#@masses = (120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600);
@masses = (120);

for($i=0; $i<($#masses+1); $i++) {
    $mass = $masses[$i];
    print "-------------------------->\n";
    print "SUBMITTING MASS SELECTION: mH = $mass ...\n\n";
    open(MASSFILE,">config/higgs/higgsMass.txt");
    print MASSFILE "HiggsMass\t$mass\n";

    print  "submitting double ele...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW DoubleElectron 10 HiggsApp 8nh $prefix 0");
    print  "done with double ele. Sleeping 600s to recover.\n";
    
    sleep 600;
    
    print  "submitting double mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW DoubleMu 10 HiggsApp 8nh $prefix 0");
    print  "done with double mu. Sleeping 600s to recover.\n";

    sleep 600;

    print  "submitting mu eg...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW MuEG 10 HiggsApp 8nh $prefix 0");
    print  "done with mu eg. Sleeping 600s to recover.\n";

    sleep 600;
    
    print  "submitting single mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW SingleMu 5 HiggsApp 8nh $prefix 0");
    print  "done with single mu. Sleeping 600s to recover.\n";

    sleep 600;

    print  "submitting single ele...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW SingleElectron 10 HiggsApp 8nh $prefix 0");
    print  "done with single electron. Sleeping 600s to recover.\n";

    sleep 600;

    print  "submitting double ele...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW DoubleElectron_Run2011B 10 HiggsApp 8nh $prefix 0");
    print  "done with double ele. Sleeping 600s to recover.\n";
    
    sleep 600;
    
    print  "submitting double mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW DoubleMu_Run2011B 10 HiggsApp 8nh $prefix 0");
    print  "done with double mu. Sleeping 600s to recover.\n";

    sleep 600;

    print  "submitting mu eg...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW MuEG_Run2011B 10 HiggsApp 8nh $prefix 0");
    print  "done with mu eg. Sleeping 600s to recover.\n";

    sleep 600;
    
    print  "submitting single mu...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW SingleMu_Run2011B 5 HiggsApp 8nh $prefix 0");
    print  "done with single mu. Sleeping 600s to recover.\n";

    sleep 600;

    print  "submitting single ele...\n";
    system("python cmst3_submit_manyfilesperjob.py Data7TeVHWW SingleElectron_Run2011B 10 HiggsApp 8nh $prefix 0");
    print  "done with single electron.\n";
    
    print "\nDONE\n";
    print "<--------------------------";

}
