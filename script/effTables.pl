#! /usr/bin/perl

use Getopt::Std;
getopts('m:');

$lumi = 100;
$HiggsMass = $opt_m;
print "creating tex table for Higgs mass $HiggsMass GeV/c^2 and lumi = $lumi pb-1....... \n";

$HiggsLogFile = "./h".$HiggsMass."-kFactor-Preselection.log";
$WWLogFile    = "./ww".$HiggsMass."-newkFactor-Preselection.log";
$ttbarLogFile = "./ttbar".$HiggsMass."-Preselection.log";
$DYLogFile    = "./DY".$HiggsMass."-Preselection.log";
$ZZLogFile    = "./ZZ".$HiggsMass."-Preselection.log";


$HiggsXsec = 0;                                    # mtop = 175; BR(W->lnu) = (0.1075+0.1057+0.1125) included 
if($HiggsMass==120) { $HiggsXsec=0.57846826; }
if($HiggsMass==130) { $HiggsXsec=1.0979959; }      # pb
if($HiggsMass==140) { $HiggsXsec=1.6286752; } 
if($HiggsMass==150) { $HiggsXsec=2.0408732; } 
if($HiggsMass==160) { $HiggsXsec=2.4120955; } 
if($HiggsMass==170) { $HiggsXsec=2.3339263; } 
if($HiggsMass==180) { $HiggsXsec=2.0502736; } 
if($HiggsMass==190) { $HiggsXsec=1.5568427; } 
if($HiggsMass==200) { $HiggsXsec=1.350654; }

$WWXsec    = 12.125;                               # BR(W->lnu) included
$ttbarXsec = 89.107612;                            # BR(W->lnu) included
$dyXsec    = 773.500;                              # inclusive
$ZZXsec    = 1.52;                                 # ll, lnu, nunu 

#init table
$textable = "eff".$HiggsMass.".tex";
open(TEXFILE,">$textable");
print TEXFILE "\\begin\{sidewaystable\}\n";
print TEXFILE "\\begin\{center\}\n";
print TEXFILE "\\begin\{tabular\}\[t\]\{|l|c|c|c|c|c|\}\n";
print TEXFILE "\\hline\n";
print TEXFILE "\\hline\n";
print TEXFILE "& & & & & \\\\ \n";
print TEXFILE "& \n";
print TEXFILE "\$H \\rightarrow WW \\rightarrow l\^\+l\^\-n\\nu\$  \&\n";
print TEXFILE "\$WW \\rightarrow l\^\+l\^\-n\\nu\$  \&\n";
print TEXFILE "\$t \\bar\{t\} \\rightarrow b\\bar b l\^\+l\^\-n\\nu \$ \&\n";
print TEXFILE "Drell-Yan (inclusive) \&\n";
print TEXFILE "\$ZZ \\rightarrow\$ 4 leptons\n";
print TEXFILE "\\\\\n";
print TEXFILE "\\hline\n";
print TEXFILE "\\hline\n";

$Higgs_MCtruth=0;
$Higgs_trigger = 0;
$Higgs_nRecoEle = 0;
$Higgs_twoGoodRec=0;
$Higgs_jetVeto-0;
$Higgs_eleID=0;
$Higgs_trackerIsol=0;
$Higgs_hcalIsol=0;
$Higgs_presel=0;
$Higgs_MET=0;
$Higgs_deltaPhi=0;
$Higgs_eleInvMass=0;
$Higgs_maxPtEle=0;
$Higgs_minPtEle=0;
$Higgs_detaLeptons=0;
$Higgs_final=0;

open(HIGGSLOGFILE,$HiggsLogFile) or die "cannot open $HiggsLogFile !\n";
while($row=<HIGGSLOGFILE>) {
    # normalization
    if($row=~/\*\s+MCtruth\:\s+(\S+)$/) {
	if($1!=0) {$Higgs_MCtruth=$1;}
    }

    if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	if($1!=0) {$Higgs_trigger=$1;}
    }

    if($row=~/\*\s+nRecoEle\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_nRecoEle=$1;}
    }

    if($row=~/\*\s+twoGoodRec\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_twoGoodRec=$1;}
    }

    if($row=~/\*\s+jetVeto\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_jetVeto=$1;}
    }

    if($row=~/\*\s+eleID\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_eleID=$1;}
    }

    if($row=~/\*\s+trackerIsol\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_trackerIsol=$1;}
    }

    if($row=~/\*\s+hcalIsol\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_hcalIsol=$1;}
    }

    if($row=~/\*\s+presel\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_presel=$1;}
    }

    if($row=~/\*\s+MET\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_MET=$1;}
    }

    if($row=~/\*\s+deltaPhi\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_deltaPhi=$1;}
    }

    if($row=~/\*\s+eleInvMass\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_eleInvMass=$1;}
    }

    if($row=~/\*\s+maxPtEle\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_maxPtEle=$1;}
    }

    if($row=~/\*\s+minPtEle\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_minPtEle=$1;}
    }

    if($row=~/\*\s+detaLeptons\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_detaLeptons=$1;}
    }

    if($row=~/\*\s+final\:\s+(\S+)$/) {
        if($1!=0) {$Higgs_final=$1;}
    }
	
}


$WW_MCtruth=0;
$WW_trigger = 0;
$WW_nRecoEle = 0;
$WW_twoGoodRec=0;
$WW_jetVeto-0;
$WW_eleID=0;
$WW_trackerIsol=0;
$WW_hcalIsol=0;
$WW_presel=0;
$WW_MET=0;
$WW_deltaPhi=0;
$WW_eleInvMass=0;
$WW_maxPtEle=0;
$WW_minPtEle=0;
$WW_detaLeptons=0;
$WW_final=0;

open(WWLOGFILE,$WWLogFile) or die "cannot open $WWLogFile !\n";
while($row=<WWLOGFILE>) {
    # normalization
    if($row=~/\*\s+MCtruth\:\s+(\S+)$/) {
	if($1!=0) {$WW_MCtruth=$1;}
    }

    if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	if($1!=0) {$WW_trigger=$1;}
    }

    if($row=~/\*\s+nRecoEle\:\s+(\S+)$/) {
        if($1!=0) {$WW_nRecoEle=$1;}
    }

    if($row=~/\*\s+twoGoodRec\:\s+(\S+)$/) {
        if($1!=0) {$WW_twoGoodRec=$1;}
    }

    if($row=~/\*\s+jetVeto\:\s+(\S+)$/) {
        if($1!=0) {$WW_jetVeto=$1;}
    }

    if($row=~/\*\s+eleID\:\s+(\S+)$/) {
        if($1!=0) {$WW_eleID=$1;}
    }

    if($row=~/\*\s+trackerIsol\:\s+(\S+)$/) {
        if($1!=0) {$WW_trackerIsol=$1;}
    }

    if($row=~/\*\s+hcalIsol\:\s+(\S+)$/) {
        if($1!=0) {$WW_hcalIsol=$1;}
    }

    if($row=~/\*\s+presel\:\s+(\S+)$/) {
        if($1!=0) {$WW_presel=$1;}
    }

    if($row=~/\*\s+MET\:\s+(\S+)$/) {
        if($1!=0) {$WW_MET=$1;}
    }

    if($row=~/\*\s+deltaPhi\:\s+(\S+)$/) {
        if($1!=0) {$WW_deltaPhi=$1;}
    }

    if($row=~/\*\s+eleInvMass\:\s+(\S+)$/) {
        if($1!=0) {$WW_eleInvMass=$1;}
    }

    if($row=~/\*\s+maxPtEle\:\s+(\S+)$/) {
        if($1!=0) {$WW_maxPtEle=$1;}
    }

    if($row=~/\*\s+minPtEle\:\s+(\S+)$/) {
        if($1!=0) {$WW_minPtEle=$1;}
    }

    if($row=~/\*\s+detaLeptons\:\s+(\S+)$/) {
        if($1!=0) {$WW_detaLeptons=$1;}
    }

    if($row=~/\*\s+final\:\s+(\S+)$/) {
        if($1!=0) {$WW_final=$1;}
    }
	
}

$ttbar_MCtruth=0;
$ttbar_trigger = 0;
$ttbar_nRecoEle = 0;
$ttbar_twoGoodRec=0;
$ttbar_jetVeto-0;
$ttbar_eleID=0;
$ttbar_trackerIsol=0;
$ttbar_hcalIsol=0;
$ttbar_presel=0;
$ttbar_MET=0;
$ttbar_deltaPhi=0;
$ttbar_eleInvMass=0;
$ttbar_maxPtEle=0;
$ttbar_minPtEle=0;
$ttbar_detaLeptons=0;
$ttbar_final=0;

open(TTBARLOGFILE,$ttbarLogFile) or die "cannot open $ttbarLogFile !\n";
while($row=<TTBARLOGFILE>) {
    # normalization
    if($row=~/\*\s+MCtruth\:\s+(\S+)$/) {
	if($1!=0) {$ttbar_MCtruth=$1;}
    }

    if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	if($1!=0) {$ttbar_trigger=$1;}
    }

    if($row=~/\*\s+nRecoEle\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_nRecoEle=$1;}
    }

    if($row=~/\*\s+twoGoodRec\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_twoGoodRec=$1;}
    }

    if($row=~/\*\s+jetVeto\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_jetVeto=$1;}
    }

    if($row=~/\*\s+eleID\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_eleID=$1;}
    }

    if($row=~/\*\s+trackerIsol\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_trackerIsol=$1;}
    }

    if($row=~/\*\s+hcalIsol\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_hcalIsol=$1;}
    }

    if($row=~/\*\s+presel\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_presel=$1;}
    }

    if($row=~/\*\s+MET\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_MET=$1;}
    }

    if($row=~/\*\s+deltaPhi\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_deltaPhi=$1;}
    }

    if($row=~/\*\s+eleInvMass\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_eleInvMass=$1;}
    }

    if($row=~/\*\s+maxPtEle\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_maxPtEle=$1;}
    }

    if($row=~/\*\s+minPtEle\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_minPtEle=$1;}
    }

    if($row=~/\*\s+detaLeptons\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_detaLeptons=$1;}
    }

    if($row=~/\*\s+final\:\s+(\S+)$/) {
        if($1!=0) {$ttbar_final=$1;}
    }
	
}




$dy_MCtruth=0;
$dy_trigger = 0;
$dy_nRecoEle = 0;
$dy_twoGoodRec=0;
$dy_jetVeto-0;
$dy_eleID=0;
$dy_trackerIsol=0;
$dy_hcalIsol=0;
$dy_presel=0;
$dy_MET=0;
$dy_deltaPhi=0;
$dy_eleInvMass=0;
$dy_maxPtEle=0;
$dy_minPtEle=0;
$dy_detaLeptons=0;
$dy_final=0;

open(DYLOGFILE,$DYLogFile) or die "cannot open $DYLogFile !\n";
while($row=<DYLOGFILE>) {
    # normalization
    if($row=~/\*\s+MCtruth\:\s+(\S+)$/) {
	if($1!=0) {$dy_MCtruth=$1;}
    }

    if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	if($1!=0) {$dy_trigger=$1;}
    }

    if($row=~/\*\s+nRecoEle\:\s+(\S+)$/) {
        if($1!=0) {$dy_nRecoEle=$1;}
    }

    if($row=~/\*\s+twoGoodRec\:\s+(\S+)$/) {
        if($1!=0) {$dy_twoGoodRec=$1;}
    }

    if($row=~/\*\s+jetVeto\:\s+(\S+)$/) {
        if($1!=0) {$dy_jetVeto=$1;}
    }

    if($row=~/\*\s+eleID\:\s+(\S+)$/) {
        if($1!=0) {$dy_eleID=$1;}
    }

    if($row=~/\*\s+trackerIsol\:\s+(\S+)$/) {
        if($1!=0) {$dy_trackerIsol=$1;}
    }

    if($row=~/\*\s+hcalIsol\:\s+(\S+)$/) {
        if($1!=0) {$dy_hcalIsol=$1;}
    }

    if($row=~/\*\s+presel\:\s+(\S+)$/) {
        if($1!=0) {$dy_presel=$1;}
    }

    if($row=~/\*\s+MET\:\s+(\S+)$/) {
        if($1!=0) {$dy_MET=$1;}
    }

    if($row=~/\*\s+deltaPhi\:\s+(\S+)$/) {
        if($1!=0) {$dy_deltaPhi=$1;}
    }

    if($row=~/\*\s+eleInvMass\:\s+(\S+)$/) {
        if($1!=0) {$dy_eleInvMass=$1;}
    }

    if($row=~/\*\s+maxPtEle\:\s+(\S+)$/) {
        if($1!=0) {$dy_maxPtEle=$1;}
    }

    if($row=~/\*\s+minPtEle\:\s+(\S+)$/) {
        if($1!=0) {$dy_minPtEle=$1;}
    }

    if($row=~/\*\s+detaLeptons\:\s+(\S+)$/) {
        if($1!=0) {$dy_detaLeptons=$1;}
    }

    if($row=~/\*\s+final\:\s+(\S+)$/) {
        if($1!=0) {$dy_final=$1;}
    }
	
}












$ZZ_MCtruth=0;
$ZZ_trigger = 0;
$ZZ_nRecoEle = 0;
$ZZ_twoGoodRec=0;
$ZZ_jetVeto-0;
$ZZ_eleID=0;
$ZZ_trackerIsol=0;
$ZZ_hcalIsol=0;
$ZZ_presel=0;
$ZZ_MET=0;
$ZZ_deltaPhi=0;
$ZZ_eleInvMass=0;
$ZZ_maxPtEle=0;
$ZZ_minPtEle=0;
$ZZ_detaLeptons=0;
$ZZ_final=0;

open(ZZLOGFILE,$ZZLogFile) or die "cannot open $ZZLogFile !\n";
while($row=<ZZLOGFILE>) {
    # normalization
    if($row=~/\*\s+MCtruth\:\s+(\S+)$/) {
	if($1!=0) {$ZZ_MCtruth=$1;}
    }

    if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	if($1!=0) {$ZZ_trigger=$1;}
    }

    if($row=~/\*\s+nRecoEle\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_nRecoEle=$1;}
    }

    if($row=~/\*\s+twoGoodRec\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_twoGoodRec=$1;}
    }

    if($row=~/\*\s+jetVeto\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_jetVeto=$1;}
    }

    if($row=~/\*\s+eleID\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_eleID=$1;}
    }

    if($row=~/\*\s+trackerIsol\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_trackerIsol=$1;}
    }

    if($row=~/\*\s+hcalIsol\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_hcalIsol=$1;}
    }

    if($row=~/\*\s+presel\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_presel=$1;}
    }

    if($row=~/\*\s+MET\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_MET=$1;}
    }

    if($row=~/\*\s+deltaPhi\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_deltaPhi=$1;}
    }

    if($row=~/\*\s+eleInvMass\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_eleInvMass=$1;}
    }

    if($row=~/\*\s+maxPtEle\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_maxPtEle=$1;}
    }

    if($row=~/\*\s+minPtEle\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_minPtEle=$1;}
    }

    if($row=~/\*\s+detaLeptons\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_detaLeptons=$1;}
    }

    if($row=~/\*\s+final\:\s+(\S+)$/) {
        if($1!=0) {$ZZ_final=$1;}
    }
	
}




# --- HLT ---   
print TEXFILE "HLT &\n";
$decimals = 0;
$n_Higgs_trigger= sprintf("%.4g", $Higgs_trigger/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_trigger = sprintf("%.0f", 100 * $Higgs_trigger/$Higgs_MCtruth);

$n_WW_trigger= sprintf("%.4g", $WW_trigger/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_trigger = sprintf("%.0f", 100 * $WW_trigger/$WW_MCtruth);

$n_ttbar_trigger= sprintf("%.4g", $ttbar_trigger/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_trigger = sprintf("%.0f", 100 * $ttbar_trigger/$ttbar_MCtruth);

$n_dy_trigger= sprintf("%.4g", $dy_trigger/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_trigger = sprintf("%.0f", 100 * $dy_trigger/$dy_MCtruth);

$n_ZZ_trigger= sprintf("%.4g", $ZZ_trigger/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_trigger = sprintf("%.0f", 100 * $ZZ_trigger/$ZZ_MCtruth);

print TEXFILE "$n_Higgs_trigger ($eff_Higgs_trigger \\%)&\n";
print TEXFILE "$n_WW_trigger ($eff_WW_trigger \\%)&\n";
print TEXFILE "$n_ttbar_trigger ($eff_ttbar_trigger \\%)&\n";
print TEXFILE "$n_dy_trigger ($eff_dy_trigger \\%)&\n";
print TEXFILE "$n_ZZ_trigger ($eff_ZZ_trigger \\%)\\\\\n";




# --- exactly 2 reco electrons ---
print TEXFILE "\$e\^+e\^-\$ requirement &\n";
$n_Higgs_nRecoEle= sprintf("%.4g", $Higgs_nRecoEle/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_nRecoEle = sprintf("%.0f", 100 * $Higgs_nRecoEle/$Higgs_trigger);

$n_WW_nRecoEle= sprintf("%.4g", $WW_nRecoEle/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_nRecoEle = sprintf("%.0f", 100 * $WW_nRecoEle/$WW_trigger);

$n_ttbar_nRecoEle= sprintf("%.4g", $ttbar_nRecoEle/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_nRecoEle = sprintf("%.0f", 100 * $ttbar_nRecoEle/$ttbar_trigger);

$n_dy_nRecoEle= sprintf("%.4g", $dy_nRecoEle/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_nRecoEle = sprintf("%.0f", 100 * $dy_nRecoEle/$dy_trigger);

$n_ZZ_nRecoEle= sprintf("%.4g", $ZZ_nRecoEle/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_nRecoEle = sprintf("%.0f", 100 * $ZZ_nRecoEle/$ZZ_trigger);

print TEXFILE "$n_Higgs_nRecoEle ($eff_Higgs_nRecoEle \\%)&\n";
print TEXFILE "$n_WW_nRecoEle ($eff_WW_nRecoEle \\%)&\n";
print TEXFILE "$n_ttbar_nRecoEle ($eff_ttbar_nRecoEle \\%)&\n";
print TEXFILE "$n_dy_nRecoEle ($eff_dy_nRecoEle \\%)&\n";
print TEXFILE "$n_ZZ_nRecoEle ($eff_ZZ_nRecoEle \\%)\\\\\n";





# --- electron acceptance ---
print TEXFILE "electron \$|\\eta|<2.5\,\~p_\{T\}>\$ 10 GeV/c &\n";
$n_Higgs_twoGoodRec= sprintf("%.4g", $Higgs_twoGoodRec/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_twoGoodRec = sprintf("%.0f", 100 * $Higgs_twoGoodRec/$Higgs_nRecoEle);

$n_WW_twoGoodRec= sprintf("%.4g", $WW_twoGoodRec/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_twoGoodRec = sprintf("%.0f", 100 * $WW_twoGoodRec/$WW_nRecoEle);

$n_ttbar_twoGoodRec= sprintf("%.4g", $ttbar_twoGoodRec/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_twoGoodRec = sprintf("%.0f", 100 * $ttbar_twoGoodRec/$ttbar_nRecoEle);

$n_dy_twoGoodRec= sprintf("%.4g", $dy_twoGoodRec/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_twoGoodRec = sprintf("%.0f", 100 * $dy_twoGoodRec/$dy_nRecoEle);

$n_ZZ_twoGoodRec= sprintf("%.4g", $ZZ_twoGoodRec/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_twoGoodRec = sprintf("%.0f", 100 * $ZZ_twoGoodRec/$ZZ_nRecoEle);

print TEXFILE "$n_Higgs_twoGoodRec ($eff_Higgs_twoGoodRec \\%)&\n";
print TEXFILE "$n_WW_twoGoodRec ($eff_WW_twoGoodRec \\%)&\n";
print TEXFILE "$n_ttbar_twoGoodRec ($eff_ttbar_twoGoodRec \\%)&\n";
print TEXFILE "$n_dy_twoGoodRec ($eff_dy_twoGoodRec \\%)&\n";
print TEXFILE "$n_ZZ_twoGoodRec ($eff_ZZ_twoGoodRec \\%)\\\\\n";








# --- electron ID ---
print TEXFILE "electron ID &\n";
$n_Higgs_eleID= sprintf("%.4g", $Higgs_eleID/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_eleID = sprintf("%.0f", 100 * $Higgs_eleID/$Higgs_twoGoodRec);

$n_WW_eleID= sprintf("%.4g", $WW_eleID/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_eleID = sprintf("%.0f", 100 * $WW_eleID/$WW_twoGoodRec);

$n_ttbar_eleID= sprintf("%.4g", $ttbar_eleID/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_eleID = sprintf("%.0f", 100 * $ttbar_eleID/$ttbar_twoGoodRec);

$n_dy_eleID= sprintf("%.4g", $dy_eleID/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_eleID = sprintf("%.0f", 100 * $dy_eleID/$dy_twoGoodRec);

$n_ZZ_eleID= sprintf("%.4g", $ZZ_eleID/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_eleID = sprintf("%.0f", 100 * $ZZ_eleID/$ZZ_twoGoodRec);

print TEXFILE "$n_Higgs_eleID ($eff_Higgs_eleID \\%)&\n";
print TEXFILE "$n_WW_eleID ($eff_WW_eleID \\%)&\n";
print TEXFILE "$n_ttbar_eleID ($eff_ttbar_eleID \\%)&\n";
print TEXFILE "$n_dy_eleID ($eff_dy_eleID \\%)&\n";
print TEXFILE "$n_ZZ_eleID ($eff_ZZ_eleID \\%)\\\\\n";









# --- tracker isolation ---
print TEXFILE "tracker isolation &\n";
$n_Higgs_trackerIsol= sprintf("%.4g", $Higgs_trackerIsol/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_trackerIsol = sprintf("%.0f", 100 * $Higgs_trackerIsol/$Higgs_eleID);

$n_WW_trackerIsol= sprintf("%.4g", $WW_trackerIsol/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_trackerIsol = sprintf("%.0f", 100 * $WW_trackerIsol/$WW_eleID);

$n_ttbar_trackerIsol= sprintf("%.4g", $ttbar_trackerIsol/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_trackerIsol = sprintf("%.0f", 100 * $ttbar_trackerIsol/$ttbar_eleID);

$n_dy_trackerIsol= sprintf("%.4g", $dy_trackerIsol/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_trackerIsol = sprintf("%.0f", 100 * $dy_trackerIsol/$dy_eleID);

$n_ZZ_trackerIsol= sprintf("%.4g", $ZZ_trackerIsol/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_trackerIsol = sprintf("%.0f", 100 * $ZZ_trackerIsol/$ZZ_eleID);

print TEXFILE "$n_Higgs_trackerIsol ($eff_Higgs_trackerIsol \\%)&\n";
print TEXFILE "$n_WW_trackerIsol ($eff_WW_trackerIsol \\%)&\n";
print TEXFILE "$n_ttbar_trackerIsol ($eff_ttbar_trackerIsol \\%)&\n";
print TEXFILE "$n_dy_trackerIsol ($eff_dy_trackerIsol \\%)&\n";
print TEXFILE "$n_ZZ_trackerIsol ($eff_ZZ_trackerIsol \\%)\\\\\n";






# --- HCAL isolation ---
print TEXFILE "HCAL isolation &\n";
$n_Higgs_hcalIsol= sprintf("%.4g", $Higgs_hcalIsol/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_hcalIsol = sprintf("%.0f", 100 * $Higgs_hcalIsol/$Higgs_trackerIsol);

$n_WW_hcalIsol= sprintf("%.4g", $WW_hcalIsol/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_hcalIsol = sprintf("%.0f", 100 * $WW_hcalIsol/$WW_trackerIsol);

$n_ttbar_hcalIsol= sprintf("%.4g", $ttbar_hcalIsol/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_hcalIsol = sprintf("%.0f", 100 * $ttbar_hcalIsol/$ttbar_trackerIsol);

$n_dy_hcalIsol= sprintf("%.4g", $dy_hcalIsol/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_hcalIsol = sprintf("%.0f", 100 * $dy_hcalIsol/$dy_trackerIsol);

$n_ZZ_hcalIsol= sprintf("%.4g", $ZZ_hcalIsol/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_hcalIsol = sprintf("%.0f", 100 * $ZZ_hcalIsol/$ZZ_trackerIsol);

print TEXFILE "$n_Higgs_hcalIsol ($eff_Higgs_hcalIsol \\%)&\n";
print TEXFILE "$n_WW_hcalIsol ($eff_WW_hcalIsol \\%)&\n";
print TEXFILE "$n_ttbar_hcalIsol ($eff_ttbar_hcalIsol \\%)&\n";
print TEXFILE "$n_dy_hcalIsol ($eff_dy_hcalIsol \\%)&\n";
print TEXFILE "$n_ZZ_hcalIsol ($eff_ZZ_hcalIsol \\%)\\\\\n";





# --- PRESELECTION ---
print TEXFILE "Pre-selection &\n";
$n_Higgs_presel= sprintf("%.4g", $Higgs_presel/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_presel = sprintf("%.0f", 100 * $Higgs_presel/$Higgs_hcalIsol);

$n_WW_presel= sprintf("%.4g", $WW_presel/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_presel = sprintf("%.0f", 100 * $WW_presel/$WW_hcalIsol);

$n_ttbar_presel= sprintf("%.4g", $ttbar_presel/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_presel = sprintf("%.0f", 100 * $ttbar_presel/$ttbar_hcalIsol);

$n_dy_presel= sprintf("%.4g", $dy_presel/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_presel = sprintf("%.0f", 100 * $dy_presel/$dy_hcalIsol);

$n_ZZ_presel= sprintf("%.4g", $ZZ_presel/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_presel = sprintf("%.0f", 100 * $ZZ_presel/$ZZ_hcalIsol);

print TEXFILE "$n_Higgs_presel ($eff_Higgs_presel \\%)&\n";
print TEXFILE "$n_WW_presel ($eff_WW_presel \\%)&\n";
print TEXFILE "$n_ttbar_presel ($eff_ttbar_presel \\%)&\n";
print TEXFILE "$n_dy_presel ($eff_dy_presel \\%)&\n";
print TEXFILE "$n_ZZ_presel ($eff_ZZ_presel \\%)\\\\\n";





# --- Central Jet Veto ---
print TEXFILE "Central Jet Veto &\n";
$n_Higgs_jetVeto= sprintf("%.4g", $Higgs_jetVeto/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_jetVeto = sprintf("%.0f", 100 * $Higgs_jetVeto/$Higgs_presel);

$n_WW_jetVeto= sprintf("%.4g", $WW_jetVeto/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_jetVeto = sprintf("%.0f", 100 * $WW_jetVeto/$WW_presel);

$n_ttbar_jetVeto= sprintf("%.4g", $ttbar_jetVeto/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_jetVeto = sprintf("%.0f", 100 * $ttbar_jetVeto/$ttbar_presel);

$n_dy_jetVeto= sprintf("%.4g", $dy_jetVeto/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_jetVeto = sprintf("%.0f", 100 * $dy_jetVeto/$dy_presel);

$n_ZZ_jetVeto= sprintf("%.4g", $ZZ_jetVeto/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_jetVeto = sprintf("%.0f", 100 * $ZZ_jetVeto/$ZZ_presel);

print TEXFILE "$n_Higgs_jetVeto ($eff_Higgs_jetVeto \\%)&\n";
print TEXFILE "$n_WW_jetVeto ($eff_WW_jetVeto \\%)&\n";
print TEXFILE "$n_ttbar_jetVeto ($eff_ttbar_jetVeto \\%)&\n";
print TEXFILE "$n_dy_jetVeto ($eff_dy_jetVeto \\%)&\n";
print TEXFILE "$n_ZZ_jetVeto ($eff_ZZ_jetVeto \\%)\\\\\n";








# --- MET ---
print TEXFILE "\$\\met \$ &\n";
$n_Higgs_MET= sprintf("%.4g", $Higgs_MET/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_MET = sprintf("%.0f", 100 * $Higgs_MET/$Higgs_jetVeto);

$n_WW_MET= sprintf("%.4g", $WW_MET/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_MET = sprintf("%.0f", 100 * $WW_MET/$WW_jetVeto);

$n_ttbar_MET= sprintf("%.4g", $ttbar_MET/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_MET = sprintf("%.0f", 100 * $ttbar_MET/$ttbar_jetVeto);

$n_dy_MET= sprintf("%.4g", $dy_MET/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_MET = sprintf("%.0f", 100 * $dy_MET/$dy_jetVeto);

$n_ZZ_MET= sprintf("%.4g", $ZZ_MET/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_MET = sprintf("%.0f", 100 * $ZZ_MET/$ZZ_jetVeto);

print TEXFILE "$n_Higgs_MET ($eff_Higgs_MET \\%)&\n";
print TEXFILE "$n_WW_MET ($eff_WW_MET \\%)&\n";
print TEXFILE "$n_ttbar_MET ($eff_ttbar_MET \\%)&\n";
print TEXFILE "$n_dy_MET ($eff_dy_MET \\%)&\n";
print TEXFILE "$n_ZZ_MET ($eff_ZZ_MET \\%)\\\\\n";






# --- delta phi ---
print TEXFILE "\$\\delphill\$ &\n";
$n_Higgs_deltaPhi= sprintf("%.4g", $Higgs_deltaPhi/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_deltaPhi = sprintf("%.0f", 100 * $Higgs_deltaPhi/$Higgs_MET);

$n_WW_deltaPhi= sprintf("%.4g", $WW_deltaPhi/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_deltaPhi = sprintf("%.0f", 100 * $WW_deltaPhi/$WW_MET);

$n_ttbar_deltaPhi= sprintf("%.4g", $ttbar_deltaPhi/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_deltaPhi = sprintf("%.0f", 100 * $ttbar_deltaPhi/$ttbar_MET);

$n_dy_deltaPhi= sprintf("%.4g", $dy_deltaPhi/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_deltaPhi = sprintf("%.0f", 100 * $dy_deltaPhi/$dy_MET);

$n_ZZ_deltaPhi= sprintf("%.4g", $ZZ_deltaPhi/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_deltaPhi = sprintf("%.0f", 100 * $ZZ_deltaPhi/$ZZ_MET);

print TEXFILE "$n_Higgs_deltaPhi ($eff_Higgs_deltaPhi \\%)&\n";
print TEXFILE "$n_WW_deltaPhi ($eff_WW_deltaPhi \\%)&\n";
print TEXFILE "$n_ttbar_deltaPhi ($eff_ttbar_deltaPhi \\%)&\n";
print TEXFILE "$n_dy_deltaPhi ($eff_dy_deltaPhi \\%)&\n";
print TEXFILE "$n_ZZ_deltaPhi ($eff_ZZ_deltaPhi \\%)\\\\\n";








# --- e+e- Invariant Mass ---
print TEXFILE "\$ \\mll \$ &\n";
$n_Higgs_eleInvMass= sprintf("%.4g", $Higgs_eleInvMass/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_eleInvMass = sprintf("%.0f", 100 * $Higgs_eleInvMass/$Higgs_deltaPhi);

$n_WW_eleInvMass= sprintf("%.4g", $WW_eleInvMass/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_eleInvMass = sprintf("%.0f", 100 * $WW_eleInvMass/$WW_deltaPhi);

$n_ttbar_eleInvMass= sprintf("%.4g", $ttbar_eleInvMass/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_eleInvMass = sprintf("%.0f", 100 * $ttbar_eleInvMass/$ttbar_deltaPhi);

$n_dy_eleInvMass= sprintf("%.4g", $dy_eleInvMass/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_eleInvMass = sprintf("%.0f", 100 * $dy_eleInvMass/$dy_deltaPhi);

$n_ZZ_eleInvMass= sprintf("%.4g", $ZZ_eleInvMass/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_eleInvMass = sprintf("%.0f", 100 * $ZZ_eleInvMass/$ZZ_deltaPhi);

print TEXFILE "$n_Higgs_eleInvMass ($eff_Higgs_eleInvMass \\%)&\n";
print TEXFILE "$n_WW_eleInvMass ($eff_WW_eleInvMass \\%)&\n";
print TEXFILE "$n_ttbar_eleInvMass ($eff_ttbar_eleInvMass \\%)&\n";
print TEXFILE "$n_dy_eleInvMass ($eff_dy_eleInvMass \\%)&\n";
print TEXFILE "$n_ZZ_eleInvMass ($eff_ZZ_eleInvMass \\%)\\\\\n";











# --- max Pt Electron Pt cut ---
print TEXFILE "\$\\ptlmax \$ &\n";
$n_Higgs_maxPtEle= sprintf("%.4g", $Higgs_maxPtEle/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_maxPtEle = sprintf("%.0f", 100 * $Higgs_maxPtEle/$Higgs_eleInvMass);

$n_WW_maxPtEle= sprintf("%.4g", $WW_maxPtEle/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_maxPtEle = sprintf("%.0f", 100 * $WW_maxPtEle/$WW_eleInvMass);

$n_ttbar_maxPtEle= sprintf("%.4g", $ttbar_maxPtEle/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_maxPtEle = sprintf("%.0f", 100 * $ttbar_maxPtEle/$ttbar_eleInvMass);

$n_dy_maxPtEle= sprintf("%.4g", $dy_maxPtEle/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_maxPtEle = sprintf("%.0f", 100 * $dy_maxPtEle/$dy_eleInvMass);

if($ZZ_eleInvMass==0) {$ZZ_eleInvMass=1;}
$n_ZZ_maxPtEle= sprintf("%.4g", $ZZ_maxPtEle/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_maxPtEle = sprintf("%.0f", 100 * $ZZ_maxPtEle/$ZZ_eleInvMass);

print TEXFILE "$n_Higgs_maxPtEle ($eff_Higgs_maxPtEle \\%)&\n";
print TEXFILE "$n_WW_maxPtEle ($eff_WW_maxPtEle \\%)&\n";
print TEXFILE "$n_ttbar_maxPtEle ($eff_ttbar_maxPtEle \\%)&\n";
print TEXFILE "$n_dy_maxPtEle ($eff_dy_maxPtEle \\%)&\n";
print TEXFILE "$n_ZZ_maxPtEle ($eff_ZZ_maxPtEle \\%)\\\\\n";
















# --- min Pt Electron Pt cut ---
print TEXFILE "\$\\ptlmin \$ &\n";
$n_Higgs_minPtEle= sprintf("%.4g", $Higgs_minPtEle/$Higgs_MCtruth * $HiggsXsec * $lumi);
$eff_Higgs_minPtEle = sprintf("%.0f", 100 * $Higgs_minPtEle/$Higgs_maxPtEle);

$n_WW_minPtEle= sprintf("%.4g", $WW_minPtEle/$WW_MCtruth * $WWXsec * $lumi);
$eff_WW_minPtEle = sprintf("%.0f", 100 * $WW_minPtEle/$WW_maxPtEle);

$n_ttbar_minPtEle= sprintf("%.4g", $ttbar_minPtEle/$ttbar_MCtruth * $ttbarXsec * $lumi);
$eff_ttbar_minPtEle = sprintf("%.0f", 100 * $ttbar_minPtEle/$ttbar_maxPtEle);

$n_dy_minPtEle= sprintf("%.4g", $dy_minPtEle/$dy_MCtruth * $dyXsec * $lumi);
$eff_dy_minPtEle = sprintf("%.0f", 100 * $dy_minPtEle/$dy_maxPtEle);

if($ZZ_maxPtEle==0) {$ZZ_maxPtEle=1;}
$n_ZZ_minPtEle= sprintf("%.4g", $ZZ_minPtEle/$ZZ_MCtruth * $ZZXsec * $lumi);
$eff_ZZ_minPtEle = sprintf("%.0f", 100 * $ZZ_minPtEle/$ZZ_maxPtEle);

print TEXFILE "$n_Higgs_minPtEle ($eff_Higgs_minPtEle \\%)&\n";
print TEXFILE "$n_WW_minPtEle ($eff_WW_minPtEle \\%)&\n";
print TEXFILE "$n_ttbar_minPtEle ($eff_ttbar_minPtEle \\%)&\n";
print TEXFILE "$n_dy_minPtEle ($eff_dy_minPtEle \\%)&\n";
print TEXFILE "$n_ZZ_minPtEle ($eff_ZZ_minPtEle \\%)\\\\\n";












# --- final efficiency ---
print TEXFILE "& & & & & \\\\ \n";
print TEXFILE "\\hline \n";
print TEXFILE "\\hline \n";
print TEXFILE "Total selection efficiency (\\%) &\n";
$eff_Higgs_final_u = $Higgs_final/$Higgs_MCtruth;
$erreff_Higgs_final_u = sqrt($eff_Higgs_final_u * (1-$eff_Higgs_final_u)/$Higgs_MCtruth);
$eff_Higgs_final = sprintf("%.3f", 100 * $eff_Higgs_final_u);
$erreff_Higgs_final = sprintf("%.3f", 100 * $erreff_Higgs_final_u);

$eff_WW_final_u = $WW_final/$WW_MCtruth;
$erreff_WW_final_u = sqrt($eff_WW_final_u * (1-$eff_WW_final_u)/$WW_MCtruth);
$eff_WW_final = sprintf("%.3f", 100 * $eff_WW_final_u);
$erreff_WW_final = sprintf("%.3f", 100 * $erreff_WW_final_u);

$eff_ttbar_final_u = $ttbar_final/$ttbar_MCtruth;
$erreff_ttbar_final_u = sqrt($eff_ttbar_final_u * (1-$eff_ttbar_final_u)/$ttbar_MCtruth);
$eff_ttbar_final = sprintf("%.3f", 100 * $eff_ttbar_final_u);
$erreff_ttbar_final = sprintf("%.3f", 100 * $erreff_ttbar_final_u);

$eff_dy_final_u = $dy_final/$dy_MCtruth;
$erreff_dy_final_u = sqrt($eff_dy_final_u * (1-$eff_dy_final_u)/$dy_MCtruth);
$eff_dy_final = sprintf("%.3f", 100 * $eff_dy_final_u);
$erreff_dy_final = sprintf("%.3f", 100 * $erreff_dy_final_u);

$eff_ZZ_final_u = $ZZ_final/$ZZ_MCtruth;
$erreff_ZZ_final_u = sqrt($eff_ZZ_final_u * (1-$eff_ZZ_final_u)/$ZZ_MCtruth);
$eff_ZZ_final = sprintf("%.3f", 100 * $eff_ZZ_final_u);
$erreff_ZZ_final = sprintf("%.3f", 100 * $erreff_ZZ_final_u);

print TEXFILE "$eff_Higgs_final \$\\pm\$ $erreff_Higgs_final &\n";
print TEXFILE "$eff_WW_final \$\\pm\$ $erreff_WW_final &\n";
print TEXFILE "$eff_ttbar_final \$\\pm\$ $erreff_ttbar_final &\n";
print TEXFILE "$eff_dy_final \$\\pm\$ $erreff_dy_final &\n";
print TEXFILE "$eff_ZZ_final \$\\pm\$ $erreff_ZZ_final \\\\ \n";




# --- expected events ---
print TEXFILE "& & & & & \\\\ \n";
print TEXFILE "\\hline \n";
print TEXFILE "Expected events &\n";

$n_Higgs_final = sprintf("%.3f", $eff_Higgs_final_u * $HiggsXsec * $lumi);
$nerr_Higgs_final = sprintf("%.3f", $erreff_Higgs_final_u * $HiggsXsec * $lumi);
$n_WW_final = sprintf("%.3f", $eff_WW_final_u * $WWXsec * $lumi);
$nerr_WW_final = sprintf("%.3f", $erreff_WW_final_u * $WWXsec * $lumi);
$n_ttbar_final = sprintf("%.3f", $eff_ttbar_final_u * $ttbarXsec * $lumi);
$nerr_ttbar_final = sprintf("%.3f", $erreff_ttbar_final_u * $ttbarXsec * $lumi);
$n_dy_final = sprintf("%.3f", $eff_dy_final_u * $dyXsec * $lumi);
$nerr_dy_final = sprintf("%.3f", $erreff_dy_final_u * $dyXsec * $lumi);
$n_ZZ_final = sprintf("%.3f", $eff_ZZ_final_u * $ZZXsec * $lumi);
$nerr_ZZ_final = sprintf("%.3f", $erreff_ZZ_final_u * $ZZXsec * $lumi);

print TEXFILE "$n_Higgs_final \$\\pm\$ $nerr_Higgs_final &\n";
print TEXFILE "$n_WW_final \$\\pm\$ $nerr_WW_final &\n";
print TEXFILE "$n_ttbar_final \$\\pm\$ $nerr_ttbar_final &\n";
print TEXFILE "$n_dy_final \$\\pm\$ $nerr_dy_final &\n";
print TEXFILE "$n_ZZ_final \$\\pm\$ $nerr_ZZ_final \n";









# final table

print TEXFILE "\\\\ \n";
print TEXFILE "& & & & & \\\\ \n";
print TEXFILE "\\hline \n";
print TEXFILE "\\hline \n";
print TEXFILE "\\end\{tabular\}\n";
print TEXFILE "\\caption\{\\emph\{Mass hypothesis: $HiggsMass GeV\/c\$^2\$. Number of expected events for an integrated luminosity of $lumi pb\$^\{-1\}\$ after each selection. The relative efficiencies with respect to the previous cut are given in percent within the brackets. The quoted errors include the statistical uncertainty only.\}\}\n";
print TEXFILE "\\label\{tab:eff$HiggsMass\}\n";
print TEXFILE "\\end\{center\}\n";
print TEXFILE "\\end\{sidewaystable\}\n";

print "Done. Include in your tex and compile.\n";


