#!/bin/bash
# set -x

source "$CMSSW_BASE/src/HWWAnalysis/ShapeAnalysis/test/env.sh"

tailsf="_shape.txt"
tailof="_shape.txt"
cwd=$PWD
masses="110 115 120 125 130 135 140 145 150 155 160 170 180"
lumi=12.1

function usage() {
	echo "$( basename $0) -p <prefix> -t (to use 2D OF)"
}

prefix=
suffix2D=
while getopts "hp:t" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        p)
            prefix=$OPTARG
            ;;
        t)
            tailof="_shape_2D.txt"
            ;;
    esac
done

echo `shape-config.py`
eval `shape-config.py $lumi`

if [[ $prefix ]]; then
    echo "Running in $prefix"
    cd $prefix
fi

echo "Lumi $lumi"
head="hww-${lumi}fb.mH"
head2j="hww-${lumi}fb.mH"
hline=$(printf '%.0s-' {1..80})
combCmd='combineCards.py -S'
echo $PWD
for mass in $masses
do

	for jet in 0j 1j
	do
        echo $hline
        echo " Combining cards $mass $jet"
        echo $hline
		if [ -e $head$mass".of_"$jet$tailof ] && [ -e $head$mass".sf_"$jet$tailsf ] 
		then
			echo " ${mass}: of_${jet} sf_${jet} => comb_${jet}"
			$combCmd "of_"$jet"="$head$mass".of_"$jet$tailof  "sf_"$jet"="$head$mass".sf_"$jet$tailsf > $head$mass".comb_"$jet$tailof 
		fi
	done

	for fl in of sf
	do
                tailfl=
                if [[ $fl == "of" ]]; then
                    tailfl=$tailof
                else
                    tailfl=$tailsf
                fi
		file0="${head}${mass}.${fl}_0j${tailfl}"
		file1="${head}${mass}.${fl}_1j${tailfl}"
		cmb="${head}${mass}.comb_${fl}${tailfl}"

        echo $hline
		echo " Combining cards $mass $fl ( ${file0} + ${file1} )"
        echo $hline
		if [ -e "${file0}" ] && [ -e "${file1}" ] 
		then
			echo " ${mass}: ${fl}_0j ${fl}_1j => comb_${fl}"
			$combCmd "${fl}_0j=${file0}"  "${fl}_1j=${file1}"  > "${cmb}"
		fi
	done

	if [ -e $head$mass".of_0j"$tailof ] && [ -e $head$mass".sf_0j"$tailsf ] && [ -e $head$mass".of_1j"$tailof ] && [ -e $head$mass".sf_1j"$tailsf ]
	then
		echo " $mass: sf_0j of_0j of_1j sf_1j => comb_0j1j "
		$combCmd "of_0j="$head$mass".of_0j"$tailof  "sf_0j="$head$mass".sf_0j"$tailsf "of_1j="$head$mass".of_1j"$tailof  "sf_1j="$head$mass".sf_1j"$tailsf > $head$mass".comb_0j1j"$tailof
	fi

	if [ -e $head$mass".comb_0j"$tailof ] && [ -e $head$mass".comb_1j"$tailof ] && [ -e $head$mass".of_2j"$tailsf ] && [ -e $head$mass".sf_2j"$tailsf ]
	then
		echo " $mass: comb_0j comb_1j of_2j sf_2j => comb_0j1j2j "
		$combCmd HWW_0j=$head$mass".comb_0j"$tailof HWW_1j=$head$mass".comb_1j"$tailof HWW_2j_of=$head$mass".of_2j"$tailsf HWW_2j_sf=$head$mass".sf_2j"$tailsf > $head$mass.comb_0j1j2j$tailof
fi


done
echo "...Done"


