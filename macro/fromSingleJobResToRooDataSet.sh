#! /bin/sh 
# ./fromSingleJobResToRooDataSet.sh:
# it runs several scripts that bring from the results of the selection stored in the directory "results/"
# to the merged RooDataSets of the species used for the fit

lumi=$1

echo "MERGING TREES STEP..."
chmod u+x ./mergeTrees.sh
./mergeTrees.sh
echo "MERGING TREES STEP DONE."

echo "WEIGHTING TREES STEP..."
chmod u+x ./weightTrees.sh
./weightTrees.sh $lumi $lumi $lumi
echo "MERGING TREES STEP DONE."

echo "MERGING WEIGHTED TREES ACCORDING FIT SPECIES DEFINITION..."
chmod u+x ./mergeMultiSamples.sh
./mergeMultiSamples.sh
echo "MERGING WEIGHTED TREES ACCORDING FIT SPECIES DEFINITION DONE."

echo "MERGING FINAL STATES..."
chmod u+x ./mergeFinalStates.sh
./mergeFinalStates.sh
echo "MERGING FINAL STATES DONE."

#echo "CREATING FIT ROODATASETS STEP..."
#chmod u+x ./createFitDatasets.sh
#./createFitDatasets.sh $Hmass
#echo "CREATING FIT ROODATASETS DONE."
