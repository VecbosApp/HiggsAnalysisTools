#! /bin/sh 
# ./fromSingleJobResToRooDataSetData.sh:

echo "===> STARTING... <==="

echo "MERGING TREES STEP..."
chmod u+x ./mergeTreesData.sh
./mergeTreesData.sh
echo "MERGING TREES STEP DONE."

echo "WEIGHTING TREES STEP..."
chmod u+x ./weightTreesData.sh
./weightTreesData.sh
echo "MERGING TREES STEP DONE."

echo "MERGING MULTIPLE DATASETS IN THE SAME FINAL STATE..."
chmod u+x ./mergeMultiSamplesData.sh
./mergeMultiSamplesData.sh
./mergeMultiSamplesDataSkim.sh
echo "MERGING MULTIPLE DATASETS IN THE SAME FINAL STATE DONE."

echo "MERGING FINAL STATES..."
chmod u+x ./mergeFinalStates.sh
./mergeFinalStates.sh
./mergeFinalStatesSkim.sh
echo "MERGING FINAL STATES DONE."

#echo "CREATING ROODATASETS FOR THE FIT..."
#chmod u+x createFitDatasetsData.sh
#./createFitDatasetsData.sh
#echo "CREATING ROODATASETS FOR THE FIT DONE."

echo "===> VERY DONE. <==="
