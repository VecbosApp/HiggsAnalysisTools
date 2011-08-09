#! /bin/sh 
# ./fromSingleJobResToRooDataSetData.sh:

echo "===> STARTING... <==="

echo "MERGING TREES STEP..."
./mergeTreesData.sh
echo "MERGING TREES STEP DONE."

echo "WEIGHTING TREES STEP..."
./weightTreesData.sh
echo "MERGING TREES STEP DONE."

echo "MERGING MULTIPLE DATASETS IN THE SAME FINAL STATE..."
./mergeMultiSamplesData.sh
./mergeMultiSamplesDataSkim.sh
echo "MERGING MULTIPLE DATASETS IN THE SAME FINAL STATE DONE."

echo "MERGING FINAL STATES..."
./mergeFinalStatesData.sh
./mergeFinalStatesSkimData.sh
echo "MERGING FINAL STATES DONE."

#echo "CREATING ROODATASETS FOR THE FIT..."
#chmod u+x createFitDatasetsData.sh
#./createFitDatasetsData.sh
#echo "CREATING ROODATASETS FOR THE FIT DONE."

echo "===> VERY DONE. <==="
