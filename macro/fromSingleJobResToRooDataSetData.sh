#! /bin/sh 
# ./fromSingleJobResToRooDataSetData.sh:

echo "===> STARTING... <==="

echo "MERGING TREES STEP..."
chmod u+x ./mergeTreesData.sh
./mergeTreesData.sh eg
./mergeTreesData.sh mu
echo "MERGING TREES STEP DONE."

echo "WEIGHTING TREES STEP..."
chmod u+x ./weightTrees.sh
./weightTreesData.sh
echo "MERGING TREES STEP DONE."

echo "DECOUPLE EG/MU DATASETS FOR EMU..."
chmod u+x ./decouplePDs.sh
./decouplePDs.sh
echo "DECOUPLING EG/MU DATASETS FOR EMU DONE."

echo "MERGING FINAL STATES..."
chmod u+x ./mergeFinalStates.sh
./mergeFinalStates.sh
echo "MERGING FINAL STATES DONE."

echo "CREATING ROODATASETS FOR THE FIT..."
chmod u+x createFitDatasetsData.sh
./createFitDatasetsData.sh
echo "CREATING ROODATASETS FOR THE FIT DONE."

echo "===> VERY DONE. <==="
