#!/bin/bash

MODEL_PATH=$1

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
STEERING_FILE=Tracking/test/testTrackFinder/SteeringFile_IDEA_o1_v03.py
TBETA=0.6
TD=0.3

src_file="$FCCCONFIG/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py"
cp "$src_file" "$STEERING_FILE"
sed -i 's/simulateCalo *= *True/simulateCalo = False/' "$STEERING_FILE"



ddsim --steeringFile $STEERING_FILE \
      --compactFile  $XML_FILE \
      -G --gun.distribution uniform --gun.particle mu- \
      --random.seed 42 \
      --numberOfEvents 1 \
      --outputFile Tracking/test/testTrackFinder/out_sim_edm4hep.root

k4run Tracking/test/testTrackFinder/runTestTrackFinder.py --inputFile Tracking/test/testTrackFinder/out_sim_edm4hep.root --outputFile Tracking/test/testTrackFinder/out_tracks.root --modelPath $MODEL_PATH --tbeta $TBETA --td $TD
