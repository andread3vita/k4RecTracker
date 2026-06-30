#!/bin/bash

# clean up previous output files
rm -f out_sim_edm4hep.root
rm -f out_vertexing.root

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml

wget https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

# ddsim --steeringFile SteeringFile_IDEA_o1_v03.py \
#       --compactFile  $XML_FILE \
#       -G --gun.distribution uniform --gun.particle mu- --gun.multiplicity 3 \
#       --random.seed 10 \
#       --numberOfEvents 1 \
#       --outputFile out_sim_edm4hep.root

k4run runVertexFitter.py --inputFile out_sim_edm4hep.root --outputFile out_vertexing.root