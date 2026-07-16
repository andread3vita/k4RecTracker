#!/bin/bash

set -e  # stop if anything fails

# Never leave a successful-looking artifact behind after a failed rerun.
rm -f out_sim_edm4hep.root out_vertices.root

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml

wget -q https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

ddsim --steeringFile SteeringFile_IDEA_o1_v03.py \
      --compactFile  $XML_FILE \
      -G --gun.distribution uniform --gun.particle mu- --gun.multiplicity 4 \
      --random.seed 42 \
      --numberOfEvents 1 \
      --outputFile out_sim_edm4hep.root \
      --part.minimalKineticEnergy "0.00*MeV"

k4run runVertexFinder.py \
      --inputFile out_sim_edm4hep.root \
      --outputFile out_vertices.root
