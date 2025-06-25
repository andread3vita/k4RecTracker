#!/bin/bash
# file: test_GenfitFitter.sh
# author: Andrea De Vita, CERN 2025
# to run: sh + test_GenfitFitter.sh

ddsim   --enableGun \
        --gun.distribution uniform \
        --gun.momentumMin "0.5*GeV" \
        --gun.momentumMax "20*GeV" \
        --gun.particle pi+ \
        --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml \
        --outputFile out_sim_edm4hep.root \
        --numberOfEvents 50 \
        --random.seed 42 \
        --steeringFile  sim_steering.py

# download file for cluster counting technique
ifilename="https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root"
# -nc, --no-clobber: skip downloads that would download to existing files
wget --no-clobber $ifilename

# Check if the input file exists
if [[ ! -f "$(basename $ifilename)" ]]; then
    echo "Error: File '$(basename $ifilename)' not found."
    exit 1
fi

# run digitizer for position smearing and cluster counting calculation
k4run runGenfitFitterPipeline.py
