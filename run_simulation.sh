#!/bin/bash

CURRPATH=$(pwd)
ORIG_PARAMS=("$@")
set --
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2025-02-12
set -- "${ORIG_PARAMS[@]}"

WORKDIR="/eos/user/a/adevita/test/simulation"
PATH_TO_K4GEO="/eos/user/a/adevita/workDir/k4geo"
K4RECTRACKER_dir="/eos/user/a/adevita/workDir/k4RecTracker"

ORIG_PARAMS=("$@")
set --
cd $K4RECTRACKER_dir
source setup.sh
k4_local_repo
set -- "${ORIG_PARAMS[@]}"
echo ""

# ORIG_PARAMS=("$@")
# set --
# cd $PATH_TO_K4GEO
# k4_local_repo
# set -- "${ORIG_PARAMS[@]}"
# echo ""

cd $WORKDIR
k4run pythia.py -n 10 --Dumper.Filename out.hepmc --Pythia8.PythiaInterface.pythiacard Zcard.cmd

ddsim --compactFile $PATH_TO_K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml \
      --outputFile out_sim_edm4hep.root \
      --inputFiles out.hepmc \
      --numberOfEvents 10 \
      --random.seed 42 \
      --part.minimalKineticEnergy "0.001*MeV"
