#!/bin/bash

k4run runTracksToMCParticlesLinker.py \
  --input ../testTrackFinder/out_tracks.root \
  --output tracks_mc_particles_links.root
