#!/bin/bash

set -e  # stop if anything fails

k4run runVertexFitter.py \
      --inputFile ../testDeterministicAnnealingVertexFinder/out_vertices.root \
      --outputFile out_vertexing.root