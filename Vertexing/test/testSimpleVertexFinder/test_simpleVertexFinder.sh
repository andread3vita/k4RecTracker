#!/bin/bash

set -e

finder_input=../testLCFIPlusVertexFinder/out_sim_edm4hep.root
finder_output=out_simple_vertices.root

test -f "${finder_input}"
rm -f "${finder_output}"

k4run runSimpleVertexFinder.py \
      --inputFile "${finder_input}" \
      --outputFile "${finder_output}"
