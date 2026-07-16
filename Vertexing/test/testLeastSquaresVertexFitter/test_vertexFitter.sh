#!/bin/bash

set -e

finder_output=../testLCFIPlusVertexFinder/out_vertices.root
fitter_output=out_fitted_vertices.root

test -f "${finder_output}"
rm -f "${fitter_output}"

k4run runVertexFitter.py \
      --inputFile "${finder_output}" \
      --outputFile "${fitter_output}"

python checkVertexFitterOutput.py "${finder_output}" "${fitter_output}"
