#!/bin/sh
set -eu

python create_primary_vertex_input.py
k4run runPrimaryVertexFinder.py
python checkPrimaryVertexOutput.py
