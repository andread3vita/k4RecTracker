#!/bin/sh

set -eu

rm -f two_vertex_input.root two_vertex_output.root
python3 create_two_vertex_input.py
k4run runTwoVertexFinder.py
python3 checkTwoVertexFinderOutput.py
