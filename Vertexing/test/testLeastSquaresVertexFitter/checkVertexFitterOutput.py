import math
import sys

from podio import root_io


input_file, output_file = sys.argv[1:3]
input_frame = next(iter(root_io.Reader(input_file).get("events")))
output_frame = next(iter(root_io.Reader(output_file).get("events")))

input_vertices = input_frame.get("VertexCandidates")
fitted_vertices = output_frame.get("FittedVertices")

assert len(input_vertices) > 0, "the finder produced no input vertex candidates"
assert len(fitted_vertices) == len(input_vertices), (
    f"expected {len(input_vertices)} fitted vertices, got {len(fitted_vertices)}"
)

for candidate, fitted in zip(input_vertices, fitted_vertices):
    assert len(fitted.getTracks()) == len(candidate.getTracks())
    assert len(fitted.getTracks()) >= 2
    assert math.isfinite(fitted.getChi2())
    assert fitted.getNdf() > 0
    assert all(math.isfinite(coordinate) for coordinate in fitted.getPosition())
    if candidate.isPrimary():
        assert fitted.isPrimary(), "the fitter did not preserve the primary-vertex label"
