import math

from podio import root_io


frame = next(iter(root_io.Reader("primary_vertex_output.root").get("events")))
vertices = frame.get("PrimaryVertex")

assert len(vertices) == 1, f"expected exactly one primary vertex, got {len(vertices)}"
vertex = vertices[0]
position = vertex.getPosition()

assert vertex.isPrimary()
assert not vertex.isSecondary()
assert len(vertex.getTracks()) == 3, "only the selected prompt tracks should be related"
assert all(math.isclose(track.getTrackStates()[0].Z0, 0.0) for track in vertex.getTracks())
assert math.isclose(position[0], 0.0, abs_tol=0.1)
assert math.isclose(position[1], 0.0, abs_tol=0.1)
assert math.isclose(position[2], 0.0, abs_tol=0.1)
assert math.isfinite(vertex.getChi2())
assert vertex.getNdf() > 0
assert vertex.getAlgorithmType() == 1
