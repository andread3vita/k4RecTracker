import math

from podio import root_io


frame = next(iter(root_io.Reader("two_vertex_output.root").get("events")))
vertices = frame.get("VertexCandidates")

assert len(vertices) == 2, f"expected one primary and one secondary vertex, got {len(vertices)}"
primary = [vertex for vertex in vertices if vertex.isPrimary()]
secondary = [vertex for vertex in vertices if vertex.isSecondary()]
assert len(primary) == 1, f"expected one primary vertex, got {len(primary)}"
assert len(secondary) == 1, f"expected one secondary vertex, got {len(secondary)}"

for vertex, expected_z in ((primary[0], 0.0), (secondary[0], 10.0)):
    position = vertex.getPosition()
    assert len(vertex.getTracks()) == 3
    assert math.isclose(position[0], 0.0, abs_tol=0.1)
    assert math.isclose(position[1], 0.0, abs_tol=0.1)
    assert math.isclose(position[2], expected_z, abs_tol=0.1), position
    assert math.isfinite(vertex.getChi2())
    assert vertex.getNdf() > 0
    assert vertex.getAlgorithmType() == 1
