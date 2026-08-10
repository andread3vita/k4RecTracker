import math

from podio import root_io


frame = next(iter(root_io.Reader("two_vertex_output.root").get("events")))
vertices = frame.get("VertexCandidates")
v0_vertices = frame.get("V0Vertices")

assert len(vertices) == 3, f"expected primary, secondary and V0 vertices, got {len(vertices)}"
primary = [vertex for vertex in vertices if vertex.isPrimary()]
secondary = [
    vertex
    for vertex in vertices
    if vertex.isSecondary() and len(vertex.getParameters()) == 0
]
v0_in_candidates = [
    vertex
    for vertex in vertices
    if vertex.isSecondary() and len(vertex.getParameters()) == 2
]
assert len(primary) == 1, f"expected one primary vertex, got {len(primary)}"
assert len(secondary) == 1, f"expected one secondary vertex, got {len(secondary)}"
assert len(v0_in_candidates) == 1, f"expected one V0 in VertexCandidates, got {len(v0_in_candidates)}"

for vertex, expected_z in ((primary[0], 0.0), (secondary[0], 10.0)):
    position = vertex.getPosition()
    assert len(vertex.getTracks()) == 3
    assert math.isclose(position[0], 0.0, abs_tol=0.1)
    assert math.isclose(position[1], 0.0, abs_tol=0.1)
    assert math.isclose(position[2], expected_z, abs_tol=0.1), position
    assert math.isfinite(vertex.getChi2())
    assert vertex.getNdf() > 0
    assert vertex.getAlgorithmType() == 1

assert len(v0_vertices) == 1, f"expected one fitted V0 vertex, got {len(v0_vertices)}"
v0 = v0_vertices[0]
assert v0.isSecondary()
assert len(v0.getTracks()) == 2
assert len(v0.getParameters()) == 2
assert int(v0.getParameters(0)) == 310
assert math.isclose(v0.getParameters(1), 0.497611, abs_tol=1.0e-5)
assert math.isclose(v0.getPosition()[2], 20.0, abs_tol=0.1)
assert math.isfinite(v0.getChi2())
assert v0.getNdf() == 1

merged_v0 = v0_in_candidates[0]
assert int(merged_v0.getParameters(0)) == 310
assert math.isclose(merged_v0.getParameters(1), 0.497611, abs_tol=1.0e-5)
assert math.isclose(merged_v0.getPosition()[2], 20.0, abs_tol=0.1)
