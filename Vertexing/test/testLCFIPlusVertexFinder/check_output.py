import math
from podio import root_io

frame = next(iter(root_io.Reader("output.root").get("events")))
vertices = frame.get("VertexCandidates")

assert len(vertices) == 2, f"expected one primary and one secondary vertex, got {len(vertices)}"
assert sorted(len(vertex.getTracks()) for vertex in vertices) == [3, 3]

primary = [vertex for vertex in vertices if vertex.isPrimary()]
secondary = [vertex for vertex in vertices if vertex.isSecondary()]
assert len(primary) == 1, f"expected one primary vertex, got {len(primary)}"
assert len(secondary) == 1, f"expected one secondary vertex, got {len(secondary)}"

positions = sorted(
    ([vertex.getPosition()[axis] for axis in range(3)] for vertex in vertices),
    key=lambda position: position[2],
)
for actual, expected in zip(positions, ((0.0, 0.0, 0.0), (0.0, 0.0, 10.0))):
    assert all(math.isclose(a, e, abs_tol=0.1) for a, e in zip(actual, expected)), (actual, expected)

for vertex in vertices:
    assert math.isfinite(vertex.getChi2())
    assert vertex.getNdf() > 0
    assert vertex.getAlgorithmType() == 1
