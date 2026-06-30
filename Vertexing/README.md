# Vertexing

Vertex reconstruction for Key4hep, usable with any detector (the fit only needs
the helix track parameters, so it does not depend on which subdetector produced
the tracks).

## What it does

`LeastSquaresVertexFitter` is a Gaudi transformer that takes an `edm4hep::TrackCollection`
and fits the tracks to a common vertex, producing an `edm4hep::VertexCollection`
with the vertex position, covariance, chi², and number of degrees of freedom.

The numerical engine is `LinearizedHelixVertexFitter`, a native EDM4hep adaptation of
Franco Bedeschi's iterated least-squares helix vertex fitter from the Delphes
fast simulation (`TrackCovariance/VertexFit`). Unlike the FCCAnalyses wrapper,
this version works **directly on EDM4hep track parameters** — there is no
intermediate Delphes track object and no runtime parameter/covariance
conversion.

## How the EDM4hep ↔ helix convention is handled

The fit uses the standard EDM4hep perigee parameters
`(d0, phi, omega, z0, tanLambda)` in EDM4hep units (mm, 1/mm). The closed-form
helix equations internally use the signed half-curvature

```
signedHalfCurvature = -0.5 * omega
```

which reproduces the EDM4hep↔Delphes mapping validated in FCCAnalyses. Because
the derivative matrix is built with respect to the EDM4hep parameters directly,
the input track covariance from `edm4hep::TrackState::getCovMatrix()` is used
as-is: the `-1/2` Jacobian factor cancels exactly between the derivative matrix
and the covariance, so **no covariance conversion is needed**.

The fit is purely geometric and does **not** require the magnetic field (the
curvature is already encoded in `omega`).

## Usage

A single vertex is fitted from all input tracks (the typical primary-vertex
case). To reconstruct secondary/displaced vertices, run the transformer on a
pre-selected sub-collection of tracks — track-to-vertex finding/association is
left to upstream algorithms.

See [`test/runVertexFitter.py`](test/runVertexFitter.py) for a runnable example.

## Configuration (Gaudi properties)

| Property | Default | Meaning |
|---|---|---|
| `InputTracks` | `Tracks` | input `edm4hep::TrackCollection` |
| `OutputVertices` | `Vertices` | output `edm4hep::VertexCollection` |
| `TrackStateLocation` | `1` (AtIP) | which track state to use |
| `MinimumNumberOfTracks` | `2` | minimum usable tracks to attempt a fit |
| `MaxIterations` | `100` | maximum relinearisation iterations |
| `ConvergenceThreshold` | `1e-12` | stop when the covariance-weighted squared vertex shift is below this |
| `SeedStartRadius` | `-1.0` | seed phase start radius in mm; negative = perigee |
| `UseBeamSpotConstraint` | `false` | add a Gaussian beam-spot prior |
| `BeamSpotPosition` | `[0,0,0]` | beam-spot centre (mm) |
| `BeamSpotSize` | `[0.01,0.01,0.1]` | beam-spot widths σx, σy, σz (mm) |
| `MarkAsPrimaryVertex` | `true` | set the primary-vertex bit on the output |
| `AlgorithmType` | `0` | value written to `Vertex::algorithmType` |

## Files

- `include/LinearizedHelixVertexFitter.h`, `src/LinearizedHelixVertexFitter.cpp` — framework-free
  fit kernel (depends only on ROOT linear algebra and `edm4hep::TrackState`).
- `components/LeastSquaresVertexFitter.cpp` — the Gaudi transformer.

## Credit

The vertex-fit algorithm is due to Franco Bedeschi (Delphes
`TrackCovariance/VertexFit`). This module re-expresses it natively in EDM4hep.
