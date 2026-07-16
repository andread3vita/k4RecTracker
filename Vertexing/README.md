# Vertexing components

The `Vertexing` package provides event-level primary and secondary vertex reconstruction and a reusable least-squares vertex-refitting stage. Both Gaudi components operate on the package's `extension::VertexCollection`, whose vertices can store track relations, fit results, and primary or secondary flags.

Collection names shown below are the configurable Gaudi property names. The names in parentheses are their defaults.

## Component overview

| Component | Purpose | Input | Output |
| --- | --- | --- | --- |
| `LCFIPlusVertexFinder` | Find and fit primary and secondary vertices from all tracks in an event | Fitted tracks | Primary and secondary vertices |
| `LeastSquaresVertexFitter` | Independently refit existing vertex candidates | Vertex candidates with track relations | Refitted vertices |

## `LCFIPlusVertexFinder`

`LCFIPlusVertexFinder` is a complete event-level finder. It first fits a primary-vertex hypothesis using all usable tracks and iteratively rejects the track with the largest incompatible chi-squared contribution. The rejected tracks form the secondary-vertex input. The secondary stage follows the FCCAnalyses LCFIPlus workflow: V0 rejection, best two-track seed selection, greedy track addition, and candidate selection using fit quality, invariant mass, energy, and momentum pointing.

All trial and final combinations are fitted internally with `LinearizedHelixVertexFitter`; a separate fitter component is therefore not required after this finder.

Input:

- `InputFittedTracks` (`InputFittedTracks`): `edm4hep::TrackCollection`. Tracks must contain a usable `edm4hep::TrackState`, including helix parameters and covariance.

Output:

- `OutputVerticesCandidates` (`OutputVerticesCandidates`): `extension::VertexCollection`. The primary vertex is stored first and marked with `isPrimary()`. Accepted secondary vertices are marked with `isSecondary()`. Each vertex contains its associated tracks, fitted position and covariance, chi-squared, NDF, and algorithm type.

Important properties:

| Property | Default | Meaning |
| --- | ---: | --- |
| `TrackStateLocation` | `edm4hep::TrackState::AtIP` | Preferred track state for fitting |
| `FallbackToFirstTrackState` | `true` | Use the first state if the requested location is absent |
| `PrimaryTrackChi2Cut` | 25 | Maximum individual contribution retained in the primary fit |
| `UseBeamSpotConstraint` | `true` | Apply a Gaussian prior to the primary vertex |
| `BeamSpotPosition` | `[0, 0, 0]` mm | Beam-spot centre |
| `BeamSpotSize` | `[0.01, 0.01, 0.1]` mm | Beam-spot standard deviations |
| `RejectV0s` | `true` | Remove K-short, Lambda, and photon-conversion candidates before secondary finding |
| `Chi2Cut` | 9 | Maximum total candidate chi-squared |
| `InvariantMassCut` | 10 GeV | Maximum pion-hypothesis candidate mass |
| `AddedTrackChi2Cut` | 5 | Maximum contribution from a newly added track |
| `MagneticFieldZ` | 2 T | Constant longitudinal field used to derive track momenta |
| `MaxIterations` | 100 | Maximum fitter iterations |
| `ConvergenceThreshold` | `1e-12` | Fitter convergence criterion |
| `SecondarySeedStartRadius` | -1 mm | Initial helix radius; negative starts at the perigee |
| `AlgorithmType` | 1 | Identifier written to output vertices |

Example configuration:

```python
from Configurables import LCFIPlusVertexFinder

vertexFinder = LCFIPlusVertexFinder(
    "LCFIPlusVertexFinder",
    InputFittedTracks=["FittedTracks"],
    OutputVerticesCandidates=["VertexCandidates"],
    MagneticFieldZ=2.0,
    RejectV0s=True,
)
```

See [`test/testLCFIPlusVertexFinder/runVertexFinder.py`](test/testLCFIPlusVertexFinder/runVertexFinder.py) for its use in a complete digitization and tracking chain.

## `LeastSquaresVertexFitter`

`LeastSquaresVertexFitter` refits each input vertex independently using the `edm4hep::Track` relations already stored in that vertex. The underlying `LinearizedHelixVertexFitter` performs an iterative, covariance-weighted common-point fit directly from EDM4hep helix parameters. A beam-spot constraint can optionally be applied to candidates marked as primary.

This component is useful when candidates were created by a finder that only associates tracks or when a common fitting configuration should be applied to an existing vertex collection. Using it after `LCFIPlusVertexFinder` is optional because that finder already fits its accepted vertices.

Input:

- `InputVerticesCandidates` (`InputVerticesCandidates`): `extension::VertexCollection`. Each candidate must relate at least `MinimumNumberOfTracks` tracks containing the requested track state.

Output:

- `OutputFittedVertices` (`OutputFittedVertices`): `extension::VertexCollection` containing the refitted position, covariance, chi-squared, NDF, algorithm type, and input track relations. Primary labels are preserved and control whether the beam-spot constraint is applied.

Important properties:

| Property | Default | Meaning |
| --- | ---: | --- |
| `TrackStateLocation` | `edm4hep::TrackState::AtIP` | Track state supplied to the fitter |
| `MinimumNumberOfTracks` | 2 | Minimum number of related tracks required |
| `MaxIterations` | 100 | Maximum relinearization iterations |
| `ConvergenceThreshold` | `1e-12` | Covariance-weighted displacement threshold |
| `SeedStartRadius` | -1 mm | Initial helix radius; negative starts at the perigee |
| `UseBeamSpotConstraint` | `false` | Constrain vertices marked as primary |
| `BeamSpotPosition` | `[0, 0, 0]` mm | Beam-spot centre |
| `BeamSpotSize` | `[0.01, 0.01, 0.1]` mm | Beam-spot standard deviations |
| `AlgorithmType` | 0 | Identifier written to fitted vertices |

Example configuration:

```python
from Configurables import LeastSquaresVertexFitter

vertexFitter = LeastSquaresVertexFitter(
    "LeastSquaresVertexFitter",
    InputVerticesCandidates=["VertexCandidates"],
    OutputFittedVertices=["FittedVertices"],
    UseBeamSpotConstraint=True,
)
```

The CI example in [`test/testLeastSquaresVertexFitter/runVertexFitter.py`](test/testLeastSquaresVertexFitter/runVertexFitter.py) reads the output of the LCFIPlus finder and refits its candidates.

## Data flow

The standard reconstruction path is:

```text
fitted edm4hep tracks
    -> LCFIPlusVertexFinder
    -> extension::VertexCollection
```

For an explicit refitting stage:

```text
vertex candidates with track relations
    -> LeastSquaresVertexFitter
    -> refitted extension::VertexCollection
```
