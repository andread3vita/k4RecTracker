# Vertexing components

The `Vertexing` package provides event-level primary and secondary vertex reconstruction and a reusable least-squares vertex-refitting stage. Its Gaudi components operate on the package's `extension::VertexCollection`, whose vertices can store track relations, fit results, and primary or secondary flags.

Collection names shown below are the configurable Gaudi property names. The names in parentheses are their defaults.

## Fitted-track event display

[`scripts/display_fitted_tracks.py`](scripts/display_fitted_tracks.py) reads an EDM4hep ROOT file and draws the fitted tracks and reconstructed vertices in side-by-side XY and YZ projections. The positional argument is the `EventHeader` event number:

```bash
python scripts/display_fitted_tracks.py 0
```

This reads `/afs/cern.ch/work/a/adevita/public/workDir/vertexingMeeting/out_vertex_lcfiplus.root` and writes `event_display_0.png`. Use `--entry` to select by zero-based frame entry instead, or override the defaults with `--input`, `--tracks`, `--vertices`, and `--output`:

```bash
python scripts/display_fitted_tracks.py 5 --entry \
  --input events.root \
  --output event_5.png
```

## Component overview

| Component | Purpose | Input | Output |
| --- | --- | --- | --- |
| `PrimaryVertexFinderLCFIPlus` | Select prompt tracks and reconstruct one primary vertex | All fitted tracks | One primary vertex with selected-track relations |
| `LCFIPlusVertexFinder` | Find and fit primary and secondary vertices from all tracks in an event | Fitted tracks | Primary and secondary vertices |
| `LeastSquaresVertexFitter` | Independently refit existing vertex candidates | Vertex candidates with track relations | Refitted vertices |

## `PrimaryVertexFinderLCFIPlus`

`PrimaryVertexFinderLCFIPlus` implements the primary-vertex sequence used in the FCCAnalyses vertex example. It starts with every input track that has a state at IP (`location == 1`), fits a common vertex, and iteratively discards the track with the largest individual chi-squared contribution while that contribution is at least the configured cut. Tracks without a state at IP are discarded. The output collection is empty when fewer than two compatible tracks remain; otherwise it contains exactly one vertex, marked primary, whose track relations contain only the selected tracks.

Input and output:

- `InputTracks` (`InputTracks`): all fitted `edm4hep::Track` objects in the event.
- `OutputPrimaryVertex` (`PrimaryVertex`): an `extension::VertexCollection` containing zero or one primary vertex.

Important properties:

| Property | Default | Meaning |
| --- | ---: | --- |
| `TrackChi2Cut` | 25 | Reject the largest contributor while its individual chi-squared is at least this value |
| `UseBeamSpotConstraint` | `true` | Apply a Gaussian beam-spot prior |
| `BeamSpotPosition` | `[0, 0, 0]` mm | Beam-spot centre |
| `BeamSpotSize` | `[0.0045, 0.00002, 0.3]` mm | FCC-ee Z-pole beam-spot widths used by the reference example |
| `MaxIterations` | 100 | Maximum fitter iterations |
| `ConvergenceThreshold` | `1e-12` | Fitter convergence criterion |
| `AlgorithmType` | 1 | Identifier written to the output vertex |

```python
from Configurables import PrimaryVertexFinderLCFIPlus

primary_vertex_finder = PrimaryVertexFinderLCFIPlus(
    "PrimaryVertexFinderLCFIPlus",
    InputTracks=["FittedTracks"],
    OutputPrimaryVertex=["PrimaryVertex"],
)
```

## `LCFIPlusVertexFinder`

`LCFIPlusVertexFinder` is a complete event-level finder. Its primary stage follows the same selection and fitting strategy as `PrimaryVertexFinderLCFIPlus`: only tracks with the requested state at IP are used by default, all selected tracks enter the first fit, and the largest incompatible chi-squared contributor is removed before refitting. The rejected tracks can optionally form the secondary-vertex input. When enabled, the secondary stage follows the FCCAnalyses LCFIPlus workflow: tight event-level V0 rejection, loose V0 rejection for seed pairs, best two-track seed selection, greedy track addition, and candidate selection using fit quality, invariant mass, energy, and momentum pointing. The finder can independently run the FCCAnalyses `get_V0s` reconstruction, including the two-track vertex fit and K-short, Lambda, and photon-conversion classification.

All trial and final combinations are fitted internally with `LinearizedHelixVertexFitter`; a separate fitter component is therefore not required after this finder.

Input:

- `InputFittedTracks` (`InputFittedTracks`): `edm4hep::TrackCollection`. Tracks must contain a usable `edm4hep::TrackState`, including helix parameters and covariance.

Output:

- `OutputVerticesCandidates` (`OutputVerticesCandidates`): `extension::VertexCollection`. The primary vertex is stored first and marked with `isPrimary()`. Accepted secondary vertices are marked with `isSecondary()`. When `IncludeV0InVertexCandidates=true`, fitted V0 vertices are included as secondary vertices in this collection as well. Each vertex contains its associated tracks, fitted position and covariance, chi-squared, NDF, and algorithm type.
- `OutputV0Vertices` (`V0Vertices`): `extension::VertexCollection`. Fitted V0 candidates are marked with `isSecondary()`. Because the extension EDM has no `FCCAnalysesV0` type, `parameters[0]` stores the absolute PDG code (`310`, `3122`, or `22`) and `parameters[1]` stores the selected invariant mass in GeV.

Important properties:

| Property | Default | Meaning |
| --- | ---: | --- |
| `TrackStateLocation` | `edm4hep::TrackState::AtIP` | Preferred track state for fitting |
| `FallbackToFirstTrackState` | `false` | Use the first state if the requested location is absent |
| `PrimaryTrackChi2Cut` | 25 | Maximum individual contribution retained in the primary fit |
| `UseBeamSpotConstraint` | `true` | Apply a Gaussian prior to the primary vertex |
| `BeamSpotPosition` | `[0, 0, 0]` mm | Beam-spot centre |
| `BeamSpotSize` | `[0.0045, 0.00002, 0.3]` mm | FCC-ee Z-pole beam-spot standard deviations |
| `FindSecondaryVertices` | `true` | Also find and fit secondary vertices from tracks rejected by the primary fit |
| `ReconstructV0Vertices` | `true` | Independently fit V0 candidates and write `OutputV0Vertices` |
| `IncludeV0InVertexCandidates` | `false` | Also include fitted V0 candidates in `OutputVerticesCandidates` |
| `V0UseTightConstraints` | `true` | Use the FCCAnalyses tight V0 windows; `false` selects the loose windows |
| `V0Chi2Cut` | 9 | Maximum chi-squared for a reconstructed V0 candidate |
| `V0KsConstraints` | `[]` | Optional custom `[massLow, massHigh, minimumDistance, minimumPointingCosine]` for K-short reconstruction |
| `V0LambdaConstraints` | `[]` | Optional custom window in the same order for Lambda reconstruction |
| `V0GammaConstraints` | `[]` | Optional custom window in the same order for photon-conversion reconstruction |
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
    OutputV0Vertices=["V0Vertices"],
    FindSecondaryVertices=True,
    ReconstructV0Vertices=True,
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
