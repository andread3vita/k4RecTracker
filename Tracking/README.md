# Tracking components

The `Tracking` package provides Gaudi components for truth-assisted and machine-learning-based track finding, track fitting, generator-level track construction, tracking validation, and drift-chamber particle-identification observables.

Collection names shown below are the configurable Gaudi property names. The names in parentheses are their defaults.

## Component overview

| Component | Purpose | Main input | Main output |
| --- | --- | --- | --- |
| `GGTFTrackFinder` | Detector-agnostic tracking pattern recognition with an ONNX graph model | Planar and wire tracker hits | Reconstructed tracks |
| `PerfectTrackFinder` | MC-truth-based hit grouping for validation | Digitized-to-simulated hit links and MC particles | Truth-assisted tracks |
| `GenfitTrackFitter` | Fit reconstructed tracks with GENFIT | Tracks containing tracker hits | Fitted tracks and filtered hits |
| `TracksFromGenParticles` | Build idealized helix tracks directly from generator particles | MC particles and simulated tracker hits | Tracks and track-to-particle links |
| `PlotTrackHitDistances` | Fill track-to-simulated-hit residual histograms | Simulated hits and track-to-particle links | Gaudi histogram |
| `TrackdNdxDelphesBased` | Smear the expected drift-chamber cluster density | Track-to-particle links and event header | Reconstructed dN/dx quantities |

## `GGTFTrackFinder`

`GGTFTrackFinder` performs tracking pattern recognition directly on digitized silicon and drift-chamber hits. It converts every hit into a seven-component feature vector, evaluates an ONNX Geometric Graph Track Finding model, clusters the learned embedding, and creates one `edm4hep::Track` per non-noise cluster. The produced tracks contain relations to their assigned input hits but are not fitted.

Inputs:

- `InputPlanarHitCollections`: vector of `edm4hep::TrackerHitPlaneCollection`, normally vertex-detector and silicon-wrapper hits.
- `InputWireHitCollections`: vector of `edm4hep::SenseWireHitCollection`, normally drift-chamber hits.

Output:

- `OutputTracksGGTF` (`OutputTracksGGTF`): `edm4hep::TrackCollection` containing the hit clusters reconstructed as tracks.

Important properties:

- `ModelPath`: path to the required ONNX model.
- `Tbeta` (0.6): threshold used to identify cluster cores.
- `Td` (0.3): radius used to associate embedded points with a cluster core.

An example steering file is available in [`test/testTrackFinder/runTestTrackFinder.py`](test/testTrackFinder/runTestTrackFinder.py).

## `PerfectTrackFinder`

`PerfectTrackFinder` groups digitized hits using their Monte Carlo associations. For each stable MC particle, it collects linked planar and wire hits, orders them by simulated-hit time, and creates a track. It is intended for reconstruction validation and for studying later stages without pattern-recognition inefficiency.

Inputs:

- `InputPlanarHitCollections`: vector of `edm4hep::TrackerHitSimTrackerHitLinkCollection` for planar detectors.
- `InputWireHitCollections`: vector of `edm4hep::TrackerHitSimTrackerHitLinkCollection` for wire detectors.
- `InputMCParticles`: `edm4hep::MCParticleCollection`.

Output:

- `OutputPerfectTracks`: `edm4hep::TrackCollection` with truth-associated tracker hits.

These tracks still need a track fitter before algorithms requiring fitted `TrackState` parameters are run.

## `GenfitTrackFitter`

`GenfitTrackFitter` converts EDM4hep tracks and their planar or wire measurements into GENFIT objects, performs the configured Kalman or deterministic-annealing fit, and converts the result back to EDM4hep. Magnetic-field propagation, detector material, multiple scattering, energy loss, left/right drift-chamber ambiguity, and optional calorimeter extrapolation are supported.

Input:

- `InputTracks` (`InputTracks`): `edm4hep::TrackCollection` containing tracker-hit relations.

Outputs:

- `OutputFittedTracks` (`Fitted_tracks`): `edm4hep::TrackCollection` with fitted states and fit-quality information.
- `OutputFittedTracksWithFilteredHits` (`Fitted_tracks_with_filtered_hits`): fitted tracks whose hit relations reflect hit filtering and resolved drift-chamber ambiguities.
- `OutputFittedHits` (`Fitted_hits`): `edm4hep::TrackerHitPlaneCollection` containing the filtered or reconstructed measurement positions.

Important properties include `FitterType`, `ParticleHypothesisList`, `TrackStateLocation`, `UseBrems`, `FilterTrackHits`, `RunSingleEvaluation`, `SkipTrackOrdering`, `BetaInit`, `BetaFinal`, and `BetaSteps`. The component requires `GeoSvc`; calorimeter extrapolation additionally requires suitable calorimeter geometry extensions.

An example is available in [`test/testTrackFitter/runTestTrackFitter.py`](test/testTrackFitter/runTestTrackFitter.py).

## `TracksFromGenParticles`

`TracksFromGenParticles` constructs idealized helix tracks from charged generator particles. It derives helix parameters from particle position, momentum, charge, and the DD4hep magnetic field, then creates states at the interaction point and, when available, at the first and last simulated tracker hits. Optional extrapolation creates a state at the electromagnetic calorimeter.

Inputs:

- `InputGenParticles` (`MCParticles`): `edm4hep::MCParticleCollection`.
- `InputSimTrackerHits` (`SimTrackerHits`): vector of `edm4hep::SimTrackerHitCollection`.

Outputs:

- `OutputTracks` (`TracksFromGenParticles`): generated `edm4hep::TrackCollection`.
- `OutputMCRecoTrackParticleAssociation` (`TracksFromGenParticlesAssociation`): `edm4hep::TrackMCParticleLinkCollection` connecting each output track to its source particle.

Useful properties include `MinimumParticleMomentum`, `TrackerIDs`, `ExtrapolateToECal`, and `KeepOnlyBestExtrapolation`.

## `PlotTrackHitDistances`

`PlotTrackHitDistances` is a validation consumer. For every track-to-particle association, it builds a helix from the track state at the interaction point and fills the three-dimensional closest-approach distance to simulated hits produced by the same MC particle.

Inputs:

- `InputSimTrackerHits` (`DCHCollection`): `edm4hep::SimTrackerHitCollection`.
- `InputTracksFromGenParticlesAssociation` (`TracksFromGenParticlesAssociation`): `edm4hep::TrackMCParticleLinkCollection`.

Output:

- Gaudi histogram `track_hits_distance_closest_approach`; no event-data collection is produced.

The `Bz` property (2 T by default) sets the constant longitudinal magnetic field used by the validation helix.

## `TrackdNdxDelphesBased`

`TrackdNdxDelphesBased` estimates the drift-chamber cluster density using the Delphes parametrization. It obtains particle kinematics from the MC association, calculates the expected number of clusters along the track within geometry boundaries, applies statistical fluctuations and detector fill factor, and stores a reconstructed dN/dx measurement.

Inputs:

- `InputLinkCollection` (`TrackMCParticleLinks`): `edm4hep::TrackMCParticleLinkCollection`.
- `HeaderName` (`EventHeader`): `edm4hep::EventHeaderCollection`, used to seed reproducible random fluctuations.

Output:

- `OutputCollection` (`RecDqdxCollection`): `edm4hep::RecDqdxCollection` associated with the input tracks.

The component requires `GeoSvc` and `UniqueIDGenSvc`. Its gas-mixture, fill-factor, and detector-boundary property names must match the detector description.

## Typical reconstruction sequences

Data-like reconstruction:

```text
digitized tracker hits
    -> GGTFTrackFinder
    -> GenfitTrackFitter
    -> fitted tracks
```

Truth-assisted validation:

```text
digitized-to-simulated hit links + MC particles
    -> PerfectTrackFinder
    -> GenfitTrackFitter
```

Generator-level performance studies:

```text
MC particles + simulated hits
    -> TracksFromGenParticles
    -> PlotTrackHitDistances
```
