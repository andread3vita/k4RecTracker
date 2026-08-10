#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "LinearizedHelixVertexFitter.h"
#include "edm4hep/TrackCollection.h"
#include "extension/MutableVertex.h"
#include "extension/VertexCollection.h"

#include <TMatrixDSym.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/** @class LCFIPlusVertexFinder
 *
 *  Gaudi transformer that reconstructs primary and secondary vertices from a
 *  collection of fitted edm4hep tracks following the LCFIPlus vertex-finding
 *  workflow implemented in FCCAnalyses.
 *
 *  The primary vertex is reconstructed first from all usable input tracks.
 *  Tracks incompatible with the primary hypothesis are iteratively removed
 *  according to their individual chi-squared contribution. An optional
 *  Gaussian beam-spot constraint can be included in this primary fit.
 *
 *  The rejected tracks are subsequently treated as non-primary tracks. The
 *  secondary-vertex stage applies tight event-level V0 rejection, selects the
 *  best two-track seed, and greedily adds compatible tracks. Candidate
 *  selection uses the fitted chi-squared, individual track contribution,
 *  invariant mass, track-energy sum and momentum-pointing constraints.
 *
 *  Every temporary and final vertex hypothesis is fitted directly with
 *  LinearizedHelixVertexFitter. Track momenta used by the mass and pointing
 *  constraints are derived from the EDM4hep helix parameters and the
 *  configured longitudinal magnetic field.
 *
 *  OutputVerticesCandidates contains the fitted primary vertex followed by all
 *  accepted secondary vertices. OutputV0Vertices independently contains the
 *  fitted K-short, Lambda and photon-conversion candidates reconstructed by
 *  the FCCAnalyses get_V0s strategy. Each V0 stores abs(PDG) and invariant
 *  mass [GeV] as parameters[0] and parameters[1], respectively. V0 candidates
 *  can optionally also be included in OutputVerticesCandidates.
 *
 *  @author Andrea De Vita (adapted from FCCAnalyses VertexFinderLCFIPlus and
 *  F. Bedeschi's Delphes VertexFit)
 */
struct LCFIPlusVertexFinder final
    : k4FWCore::MultiTransformer<std::tuple<extension::VertexCollection, extension::VertexCollection>(
          const edm4hep::TrackCollection&)> {

  LCFIPlusVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : MultiTransformer(name, serviceLocator, {KeyValues("InputFittedTracks", {"InputFittedTracks"})},
                         {KeyValues("OutputVerticesCandidates", {"OutputVerticesCandidates"}),
                          KeyValues("OutputV0Vertices", {"V0Vertices"})}) {}

  StatusCode initialize() override {
    const auto status = MultiTransformer::initialize();
    if (!status.isSuccess())
      return status;
    if (m_beamSpotPosition.value().size() != 3U || m_beamSpotSize.value().size() != 3U) {
      error() << "BeamSpotPosition and BeamSpotSize must each contain three values" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_primaryTrackChi2Cut <= 0.) {
      error() << "PrimaryTrackChi2Cut must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_maxIterations <= 0 || m_convergenceThreshold < 0.) {
      error() << "MaxIterations must be positive and ConvergenceThreshold "
                 "non-negative"
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_useBeamSpotConstraint) {
      for (const auto size : m_beamSpotSize.value()) {
        if (!std::isfinite(size) || size <= 0.) {
          error() << "BeamSpotSize values must be finite and positive when the "
                     "constraint is enabled"
                  << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
    if (m_findSecondaryVertices && (m_chi2Cut <= 0. || m_invariantMassCut <= 0. || m_addedTrackChi2Cut <= 0.)) {
      error() << "Secondary-vertex cuts must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if ((m_findSecondaryVertices || m_reconstructV0Vertices) && m_magneticFieldZ == 0.) {
      error() << "MagneticFieldZ must be non-zero when secondary or V0 reconstruction is enabled" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_reconstructV0Vertices && m_v0Chi2Cut <= 0.) {
      error() << "V0Chi2Cut must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    const auto customV0WindowCount = static_cast<unsigned int>(!m_v0KsConstraints.value().empty()) +
                                     static_cast<unsigned int>(!m_v0LambdaConstraints.value().empty()) +
                                     static_cast<unsigned int>(!m_v0GammaConstraints.value().empty());
    if (customV0WindowCount != 0U &&
        (customV0WindowCount != 3U || m_v0KsConstraints.value().size() != 4U ||
         m_v0LambdaConstraints.value().size() != 4U || m_v0GammaConstraints.value().size() != 4U)) {
      error() << "V0KsConstraints, V0LambdaConstraints and V0GammaConstraints must either all be empty or each "
                 "contain {massLow, massHigh, minimumDistance, minimumPointingCosine}"
              << endmsg;
      return StatusCode::FAILURE;
    }
    if (customV0WindowCount == 3U) {
      const auto validWindow = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); }) &&
               values[0] < values[1] && values[2] >= 0. && values[3] >= -1. && values[3] <= 1.;
      };
      if (!validWindow(m_v0KsConstraints.value()) || !validWindow(m_v0LambdaConstraints.value()) ||
          !validWindow(m_v0GammaConstraints.value())) {
        error() << "Custom V0 constraints must be finite, have massLow < massHigh, non-negative distance and a "
                   "pointing cosine in [-1, 1]"
                << endmsg;
        return StatusCode::FAILURE;
      }
    }
    debug() << "Configured primary selection: TrackStateLocation = " << m_trackStateLocation.value()
            << ", FallbackToFirstTrackState = " << m_fallbackToFirstTrackState.value()
            << ", PrimaryTrackChi2Cut = " << m_primaryTrackChi2Cut.value()
            << ", UseBeamSpotConstraint = " << m_useBeamSpotConstraint.value()
            << ", FindSecondaryVertices = " << m_findSecondaryVertices.value()
            << ", ReconstructV0Vertices = " << m_reconstructV0Vertices.value()
            << ", V0UseTightConstraints = " << m_v0UseTightConstraints.value()
            << ", custom V0 constraints = " << (customV0WindowCount == 3U) << ", V0Chi2Cut = " << m_v0Chi2Cut.value()
            << endmsg;
    return StatusCode::SUCCESS;
  }

  std::tuple<extension::VertexCollection, extension::VertexCollection>
  operator()(const edm4hep::TrackCollection& inputTracks) const override {
    extension::VertexCollection outputVertices;
    extension::VertexCollection outputV0Vertices;

    debug() << "Starting vertex reconstruction with " << inputTracks.size() << " input tracks" << endmsg;

    // Keep the EDM track together with the state used by all trial fits.
    const auto tracks = makeTrackEntries(inputTracks);
    debug() << "Selected " << tracks.size() << " tracks with a usable TrackState" << endmsg;
    if (tracks.size() < 2) {
      debug() << "No vertex reconstruction attempted: at least two usable "
                 "tracks are required"
              << endmsg;
      return {std::move(outputVertices), std::move(outputV0Vertices)};
    }

    // Start inclusively: the primary fit itself decides which tracks do not
    // belong to the interaction point.
    auto primaryIndices = allIndices(tracks.size());
    std::vector<std::size_t> nonPrimaryIndices;
    auto primaryFit = selectPrimary(tracks, primaryIndices, nonPrimaryIndices);

    if (primaryFit.valid) {
      appendVertex(outputVertices, tracks, primaryIndices, primaryFit, true);
    } else {
      warning() << "No valid primary vertex could be reconstructed" << endmsg;
      return {std::move(outputVertices), std::move(outputV0Vertices)};
    }

    // FCCAnalyses reconstructs V0s independently of their rejection from the
    // secondary-vertex input. Both operations therefore see the full set of
    // non-primary tracks.
    auto nonPrimaryTracks = gather(tracks, nonPrimaryIndices);
    if (m_reconstructV0Vertices) {
      reconstructV0s(outputV0Vertices, m_includeV0InVertexCandidates ? &outputVertices : nullptr, nonPrimaryTracks,
                     primaryFit.position, m_v0UseTightConstraints.value());
      debug() << "V0 reconstruction produced " << outputV0Vertices.size() << " fitted candidates" << endmsg;
    }

    if (!m_findSecondaryVertices) {
      debug() << "Secondary-vertex finding and fitting disabled; returning the primary vertex and "
              << outputV0Vertices.size() << " independently reconstructed V0 vertices"
              << (m_includeV0InVertexCandidates ? " (also included in OutputVerticesCandidates)" : "") << endmsg;
      return {std::move(outputVertices), std::move(outputV0Vertices)};
    }

    // Secondary finding only sees tracks rejected by the primary hypothesis.
    auto remaining = std::move(nonPrimaryTracks);
    debug() << "Starting secondary-vertex reconstruction with " << remaining.size()
            << " tracks rejected by the primary fit" << endmsg;
    if (m_rejectV0s) {
      const auto beforeV0Rejection = remaining.size();
      remaining = removeV0Tracks(remaining, primaryFit.position, true);
      debug() << "V0 rejection removed " << beforeV0Rejection - remaining.size() << " tracks; " << remaining.size()
              << " remain for secondary vertices" << endmsg;
    }

    // Build one secondary vertex at a time, remove its tracks, and repeat.
    std::size_t secondaryIndex = 0;
    while (remaining.size() > 1) {
      debug() << "Searching for secondary vertex " << secondaryIndex << " using " << remaining.size()
              << " available tracks" << endmsg;
      const auto seed = bestSeed(remaining, primaryFit.position);
      if (!seed.has_value()) {
        debug() << "No valid two-track seed found; secondary reconstruction is "
                   "complete"
                << endmsg;
        break;
      }

      auto selected = seed->indices;
      debug() << "Selected secondary seed " << formatIndices(remaining, selected) << ", chi2 = " << seed->fit.chi2
              << ", ndf = " << seed->fit.ndf << endmsg;
      // addBestTrack returns the same list when no compatible track remains.
      while (true) {
        const auto enlarged = addBestTrack(remaining, selected, primaryFit.position);
        if (enlarged.size() == selected.size())
          break;
        selected = enlarged;
      }

      const auto finalFit = fit(remaining, selected, false);
      if (finalFit.valid) {
        appendVertex(outputVertices, remaining, selected, finalFit, false);
        ++secondaryIndex;
      } else {
        warning() << "Final fit failed for secondary candidate " << formatIndices(remaining, selected) << endmsg;
      }
      remaining = removeSelected(remaining, selected);
    }
    debug() << "Vertex reconstruction produced " << outputVertices.size() << " main-output vertices (1 primary, "
            << secondaryIndex << " ordinary secondary, and "
            << (m_includeV0InVertexCandidates ? outputV0Vertices.size() : 0U) << " included V0) plus "
            << outputV0Vertices.size() << " vertices in the dedicated V0 output" << endmsg;
    return {std::move(outputVertices), std::move(outputV0Vertices)};
  }

private:
  static constexpr double pionMass = 0.13957039;
  static constexpr double protonMass = 0.93827208;
  static constexpr double electronMass = 0.00051099;
  static constexpr double ptConversion = 0.000299792458; // pT[GeV] = factor B[T] R[mm]

  struct TrackEntry {
    edm4hep::Track track;
    edm4hep::TrackState state;
    std::size_t inputIndex{};
    std::size_t stateIndex{};
  };

  struct SelectedTrackState {
    edm4hep::TrackState state;
    std::size_t stateIndex{};
    bool usedFallback{false};
  };

  struct FitResult {
    bool valid{false};
    TVector3 position;
    TMatrixDSym covariance{3};
    double chi2{std::numeric_limits<double>::infinity()};
    int ndf{};
    std::vector<double> trackChi2;
  };

  struct Candidate {
    FitResult fit;
    std::vector<std::size_t> indices;
  };

  struct V0Window {
    double massLow{};
    double massHigh{};
    double minimumDistance{};
    double minimumPointingCosine{};
  };

  struct V0Candidate {
    FitResult fit;
    int absolutePdg{};
    double invariantMass{};
    double distance{};
    double pointing{};
  };

  static V0Window ksWindow(bool tight) {
    return tight ? V0Window{0.493, 0.503, 0.5, 0.999} : V0Window{0.488, 0.508, 0.3, 0.999};
  }

  static V0Window lambdaWindow(bool tight) {
    return tight ? V0Window{1.111, 1.121, 0.5, 0.99995} : V0Window{1.106, 1.126, 0.3, 0.999};
  }

  static V0Window gammaWindow(bool tight) {
    return tight ? V0Window{0., 0.005, 9., 0.99995} : V0Window{0., 0.010, 9., 0.999};
  }

  static bool insideWindow(double mass, double distance, double pointing, const V0Window& window) {
    return mass > window.massLow && mass < window.massHigh && distance > window.minimumDistance &&
           pointing > window.minimumPointingCosine;
  }

  static V0Window makeWindow(const std::vector<double>& values) { return {values[0], values[1], values[2], values[3]}; }

  bool useCustomV0Constraints() const { return !m_v0KsConstraints.value().empty(); }

  std::optional<SelectedTrackState> selectState(const edm4hep::Track& track) const {
    std::optional<SelectedTrackState> first;
    std::size_t stateIndex = 0;
    for (const auto& state : track.getTrackStates()) {
      if (!first)
        first = SelectedTrackState{state, stateIndex, true};
      if (state.location == m_trackStateLocation)
        return SelectedTrackState{state, stateIndex, false};
      ++stateIndex;
    }
    return m_fallbackToFirstTrackState ? first : std::nullopt;
  }

  std::vector<TrackEntry> makeTrackEntries(const edm4hep::TrackCollection& tracks) const {
    std::vector<TrackEntry> result;
    result.reserve(tracks.size());
    std::size_t inputIndex = 0;
    for (const auto& track : tracks) {
      const auto selected = selectState(track);
      if (selected) {
        result.push_back({track, selected->state, inputIndex, selected->stateIndex});
        debug() << "Selected Track[" << inputIndex << "] TrackState[" << selected->stateIndex << "] at location "
                << selected->state.location << (selected->usedFallback ? " (fallback)" : "") << endmsg;
      } else {
        debug() << "Discarding Track[" << inputIndex << "]: no TrackState at requested location "
                << m_trackStateLocation.value() << endmsg;
      }
      ++inputIndex;
    }
    return result;
  }

  static std::string formatIndices(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices) {
    std::ostringstream stream;
    stream << "Track[";
    for (std::size_t offset = 0; offset < indices.size(); ++offset) {
      if (offset != 0)
        stream << ", ";
      stream << tracks[indices[offset]].inputIndex;
    }
    stream << "]";
    return stream.str();
  }

  static std::vector<std::size_t> allIndices(std::size_t size) {
    std::vector<std::size_t> result(size);
    for (std::size_t index = 0; index < size; ++index)
      result[index] = index;
    return result;
  }

  static std::vector<TrackEntry> gather(const std::vector<TrackEntry>& tracks,
                                        const std::vector<std::size_t>& indices) {
    std::vector<TrackEntry> result;
    result.reserve(indices.size());
    for (const auto index : indices)
      result.push_back(tracks[index]);
    return result;
  }

  /**
   * @brief Fit one track combination to a common vertex.
   *
   * The helper configures a fresh LinearizedHelixVertexFitter for every trial.
   * Only primary hypotheses may receive the beam-spot constraint; secondary
   * hypotheses use the displaced-vertex seed-radius configuration.
   *
   * @param tracks Full track list from which the combination is selected.
   * @param indices Indices of the tracks participating in this fit.
   * @param primary Whether the combination is a primary-vertex hypothesis.
   * @return Fit position, covariance, total chi2, NDF and per-track chi2. The
   *         returned result is invalid when the numerical fit fails.
   */
  FitResult fit(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices, bool primary) const {
    FitResult result;
    debug() << "Fitting " << (primary ? "primary" : "secondary") << " candidate " << formatIndices(tracks, indices)
            << endmsg;
    LinearizedHelixVertexFitter fitter;
    fitter.setMaxIterations(m_maxIterations);
    fitter.setConvergenceThreshold(m_convergenceThreshold);
    fitter.setSeedStartRadius(primary ? -1. : m_secondarySeedStartRadius.value());

    // The beam spot is a prior for the primary vertex only. Displaced vertices
    // must remain free to move away from the interaction point.
    if (primary && m_useBeamSpotConstraint) {
      TMatrixDSym covariance(3);
      covariance.Zero();
      for (int axis = 0; axis < 3; ++axis)
        covariance(axis, axis) = m_beamSpotSize.value()[axis] * m_beamSpotSize.value()[axis];
      fitter.setBeamSpotConstraint(
          TVector3(m_beamSpotPosition.value()[0], m_beamSpotPosition.value()[1], m_beamSpotPosition.value()[2]),
          covariance);
    }
    for (const auto index : indices)
      fitter.addTrack(tracks[index].state);
    if (!fitter.fit() || !std::isfinite(fitter.chiSquared())) {
      debug() << "Fit failed for " << formatIndices(tracks, indices) << endmsg;
      return result;
    }

    result.position = fitter.vertexPosition();
    result.covariance = fitter.vertexCovariance();
    result.chi2 = fitter.chiSquared();
    result.ndf = fitter.numberOfDegreesOfFreedom();
    result.trackChi2.reserve(indices.size());
    for (std::size_t index = 0; index < indices.size(); ++index) {
      const auto contribution = fitter.trackChiSquared(index);
      if (!std::isfinite(contribution)) {
        debug() << "Fit produced a non-finite chi2 contribution for Track[" << tracks[indices[index]].inputIndex << "]"
                << endmsg;
        return {};
      }
      result.trackChi2.push_back(contribution);
      debug() << "Track[" << tracks[indices[index]].inputIndex << "] chi2 contribution = " << contribution << endmsg;
    }
    result.valid = true;
    debug() << "Fit succeeded at (" << result.position.X() << ", " << result.position.Y() << ", " << result.position.Z()
            << "), chi2 = " << result.chi2 << ", ndf = " << result.ndf << endmsg;
    return result;
  }

  /**
   * @brief Select the tracks compatible with the primary vertex.
   *
   * All tracks initially participate in the primary fit. If the largest
   * individual chi2 exceeds PrimaryTrackChi2Cut, that track is moved to the
   * rejected list and the reduced combination is fitted again.
   *
   * @param tracks All usable event tracks.
   * @param selected In/out indices of tracks retained in the primary vertex.
   * @param rejected Output indices passed to secondary-vertex reconstruction.
   * @return The accepted primary fit, or an invalid result if fewer than two
   *         compatible tracks remain or a fit fails.
   */
  FitResult selectPrimary(const std::vector<TrackEntry>& tracks, std::vector<std::size_t>& selected,
                          std::vector<std::size_t>& rejected) const {
    // Refit after every removal because all per-track chi2 contributions
    // change when the common vertex moves.
    while (selected.size() >= 2U) {
      debug() << "Attempting primary fit with " << formatIndices(tracks, selected) << endmsg;
      const auto candidate = fit(tracks, selected, true);
      if (!candidate.valid) {
        warning() << "Primary fit failed for " << formatIndices(tracks, selected) << endmsg;
        return {};
      }
      const auto worst = std::max_element(candidate.trackChi2.begin(), candidate.trackChi2.end());
      if (worst == candidate.trackChi2.end() || *worst < m_primaryTrackChi2Cut) {
        debug() << "Accepted primary candidate " << formatIndices(tracks, selected)
                << "; largest track chi2 = " << (worst == candidate.trackChi2.end() ? 0. : *worst) << endmsg;
        return candidate;
      }
      const auto offset = static_cast<std::size_t>(std::distance(candidate.trackChi2.begin(), worst));
      debug() << "Rejecting Track[" << tracks[selected[offset]].inputIndex
              << "] from the primary candidate: chi2 = " << *worst
              << " is not below PrimaryTrackChi2Cut = " << m_primaryTrackChi2Cut.value() << endmsg;
      rejected.push_back(selected[offset]);
      selected.erase(selected.begin() + static_cast<std::ptrdiff_t>(offset));
    }
    rejected.insert(rejected.end(), selected.begin(), selected.end());
    selected.clear();
    debug() << "Primary selection stopped because fewer than two tracks survived" << endmsg;
    return {};
  }

  /**
   * @brief Recover a track momentum from its EDM4hep helix parameters.
   *
   * @param track Track and selected TrackState.
   * @return Momentum vector in GeV for the configured constant Bz field.
   */
  TVector3 momentum(const TrackEntry& track) const {
    // EDM4hep stores curvature in 1/mm, hence the mm-to-GeV conversion in
    // ptConversion. The input phi and tanLambda then set the 3D direction.
    const auto transverseMomentum = ptConversion * std::abs(m_magneticFieldZ.value() / track.state.omega);
    return {transverseMomentum * std::cos(track.state.phi), transverseMomentum * std::sin(track.state.phi),
            transverseMomentum * track.state.tanLambda};
  }

  /**
   * @brief Compute the invariant mass of a selected track combination.
   *
   * @param tracks Full track list.
   * @param indices Selected track indices.
   * @param masses Mass hypothesis in GeV for each selected track, in the same
   * order.
   * @return Invariant mass in GeV.
   */
  double invariantMass(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices,
                       const std::vector<double>& masses) const {
    TVector3 momentumSum;
    double energySum = 0.;
    for (std::size_t entry = 0; entry < indices.size(); ++entry) {
      const auto p = momentum(tracks[indices[entry]]);
      momentumSum += p;
      energySum += std::sqrt(p.Mag2() + masses[entry] * masses[entry]);
    }
    return std::sqrt(std::max(0., energySum * energySum - momentumSum.Mag2()));
  }

  double energySum(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices) const {
    double result = 0.;
    for (const auto index : indices) {
      const auto p = momentum(tracks[index]);
      result += std::sqrt(p.Mag2() + pionMass * pionMass);
    }
    return result;
  }

  /**
   * @brief Measure whether a candidate momentum points away from the primary
   * vertex.
   *
   * @param tracks Full track list.
   * @param indices Tracks forming the candidate.
   * @param position Fitted candidate position.
   * @param primaryPosition Fitted primary-vertex position.
   * @return Cosine between the summed momentum and PV-to-candidate
   * displacement, or -1 when either vector has zero length.
   */
  double pointingCosine(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices,
                        const TVector3& position, const TVector3& primaryPosition) const {
    TVector3 momentumSum;
    for (const auto index : indices)
      momentumSum += momentum(tracks[index]);
    const auto displacement = position - primaryPosition;
    if (momentumSum.Mag2() == 0. || displacement.Mag2() == 0.)
      return -1.;
    return momentumSum.Dot(displacement) / (momentumSum.Mag() * displacement.Mag());
  }

  /**
   * @brief Apply the FCCAnalyses secondary-vertex candidate cuts.
   *
   * Every candidate must pass the total chi2, pion-mass, energy and pointing
   * requirements. When growing an existing seed, the last track must also pass
   * the individual AddedTrackChi2Cut.
   *
   * @param tracks Full track list.
   * @param indices Candidate track indices; the newly tested track is last.
   * @param candidate Result of fitting this exact combination.
   * @param primaryPosition Fitted primary-vertex position.
   * @param seed True for a two-track seed, false while adding a track.
   * @return True when all applicable constraints are satisfied.
   */
  bool passesConstraints(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices,
                         const FitResult& candidate, const TVector3& primaryPosition, bool seed) const {
    if (!candidate.valid) {
      debug() << "Rejecting secondary candidate " << formatIndices(tracks, indices) << ": fit failed" << endmsg;
      return false;
    }
    if (candidate.chi2 >= m_chi2Cut) {
      debug() << "Rejecting secondary candidate " << formatIndices(tracks, indices) << ": chi2 = " << candidate.chi2
              << " is not below Chi2Cut = " << m_chi2Cut.value() << endmsg;
      return false;
    }
    const std::vector<double> masses(indices.size(), pionMass);
    const auto mass = invariantMass(tracks, indices, masses);
    const auto energy = energySum(tracks, indices);
    if (!std::isfinite(mass) || mass >= m_invariantMassCut || mass >= energy) {
      debug() << "Rejecting secondary candidate " << formatIndices(tracks, indices) << ": mass = " << mass
              << " GeV, energy sum = " << energy << " GeV" << endmsg;
      return false;
    }
    const auto pointing = pointingCosine(tracks, indices, candidate.position, primaryPosition);
    if (pointing < 0.) {
      debug() << "Rejecting secondary candidate " << formatIndices(tracks, indices)
              << ": pointing cosine = " << pointing << endmsg;
      return false;
    }

    // For an enlarged candidate the newly tested track is deliberately last.
    if (!seed && (candidate.trackChi2.empty() || candidate.trackChi2.back() >= m_addedTrackChi2Cut)) {
      debug() << "Rejecting added Track[" << tracks[indices.back()].inputIndex << "]: chi2 contribution = "
              << (candidate.trackChi2.empty() ? std::numeric_limits<double>::infinity() : candidate.trackChi2.back())
              << " is not below AddedTrackChi2Cut = " << m_addedTrackChi2Cut.value() << endmsg;
      return false;
    }
    debug() << "Accepted secondary candidate " << formatIndices(tracks, indices) << ", chi2 = " << candidate.chi2
            << ", mass = " << mass << " GeV, pointing cosine = " << pointing << endmsg;
    return true;
  }

  /**
   * @brief Test whether an oppositely charged pair is compatible with a V0
   * decay.
   *
   * The fitted pair is checked against the K-short, both Lambda mass
   * assignments, and photon-conversion windows. Each window combines invariant
   * mass, displacement from the primary vertex and momentum pointing.
   *
   * @param tracks Full track list.
   * @param first Index of the first track.
   * @param second Index of the second track.
   * @param primaryPosition Fitted primary-vertex position.
   * @param tight Select the tight event-level or loose seed-level windows.
   * @return True if the pair matches at least one V0 hypothesis.
   */
  std::optional<V0Candidate> getV0Candidate(const std::vector<TrackEntry>& tracks, std::size_t first,
                                            std::size_t second, const TVector3& primaryPosition, bool tight,
                                            bool applyChi2Cut, bool customConstraints = false) const {
    // Equal curvature signs correspond to equal charges and cannot form the
    // neutral two-body candidates considered here.
    if (tracks[first].state.omega * tracks[second].state.omega > 0.)
      return std::nullopt;
    const std::vector<std::size_t> indices{first, second};
    const auto candidate = fit(tracks, indices, false);
    if (!candidate.valid)
      return std::nullopt;
    if (applyChi2Cut && candidate.chi2 >= m_v0Chi2Cut) {
      debug() << "Rejecting V0 candidate " << formatIndices(tracks, indices) << ": chi2 = " << candidate.chi2
              << " is not below V0Chi2Cut = " << m_v0Chi2Cut.value() << endmsg;
      return std::nullopt;
    }

    const auto distance = (candidate.position - primaryPosition).Mag();
    const auto pointing = pointingCosine(tracks, indices, candidate.position, primaryPosition);
    const auto useCustom = customConstraints && useCustomV0Constraints();
    const auto ks = useCustom ? makeWindow(m_v0KsConstraints.value()) : ksWindow(tight);
    const auto lambda = useCustom ? makeWindow(m_v0LambdaConstraints.value()) : lambdaWindow(tight);
    const auto gamma = useCustom ? makeWindow(m_v0GammaConstraints.value()) : gammaWindow(tight);
    const auto ksMass = invariantMass(tracks, indices, {pionMass, pionMass});
    const auto lambdaMassFirstPion = invariantMass(tracks, indices, {pionMass, protonMass});
    const auto lambdaMassFirstProton = invariantMass(tracks, indices, {protonMass, pionMass});
    const auto gammaMass = invariantMass(tracks, indices, {electronMass, electronMass});

    std::optional<V0Candidate> result;
    // Keep the hypothesis priority used by FCCAnalyses get_V0s: K-short,
    // Lambda for each mass ordering, then photon conversion.
    if (insideWindow(ksMass, distance, pointing, ks))
      result = V0Candidate{candidate, 310, ksMass, distance, pointing};
    else if (insideWindow(lambdaMassFirstPion, distance, pointing, lambda))
      result = V0Candidate{candidate, 3122, lambdaMassFirstPion, distance, pointing};
    else if (insideWindow(lambdaMassFirstProton, distance, pointing, lambda))
      result = V0Candidate{candidate, 3122, lambdaMassFirstProton, distance, pointing};
    else if (insideWindow(gammaMass, distance, pointing, gamma))
      result = V0Candidate{candidate, 22, gammaMass, distance, pointing};

    if (result)
      debug() << "Identified V0 pair " << formatIndices(tracks, indices) << " as abs(PDG) = " << result->absolutePdg
              << ", invariant mass = " << result->invariantMass << " GeV, distance = " << distance
              << " mm, pointing cosine = " << pointing << (tight ? " (tight cuts)" : " (loose cuts)") << endmsg;
    return result;
  }

  bool isV0Pair(const std::vector<TrackEntry>& tracks, std::size_t first, std::size_t second,
                const TVector3& primaryPosition, bool tight) const {
    // isV0 in FCCAnalyses performs the fit and topology/mass tests without a
    // fit-quality cut. The chi2 cut is applied only by get_V0s.
    return getV0Candidate(tracks, first, second, primaryPosition, tight, false).has_value();
  }

  /**
   * @brief Reconstruct and persist non-overlapping fitted V0 candidates.
   *
   * This is the event-level FCCAnalyses get_V0s loop. The first accepted
   * opposite-charge partner claims both tracks. The fitted vertex is stored,
   * and abs(PDG) plus invariant mass are preserved in the parameter vector.
   */
  void reconstructV0s(extension::VertexCollection& output, extension::VertexCollection* combinedOutput,
                      const std::vector<TrackEntry>& tracks, const TVector3& primaryPosition, bool tight) const {
    std::vector<bool> assigned(tracks.size(), false);
    for (std::size_t first = 0; first + 1U < tracks.size(); ++first) {
      if (assigned[first])
        continue;
      for (std::size_t second = first + 1U; second < tracks.size(); ++second) {
        if (assigned[second])
          continue;
        const auto candidate = getV0Candidate(tracks, first, second, primaryPosition, tight, true, true);
        if (!candidate)
          continue;

        const std::vector<std::size_t> indices{first, second};
        const auto appendCandidate = [&](extension::VertexCollection& destination) {
          auto vertex = appendVertex(destination, tracks, indices, candidate->fit, false, "V0");
          // extension::Vertex has no dedicated FCCAnalysesV0 payload. Preserve
          // the two parallel fields from FCCAnalysesV0 using a documented order.
          vertex.addToParameters(static_cast<float>(candidate->absolutePdg));
          vertex.addToParameters(static_cast<float>(candidate->invariantMass));
        };
        appendCandidate(output);
        if (combinedOutput)
          appendCandidate(*combinedOutput);
        assigned[first] = true;
        assigned[second] = true;
        break;
      }
    }
  }

  /**
   * @brief Remove tracks assigned to non-overlapping V0 candidates.
   *
   * Once a track is assigned to a V0 pair it is not considered in another pair,
   * matching the ordering used by the FCCAnalyses implementation.
   *
   * @param tracks Non-primary tracks to inspect.
   * @param primaryPosition Fitted primary-vertex position.
   * @param tight Select the tight or loose V0 windows.
   * @return A copy containing only tracks not assigned to a V0 pair.
   */
  std::vector<TrackEntry> removeV0Tracks(const std::vector<TrackEntry>& tracks, const TVector3& primaryPosition,
                                         bool tight) const {
    std::vector<bool> rejected(tracks.size(), false);
    for (std::size_t first = 0; first + 1U < tracks.size(); ++first) {
      if (rejected[first])
        continue;
      for (std::size_t second = first + 1U; second < tracks.size(); ++second) {
        if (!rejected[second] && isV0Pair(tracks, first, second, primaryPosition, tight)) {
          rejected[first] = true;
          rejected[second] = true;
          break;
        }
      }
    }
    std::vector<TrackEntry> result;
    for (std::size_t index = 0; index < tracks.size(); ++index)
      if (!rejected[index])
        result.push_back(tracks[index]);
    return result;
  }

  /**
   * @brief Find the best valid two-track secondary-vertex seed.
   *
   * Every pair is subjected to loose V0 rejection and the common secondary
   * constraints. Among the surviving pairs, the one with the smallest chi2/NDF
   * is selected.
   *
   * @param tracks Available non-primary tracks.
   * @param primaryPosition Fitted primary-vertex position.
   * @return The best fitted seed and its indices, or std::nullopt if none
   * passes.
   */
  std::optional<Candidate> bestSeed(const std::vector<TrackEntry>& tracks, const TVector3& primaryPosition) const {
    std::optional<Candidate> best;
    auto minimumChi2 = std::numeric_limits<double>::infinity();
    for (std::size_t first = 0; first + 1U < tracks.size(); ++first) {
      for (std::size_t second = first + 1U; second < tracks.size(); ++second) {
        // FCCAnalyses uses the looser V0 windows while choosing SV seeds.
        if (isV0Pair(tracks, first, second, primaryPosition, false))
          continue;
        const std::vector<std::size_t> indices{first, second};
        const auto candidate = fit(tracks, indices, false);
        if (!passesConstraints(tracks, indices, candidate, primaryPosition, true))
          continue;
        const auto normalizedChi2 = candidate.chi2 / candidate.ndf;
        if (normalizedChi2 < minimumChi2) {
          minimumChi2 = normalizedChi2;
          best = Candidate{candidate, indices};
          debug() << "New best secondary seed " << formatIndices(tracks, indices)
                  << " with chi2/ndf = " << normalizedChi2 << endmsg;
        }
      }
    }
    return best;
  }

  /**
   * @brief Add the best compatible remaining track to a secondary candidate.
   *
   * Each unused track is appended temporarily, fitted, and checked. The valid
   * extension with the smallest total chi2 wins.
   *
   * @param tracks Available non-primary tracks.
   * @param selected Indices currently assigned to the candidate.
   * @param primaryPosition Fitted primary-vertex position.
   * @return The enlarged index list, or the unchanged list when no track
   * passes.
   */
  std::vector<std::size_t> addBestTrack(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& selected,
                                        const TVector3& primaryPosition) const {
    std::optional<std::size_t> best;
    auto minimumChi2 = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < tracks.size(); ++index) {
      if (std::find(selected.begin(), selected.end(), index) != selected.end())
        continue;
      auto trial = selected;
      trial.push_back(index);
      const auto candidate = fit(tracks, trial, false);
      if (passesConstraints(tracks, trial, candidate, primaryPosition, false) && candidate.chi2 < minimumChi2) {
        minimumChi2 = candidate.chi2;
        best = index;
        debug() << "Track[" << tracks[index].inputIndex << "] is the current best addition to "
                << formatIndices(tracks, selected) << ", candidate chi2 = " << candidate.chi2 << endmsg;
      }
    }
    auto result = selected;
    if (best) {
      result.push_back(*best);
      debug() << "Added Track[" << tracks[*best].inputIndex << "] to secondary candidate" << endmsg;
    } else {
      debug() << "No compatible track can be added to secondary candidate " << formatIndices(tracks, selected)
              << endmsg;
    }
    return result;
  }

  static std::vector<TrackEntry> removeSelected(const std::vector<TrackEntry>& tracks,
                                                const std::vector<std::size_t>& selected) {
    std::vector<TrackEntry> result;
    for (std::size_t index = 0; index < tracks.size(); ++index)
      if (std::find(selected.begin(), selected.end(), index) == selected.end())
        result.push_back(tracks[index]);
    return result;
  }

  /**
   * @brief Store a successful fit in the extension vertex collection.
   *
   * @param output Destination collection.
   * @param tracks Track list used by the fit.
   * @param indices Tracks associated with the new vertex.
   * @param fitResult Successful fit result to persist.
   * @param primary Whether to set the primary rather than secondary flag.
   */
  extension::MutableVertex appendVertex(extension::VertexCollection& output, const std::vector<TrackEntry>& tracks,
                                        const std::vector<std::size_t>& indices, const FitResult& fitResult,
                                        bool primary, const std::string& label = "") const {
    auto vertex = output.create();
    vertex.setPosition({static_cast<float>(fitResult.position.X()), static_cast<float>(fitResult.position.Y()),
                        static_cast<float>(fitResult.position.Z())});
    // EDM4hep packs the symmetric 3x3 covariance as xx, xy, yy, xz, yz, zz.
    const std::array<float, 6> covariance{
        static_cast<float>(fitResult.covariance(0, 0)), static_cast<float>(fitResult.covariance(1, 0)),
        static_cast<float>(fitResult.covariance(1, 1)), static_cast<float>(fitResult.covariance(2, 0)),
        static_cast<float>(fitResult.covariance(2, 1)), static_cast<float>(fitResult.covariance(2, 2))};
    vertex.setCovMatrix(edm4hep::CovMatrix3f(covariance));
    vertex.setChi2(static_cast<float>(fitResult.chi2));
    vertex.setNdf(fitResult.ndf);
    vertex.setAlgorithmType(m_algorithmType);
    if (primary)
      vertex.setPrimary();
    else
      vertex.setSecondary();
    for (const auto index : indices)
      vertex.addToTracks(tracks[index].track);
    const auto vertexLabel = label.empty() ? (primary ? "primary" : "secondary") : label;
    debug() << "Created " << vertexLabel << " vertex from " << formatIndices(tracks, indices) << " at ("
            << fitResult.position.X() << ", " << fitResult.position.Y() << ", " << fitResult.position.Z()
            << "), chi2 = " << fitResult.chi2 << ", ndf = " << fitResult.ndf << endmsg;
    return vertex;
  }

  Gaudi::Property<int> m_trackStateLocation{this, "TrackStateLocation", edm4hep::TrackState::AtIP};
  Gaudi::Property<bool> m_fallbackToFirstTrackState{this, "FallbackToFirstTrackState", false,
                                                    "Use the first TrackState when the requested location is absent; "
                                                    "disabled to match primary-only selection"};
  Gaudi::Property<double> m_primaryTrackChi2Cut{this, "PrimaryTrackChi2Cut", 25.};
  Gaudi::Property<bool> m_useBeamSpotConstraint{this, "UseBeamSpotConstraint", true};
  Gaudi::Property<std::vector<double>> m_beamSpotPosition{this, "BeamSpotPosition", {0., 0., 0.}};
  Gaudi::Property<std::vector<double>> m_beamSpotSize{this, "BeamSpotSize", {0.0045, 0.00002, 0.3}};
  Gaudi::Property<bool> m_findSecondaryVertices{this, "FindSecondaryVertices", true,
                                                "Run secondary-vertex finding and fitting after the primary vertex has "
                                                "been reconstructed"};
  Gaudi::Property<bool> m_reconstructV0Vertices{
      this, "ReconstructV0Vertices", true,
      "Fit K-short, Lambda and photon-conversion candidates independently and write OutputV0Vertices"};
  Gaudi::Property<bool> m_includeV0InVertexCandidates{
      this, "IncludeV0InVertexCandidates", false,
      "Also copy fitted V0 candidates into OutputVerticesCandidates while retaining OutputV0Vertices"};
  Gaudi::Property<bool> m_v0UseTightConstraints{
      this, "V0UseTightConstraints", true,
      "Use the FCCAnalyses tight V0 mass, displacement and pointing windows; false selects the loose windows"};
  Gaudi::Property<double> m_v0Chi2Cut{this, "V0Chi2Cut", 9., "Maximum chi2 for a reconstructed V0 candidate"};
  Gaudi::Property<std::vector<double>> m_v0KsConstraints{
      this,
      "V0KsConstraints",
      {},
      "Optional custom {massLow, massHigh, minimumDistance, minimumPointingCosine} for K-short reconstruction"};
  Gaudi::Property<std::vector<double>> m_v0LambdaConstraints{
      this,
      "V0LambdaConstraints",
      {},
      "Optional custom {massLow, massHigh, minimumDistance, minimumPointingCosine} for Lambda reconstruction"};
  Gaudi::Property<std::vector<double>> m_v0GammaConstraints{
      this,
      "V0GammaConstraints",
      {},
      "Optional custom {massLow, massHigh, minimumDistance, minimumPointingCosine} for photon conversions"};
  Gaudi::Property<bool> m_rejectV0s{this, "RejectV0s", true};
  Gaudi::Property<double> m_chi2Cut{this, "Chi2Cut", 9.};
  Gaudi::Property<double> m_invariantMassCut{this, "InvariantMassCut", 10.};
  Gaudi::Property<double> m_addedTrackChi2Cut{this, "AddedTrackChi2Cut", 5.};
  Gaudi::Property<double> m_magneticFieldZ{this, "MagneticFieldZ", 2.};
  Gaudi::Property<int> m_maxIterations{this, "MaxIterations", 100};
  Gaudi::Property<double> m_convergenceThreshold{this, "ConvergenceThreshold", 1.e-12};
  Gaudi::Property<double> m_secondarySeedStartRadius{this, "SecondarySeedStartRadius", -1.};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 1};
};

DECLARE_COMPONENT(LCFIPlusVertexFinder)
