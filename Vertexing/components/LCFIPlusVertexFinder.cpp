#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "LinearizedHelixVertexFitter.h"
#include "edm4hep/TrackCollection.h"
#include "extension/VertexCollection.h"

#include <TMatrixDSym.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
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
 *  The output extension::VertexCollection contains the fitted primary vertex
 *  followed by all accepted secondary vertices. Each output vertex stores its
 *  associated tracks, fitted position and covariance, chi-squared, number of
 *  degrees of freedom and algorithm identifier. Primary and secondary
 *  vertices are labelled through the corresponding extension::Vertex flags.
 *
 *  @author Andrea De Vita (adapted from FCCAnalyses VertexFinderLCFIPlus and
 *  F. Bedeschi's Delphes VertexFit)
 */
struct LCFIPlusVertexFinder final
    : k4FWCore::Transformer<extension::VertexCollection(const edm4hep::TrackCollection&)> {

  LCFIPlusVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : Transformer(name, serviceLocator, {KeyValues("InputFittedTracks", {"InputFittedTracks"})},
                    {KeyValues("OutputVerticesCandidates", {"OutputVerticesCandidates"})}) {}

  StatusCode initialize() override {
    const auto status = Transformer::initialize();
    if (!status.isSuccess())
      return status;
    if (m_beamSpotPosition.value().size() != 3U || m_beamSpotSize.value().size() != 3U) {
      error() << "BeamSpotPosition and BeamSpotSize must each contain three values" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_primaryTrackChi2Cut <= 0. || m_chi2Cut <= 0. || m_invariantMassCut <= 0. || m_addedTrackChi2Cut <= 0. ||
        m_magneticFieldZ == 0.) {
      error() << "All selection cuts must be positive and MagneticFieldZ must be non-zero" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  extension::VertexCollection operator()(const edm4hep::TrackCollection& inputTracks) const override {
    extension::VertexCollection output;
    const auto tracks = makeTrackEntries(inputTracks);
    if (tracks.size() < 2U)
      return output;

    auto primaryIndices = allIndices(tracks.size());
    std::vector<std::size_t> nonPrimaryIndices;
    auto primaryFit = selectPrimary(tracks, primaryIndices, nonPrimaryIndices);

    if (primaryFit.valid) {
      appendVertex(output, tracks, primaryIndices, primaryFit, true);
    } else {
      warning() << "No valid primary vertex could be reconstructed" << endmsg;
      nonPrimaryIndices = allIndices(tracks.size());
    }

    auto remaining = gather(tracks, nonPrimaryIndices);
    if (m_rejectV0s)
      remaining = removeV0Tracks(remaining, primaryFit.position, true);

    while (remaining.size() > 1U) {
      const auto seed = bestSeed(remaining, primaryFit.position);
      if (!seed.has_value())
        break;

      auto selected = seed->indices;
      while (true) {
        const auto enlarged = addBestTrack(remaining, selected, primaryFit.position);
        if (enlarged.size() == selected.size())
          break;
        selected = enlarged;
      }

      const auto finalFit = fit(remaining, selected, false);
      if (finalFit.valid)
        appendVertex(output, remaining, selected, finalFit, false);
      remaining = removeSelected(remaining, selected);
    }
    return output;
  }

private:
  static constexpr double pionMass = 0.13957039;
  static constexpr double protonMass = 0.93827208;
  static constexpr double electronMass = 0.00051099;
  static constexpr double ptConversion = 0.000299792458; // pT[GeV] = factor B[T] R[mm]

  struct TrackEntry {
    edm4hep::Track track;
    edm4hep::TrackState state;
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

  std::optional<edm4hep::TrackState> selectState(const edm4hep::Track& track) const {
    std::optional<edm4hep::TrackState> first;
    for (const auto& state : track.getTrackStates()) {
      if (!first)
        first = state;
      if (state.location == m_trackStateLocation)
        return state;
    }
    return m_fallbackToFirstTrackState ? first : std::nullopt;
  }

  std::vector<TrackEntry> makeTrackEntries(const edm4hep::TrackCollection& tracks) const {
    std::vector<TrackEntry> result;
    result.reserve(tracks.size());
    for (const auto& track : tracks) {
      const auto state = selectState(track);
      if (state)
        result.push_back({track, *state});
      else
        warning() << "Skipping track without the requested TrackState" << endmsg;
    }
    return result;
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

  FitResult fit(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices, bool primary) const {
    FitResult result;
    LinearizedHelixVertexFitter fitter;
    fitter.setMaxIterations(m_maxIterations);
    fitter.setConvergenceThreshold(m_convergenceThreshold);
    fitter.setSeedStartRadius(primary ? -1. : m_secondarySeedStartRadius.value());
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
    if (!fitter.fit())
      return result;

    result.valid = std::isfinite(fitter.chiSquared());
    result.position = fitter.vertexPosition();
    result.covariance = fitter.vertexCovariance();
    result.chi2 = fitter.chiSquared();
    result.ndf = fitter.numberOfDegreesOfFreedom();
    result.trackChi2.reserve(indices.size());
    for (std::size_t index = 0; index < indices.size(); ++index)
      result.trackChi2.push_back(fitter.trackChiSquared(index));
    return result;
  }

  FitResult selectPrimary(const std::vector<TrackEntry>& tracks, std::vector<std::size_t>& selected,
                          std::vector<std::size_t>& rejected) const {
    while (selected.size() >= 2U) {
      const auto candidate = fit(tracks, selected, true);
      if (!candidate.valid)
        return {};
      const auto worst = std::max_element(candidate.trackChi2.begin(), candidate.trackChi2.end());
      if (worst == candidate.trackChi2.end() || *worst <= m_primaryTrackChi2Cut)
        return candidate;
      const auto offset = static_cast<std::size_t>(std::distance(candidate.trackChi2.begin(), worst));
      rejected.push_back(selected[offset]);
      selected.erase(selected.begin() + static_cast<std::ptrdiff_t>(offset));
    }
    rejected.insert(rejected.end(), selected.begin(), selected.end());
    selected.clear();
    return {};
  }

  TVector3 momentum(const TrackEntry& track) const {
    const auto transverseMomentum = ptConversion * std::abs(m_magneticFieldZ.value() / track.state.omega);
    return {transverseMomentum * std::cos(track.state.phi), transverseMomentum * std::sin(track.state.phi),
            transverseMomentum * track.state.tanLambda};
  }

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

  bool passesConstraints(const std::vector<TrackEntry>& tracks, const std::vector<std::size_t>& indices,
                         const FitResult& candidate, const TVector3& primaryPosition, bool seed) const {
    if (!candidate.valid || candidate.chi2 >= m_chi2Cut)
      return false;
    const std::vector<double> masses(indices.size(), pionMass);
    const auto mass = invariantMass(tracks, indices, masses);
    if (!std::isfinite(mass) || mass >= m_invariantMassCut || mass >= energySum(tracks, indices))
      return false;
    if (pointingCosine(tracks, indices, candidate.position, primaryPosition) < 0.)
      return false;
    return seed || (!candidate.trackChi2.empty() && candidate.trackChi2.back() < m_addedTrackChi2Cut);
  }

  bool isV0Pair(const std::vector<TrackEntry>& tracks, std::size_t first, std::size_t second,
                const TVector3& primaryPosition, bool tight) const {
    if (tracks[first].state.omega * tracks[second].state.omega > 0.)
      return false;
    const std::vector<std::size_t> indices{first, second};
    const auto candidate = fit(tracks, indices, false);
    if (!candidate.valid)
      return false;
    const auto distance = (candidate.position - primaryPosition).Mag();
    const auto pointing = pointingCosine(tracks, indices, candidate.position, primaryPosition);
    const auto ks = ksWindow(tight);
    const auto lambda = lambdaWindow(tight);
    const auto gamma = gammaWindow(tight);
    return insideWindow(invariantMass(tracks, indices, {pionMass, pionMass}), distance, pointing, ks) ||
           insideWindow(invariantMass(tracks, indices, {pionMass, protonMass}), distance, pointing, lambda) ||
           insideWindow(invariantMass(tracks, indices, {protonMass, pionMass}), distance, pointing, lambda) ||
           insideWindow(invariantMass(tracks, indices, {electronMass, electronMass}), distance, pointing, gamma);
  }

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

  std::optional<Candidate> bestSeed(const std::vector<TrackEntry>& tracks, const TVector3& primaryPosition) const {
    std::optional<Candidate> best;
    auto minimumChi2 = std::numeric_limits<double>::infinity();
    for (std::size_t first = 0; first + 1U < tracks.size(); ++first) {
      for (std::size_t second = first + 1U; second < tracks.size(); ++second) {
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
        }
      }
    }
    return best;
  }

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
      }
    }
    auto result = selected;
    if (best)
      result.push_back(*best);
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

  void appendVertex(extension::VertexCollection& output, const std::vector<TrackEntry>& tracks,
                    const std::vector<std::size_t>& indices, const FitResult& fitResult, bool primary) const {
    auto vertex = output.create();
    vertex.setPosition({static_cast<float>(fitResult.position.X()), static_cast<float>(fitResult.position.Y()),
                        static_cast<float>(fitResult.position.Z())});
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
  }

  Gaudi::Property<int> m_trackStateLocation{this, "TrackStateLocation", edm4hep::TrackState::AtIP};
  Gaudi::Property<bool> m_fallbackToFirstTrackState{this, "FallbackToFirstTrackState", true};
  Gaudi::Property<double> m_primaryTrackChi2Cut{this, "PrimaryTrackChi2Cut", 25.};
  Gaudi::Property<bool> m_useBeamSpotConstraint{this, "UseBeamSpotConstraint", true};
  Gaudi::Property<std::vector<double>> m_beamSpotPosition{this, "BeamSpotPosition", {0., 0., 0.}};
  Gaudi::Property<std::vector<double>> m_beamSpotSize{this, "BeamSpotSize", {0.01, 0.01, 0.1}};
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
