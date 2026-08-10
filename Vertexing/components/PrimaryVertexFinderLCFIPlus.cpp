#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "LinearizedHelixVertexFitter.h"
#include "edm4hep/TrackCollection.h"
#include "extension/VertexCollection.h"

#include <TMatrixDSym.h>
#include <TVector3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/** @class PrimaryVertexFinderLCFIPlus
 *
 * Gaudi functional implementing the primary-track selection used by the
 * FCCAnalyses vertex example. All usable input tracks are fitted to a common
 * vertex. While the largest individual track chi-squared contribution is at
 * least TrackChi2Cut, that track is discarded and the remaining tracks are
 * refitted. The one accepted primary vertex is returned together with
 * relations to only its selected edm4hep::Track objects.
 *
 * A Gaussian beam-spot constraint is enabled by default. Property units are
 * mm and mm^2, consistently with EDM4hep TrackState and extension::Vertex.
 */
struct PrimaryVertexFinderLCFIPlus final
    : k4FWCore::Transformer<extension::VertexCollection(const edm4hep::TrackCollection&)> {

  PrimaryVertexFinderLCFIPlus(const std::string& name, ISvcLocator* serviceLocator)
      : Transformer(name, serviceLocator, {KeyValues("InputTracks", {"InputTracks"})},
                    {KeyValues("OutputPrimaryVertex", {"PrimaryVertex"})}) {}

  StatusCode initialize() override {
    const auto status = Transformer::initialize();
    if (!status.isSuccess())
      return status;

    if (m_beamSpotPosition.value().size() != 3 || m_beamSpotSize.value().size() != 3) {
      error() << "BeamSpotPosition and BeamSpotSize must each contain three values" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_trackChi2Cut <= 0.) {
      error() << "TrackChi2Cut must be positive" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_maxIterations <= 0 || m_convergenceThreshold < 0.) {
      error() << "MaxIterations must be positive and ConvergenceThreshold non-negative" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_useBeamSpotConstraint) {
      for (const auto size : m_beamSpotSize.value()) {
        if (!std::isfinite(size) || size <= 0.) {
          error() << "BeamSpotSize values must be finite and positive when the constraint is enabled" << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
    return StatusCode::SUCCESS;
  }

  extension::VertexCollection operator()(const edm4hep::TrackCollection& inputTracks) const override {
    extension::VertexCollection output;

    std::vector<edm4hep::Track> tracks;
    std::vector<edm4hep::TrackState> trackStates;
    tracks.reserve(inputTracks.size());
    trackStates.reserve(inputTracks.size());

    std::size_t trackIndex = 0;
    for (const auto& track : inputTracks) {
      bool foundTrackStateAtIP = false;
      std::size_t trackStateIndex = 0;
      for (const auto& trackState : track.getTrackStates()) {
        std::ostringstream trackStateContent;
        trackStateContent << trackState;
        info() << "Track[" << trackIndex << "] TrackState[" << trackStateIndex << "]: " << trackStateContent.str()
               << endmsg;

        if (trackState.location == 1 && !foundTrackStateAtIP) {
          tracks.push_back(track);
          trackStates.push_back(trackState);
          foundTrackStateAtIP = true;
        }
        ++trackStateIndex;
      }
      if (!foundTrackStateAtIP)
        debug() << "Discarding a track without a TrackState at IP (location == 1)" << endmsg;
      ++trackIndex;
    }

    // A geometrical vertex requires at least two tracks. Tracks without a
    // TrackState at IP have already been deliberately discarded.
    while (tracks.size() >= 2) {
      LinearizedHelixVertexFitter fitter;
      fitter.setMaxIterations(m_maxIterations);
      fitter.setConvergenceThreshold(m_convergenceThreshold);
      fitter.setSeedStartRadius(-1.);

      if (m_useBeamSpotConstraint) {
        TMatrixDSym beamSpotCovariance(3);
        beamSpotCovariance.Zero();
        for (int axis = 0; axis < 3; ++axis) {
          const auto sigma = m_beamSpotSize.value()[axis];
          beamSpotCovariance(axis, axis) = sigma * sigma;
        }
        fitter.setBeamSpotConstraint(
            TVector3(m_beamSpotPosition.value()[0], m_beamSpotPosition.value()[1], m_beamSpotPosition.value()[2]),
            beamSpotCovariance);
      }

      for (const auto& trackState : trackStates)
        fitter.addTrack(trackState);

      if (!fitter.fit() || !std::isfinite(fitter.chiSquared())) {
        warning() << "Primary-vertex fit failed for " << tracks.size() << " tracks" << endmsg;
        return output;
      }

      double largestTrackChi2 = -std::numeric_limits<double>::infinity();
      std::size_t worstTrackIndex = 0;
      for (std::size_t index = 0; index < trackStates.size(); ++index) {
        const auto trackChi2 = fitter.trackChiSquared(index);
        if (!std::isfinite(trackChi2)) {
          warning() << "Primary-vertex fit produced a non-finite track chi2" << endmsg;
          return output;
        }
        if (trackChi2 > largestTrackChi2) {
          largestTrackChi2 = trackChi2;
          worstTrackIndex = index;
        }
      }

      if (largestTrackChi2 < m_trackChi2Cut) {
        auto vertex = output.create();
        const auto& position = fitter.vertexPosition();
        vertex.setPosition(
            {static_cast<float>(position.X()), static_cast<float>(position.Y()), static_cast<float>(position.Z())});

        const auto& fittedCovariance = fitter.vertexCovariance();
        const std::array<float, 6> packedCovariance{
            static_cast<float>(fittedCovariance(0, 0)), static_cast<float>(fittedCovariance(1, 0)),
            static_cast<float>(fittedCovariance(1, 1)), static_cast<float>(fittedCovariance(2, 0)),
            static_cast<float>(fittedCovariance(2, 1)), static_cast<float>(fittedCovariance(2, 2))};
        vertex.setCovMatrix(edm4hep::CovMatrix3f(packedCovariance));
        vertex.setChi2(static_cast<float>(fitter.chiSquared()));
        vertex.setNdf(fitter.numberOfDegreesOfFreedom());
        vertex.setAlgorithmType(m_algorithmType);
        vertex.setPrimary();
        for (const auto& track : tracks)
          vertex.addToTracks(track);

        return output;
      }

      // As in FCCAnalyses get_PrimaryTracks, remove the worst contributor and
      // refit because every contribution changes when the vertex moves.
      tracks.erase(tracks.begin() + static_cast<std::ptrdiff_t>(worstTrackIndex));
      trackStates.erase(trackStates.begin() + static_cast<std::ptrdiff_t>(worstTrackIndex));
    }

    debug() << "Fewer than two tracks survived primary-track selection" << endmsg;
    return output;
  }

private:
  Gaudi::Property<double> m_trackChi2Cut{this, "TrackChi2Cut", 25.,
                                         "Discard the worst track while its individual chi2 is at least this value"};
  Gaudi::Property<bool> m_useBeamSpotConstraint{this, "UseBeamSpotConstraint", true,
                                                "Apply a Gaussian beam-spot constraint to each primary fit"};
  Gaudi::Property<std::vector<double>> m_beamSpotPosition{
      this, "BeamSpotPosition", {0., 0., 0.}, "Beam-spot centre (x, y, z) in mm"};
  Gaudi::Property<std::vector<double>> m_beamSpotSize{
      this, "BeamSpotSize", {0.0045, 0.00002, 0.3}, "Beam-spot widths (sigma x, y, z) in mm"};
  Gaudi::Property<int> m_maxIterations{this, "MaxIterations", 100, "Maximum fitter iterations"};
  Gaudi::Property<double> m_convergenceThreshold{this, "ConvergenceThreshold", 1.e-12,
                                                 "Covariance-weighted fitter convergence threshold"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 1, "Identifier stored on the output vertex"};
};

DECLARE_COMPONENT(PrimaryVertexFinderLCFIPlus)
