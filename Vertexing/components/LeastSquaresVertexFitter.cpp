// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackCollection.h"
#include "edm4hep/VertexCollection.h"

// Vertexing kernel
#include "LinearizedHelixVertexFitter.h"

// C++
#include <array>
#include <vector>

/** @class LeastSquaresVertexFitter
 *
 *  Gaudi transformer that fits a set of reconstructed tracks to a common
 *  vertex and writes the result as an edm4hep::Vertex.
 *
 *  It is detector independent: it works on edm4hep::Track objects regardless of
 *  which subdetector produced them, because the fit only uses the helix track
 *  parameters stored in the track states. The numerical fit is performed by
 *  LinearizedHelixVertexFitter, which is a native EDM4hep adaptation of Franco Bedeschi's
 *  Delphes vertex fitter (no Delphes data types, no parameter conversion).
 *
 *  In this first version all input tracks are fitted to a single common vertex
 *  (the typical "primary vertex" use case). To reconstruct secondary/displaced
 *  vertices, run the transformer on a pre-selected sub-collection of tracks
 *  (track-to-vertex association/finding is left to upstream algorithms).
 *
 *  @author Mahmoud Althakeel (adapted from F. Bedeschi's Delphes VertexFit)
 */
struct LeastSquaresVertexFitter final
    : k4FWCore::Transformer<edm4hep::VertexCollection(const edm4hep::TrackCollection&)> {

  LeastSquaresVertexFitter(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc, {KeyValues("InputTracks", {"Tracks"})}, {KeyValues("OutputVertices", {"Vertices"})}) {
  }

  StatusCode initialize() override {
    if (m_beamSpotPosition.value().size() != 3) {
      error() << "BeamSpotPosition must have exactly 3 entries (x, y, z) in mm." << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_beamSpotSize.value().size() != 3) {
      error() << "BeamSpotSize must have exactly 3 entries (sigma_x, sigma_y, sigma_z) in mm." << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& inputTracks) const override {
    edm4hep::VertexCollection outputVertices;

    // ---- configure the fitter from the Gaudi properties ----
    LinearizedHelixVertexFitter fitter;
    fitter.setMaxIterations(m_maxIterations);
    fitter.setConvergenceThreshold(m_convergenceThreshold);
    fitter.setSeedStartRadius(m_seedStartRadius);

    if (m_useBeamSpotConstraint) {
      const TVector3 beamSpotPosition(m_beamSpotPosition.value()[0], m_beamSpotPosition.value()[1],
                                      m_beamSpotPosition.value()[2]);
      TMatrixDSym beamSpotCovariance(3);
      beamSpotCovariance.Zero();
      for (int axis = 0; axis < 3; ++axis) {
        const double sigma = m_beamSpotSize.value()[axis];
        beamSpotCovariance(axis, axis) = sigma * sigma;
      }
      fitter.setBeamSpotConstraint(beamSpotPosition, beamSpotCovariance);
    }

    // ---- collect the requested track state from every input track ----
    for (const edm4hep::Track& track : inputTracks) {
      bool foundRequestedState = false;
      for (const edm4hep::TrackState& trackState : track.getTrackStates()) {
        if (trackState.location == m_trackStateLocation) {
          fitter.addTrack(trackState);
          foundRequestedState = true;
          break;
        }
      }
      if (!foundRequestedState)
        debug() << "Track has no track state at the requested location " << m_trackStateLocation.value()
                << ", skipping it." << endmsg;
    }

    const std::size_t numberOfTracks = fitter.numberOfTracks();
    if (numberOfTracks < static_cast<std::size_t>(m_minimumNumberOfTracks)) {
      debug() << "Only " << numberOfTracks << " usable tracks (minimum is " << m_minimumNumberOfTracks.value()
              << "), no vertex produced." << endmsg;
      return outputVertices;
    }

    // ---- run the fit ----
    if (!fitter.fit()) {
      warning() << "Vertex fit failed for an event with " << numberOfTracks << " tracks." << endmsg;
      return outputVertices;
    }

    // ---- fill the output edm4hep::Vertex ----
    auto vertex = outputVertices.create();

    const TVector3& position = fitter.vertexPosition(); // mm
    vertex.setPosition(edm4hep::Vector3f(static_cast<float>(position.X()), static_cast<float>(position.Y()),
                                         static_cast<float>(position.Z())));

    // 3x3 vertex covariance, packed lower-triangular as EDM4hep expects:
    // (xx, xy, yy, xz, yz, zz)
    const TMatrixDSym& covariance = fitter.vertexCovariance(); // mm^2
    const std::array<float, 6> packedCovariance = {
        static_cast<float>(covariance(0, 0)), static_cast<float>(covariance(1, 0)),
        static_cast<float>(covariance(1, 1)), static_cast<float>(covariance(2, 0)),
        static_cast<float>(covariance(2, 1)), static_cast<float>(covariance(2, 2))};
    vertex.setCovMatrix(edm4hep::CovMatrix3f(packedCovariance));

    vertex.setChi2(static_cast<float>(fitter.chiSquared()));
    vertex.setNdf(fitter.numberOfDegreesOfFreedom());
    vertex.setAlgorithmType(m_algorithmType);
    if (m_markAsPrimaryVertex)
      vertex.setPrimary(true);

    debug() << "Fitted vertex at (" << position.X() << ", " << position.Y() << ", " << position.Z() << ") mm from "
            << numberOfTracks << " tracks, chi2/ndf = " << fitter.chiSquared() << "/"
            << fitter.numberOfDegreesOfFreedom() << endmsg;

    return outputVertices;
  }

private:
  // --- track selection ---
  Gaudi::Property<int> m_trackStateLocation{
      this, "TrackStateLocation", edm4hep::TrackState::AtIP,
      "Which track state to use for the fit (EDM4hep TrackState location code; 1 = AtIP)"};
  Gaudi::Property<int> m_minimumNumberOfTracks{this, "MinimumNumberOfTracks", 2,
                                               "Minimum number of usable tracks required to attempt a fit"};

  // --- fit control ---
  Gaudi::Property<int> m_maxIterations{this, "MaxIterations", 100, "Maximum number of relinearisation iterations"};
  Gaudi::Property<double> m_convergenceThreshold{
      this, "ConvergenceThreshold", 1.0e-12,
      "Stop iterating once the squared, covariance-weighted vertex shift falls below this"};
  Gaudi::Property<double> m_seedStartRadius{
      this, "SeedStartRadius", -1.0,
      "Radius (mm) at which to start the seed helix phase; negative = start at the perigee"};

  // --- beam-spot / prior-vertex constraint ---
  Gaudi::Property<bool> m_useBeamSpotConstraint{this, "UseBeamSpotConstraint", false,
                                                "Add a Gaussian beam-spot (prior vertex) constraint to the fit"};
  Gaudi::Property<std::vector<double>> m_beamSpotPosition{
      this, "BeamSpotPosition", {0.0, 0.0, 0.0}, "Beam-spot centre (x, y, z) in mm"};
  Gaudi::Property<std::vector<double>> m_beamSpotSize{
      this, "BeamSpotSize", {0.01, 0.01, 0.1}, "Beam-spot Gaussian widths (sigma_x, sigma_y, sigma_z) in mm"};

  // --- output decoration ---
  Gaudi::Property<bool> m_markAsPrimaryVertex{this, "MarkAsPrimaryVertex", true,
                                              "Set the primary-vertex bit on the produced vertex"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 0, "Value written to edm4hep::Vertex::algorithmType"};
};

DECLARE_COMPONENT(LeastSquaresVertexFitter)
