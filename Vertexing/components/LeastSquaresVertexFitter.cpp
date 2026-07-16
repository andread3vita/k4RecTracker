// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackCollection.h"
#include "edm4hep/VertexCollection.h"
#include "extension/VertexCollection.h"

// Vertexing kernel
#include "LinearizedHelixVertexFitter.h"

// C++
#include <array>
#include <vector>

/** @class LeastSquaresVertexFitter
 *
 *  Gaudi transformer that refines a collection of candidate vertices by
 *  performing a least-squares vertex fit on the associated tracks.
 *
 *  Each input extension::Vertex is treated as a vertex hypothesis containing
 *  a set of associated edm4hep::Track objects. The algorithm re-fits each
 *  vertex position using a linearised helix-based least-squares method
 *  implemented in LinearizedHelixVertexFitter.
 *
 *  The fit is detector-independent, since it relies only on the helix track
 *  parameters stored in edm4hep::TrackState objects and does not depend on
 *  any detector-specific reconstruction details.
 *
 *  Optionally, a beam-spot constraint can be applied to vertices flagged as
 *  primary, introducing a Gaussian prior on the vertex position.
 *
 *  In this implementation, each input vertex is processed independently:
 *  tracks associated to the vertex are refitted, and a new refined vertex is
 *  produced in the output collection.
 *
 *  This transformer is typically used as a refinement stage after an initial
 *  vertex finding algorithm.
 *
 *  @author Mahmoud Althakeel, Andrea De Vita (adapted from F. Bedeschi’s Delphes VertexFit)
 */
struct LeastSquaresVertexFitter final
    : k4FWCore::Transformer<extension::VertexCollection(const extension::VertexCollection&)> {

  LeastSquaresVertexFitter(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,

                    {KeyValues("InputVerticesCandidates", {"InputVerticesCandidates"})},
                    {KeyValues("OutputFittedVertices", {"OutputFittedVertices"})}) {}

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

  extension::VertexCollection operator()(const extension::VertexCollection& verticesCandidates) const override {

    extension::VertexCollection fittedVertices;
    for (const auto& vertex : verticesCandidates) {

      auto numberOfTracks = vertex.getTracks().size();
      if (numberOfTracks < static_cast<std::size_t>(m_minimumNumberOfTracks)) {
        debug() << "Vertex has only " << numberOfTracks << " tracks (minimum is " << m_minimumNumberOfTracks.value()
                << "), skipping it." << endmsg;
        continue;
      }

      // ---- configure the fitter from the Gaudi properties ----
      LinearizedHelixVertexFitter fitter;
      fitter.setMaxIterations(m_maxIterations);
      fitter.setConvergenceThreshold(m_convergenceThreshold);
      fitter.setSeedStartRadius(m_seedStartRadius);

      if (vertex.isPrimary() && m_useBeamSpotConstraint) {

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

      // ---- run the fit ----
      auto tracks = vertex.getTracks();

      // ---- collect the requested track state from every input track ----
      auto vertexFitted = fittedVertices.create();
      for (const edm4hep::Track& track : tracks) {
        bool foundRequestedState = false;
        for (const edm4hep::TrackState& trackState : track.getTrackStates()) {
          if (trackState.location == m_trackStateLocation) {
            fitter.addTrack(trackState);
            foundRequestedState = true;

            vertexFitted.addToTracks(track);
            break;
          }
        }
        if (!foundRequestedState)
          debug() << "Track has no track state at the requested location " << m_trackStateLocation.value()
                  << ", skipping it." << endmsg;
      }

      // ---- run the fit ----
      if (!fitter.fit()) {
        warning() << "Vertex fit failed for a vertex with " << numberOfTracks << " tracks." << endmsg;
      }

      const TVector3& position = fitter.vertexPosition(); // mm
      vertexFitted.setPosition(edm4hep::Vector3f(static_cast<float>(position.X()), static_cast<float>(position.Y()),
                                                 static_cast<float>(position.Z())));

      // 3x3 vertex covariance, packed lower-triangular as EDM4hep expects:
      // (xx, xy, yy, xz, yz, zz)
      const TMatrixDSym& covariance = fitter.vertexCovariance(); // mm^2
      const std::array<float, 6> packedCovariance = {
          static_cast<float>(covariance(0, 0)), static_cast<float>(covariance(1, 0)),
          static_cast<float>(covariance(1, 1)), static_cast<float>(covariance(2, 0)),
          static_cast<float>(covariance(2, 1)), static_cast<float>(covariance(2, 2))};
      vertexFitted.setCovMatrix(edm4hep::CovMatrix3f(packedCovariance));

      vertexFitted.setChi2(static_cast<float>(fitter.chiSquared()));
      vertexFitted.setNdf(fitter.numberOfDegreesOfFreedom());
      vertexFitted.setAlgorithmType(m_algorithmType);

      if (vertex.isPrimary())
        vertexFitted.setPrimary(true);

      debug() << "Fitted vertex at (" << position.X() << ", " << position.Y() << ", " << position.Z() << ") mm from "
              << numberOfTracks << " tracks, chi2/ndf = " << fitter.chiSquared() << "/"
              << fitter.numberOfDegreesOfFreedom() << endmsg;
    }

    return fittedVertices;
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

  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 0, "Value written to edm4hep::Vertex::algorithmType"};
};

DECLARE_COMPONENT(LeastSquaresVertexFitter)