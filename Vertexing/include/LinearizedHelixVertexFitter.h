#ifndef VERTEXING_LINEARIZEDHELIXVERTEXFITTER_H
#define VERTEXING_LINEARIZEDHELIXVERTEXFITTER_H

// ROOT linear algebra (only external dependency of the kernel)
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVector3.h>
#include <TVectorD.h>

// EDM4hep input is converted explicitly to the Delphes TrackCovariance convention.
#include <edm4hep/TrackState.h>

#include <vector>

/**
 *  @class LinearizedHelixVertexFitter
 *
 *  Least-squares vertex fit of several charged-particle helices to a common
 *  space point. The numerical method is the iterated, linearised vertex fit
 *  originally written by Franco Bedeschi for the Delphes fast simulation
 *  (the "TrackCovariance/VertexFit" class). EDM4hep inputs are converted at
 *  the boundary and the fitting kernel uses the original Delphes convention:
 *
 *    - EDM4hep (d0, phi, omega, z0, tanLambda) is converted to Delphes
 *      (D, phi0, C, z0, cot(theta)), with C = -omega/2;
 *    - the covariance is transformed by the same Jacobian, so every covariance
 *      involving omega gets a factor -1/2 and Var(omega) gets a factor 1/4;
 *    - TrackState::referencePoint is added to every point on the local helix;
 *    - the fitted vertex position and covariance are returned in mm and mm^2.
 *
 *  Internal 5-parameter convention (Delphes, in mm):
 *    parameters[0] = D          transverse impact parameter            [mm]
 *    parameters[1] = phi0       azimuth of the momentum at the perigee [rad]
 *    parameters[2] = C          signed half-curvature                  [1/mm]
 *    parameters[3] = z0         longitudinal impact parameter          [mm]
 *    parameters[4] = cot(theta) (= EDM4hep tanLambda)                  [-]
 *
 *  The fit is purely geometric: it does NOT need the magnetic field, because
 *  the curvature is already contained in C.
 */
class LinearizedHelixVertexFitter {
public:
  /// Index of each parameter inside the 5-vector, named for readability.
  enum ParameterIndex {
    kD0 = 0,
    kPhi = 1,
    kHalfCurvature = 2,
    kOmega = kHalfCurvature, // compatibility alias; internal value is C, not omega
    kZ0 = 3,
    kTanLambda = 4,
    kNumberOfParameters = 5
  };

  LinearizedHelixVertexFitter() = default;

  // ---------------------------------------------------------------------------
  //  Configuration (every value the algorithm uses is settable)
  // ---------------------------------------------------------------------------

  /// Maximum number of relinearisation iterations before giving up.
  void setMaxIterations(int maxIterations) { m_maxIterations = maxIterations; }

  /// Convergence threshold: the fit stops once the squared, covariance-weighted
  /// vertex displacement between two iterations drops below this value.
  void setConvergenceThreshold(double threshold) { m_convergenceThreshold = threshold; }

  /// Optional radius (mm) at which the fast seed should start the helix phase.
  /// Useful for displaced vertices that sit far from the beam line. A negative
  /// value (the default) means "start at the perigee".
  void setSeedStartRadius(double radius) { m_seedStartRadius = radius; }

  /// Add a Gaussian beam-spot / prior-vertex constraint at position
  /// "position" (mm) with covariance "covariance" (mm^2). The constraint is
  /// added as an extra measurement to the fit.
  void setBeamSpotConstraint(const TVector3& position, const TMatrixDSym& covariance);

  // ---------------------------------------------------------------------------
  //  Track input
  // ---------------------------------------------------------------------------

  /// Add one track to the fit, reading its parameters and covariance from the
  /// given EDM4hep track state.
  void addTrack(const edm4hep::TrackState& trackState);

  /// Remove all tracks and reset the fit so the object can be reused.
  void clear();

  std::size_t numberOfTracks() const { return m_trackParameters.size(); }

  // ---------------------------------------------------------------------------
  //  Run the fit
  // ---------------------------------------------------------------------------

  /// Perform the vertex fit. Returns true on success, false if there are too
  /// few tracks or the linear algebra fails.
  bool fit();

  // ---------------------------------------------------------------------------
  //  Results (valid only after a successful fit())
  // ---------------------------------------------------------------------------

  const TVector3& vertexPosition() const { return m_vertexPosition; }        ///< [mm]
  const TMatrixDSym& vertexCovariance() const { return m_vertexCovariance; } ///< [mm^2]
  double chiSquared() const { return m_chiSquared; }
  int numberOfDegreesOfFreedom() const { return m_numberOfDegreesOfFreedom; }
  double trackChiSquared(std::size_t trackIndex) const { return m_trackChiSquared.at(trackIndex); }

  // ---------------------------------------------------------------------------
  //  Helix geometry in the Delphes convention (public + static so they can be
  //  unit-tested and reused on their own)
  // ---------------------------------------------------------------------------

  /// Convert an EDM4hep track state to (D, phi0, C, z0, cot(theta)), with
  /// C = -omega/2. Lengths remain in mm.
  static TVectorD parametersFromTrackState(const edm4hep::TrackState& trackState);

  /// Convert the EDM4hep 5x5 covariance to the Delphes parameter convention.
  static TMatrixDSym covarianceFromTrackState(const edm4hep::TrackState& trackState);

  /// 3D point on the helix as a function of the helix phase (the angle swept
  /// from the perigee). Phase 0 is the perigee point.
  static TVector3 helixPointAtPhase(const TVectorD& parameters, double phase);

  /// Matrix of derivatives d(position)/d(parameter): 3 rows (x,y,z) by
  /// 5 columns (D, phi0, C, z0, cot(theta)).
  static TMatrixD positionDerivativesWrtParameters(const TVectorD& parameters, double phase);

  /// Derivative d(position)/d(phase): a 3-vector tangent to the helix.
  static TVector3 positionDerivativeWrtPhase(const TVectorD& parameters, double phase);

  /// Numerically robust, regularised symmetric-matrix inverse (recursive
  /// block inversion with row/column normalisation), ported from Delphes
  /// TrkUtil::RegInv. Used instead of a plain inverse because the intermediate
  /// matrices can be nearly singular.
  static TMatrixDSym regularizedInverse(const TMatrixDSym& inputMatrix);

private:
  /// Build a quick, non-iterative starting vertex and per-track phases.
  /// (Port of Delphes VertexFit::VtxFitNoSteer.)
  void computeInitialSeed();

  /// Small helpers to move between ROOT's TVector3 and a length-3 TVectorD.
  static TVectorD toVectorD(const TVector3& vector);
  static TVector3 toVector3(const TVectorD& vector);
  static double dotProduct(const TVectorD& first, const TVectorD& second);

  // --- configuration -------------------------------------------------------
  int m_maxIterations = 100;
  double m_convergenceThreshold = 1.0e-12;
  double m_seedStartRadius = -1.0;
  bool m_useBeamSpotConstraint = false;
  TVector3 m_beamSpotPosition;
  TMatrixDSym m_beamSpotInverseCovariance{3};

  // --- inputs --------------------------------------------------------------
  std::vector<TVectorD> m_trackParameters;     ///< one 5-vector per track
  std::vector<TMatrixDSym> m_trackCovariances; ///< one 5x5 matrix per track
  std::vector<TVector3> m_trackReferencePoints; ///< exact EDM4hep reference points [mm]

  // --- working / output ----------------------------------------------------
  std::vector<double> m_trackPhases; ///< fitted helix phase per track
  TVector3 m_vertexPosition;
  TMatrixDSym m_vertexCovariance{3};
  std::vector<double> m_trackChiSquared;
  double m_chiSquared = 0.0;
  int m_numberOfDegreesOfFreedom = 0;
  bool m_fitSucceeded = false;
};

#endif // VERTEXING_LINEARIZEDHELIXVERTEXFITTER_H
