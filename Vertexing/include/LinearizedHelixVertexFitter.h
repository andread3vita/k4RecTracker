#ifndef VERTEXING_LINEARIZEDHELIXVERTEXFITTER_H
#define VERTEXING_LINEARIZEDHELIXVERTEXFITTER_H

// ROOT linear algebra (only external dependency of the kernel)
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVector3.h>
#include <TVectorD.h>

// EDM4hep: the kernel reads track parameters directly from an edm4hep::TrackState,
// so there is no intermediate "Delphes track" object and no parameter conversion step.
#include <edm4hep/TrackState.h>

#include <vector>

/**
 *  @class LinearizedHelixVertexFitter
 *
 *  Least-squares vertex fit of several charged-particle helices to a common
 *  space point. The numerical method is the iterated, linearised vertex fit
 *  originally written by Franco Bedeschi for the Delphes fast simulation
 *  (the "TrackCovariance/VertexFit" class). It has been re-expressed here so
 *  that it works *natively* on EDM4hep track parameters:
 *
 *    - input track parameters come straight from edm4hep::TrackState
 *      (d0, phi, omega, z0, tanLambda), in EDM4hep units (mm and 1/mm);
 *    - the input covariance is the edm4hep::TrackState covariance, used as-is;
 *    - the fitted vertex position and covariance are returned in mm and mm^2.
 *
 *  The only place where the Delphes helix convention enters is inside the
 *  helix-geometry helpers below, where the signed half-curvature used by the
 *  closed-form helix equations is computed from omega as
 *
 *      signedHalfCurvature = -0.5 * omega
 *
 *  (the minus sign reproduces the EDM4hep<->Delphes mapping validated in
 *  FCCAnalyses). Because the derivative matrix is built with respect to the
 *  EDM4hep parameters directly, no covariance conversion is ever needed: the
 *  -1/2 factor cancels exactly between the derivative matrix and the
 *  covariance.
 *
 *  Conventions of the 5-parameter helix (EDM4hep perigee parametrisation):
 *    parameters[0] = d0        transverse impact parameter            [mm]
 *    parameters[1] = phi       azimuth of the momentum at the perigee [rad]
 *    parameters[2] = omega     signed curvature (sign of the charge)  [1/mm]
 *    parameters[3] = z0        longitudinal impact parameter          [mm]
 *    parameters[4] = tanLambda tangent of the dip angle               [-]
 *
 *  The fit is purely geometric: it does NOT need the magnetic field, because
 *  the curvature is already contained in omega.
 */
class LinearizedHelixVertexFitter {
public:
  /// Index of each parameter inside the 5-vector, named for readability.
  enum ParameterIndex { kD0 = 0, kPhi = 1, kOmega = 2, kZ0 = 3, kTanLambda = 4, kNumberOfParameters = 5 };

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
  //  Helix geometry in the EDM4hep convention (public + static so they can be
  //  unit-tested and reused on their own)
  // ---------------------------------------------------------------------------

  /// Copy the five helix parameters out of an EDM4hep track state into a vector.
  static TVectorD parametersFromTrackState(const edm4hep::TrackState& trackState);

  /// Copy the 5x5 parameter covariance out of an EDM4hep track state.
  static TMatrixDSym covarianceFromTrackState(const edm4hep::TrackState& trackState);

  /// 3D point on the helix as a function of the helix phase (the angle swept
  /// from the perigee). Phase 0 is the perigee point.
  static TVector3 helixPointAtPhase(const TVectorD& parameters, double phase);

  /// Matrix of derivatives d(position)/d(parameter): 3 rows (x,y,z) by
  /// 5 columns (d0, phi, omega, z0, tanLambda).
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

  /// Signed half-curvature used by the closed-form helix equations, derived
  /// from the EDM4hep curvature omega.
  static double signedHalfCurvature(double omega) { return -0.5 * omega; }

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
