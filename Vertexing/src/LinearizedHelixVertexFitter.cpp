#include "LinearizedHelixVertexFitter.h"

#include <cmath>
#include <iostream>
#include <limits>

// =============================================================================
//  Conversions between ROOT's TVector3 and a plain length-3 TVectorD
// =============================================================================

TVectorD LinearizedHelixVertexFitter::toVectorD(const TVector3& vector) {
  TVectorD result(3);
  result[0] = vector.X();
  result[1] = vector.Y();
  result[2] = vector.Z();
  return result;
}

TVector3 LinearizedHelixVertexFitter::toVector3(const TVectorD& vector) {
  return TVector3(vector[0], vector[1], vector[2]);
}

double LinearizedHelixVertexFitter::dotProduct(const TVectorD& first, const TVectorD& second) {
  double sum = 0.0;
  for (int i = 0; i < first.GetNrows(); ++i)
    sum += first[i] * second[i];
  return sum;
}

// =============================================================================
//  Reading EDM4hep track states (native: no Delphes conversion)
// =============================================================================

TVectorD LinearizedHelixVertexFitter::parametersFromTrackState(const edm4hep::TrackState& trackState) {
  TVectorD parameters(kNumberOfParameters);
  parameters[kD0] = trackState.D0;
  parameters[kPhi] = trackState.phi;
  parameters[kOmega] = trackState.omega;
  parameters[kZ0] = trackState.Z0;
  parameters[kTanLambda] = trackState.tanLambda;
  return parameters;
}

TMatrixDSym LinearizedHelixVertexFitter::covarianceFromTrackState(const edm4hep::TrackState& trackState) {
  // The EDM4hep track covariance is stored as a lower-triangular, packed array
  // (CovMatrix6f, where the 6th parameter is time). For the spatial 5x5 block
  // the packed index of element (row, col) with row >= col is
  //     row * (row + 1) / 2 + col.
  // We use the parameters as they are stored in EDM4hep (mm, 1/mm); no scaling
  // and no sign flip are applied here.
  TMatrixDSym covariance(kNumberOfParameters);
  for (int row = 0; row < kNumberOfParameters; ++row) {
    for (int col = 0; col <= row; ++col) {
      const int packedIndex = row * (row + 1) / 2 + col;
      const double value = trackState.covMatrix[packedIndex];
      covariance(row, col) = value;
      covariance(col, row) = value; // keep the matrix symmetric
    }
  }
  return covariance;
}

// =============================================================================
//  Helix geometry, expressed directly in EDM4hep parameters
//
//  The closed-form helix uses the signed half-curvature C = -omega/2. We keep
//  the algebra textually close to the original Delphes formulas (with C and the
//  factor 2*C) so it can be checked against the reference, while naming every
//  quantity explicitly.
// =============================================================================

TVector3 LinearizedHelixVertexFitter::helixPointAtPhase(const TVectorD& parameters, double phase) {
  const double impactParameter = parameters[kD0];
  const double phi = parameters[kPhi];
  const double tanLambda = parameters[kTanLambda];

  const double halfCurvature = signedHalfCurvature(parameters[kOmega]); // C = -omega/2
  const double twiceCurvature = 2.0 * halfCurvature;                    // 2C ( = -omega )

  const double x = -impactParameter * std::sin(phi) + (std::sin(phase + phi) - std::sin(phi)) / twiceCurvature;
  const double y = impactParameter * std::cos(phi) - (std::cos(phase + phi) - std::cos(phi)) / twiceCurvature;
  const double z = parameters[kZ0] + tanLambda * phase / twiceCurvature;

  return TVector3(x, y, z);
}

TMatrixD LinearizedHelixVertexFitter::positionDerivativesWrtParameters(const TVectorD& parameters, double phase) {
  const double impactParameter = parameters[kD0];
  const double phi = parameters[kPhi];
  const double tanLambda = parameters[kTanLambda];

  const double halfCurvature = signedHalfCurvature(parameters[kOmega]); // C = -omega/2
  const double twiceCurvature = 2.0 * halfCurvature;                    // 2C

  TMatrixD derivatives(3, kNumberOfParameters);

  // d(position) / d(d0)
  derivatives(0, kD0) = -std::sin(phi);
  derivatives(1, kD0) = std::cos(phi);
  derivatives(2, kD0) = 0.0;

  // d(position) / d(phi)
  derivatives(0, kPhi) = -impactParameter * std::cos(phi) + (std::cos(phase + phi) - std::cos(phi)) / twiceCurvature;
  derivatives(1, kPhi) = -impactParameter * std::sin(phi) + (std::sin(phase + phi) - std::sin(phi)) / twiceCurvature;
  derivatives(2, kPhi) = 0.0;

  // d(position) / d(omega).
  // We first write the derivative with respect to the half-curvature C (the
  // original Delphes form) and then apply the chain rule with dC/d(omega) = -1/2.
  const double dC_dOmega = -0.5;
  const double dx_dC = -(std::sin(phase + phi) - std::sin(phi)) / (2.0 * halfCurvature * halfCurvature);
  const double dy_dC = (std::cos(phase + phi) - std::cos(phi)) / (2.0 * halfCurvature * halfCurvature);
  const double dz_dC = -tanLambda * phase / (2.0 * halfCurvature * halfCurvature);
  derivatives(0, kOmega) = dC_dOmega * dx_dC;
  derivatives(1, kOmega) = dC_dOmega * dy_dC;
  derivatives(2, kOmega) = dC_dOmega * dz_dC;

  // d(position) / d(z0)
  derivatives(0, kZ0) = 0.0;
  derivatives(1, kZ0) = 0.0;
  derivatives(2, kZ0) = 1.0;

  // d(position) / d(tanLambda)
  derivatives(0, kTanLambda) = 0.0;
  derivatives(1, kTanLambda) = 0.0;
  derivatives(2, kTanLambda) = phase / twiceCurvature;

  return derivatives;
}

TVector3 LinearizedHelixVertexFitter::positionDerivativeWrtPhase(const TVectorD& parameters, double phase) {
  const double phi = parameters[kPhi];
  const double tanLambda = parameters[kTanLambda];
  const double twiceCurvature = 2.0 * signedHalfCurvature(parameters[kOmega]);

  return TVector3(std::cos(phase + phi) / twiceCurvature, std::sin(phase + phi) / twiceCurvature,
                  tanLambda / twiceCurvature);
}

// =============================================================================
//  Regularised symmetric-matrix inverse (recursive block inversion).
//  Ported from Delphes TrkUtil::RegInv with descriptive names.
// =============================================================================

TMatrixDSym LinearizedHelixVertexFitter::regularizedInverse(const TMatrixDSym& inputMatrix) {
  TMatrixDSym matrix = inputMatrix; // work on a private copy
  const int size = matrix.GetNrows();

  // --- trivial 1x1 case ---
  if (size == 1) {
    TMatrixDSym inverse(1);
    inverse(0, 0) = (matrix(0, 0) != 0.0) ? 1.0 / matrix(0, 0) : 1.0;
    return inverse;
  }

  // --- normalise rows/columns by 1/sqrt(|diagonal|) so the matrix is better
  //     conditioned before inversion ---
  TMatrixDSym normalisation(size);
  normalisation.Zero();
  for (int i = 0; i < size; ++i) {
    normalisation(i, i) = (matrix(i, i) != 0.0) ? 1.0 / std::sqrt(std::abs(matrix(i, i))) : 1.0;
  }
  TMatrixDSym normalised = matrix.Similarity(normalisation); // N * matrix * N
  TMatrixDSym normalisedInverse(size);

  if (size == 2) {
    // --- closed-form 2x2 inverse ---
    const double determinant = normalised(0, 0) * normalised(1, 1) - normalised(0, 1) * normalised(1, 0);
    if (determinant == 0.0) {
      std::cerr << "LinearizedHelixVertexFitter::regularizedInverse: null 2x2 determinant" << std::endl;
      normalisedInverse.Zero();
    } else {
      normalisedInverse(0, 0) = normalised(1, 1);
      normalisedInverse(0, 1) = -normalised(0, 1);
      normalisedInverse(1, 0) = normalisedInverse(0, 1);
      normalisedInverse(1, 1) = normalised(0, 0);
      normalisedInverse *= 1.0 / determinant;
    }
  } else {
    // --- recursive block inversion for size > 2 ---
    // Split the matrix into a top-left block, a border vector and a corner.
    TMatrixDSym topLeftBlock = normalised.GetSub(0, size - 2, 0, size - 2);
    TVectorD borderVector(size - 1);
    for (int i = 0; i < size - 1; ++i)
      borderVector(i) = normalised(size - 1, i);
    const double corner = normalised(size - 1, size - 1);

    if (std::abs(corner) > 1.0e-15) {
      // Case: non-zero corner.
      TMatrixDSym reduced(size - 1);
      reduced.Rank1Update(borderVector, -1.0 / corner); // -1/corner * border border^T
      reduced += topLeftBlock;
      TMatrixDSym reducedInverse = regularizedInverse(reduced); // recursion
      TMatrixDSub(normalisedInverse, 0, size - 2, 0, size - 2) = reducedInverse;

      TVectorD lastColumn = (-1.0 / corner) * (reducedInverse * borderVector);
      for (int i = 0; i < size - 1; ++i) {
        normalisedInverse(size - 1, i) = lastColumn(i);
        normalisedInverse(i, size - 1) = lastColumn(i);
      }
      double borderDotColumn = 0.0;
      for (int i = 0; i < size - 1; ++i)
        borderDotColumn += borderVector(i) * lastColumn(i);
      normalisedInverse(size - 1, size - 1) = (1.0 - borderDotColumn) / corner;
    } else {
      // Case: zero corner needs a different decomposition.
      TMatrixDSym topLeftInverse = regularizedInverse(topLeftBlock); // recursion
      const double scalar = topLeftInverse.Similarity(borderVector); // border^T * inv * border
      normalisedInverse(size - 1, size - 1) = -1.0 / scalar;

      TVectorD lastColumn = (1.0 / scalar) * (topLeftInverse * borderVector);
      for (int i = 0; i < size - 1; ++i) {
        normalisedInverse(size - 1, i) = lastColumn(i);
        normalisedInverse(i, size - 1) = lastColumn(i);
      }
      TMatrixDSym block(size - 1);
      block.Rank1Update(borderVector, -1.0 / scalar);
      block += topLeftBlock;
      block.Similarity(topLeftInverse);
      TMatrixDSub(normalisedInverse, 0, size - 2, 0, size - 2) = block;
    }
  }

  // undo the normalisation: inverse(matrix) = N * inverse(normalised) * N
  return normalisedInverse.Similarity(normalisation);
}

// =============================================================================
//  Configuration / input
// =============================================================================

void LinearizedHelixVertexFitter::setBeamSpotConstraint(const TVector3& position, const TMatrixDSym& covariance) {
  m_useBeamSpotConstraint = true;
  m_beamSpotPosition = position;
  m_beamSpotInverseCovariance.ResizeTo(3, 3);
  m_beamSpotInverseCovariance = regularizedInverse(covariance);
}

void LinearizedHelixVertexFitter::addTrack(const edm4hep::TrackState& trackState) {
  m_trackParameters.push_back(parametersFromTrackState(trackState));
  m_trackCovariances.push_back(covarianceFromTrackState(trackState));
  m_fitSucceeded = false;
}

void LinearizedHelixVertexFitter::clear() {
  m_trackParameters.clear();
  m_trackCovariances.clear();
  m_trackPhases.clear();
  m_trackChiSquared.clear();
  m_chiSquared = 0.0;
  m_numberOfDegreesOfFreedom = 0;
  m_fitSucceeded = false;
}

// =============================================================================
//  Fast, non-iterative starting vertex (Delphes VtxFitNoSteer)
// =============================================================================

void LinearizedHelixVertexFitter::computeInitialSeed() {
  const std::size_t numberOfTracks = m_trackParameters.size();

  // NOTE: ROOT's TVectorD/TMatrixDSym operator= requires matching dimensions
  // (it does not resize), so we pre-size every element with the right shape.
  std::vector<TVectorD> seedPoint(numberOfTracks, TVectorD(3));                 // helix point at the seed phase
  std::vector<TMatrixDSym> seedPointCovariance(numberOfTracks, TMatrixDSym(3)); // covariance of that point
  std::vector<TVectorD> weightedDirection(numberOfTracks, TVectorD(3));         // C^-1 * (helix tangent)
  std::vector<double> directionNorm(numberOfTracks);                            // tangent^T C^-1 tangent
  std::vector<double> startPhase(numberOfTracks, 0.0);

  // Per-track quantities at the starting phase.
  for (std::size_t track = 0; track < numberOfTracks; ++track) {
    const TVectorD& parameters = m_trackParameters[track];
    const TMatrixDSym& covariance = m_trackCovariances[track];

    // If a start radius is requested and reachable, start the phase there
    // instead of at the perigee (helps for displaced vertices).
    double phase = 0.0;
    if (m_seedStartRadius > std::abs(parameters[kD0])) {
      const double halfCurvature = signedHalfCurvature(parameters[kOmega]);
      const double numerator = m_seedStartRadius * m_seedStartRadius - parameters[kD0] * parameters[kD0];
      const double denominator = 1.0 + 2.0 * halfCurvature * parameters[kD0];
      phase = 2.0 * std::asin(halfCurvature * std::sqrt(numerator / denominator));
    }
    startPhase[track] = phase;

    seedPoint[track] = toVectorD(helixPointAtPhase(parameters, phase));

    const TVectorD tangent = toVectorD(positionDerivativeWrtPhase(parameters, phase));
    const TMatrixD parameterDerivatives = positionDerivativesWrtParameters(parameters, phase);

    TMatrixDSym pointCovariance = covariance;         // 5x5 copy
    pointCovariance.Similarity(parameterDerivatives); // -> A C A^T (3x3)
    seedPointCovariance[track] = pointCovariance;

    const TMatrixDSym pointInformation = regularizedInverse(pointCovariance);
    weightedDirection[track] = pointInformation * tangent;
    directionNorm[track] = pointCovariance.Similarity(weightedDirection[track]);
  }

  // Accumulate the linear system for the seed vertex.
  TMatrixDSym vertexInformation(3);
  vertexInformation.Zero();
  TVectorD vertexInformationTimesPosition(3);
  vertexInformationTimesPosition.Zero();

  for (std::size_t track = 0; track < numberOfTracks; ++track) {
    const TMatrixDSym pointInformation = regularizedInverse(seedPointCovariance[track]);

    // Project out the unmeasured direction along the helix tangent.
    TMatrixDSym tangentProjection(3);
    tangentProjection.Zero();
    tangentProjection.Rank1Update(weightedDirection[track], 1.0 / directionNorm[track]);

    TMatrixDSym contribution = pointInformation;
    contribution -= tangentProjection;

    vertexInformation += contribution;
    vertexInformationTimesPosition += contribution * seedPoint[track];
  }

  if (m_useBeamSpotConstraint) {
    vertexInformation += m_beamSpotInverseCovariance;
    vertexInformationTimesPosition += m_beamSpotInverseCovariance * toVectorD(m_beamSpotPosition);
  }

  const TVectorD seedVertex = regularizedInverse(vertexInformation) * vertexInformationTimesPosition;
  m_vertexPosition = toVector3(seedVertex);

  // Initial per-track phases consistent with the seed vertex.
  m_trackPhases.assign(numberOfTracks, 0.0);
  for (std::size_t track = 0; track < numberOfTracks; ++track) {
    const double phaseShift = dotProduct(weightedDirection[track], seedVertex - seedPoint[track]) /
                              seedPointCovariance[track].Similarity(weightedDirection[track]);
    m_trackPhases[track] = phaseShift + startPhase[track];
  }
}

// =============================================================================
//  Full iterated vertex fit (Delphes VertexFit::VertexFitter)
// =============================================================================

bool LinearizedHelixVertexFitter::fit() {
  const std::size_t numberOfTracks = m_trackParameters.size();

  // Need at least two tracks, unless a beam-spot constraint supplies the
  // missing information.
  if (numberOfTracks < 2 && !m_useBeamSpotConstraint) {
    std::cerr << "LinearizedHelixVertexFitter::fit: need at least two tracks (or a beam-spot constraint)." << std::endl;
    return false;
  }

  // "current" parameters are relinearised (smoothed) every iteration, while the
  // "measured" parameters are the original measurements that stay fixed.
  std::vector<TVectorD> measuredParameters = m_trackParameters;
  std::vector<TVectorD> currentParameters = m_trackParameters;

  m_trackChiSquared.assign(numberOfTracks, 0.0);

  computeInitialSeed();
  TVector3 previousVertex = m_vertexPosition;

  int iteration = 0;
  double squaredDisplacement = std::numeric_limits<double>::max();

  while (squaredDisplacement > m_convergenceThreshold && iteration < m_maxIterations) {
    // Per-track work quantities, recomputed at the current linearisation point.
    // Pre-size with the right shapes (ROOT operator= does not resize).
    std::vector<TMatrixD> transposedParameterDerivatives(numberOfTracks, TMatrixD(5, 3)); // A^T   (5x3)
    std::vector<TMatrixDSym> pointInverseCovariance(numberOfTracks, TMatrixDSym(3));      // A C A^T (3x3)
    std::vector<TVectorD> helixPoint(numberOfTracks, TVectorD(3));                        // x(phase) (3)
    std::vector<TVectorD> linearisationShift(numberOfTracks, TVectorD(3));                // A (current-measured) (3)
    std::vector<TMatrixDSym> pointWeight(numberOfTracks, TMatrixDSym(3));                 // (A C A^T)^-1 (3x3)
    std::vector<TVectorD> tangent(numberOfTracks, TVectorD(3));                           // dX/dphase (3)
    std::vector<double> tangentNorm(numberOfTracks);                                      // tangent^T W tangent
    std::vector<TMatrixDSym> projection(numberOfTracks, TMatrixDSym(3));                  // D matrix (3x3)

    TMatrixDSym hessian(3);
    hessian.Zero();
    TVectorD constantTerm(3);
    constantTerm.Zero();
    TMatrixDSym weightedProjectionSum(3); // for the vertex covariance
    weightedProjectionSum.Zero();

    // ---- first pass over tracks: build the linear system ----
    for (std::size_t track = 0; track < numberOfTracks; ++track) {
      const double phase = m_trackPhases[track];
      const TVectorD& parameters = currentParameters[track];
      const TMatrixDSym& covariance = m_trackCovariances[track];

      const TMatrixD parameterDerivatives = positionDerivativesWrtParameters(parameters, phase); // A
      transposedParameterDerivatives[track] = TMatrixD(TMatrixD::kTransposed, parameterDerivatives);

      TMatrixDSym inverseCovariance = covariance;         // 5x5 copy
      inverseCovariance.Similarity(parameterDerivatives); // -> A C A^T (3x3)
      pointInverseCovariance[track] = inverseCovariance;

      helixPoint[track] = toVectorD(helixPointAtPhase(parameters, phase));
      linearisationShift[track] = parameterDerivatives * (parameters - measuredParameters[track]);

      const TMatrixDSym weight = regularizedInverse(inverseCovariance); // W
      pointWeight[track] = weight;

      tangent[track] = toVectorD(positionDerivativeWrtPhase(parameters, phase));
      tangentNorm[track] = weight.Similarity(tangent[track]);

      // projection matrix D = W - (W tangent)(W tangent)^T / (tangent^T W tangent)
      TMatrixDSym tangentTerm(3);
      tangentTerm.Zero();
      tangentTerm.Rank1Update(tangent[track], -1.0 / tangentNorm[track]);
      tangentTerm.Similarity(weight);
      TMatrixDSym projectionMatrix = weight;
      projectionMatrix += tangentTerm;
      projection[track] = projectionMatrix;

      // accumulate global system
      TMatrixDSym weightedProjection = inverseCovariance; // copy, becomes D (A C A^T) D^T... see below
      weightedProjection.Similarity(projectionMatrix);    // D * (A C A^T) * D^T
      weightedProjectionSum += weightedProjection;
      hessian += projectionMatrix;
      constantTerm += projectionMatrix * (helixPoint[track] - linearisationShift[track]);
    }

    if (m_useBeamSpotConstraint) {
      hessian += m_beamSpotInverseCovariance;
      constantTerm += m_beamSpotInverseCovariance * toVectorD(m_beamSpotPosition);
      weightedProjectionSum += m_beamSpotInverseCovariance;
    }

    // ---- solve for the updated vertex ----
    const TMatrixDSym inverseHessian = regularizedInverse(hessian);
    const TVectorD updatedVertex = inverseHessian * constantTerm;
    m_vertexCovariance = weightedProjectionSum.Similarity(inverseHessian);

    // ---- second pass: chi2, phases and smoothed track parameters ----
    m_chiSquared = 0.0;
    for (std::size_t track = 0; track < numberOfTracks; ++track) {
      const TVectorD residual = projection[track] * (helixPoint[track] - updatedVertex - linearisationShift[track]);
      m_trackChiSquared[track] = pointInverseCovariance[track].Similarity(residual);
      m_chiSquared += m_trackChiSquared[track];

      // update this track's helix phase
      const TVectorD weightedResidual =
          pointWeight[track] * (updatedVertex - helixPoint[track] + linearisationShift[track]);
      m_trackPhases[track] += dotProduct(tangent[track], weightedResidual) / tangentNorm[track];

      // relinearise (smooth) the track parameters around the new vertex
      const TVectorD parameterUpdate = (m_trackCovariances[track] * transposedParameterDerivatives[track]) * residual;
      currentParameters[track] = measuredParameters[track] - parameterUpdate;
    }
    if (m_useBeamSpotConstraint)
      m_chiSquared += m_beamSpotInverseCovariance.Similarity(updatedVertex - toVectorD(m_beamSpotPosition));

    // ---- convergence test ----
    const TVectorD vertexStep = updatedVertex - toVectorD(previousVertex);
    const TMatrixDSym vertexInformation = regularizedInverse(m_vertexCovariance);
    squaredDisplacement = vertexInformation.Similarity(vertexStep);

    m_vertexPosition = toVector3(updatedVertex);
    previousVertex = m_vertexPosition;
    ++iteration;
  }

  // degrees of freedom: each track provides 2 constraints, the vertex has 3
  // free coordinates; a beam-spot constraint adds 3 more measurements.
  const int numberOfMeasurements = 2 * static_cast<int>(numberOfTracks) + (m_useBeamSpotConstraint ? 3 : 0);
  m_numberOfDegreesOfFreedom = numberOfMeasurements - 3;

  m_fitSucceeded = true;
  return true;
}
