/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

#include "GenfitTrack.h"

namespace GenfitInterface {

GenfitTrack::GenfitTrack(const edm4hep::Track& track, const bool skipTrackOrdering,
                         const dd4hep::rec::WireTracker_info_struct* wire_info,
                         const dd4hep::DDSegmentation::BitFieldCoder* decoder,
                         const GenfitInterface::GenfitField* fieldMap)
    : m_originalTrack(track), m_posInit(0., 0., 0.), m_momInit(0., 0., 0.), m_covInit(6), m_edm4hepTrack(),
      m_wireInfo(wire_info), m_dcDecoder(decoder), m_fieldMap(fieldMap)

{

  CheckInitialization();
  OrderHits(track, skipTrackOrdering);
}

GenfitTrack::~GenfitTrack() {}

/**
 * @brief Check if required Genfit components are properly initialized.
 *
 * This method verifies whether the singleton instances of `genfit::FieldManager`
 * and `genfit::MaterialEffects` have been initialized. These components are essential
 * for tracking operations in the Genfit framework.
 *
 * @note This method should be called before performing any tracking-related operations
 * that depend on magnetic field or material effects.
 */
void GenfitTrack::CheckInitialization() {

  if (!genfit::FieldManager::getInstance()->isInitialized()) {
    throw std::runtime_error("Error: FieldManager is not initialized!");
  }

  if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
    throw std::runtime_error("Error: MaterialEffects is not initialized!");
  }
}

/**
 * @brief Orders the tracker hits of a track according to their spatial progression.
 *
 * This method determines a suitable starting point among the track hits,
 * typically the one closest to the interaction region (|r| minimum).
 *
 * All hits are then sorted according to their 3D distance from the chosen
 * starting point, ensuring a consistent ordering along the track trajectory.
 * The ordered hits are stored in the internal MutableTrack object for
 * subsequent processing (e.g. track fitting).
 *
 * @param track The input track whose tracker hits must be reordered.
 * @param skipTrackOrdering If true, skip the track ordering step.
 */
void GenfitTrack::OrderHits(const edm4hep::Track& track, bool skipTrackOrdering) {

  if (skipTrackOrdering) {

    bool first = true;
    for (const auto& hit : track.getTrackerHits()) {

      m_edm4hepTrack.addToTrackerHits(hit);
      auto pos = hit.getPosition();

      if (first) {
        m_firstHitReferencePoint = TVector3(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);
        first = false;
      }

      m_lastHitReferencePoint = TVector3(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);
    }
    return;
  }

  const auto hits = track.getTrackerHits();
  std::vector<std::pair<float, std::size_t>> distIndex;

  // Track endpoints in cylindrical radius
  double rMin = std::numeric_limits<double>::max();
  double firstX = 0., firstY = 0., firstZ = 0.;

  // Loop over hits to find the hit with smallest r = sqrt(x^2 + y^2)
  for (const auto& hit : hits) {
    const auto p = hit.getPosition();
    double r = std::hypot(p.x, p.y);

    if (r < rMin) {
      rMin = r;
      firstX = p.x;
      firstY = p.y;
      firstZ = p.z;
    }
  }

  // Compute distances from the chosen starting point
  distIndex.reserve(hits.size());

  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto p = hits[i].getPosition();
    const float d = std::hypot(p.x - firstX, p.y - firstY, p.z - firstZ);

    distIndex.emplace_back(d, i);
  }

  // Sort hits along the track
  std::ranges::sort(distIndex, {}, &std::pair<float, std::size_t>::first);

  bool first = true;
  for (const auto& [_, idx] : distIndex) {
    const auto& hit = hits[idx];
    auto pos = hit.getPosition();
    m_edm4hepTrack.addToTrackerHits(hit);

    if (first) {
      m_firstHitReferencePoint = TVector3(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);
      first = false;
    }

    m_lastHitReferencePoint = TVector3(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);
  }
}

/**
 * @brief Initialize the GenfitTrack object with configurable initialization strategies.
 *
 * This method initializes the internal state of the `GenfitTrack` object by defining:
 * - the initial position (`m_posInit`)
 * - the initial momentum (`m_momInit`)
 * - the charge hypothesis (`m_chargeHypothesis`)
 * - the initial covariance matrix (`m_covInit`)
 *
 * The initialization procedure depends on the selected strategy and may rely on:
 * - tracker hits,
 * - pre-existing track states,
 * - or user-provided parameters.
 *
 * In all cases, the magnetic field at the reference position is evaluated and used
 * to ensure a consistent initialization of both kinematics and covariance.
 *
 * ---
 * Initialization workflow:
 *
 * 1. Looper handling:
 *    - The tracker hits are preprocessed with `SelectFirstArcHits(LooperPTThreshold)`.
 *    - For looping tracks only the hits of the first arc of the loop (the hits before the first
 *      peak of the sawtooth in the (z, phi_local) plane) are used to seed and fit the track.
 *
 * 2. Definition of initial track parameters (controlled by InitializationType):
 *
 *    - InitializationType == 0 (default: geometric seed from hits):
 *        - Requires at least two tracker hits.
 *        - The initial position is set to the first hit.
 *        - The initial momentum direction is defined as the normalized vector
 *          from the first to the second hit.
 *        - The charge hypothesis is estimated via `ComputeInitialParameters(Bz)`.
 *        - The covariance matrix is computed using `ComputeInitialCovarianceMatrix`.
 *
 *    - InitializationType == 1 (refined: parameter-based estimation):
 *        - Uses `ComputeInitialParameters(Bz)` to obtain a more accurate estimate
 *          of position and momentum, typically exploiting multiple hits.
 *        - Optionally updates the reference point (`m_vpReferencePoint`):
 *            - if `UseFirstHitAsReference` is true, or
 *            - if displaced tracking is enabled and the first hit lies outside
 *              a sphere of radius `RadiusForDisplacedTracking`.
 *        - The covariance matrix is computed consistently with the estimated parameters.
 *
 *    - InitializationType == 2 (track state based):
 *        - Extracts initialization from an existing `edm4hep::TrackState`.
 *        - The reference point is taken from the track state.
 *        - The momentum is reconstructed from helix parameters:
 *            - transverse momentum is derived from curvature (`omega`) and Bz
 *            - direction is defined using `phi` and `tanLambda`
 *        - The covariance matrix is:
 *            - first reconstructed in helix parameter space
 *            - then transformed to Cartesian coordinates via
 *              `InitialCovarianceMatrixHelixToCartesian`.
 *
 *    - InitializationType == 3 (custom user-defined):
 *        - Uses user-provided `Init_position` and `Init_momentum`.
 *        - Both parameters are mandatory.
 *        - The charge hypothesis is still inferred via `ComputeInitialParameters(Bz)`.
 *        - The covariance matrix is computed accordingly.
 *
 * ---
 *
 * Error handling:
 * - The method throws `std::runtime_error` in the following cases:
 *     - missing required parameters (e.g. TrackStateLocation, Init_position, ...)
 *     - insufficient number of hits for the selected strategy
 *     - requested track state not found
 *     - unknown InitializationType
 *
 * ---
 *
 * @param RadiusForDisplacedTracking If > 0, enables displaced tracking logic:
 *        tracks with first hit outside this radius may use a different reference point.
 * @param UseFirstHitAsReference If true, the first hit is always used as reference point.
 * @param InitializationType Strategy selector:
 *        0 = default (first two hits),
 *        1 = refined (computed parameters),
 *        2 = track state (from edm4hep::TrackState),
 *        3 = custom (user-defined).
 * @param TrackStateLocation Required if InitializationType == 2:
 *        0 = edm4hep::TrackState::AtOther
 *        1 = edm4hep::TrackState::AtFirstHit
 *        2 = edm4hep::TrackState::AtLastHit
 *        3 = edm4hep::TrackState::AtCalorimeter
 * @param Init_position Custom initial position (required if InitializationType == 3).
 * @param Init_momentum Custom initial momentum (required if InitializationType == 3).
 * @param LooperPTThreshold Transverse momentum threshold [GeV/c] for the looper treatment: below it the first-arc
 * hits are selected whenever the sawtooth of phi_local shows a wrap.
 * @param sigma_d0        Optional standard deviation for d0 (used in covariance construction).
 * @param sigma_phi       Optional standard deviation for phi (used in covariance construction).
 * @param omega_factor    Optional factor for omega (used in covariance construction).
 * @param z0_factor       Optional factor for z0 (used in covariance construction).
 * @param sigma_tanLambda Optional standard deviation for tanLambda (used in covariance construction).
 *
 *
 * @note The optional sigma parameters (sigma_d0, sigma_phi, omega_factor, z0_factor, sigma_tanLambda) are used to
 * construct the initial covariance matrix in the absence of a track state (InitializationType == 0,1,3).
 *
 * They are applied only in:
 * - InitializationType == 0 (default)
 * - InitializationType == 1 (refined initialization)
 * - InitializationType == 3 (custom initialization)
 *
 * If not provided, default physically motivated values are used internally.
 * These parameters are not used in InitializationType == 2, where the
 * covariance is directly taken from the input TrackState and transformed
 * to Cartesian coordinates.
 */
void GenfitTrack::InitializeTrack(double RadiusForDisplacedTracking, bool UseFirstHitAsReference,
                                  int InitializationType, std::optional<int> TrackStateLocation,
                                  std::optional<TVector3> Init_position, std::optional<TVector3> Init_momentum,
                                  double LooperPTThreshold, double MomentumSeedOverride, std::optional<double> sigma_d0,
                                  std::optional<double> sigma_phi, std::optional<double> omega_factor,
                                  std::optional<double> z0_factor, std::optional<double> sigma_tanLambda) {


  // Looping tracks are seeded and fitted with the hits of the first arc only: the looper
  // identification (pT estimate of the xy circle fit vs LooperPTThreshold) is made here, once,
  // and the resulting selection feeds every fit attempt of the track.
  SelectFirstArcHits(LooperPTThreshold);

  // =========================
  // Initialization strategies
  // =========================

  if (InitializationType == 0) // --- Default: first two hits
  {
    auto hits = m_edm4hepTrack.getTrackerHits();

    if (m_fitHits.size() < 2) {
      throw std::runtime_error("InitializeTrack: Not enough hits for InitializationType == 0");
    }

    auto pos1 = hits[m_fitHits[0]].getPosition();
    auto pos2 = hits[m_fitHits[1]].getPosition();

    TVector3 p1(pos1.x * dd4hep::mm, pos1.y * dd4hep::mm, pos1.z * dd4hep::mm);
    TVector3 p2(pos2.x * dd4hep::mm, pos2.y * dd4hep::mm, pos2.z * dd4hep::mm);

    // Seed definition
    m_posInit = p1;
    m_momInit = (p2 - p1).Unit();

    double Bz = m_fieldMap->getBz(m_posInit) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    HelperInitialization initInfo = ComputeInitialParameters(Bz);
    m_chargeHypothesis = initInfo.Charge;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_chargeHypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else if (InitializationType == 1) // --- Refined
  {
    auto hits = m_edm4hepTrack.getTrackerHits();

    if (hits.empty()) {
      throw std::runtime_error("InitializeTrack: No hits available for InitializationType == 1");
    }

    auto firstHit = hits[m_fitHits[0]].getPosition();

    TVector3 firstHitVec(firstHit.x * dd4hep::mm, firstHit.y * dd4hep::mm, firstHit.z * dd4hep::mm);

    double Bz = m_fieldMap->getBz(firstHitVec) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    // Optional reference point update
    if (UseFirstHitAsReference || (RadiusForDisplacedTracking > 0. && firstHitVec.Mag() > RadiusForDisplacedTracking)) {
      m_vpReferencePoint = firstHitVec;
    }

    HelperInitialization initInfo = ComputeInitialParameters(Bz);

    m_chargeHypothesis = initInfo.Charge;
    m_posInit = initInfo.Position;
    m_momInit = initInfo.Momentum;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_chargeHypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else if (InitializationType == 2) // --- From TrackState
  {
    if (!TrackStateLocation.has_value()) {
      throw std::runtime_error("InitializeTrack: TrackStateLocation is required for InitializationType == 2");
    }

    auto trackStates = m_originalTrack.getTrackStates();
    bool found = false;

    for (auto ts : trackStates) {
      if (ts.location == TrackStateLocation.value()) {
        found = true;

        m_vpReferencePoint = TVector3(ts.referencePoint.x * dd4hep::mm, ts.referencePoint.y * dd4hep::mm,
                                      ts.referencePoint.z * dd4hep::mm); // in cm

        double Bz =
            m_fieldMap->getBz(m_vpReferencePoint) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

        // Reconstruct the PCA described by this exact TrackState
        m_posInit = TVector3(m_vpReferencePoint.X() - ts.D0 * std::sin(ts.phi) * dd4hep::mm,
                             m_vpReferencePoint.Y() + ts.D0 * std::cos(ts.phi) * dd4hep::mm,
                             m_vpReferencePoint.Z() + ts.Z0 * dd4hep::mm);
        m_chargeHypothesis = (Bz / ts.omega >= 0.) ? 1 : -1;

        // Helix -> momentum conversion
        double pT = ConversionUnits::a_lcio * std::abs(Bz) / std::abs(ts.omega);
        double px = pT * std::cos(ts.phi);
        double py = pT * std::sin(ts.phi);
        double pz = pT * ts.tanLambda;

        m_momInit = TVector3(px, py, pz);

        // Build symmetric covariance matrix (helix parameters)
        TMatrixDSym C_helix(5);
        C_helix.Zero();

        int index = 0;
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j <= i; ++j) {
            C_helix(i, j) = ts.covMatrix[index++];
          }
        }

        for (int i = 0; i < 5; ++i) {
          for (int j = i + 1; j < 5; ++j) {
            C_helix(i, j) = C_helix(j, i);
          }
        }

        m_covInit = InitialCovarianceMatrixHelixToCartesian(C_helix, m_posInit, m_momInit, m_vpReferencePoint,
                                                            m_chargeHypothesis, Bz);
      }
    }

    if (!found) {
      throw std::runtime_error("InitializeTrack: Requested TrackStateLocation not found");
    }
  }

  else if (InitializationType == 3) // --- Custom
  {
    if (!(Init_position.has_value() && Init_momentum.has_value())) {
      throw std::runtime_error(
          "InitializeTrack: Init_position and Init_momentum are required for InitializationType == 3");
    }

    m_posInit = Init_position.value();
    m_momInit = Init_momentum.value();

    double Bz = m_fieldMap->getBz(m_posInit) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    HelperInitialization initInfo = ComputeInitialParameters(Bz);
    m_chargeHypothesis = initInfo.Charge;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_chargeHypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else {
    throw std::runtime_error("InitializeTrack: Unknown InitializationType");
  }

  // Optional seed override (recovery attempts): keep the seeded direction but impose the
  // momentum magnitude, with a consistent initial covariance
  if (MomentumSeedOverride > 0. && m_momInit.Mag() > 0.) {
    m_momInit = m_momInit.Unit() * MomentumSeedOverride;
    const double Bz = m_fieldMap->getBz(m_posInit) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla
    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_chargeHypothesis, sigma_d0, sigma_phi, omega_factor, z0_factor,
                                               sigma_tanLambda);
  }
}

/**
 * @brief Identify the first arc of a looping track and select the hits used for the fit.
 *
 * Below a few hundred MeV/c charged particles curl up in the magnetic field and the
 * reconstructed track contains several loops of the same particle: the loops cannot be
 * fitted with a single helix, only the hits of the first arc can.
 *
 * The first arc is the longest monotone-phi chain rooted at the trajectory start
 * (FirstArcPhiChain):
 *   1. a circle is fitted to the xy projection of all the hits (FastCircleFit);
 *   2. the anchor, i.e. the start of the trajectory, is the hit closest to the reference
 *      point (the IP for prompt tracks);
 *   3. phi_local, the angle around the center of the fitted circle, and the radius from that
 *      center order the hits along one turn (the turns of a spiraling track differ in radius
 *      by the energy loss per turn); the first arc is the first 2*pi of the sequence rooted
 *      at the anchor.
 * The construction does not use z, so it works also at theta ~ 90 deg, where the z ordering
 * degenerates.
 *
 * The pT estimated from the fitted circle (pT = 0.3 * Bz * rho) selects which tracks can
 * loop at all: the first-arc hits are used only for pT < pTThresholdGeV. If the chain cannot
 * be identified (fewer than 3 hits) all the hits are kept.
 *
 * m_edm4hepTrack keeps all its hits: only m_fitHits (the hits used for seeding and fitting)
 * is restricted to the first arc.
 *
 * @param pTThresholdGeV Transverse momentum threshold for the looper treatment [GeV/c].
 */
void GenfitTrack::SelectFirstArcHits(double pTThresholdGeV) {

  auto hits = m_edm4hepTrack.getTrackerHits();
  const std::size_t N = hits.size();

  // Default: seed and fit with all the hits
  m_fitHits.resize(N);
  std::iota(m_fitHits.begin(), m_fitHits.end(), std::size_t{0});

  if (N < 3) {
    return;
  }

  // Hit positions [mm]
  std::vector<Point2D_xy> points_xy;
  points_xy.reserve(N);

  for (const auto& hit : hits) {
    const auto pos = hit.getPosition();
    points_xy.emplace_back(pos.x, pos.y);
  }

  // The circle fit depends only on the hit geometry: it is computed once per track and
  // re-used across fit attempts and particle hypotheses through the per-track cache
  // (nothing is shared between tracks or events).
  TrackFitCache::Geometry localGeometry;
  TrackFitCache::Geometry& geo = m_cache ? m_cache->geometry : localGeometry;

  if (!geo.valid) {

    // Circle fit in the xy plane
    FastCircleFit circle(points_xy);
    geo.x0 = circle.x0();
    geo.y0 = circle.y0();
    geo.rho = circle.rho();

    // Looper gate: pT [GeV/c] = 0.3 * Bz [T] * rho [m], with rho from the fitted circle
    const auto refPos = hits[0].getPosition();
    const double Bz = m_fieldMap->getBz(TVector3(refPos.x * dd4hep::mm, refPos.y * dd4hep::mm, refPos.z * dd4hep::mm)) /
                      (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla
    geo.pT_est = std::abs(0.3 * Bz * geo.rho) / 1000.0;    // 1000 is used to convert from m to mm
    geo.valid = true;
  }

  if (geo.pT_est >= pTThresholdGeV) {
    return;
  }

  // Anchor: the start of the trajectory, i.e. the hit closest to the reference point (the IP
  // for prompt tracks). phi_local and radius around the fitted circle center: the turns of a
  // spiraling track differ in radius by the energy loss per turn, so the radius separates the 
  // turns that z cannot.
  std::size_t anchor = 0;
  double bestDist2 = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < N; ++i) {
    const auto p = hits[i].getPosition();
    const double dx = p.x * dd4hep::mm - m_vpReferencePoint.X();
    const double dy = p.y * dd4hep::mm - m_vpReferencePoint.Y();
    const double dz = p.z * dd4hep::mm - m_vpReferencePoint.Z();
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < bestDist2) {
      bestDist2 = d2;
      anchor = i;
    }
  }

  std::vector<double> phi(N);
  std::vector<double> radial(N);
  for (std::size_t i = 0; i < N; ++i) {
    phi[i] = std::atan2(points_xy[i].y - geo.y0, points_xy[i].x - geo.x0);
    radial[i] = std::hypot(points_xy[i].x - geo.x0, points_xy[i].y - geo.y0);
  }

  auto arc = FirstArcPhiChain(phi, radial, anchor);
  if (arc.size() < 3) {
    return; // the arc is not identifiable: keep all the hits
  }

  m_fitHits = std::move(arc);
  const auto arcFirst = hits[m_fitHits.front()].getPosition();
  const auto arcLast = hits[m_fitHits.back()].getPosition();
  m_firstHitReferencePoint = TVector3(arcFirst.x * dd4hep::mm, arcFirst.y * dd4hep::mm, arcFirst.z * dd4hep::mm);
  m_lastHitReferencePoint = TVector3(arcLast.x * dd4hep::mm, arcLast.y * dd4hep::mm, arcLast.z * dd4hep::mm);
}

namespace {

  // Mean radius of the last third minus the first third of a sequence: negative when the
  // sequence follows the direction of travel of a track that spirals inward (energy loss).
  double ArcRadialTrend(const std::vector<std::size_t>& seq, const std::vector<double>& radial) {
    const std::size_t n = seq.size();
    const std::size_t third = std::max<std::size_t>(1, n / 3);
    double first = 0., last = 0.;
    for (std::size_t i = 0; i < third; ++i) {
      first += radial[seq[i]];
      last += radial[seq[n - 1 - i]];
    }
    return (last - first) / third;
  }

} // namespace

/**
 * @brief First arc as the longest monotone-phi chain from the anchor hit.
 *
 * phi_local is measured from the anchor along each direction of travel (the wrap of the
 * circle is absorbed by the 2*pi-periodic coordinate) and the longest chain with increasing
 * phi_local, bounded steps and radially continuous neighbours is selected (a longest-
 * increasing-subsequence with a step constraint). Hits that cannot be chained are left out,
 * so stray hits cannot derail the sequence the way they can a greedy nearest-neighbour walk.
 */
std::vector<std::size_t> GenfitTrack::FirstArcPhiChain(const std::vector<double>& phi,
                                                       const std::vector<double>& radial, std::size_t anchor) const {
  constexpr double kTwoPi = 6.28318530717958647692;
  constexpr double kMaxStep = 0.5; // [rad]
  constexpr double kRTol = 5.0;    // [mm] radial continuity within one turn

  const std::size_t N = phi.size();
  std::vector<std::size_t> best;
  double bestTrend = 0.;

  for (double dir : {1., -1.}) {
    std::vector<double> u(N);
    for (std::size_t i = 0; i < N; ++i) {
      double d = (phi[i] - phi[anchor]) * dir;
      d -= kTwoPi * std::floor(d / kTwoPi);
      u[i] = d;
    }
    std::vector<std::size_t> order(N);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(order, [&u](std::size_t i, std::size_t j) { return u[i] < u[j]; });

    // longest chain with increasing u, bounded step and radial continuity, rooted at the anchor
    std::vector<int> length(N, 0);
    std::vector<std::size_t> prev(N, N);
    length[anchor] = 1;
    for (std::size_t a = 0; a < N; ++a) {
      const std::size_t i = order[a];
      if (length[i] == 0)
        continue;
      for (std::size_t b = a + 1; b < N; ++b) {
        const std::size_t j = order[b];
        const double step = u[j] - u[i];
        if (step > kMaxStep)
          break; // u is sorted: no later hit can be closer
        if (std::abs(radial[j] - radial[i]) > kRTol)
          continue;
        if (length[j] < length[i] + 1) {
          length[j] = length[i] + 1;
          prev[j] = i;
        }
      }
    }

    std::size_t end = anchor;
    for (std::size_t i = 0; i < N; ++i) {
      if (length[i] > length[end])
        end = i;
    }
    std::vector<std::size_t> seq;
    for (std::size_t i = end; i != N; i = prev[i]) {
      seq.push_back(i);
    }
    std::ranges::reverse(seq);

    const double trend = ArcRadialTrend(seq, radial);
    if (best.empty() || trend < bestTrend) {
      best = std::move(seq);
      bestTrend = trend;
    }
  }
  return best;
}

/**
 * @brief Restrict the hits used for seeding and fitting to the silicon (planar) hits.
 *
 * Last-resort fallback when the fit with all the hits, or with the first-arc hits of a looper,
 * fails: the silicon hits are much more precise than the drift chamber ones (whose z resolution
 * is of the order of cm) and a fit with them alone often converges where the full one does not.
 * The output track keeps all its hits as before.
 *
 * @return true if at least 3 silicon hits are available for the fit.
 */
bool GenfitTrack::RestrictToSiliconHits() {

  auto hits = m_edm4hepTrack.getTrackerHits();

  std::vector<std::size_t> silicon;
  silicon.reserve(m_fitHits.size());
  for (const std::size_t idx : m_fitHits) {
    if (hits[idx].isA<edm4hep::TrackerHitPlane>()) {
      silicon.push_back(idx);
    }
  }

  if (silicon.size() < 3) {
    return false;
  }

  m_fitHits = std::move(silicon);

  const auto firstPos = hits[m_fitHits.front()].getPosition();
  const auto lastPos = hits[m_fitHits.back()].getPosition();
  m_firstHitReferencePoint = TVector3(firstPos.x * dd4hep::mm, firstPos.y * dd4hep::mm, firstPos.z * dd4hep::mm);
  m_lastHitReferencePoint = TVector3(lastPos.x * dd4hep::mm, lastPos.y * dd4hep::mm, lastPos.z * dd4hep::mm);

  return true;
}

/**
 * @brief Thin the hits used for seeding and fitting to at most `target` evenly spaced hits.
 *
 * Recovery fallback for tracks whose fit does not converge even on the first arc: a strongly
 * thinned hit set accumulates less multiple scattering and measurement noise, and the fit often
 * converges where the full set does not. The output track keeps all its hits as before.
 *
 * @param target Maximum number of hits kept; must be >= 3, otherwise the hit set is unchanged.
 */
void GenfitTrack::RestrictToEvenlySpacedHits(std::size_t target) {

  auto hits = m_edm4hepTrack.getTrackerHits();
  const std::size_t n = m_fitHits.size();

  if (target < 3 || n <= target) {
    return;
  }

  std::vector<std::size_t> thinned;
  thinned.reserve(target);
  const std::size_t stride = (n + target - 1) / target;
  for (std::size_t k = 0; k < n && thinned.size() < target; k += stride) {
    thinned.push_back(m_fitHits[k]);
  }
  thinned.back() = m_fitHits.back(); // the last hit of the arc is always kept
  m_fitHits = std::move(thinned);

  const auto firstPos = hits[m_fitHits.front()].getPosition();
  const auto lastPos = hits[m_fitHits.back()].getPosition();
  m_firstHitReferencePoint = TVector3(firstPos.x * dd4hep::mm, firstPos.y * dd4hep::mm, firstPos.z * dd4hep::mm);
  m_lastHitReferencePoint = TVector3(lastPos.x * dd4hep::mm, lastPos.y * dd4hep::mm, lastPos.z * dd4hep::mm);
}

/**
 * @brief Restrict the hits used for the fit to the inliers of a robust (RANSAC) circle fit,
 *        ordered along that circle.
 *
 * Recovery fallback for impure tracks (e.g. two particles merged into one hit list): the global
 * circle fit lands between the particles and every fit attempt fails or returns garbage. The
 * inliers of the dominant circle follow one turn of one particle, and ordering them by their
 * angle around that circle (sequence starting after the largest angular gap) gives GENFIT the
 * path order it needs. The output track keeps all its hits as before.
 *
 * @param toleranceMm Maximum distance from the fitted circle for a hit to count as an inlier.
 * @return true if at least 8 inliers were found.
 */
bool GenfitTrack::RestrictToCircleInliers(double toleranceMm) {

  auto hits = m_edm4hepTrack.getTrackerHits();
  const std::size_t n = m_fitHits.size();
  if (n < 8) {
    return false;
  }

  std::vector<Point2D_xy> pts;
  pts.reserve(n);
  for (const std::size_t idx : m_fitHits) {
    const auto p = hits[idx].getPosition();
    pts.emplace_back(p.x, p.y);
  }

  // RANSAC: circumcircles of deterministic hit triplets, keep the largest inlier set
  std::size_t bestCount = 0;
  double bestX = 0., bestY = 0., bestRho = 0.;
  std::vector<std::size_t> bestInliers;

  for (std::size_t k = 0; k < 80; ++k) {
    const std::size_t a = (k * 13) % n;
    const std::size_t b = (a + n / 3) % n;
    const std::size_t c = (b + n / 3) % n;
    const Point2D_xy& p1 = pts[a];
    const Point2D_xy& p2 = pts[b];
    const Point2D_xy& p3 = pts[c];

    const double det = 2. * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
    if (std::abs(det) < 1e-9) {
      continue;
    }
    const double s1 = p1.x * p1.x + p1.y * p1.y;
    const double s2 = p2.x * p2.x + p2.y * p2.y;
    const double s3 = p3.x * p3.x + p3.y * p3.y;
    const double cx = (s1 * (p2.y - p3.y) + s2 * (p3.y - p1.y) + s3 * (p1.y - p2.y)) / det;
    const double cy = (s1 * (p3.x - p2.x) + s2 * (p1.x - p3.x) + s3 * (p2.x - p1.x)) / det;
    const double rho = std::hypot(p1.x - cx, p1.y - cy);
    if (!(rho > 0.) || rho > 1e6) {
      continue;
    }

    const double tol = std::max(toleranceMm, 0.01 * rho);
    std::vector<std::size_t> inliers;
    inliers.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      if (std::abs(std::hypot(pts[i].x - cx, pts[i].y - cy) - rho) <= tol) {
        inliers.push_back(i);
      }
    }
    if (inliers.size() > bestCount) {
      bestCount = inliers.size();
      bestX = cx;
      bestY = cy;
      bestRho = rho;
      bestInliers = std::move(inliers);
    }
  }

  if (bestCount < 8) {
    return false;
  }

  // Order the inliers by their angle around the dominant circle, sequence starting after the
  // largest angular gap (the seam of the turn) so that consecutive hits stay neighbours
  std::vector<std::pair<double, std::size_t>> angleIndex;
  angleIndex.reserve(bestCount);
  for (const std::size_t i : bestInliers) {
    angleIndex.emplace_back(std::atan2(pts[i].y - bestY, pts[i].x - bestX), i);
  }
  std::ranges::sort(angleIndex);

  std::size_t seam = 0;
  double gapMax = -1.;
  constexpr double kTwoPi = 6.28318530717958647692;
  for (std::size_t k = 0; k < angleIndex.size(); ++k) {
    const double gap = (k + 1 < angleIndex.size()) ? angleIndex[k + 1].first - angleIndex[k].first
                                                   : angleIndex[0].first + kTwoPi - angleIndex[k].first;
    if (gap > gapMax) {
      gapMax = gap;
      seam = k;
    }
  }

  std::vector<std::size_t> ordered;
  ordered.reserve(bestCount);
  for (std::size_t k = 0; k < angleIndex.size(); ++k) {
    ordered.push_back(angleIndex[(seam + 1 + k) % angleIndex.size()].second);
  }

  std::vector<std::size_t> fitOrdered;
  fitOrdered.reserve(bestCount);
  for (const std::size_t i : ordered) {
    fitOrdered.push_back(m_fitHits[i]);
  }
  m_fitHits = std::move(fitOrdered);

  const auto firstPos = hits[m_fitHits.front()].getPosition();
  const auto lastPos = hits[m_fitHits.back()].getPosition();
  m_firstHitReferencePoint = TVector3(firstPos.x * dd4hep::mm, firstPos.y * dd4hep::mm, firstPos.z * dd4hep::mm);
  m_lastHitReferencePoint = TVector3(lastPos.x * dd4hep::mm, lastPos.y * dd4hep::mm, lastPos.z * dd4hep::mm);

  return true;
}

/**
 * @brief Computes the initial covariance matrix of the track parameters
 *        in Cartesian coordinates starting from helix parameter uncertainties.
 *
 * This method builds an initial estimate of the track parameter covariance
 * in helix parameter space and propagates it to Cartesian coordinates using
 * a transformation that depends on the magnetic field and the initial state
 * (position and momentum).
 *
 * The procedure is as follows:
 *
 * - Initializes a 5x5 diagonal covariance matrix in helix parameter space
 *   (d0, phi, omega, z0, tanLambda), assuming no initial correlations.
 *
 * - d0 (transverse impact parameter):
 *     Uses sigma_d0 if provided, otherwise a default value of 0.05.
 *     Sets the variance as sigma_d0 squared.
 *
 * - phi (azimuthal angle):
 *     Uses sigma_phi if provided, otherwise a default value of 0.1.
 *     Sets the variance as sigma_phi squared.
 *
 * - omega (curvature):
 *     Computes curvature from the initial momentum magnitude in the transverse plane
 *     and the magnetic field Bz, using a conversion constant.
 *     The uncertainty is scaled by omega_factor (default 0.5).
 *     The variance is (omega_factor * omega)^2.
 *
 * - z0 (longitudinal impact parameter):
 *     Uses the absolute value of the initial z position as a scale.
 *     Scales it using z0_factor (default 0.1).
 *     The variance is (z0_factor * |z_init_z|)^2.
 *
 * - tanLambda (dip angle):
 *     Uses sigma_tanLambda if provided, otherwise a default value of 0.1.
 *     Sets the variance as sigma_tanLambda squared.
 *
 * - Assumes no correlations between helix parameters (diagonal covariance only).
 *
 * - Converts the helix covariance matrix to Cartesian coordinates:
 *     The transformation depends on:
 *       - Initial position
 *       - Initial momentum
 *       - Reference point
 *       - Magnetic field Bz
 *
 * - Returns the covariance matrix in Cartesian state representation.
 *
 * @param Bz Magnetic field component along the z axis.
 * @param Charge Particle charge hypothesis (+1 or -1).
 * @param sigma_d0 Optional uncertainty on d0 [mm] (default 0.05).
 * @param sigma_phi Optional uncertainty on phi [rad] (default 0.1).
 * @param omega_factor Scaling factor for curvature uncertainty (default 0.5).
 * @param z0_factor Scaling factor for z0 uncertainty (default 0.1).
 * @param sigma_tanLambda Optional uncertainty on tanLambda (default 0.1).
 *
 * @return Covariance matrix in Cartesian coordinates (6x6).
 */
TMatrixDSym GenfitTrack::ComputeInitialCovarianceMatrix(double Bz, int Charge, std::optional<double> sigma_d0,
                                                        std::optional<double> sigma_phi,
                                                        std::optional<double> omega_factor,
                                                        std::optional<double> z0_factor,
                                                        std::optional<double> sigma_tanLambda) {
  TMatrixDSym C_helix(5);
  C_helix.Zero();

  // d0 (in mm)
  double sd0 = sigma_d0.value_or(0.05);
  C_helix(0, 0) = sd0 * sd0;

  // phi (in rad)
  double sphi = sigma_phi.value_or(0.1);
  C_helix(1, 1) = sphi * sphi;

  // omega (in 1/mm)
  double pt = m_momInit.Perp();
  double omega = std::abs(ConversionUnits::a_lcio * Bz / pt);
  C_helix(2, 2) = std::pow(omega_factor.value_or(0.5) * omega, 2);

  // z0 relative to the track-parameter reference point, in mm. The scaled term vanishes when the
  // reference point lies on the track (e.g. the first hit), so it is floored to keep the covariance
  // positive definite.
  constexpr double kMinSigmaZ0 = 1e-3; // [mm]
  const double z0_scale = std::abs(m_posInit.Z() - m_vpReferencePoint.Z()) / dd4hep::mm;
  const double sigma_z0 = std::max(z0_factor.value_or(0.1) * z0_scale, kMinSigmaZ0);
  C_helix(3, 3) = sigma_z0 * sigma_z0;

  // tanLambda
  double stl = sigma_tanLambda.value_or(0.1);
  C_helix(4, 4) = stl * stl;

  TMatrixDSym covState =
      InitialCovarianceMatrixHelixToCartesian(C_helix, m_posInit, m_momInit, m_vpReferencePoint, Charge, Bz);

  return covState;
}

/**
 * @brief Computes initial track parameters (position, momentum, charge)
 *        from tracker hits using geometric fits.
 *
 * The method estimates a first approximation of the particle trajectory
 * starting from a collection of tracker hits. It reconstructs the transverse
 * motion via a circular fit in the XY plane and the longitudinal motion
 * via a linear fit in the (z, phi_local) plane.
 *
 * The procedure is as follows:
 *
 * - Extracts 3D hit positions (x, y, z) from the input track.
 *
 * - Projects the hits onto the XY plane and performs a circular fit:
 *     - Uses FastCircleFit to determine the best-fit circle.
 *     - Computes the point of closest approach (PCA) to a reference point.
 *     - Evaluates the tangent direction at the PCA.
 *
 * - Estimates the transverse momentum (pT):
 *     - Uses the relation:
 *           pT = |rho * 0.3 * Bz| / 1000
 *       where rho is the circle radius and Bz is the magnetic field.
 *     - Builds the initial momentum vector in the transverse plane
 *       from the tangent direction.
 *
 * - Reconstructs the longitudinal momentum (pZ):
 *     - Computes phi_local, the unwrapped angle around the fitted circle
 *       center along the hit order.
 *     - Performs a linear regression in the (z, phi_local) plane:
 *           z = b + a * phi_local
 *       (z is linear in phi_local for a helix).
 *     - Computes pZ from the exact helix relation
 *           dz/dphi = rho * pz / p.
 *     - Handles degenerate cases by setting pZ = 0.
 *
 * - Determines the z-coordinate at the PCA:
 *     - Evaluates the fitted z(phi_local) line at the PCA angle.
 *
 * - Estimates the particle charge:
 *     - Computes the curvature direction from the circle center.
 *     - Computes the Lorentz force direction: p × B.
 *     - Compares the two directions via a dot product:
 *           charge = +1 if aligned, -1 otherwise.
 *
 * - Fills and returns a HelperInitialization structure containing:
 *     - Position  : PCA in 3D (converted to cm)
 *     - Momentum  : (px, py, pz)
 *     - Charge    : +-1
 *
 * @param Bz Magnetic field component along the z-axis.
 *
 * @return HelperInitialization structure with initial track parameters.
 */
GenfitTrack::HelperInitialization GenfitTrack::ComputeInitialParameters(double Bz) {

  // The seed depends only on the selected hits, the field at the seed and the reference
  // point: identical inputs are re-used across fit attempts and particle hypotheses of the
  // same track (per-track cache; the cached value is the exact output of the computation).
  const TrackFitCache::SeedKey key{m_fitHits, Bz, m_vpReferencePoint.X(), m_vpReferencePoint.Y(),
                                   m_vpReferencePoint.Z()};
  if (m_cache) {
    for (const auto& entry : m_cache->seeds) {
      if (entry.first == key) {
        return entry.second;
      }
    }
  }

  HelperInitialization result = ComputeInitialParametersUncached(Bz);

  if (m_cache) {
    m_cache->seeds.emplace_back(key, result);
  }
  return result;
}

GenfitTrack::HelperInitialization GenfitTrack::ComputeInitialParametersUncached(double Bz) {

  Point2D_xy referencePoint_xy(m_vpReferencePoint.X() / dd4hep::mm, m_vpReferencePoint.Y() / dd4hep::mm); // in mm

  auto hits = m_edm4hepTrack.getTrackerHits();
  const size_t N = m_fitHits.size();

  if (N < 2) {
    throw std::runtime_error("Not enough hits to compute initial parameters");
  }

  // XY projection + z values (the z fit is done in the (z, phi_local) plane below)
  std::vector<Point2D_xy> points_xy;
  std::vector<double> z_val;
  points_xy.reserve(N);
  z_val.reserve(N);

  for (const std::size_t idx : m_fitHits) {

    const auto p = hits[idx].getPosition();

    // XY for circle fit
    points_xy.emplace_back(p.x, p.y);
    z_val.push_back(p.z);
  }

  // Circle fit (XY)
  FastCircleFit circle(points_xy);

  Point2D_xy closestPoint = circle.closestPointTo(referencePoint_xy);
  Point2D_xy tangent_xy = circle.tangentAtPCA(closestPoint, points_xy[1]);

  double rho = circle.rho();
  double init_pT = std::abs(rho * 0.3 * Bz) / 1000.0; // 1000 is used to convert from m to mm
  TVector3 init_mom(tangent_xy.x * init_pT, tangent_xy.y * init_pT, 0.0);

  // Z fit in the (z, phi_local): for a helix z is linear in the
  // unwrapped angle around the fitted circle center, z = z0 + (dz/dphi) * phi_local, so the
  // slope and intercept give z at the PCA and the longitudinal momentum.
  constexpr double kTwoPi = 6.28318530717958647692;
  const double x0c = circle.x0();
  const double y0c = circle.y0();
  std::vector<double> phi_val(N);
  for (std::size_t i = 0; i < N; ++i) {
    phi_val[i] = std::atan2(points_xy[i].y - y0c, points_xy[i].x - x0c);
    if (i > 0) {
      double dphi = phi_val[i] - phi_val[i - 1];
      dphi -= kTwoPi * std::round(dphi / kTwoPi); // wrap to (-pi, pi]
      phi_val[i] = phi_val[i - 1] + dphi;
    }
  }

  double sumP = 0., sumZ = 0., sumPP = 0., sumPZ = 0.;
  for (std::size_t i = 0; i < N; ++i) {
    sumP += phi_val[i];
    sumZ += z_val[i];
    sumPP += phi_val[i] * phi_val[i];
    sumPZ += phi_val[i] * z_val[i];
  }

  double slope = 0.0; // dz/dphi_local [mm/rad]
  double intercept = sumZ / N;
  const double denominator = N * sumPP - sumP * sumP;
  if (std::abs(denominator) > 1e-12 * (N * sumPP)) {
    slope = (N * sumPZ - sumP * sumZ) / denominator;
    intercept = (sumZ - slope * sumP) / N;
  }

  // z at the PCA: evaluate the fitted line at the PCA angle, mapped into the unwrapped frame
  double phi_pca = std::atan2(closestPoint.y - y0c, closestPoint.x - x0c);
  const double phi_mid = 0.5 * (phi_val.front() + phi_val.back());
  phi_pca += kTwoPi * std::round((phi_mid - phi_pca) / kTwoPi);
  double z_PCA = intercept + slope * phi_pca;

  // Exact helix relation, with the sense of travel along the hit order:
  // dz/dphi = sense * rho * pz / p  ->  pz = sense * slope * pT / sqrt(rho^2 - slope^2)
  const double radX = closestPoint.x - x0c;
  const double radY = closestPoint.y - y0c;
  const double sense = (tangent_xy.x * (-radY) + tangent_xy.y * radX >= 0.) ? 1. : -1.;
  const double transverse2 = rho * rho - slope * slope;
  if (transverse2 > 0.) {
    init_mom.SetZ(sense * slope * init_pT / std::sqrt(transverse2));
  } else {
    init_mom.SetZ(0.); // degenerate fit (|dz/dphi| >= rho): no longitudinal seed
  }

  // Charge
  TVector3 pos(closestPoint.x, closestPoint.y, z_PCA);
  TVector3 center(circle.x0(), circle.y0(), 0.0);

  const double rx = center.X() - pos.X();
  const double ry = center.Y() - pos.Y();

  const int sign = (init_mom.X() * ry - init_mom.Y() * rx > 0) ? -1 : 1;

  const int bzSign = (Bz > 0) ? 1 : -1;
  int charge = sign * bzSign;

  // Output
  HelperInitialization helper;
  helper.Position = TVector3(closestPoint.x * dd4hep::mm, closestPoint.y * dd4hep::mm, z_PCA * dd4hep::mm); // in cm

  helper.Momentum = init_mom;
  helper.Charge = charge;

  // Sanity check: nearly degenerate hit geometries (e.g. almost collinear points) make the
  // circle fit explode and give a meaningless position and momentum seed, which then
  // degenerates the initial covariance and the fit. Reseed along the hit direction with a
  // moderate momentum (same 1 GeV/c convention as InitializationType == 0) and at the mean
  // hit position. The bound is only meant to catch fit explosions (1e3+ GeV), not real
  // high-momentum tracks.
  const double pMag = init_mom.Mag();
  const bool badSeed = !std::isfinite(pMag) || pMag < 0.01 || pMag > 1e3 || !std::isfinite(helper.Position.Mag()) ||
                       helper.Position.Mag() > 1e5; // > 1000 m from the origin cannot be a real seed
  if (badSeed) {
    double meanX = 0., meanY = 0., meanZ = 0.;
    for (const std::size_t idx : m_fitHits) {
      const auto p = hits[idx].getPosition();
      meanX += p.x;
      meanY += p.y;
      meanZ += p.z;
    }
    meanX /= m_fitHits.size();
    meanY /= m_fitHits.size();
    meanZ /= m_fitHits.size();

    const auto& h0 = hits[m_fitHits.front()].getPosition();
    const auto& h1 = hits[m_fitHits.back()].getPosition();
    TVector3 dir(h1.x - h0.x, h1.y - h0.y, h1.z - h0.z);
    if (dir.Mag() < 1e-6) {
      dir.SetXYZ(1., 0., 0.);
    }

    helper.Position = TVector3(meanX * dd4hep::mm, meanY * dd4hep::mm, meanZ * dd4hep::mm);
    helper.Momentum = dir.Unit(); // 1 GeV/c seed
  }

  return helper;
}

/**
 * @brief Builds the genfit::Track from the initialized state and tracker hits.
 *
 * The method performs the following steps:
 *
 * - Sets the signed particle hypothesis (including charge sign handling
 *   for leptons).
 * - Deletes any existing genfit::TrackRep and genfit::Track objects.
 * - Constructs the 6D state vector (x, y, z, px, py, pz) from m_posInit
 *   and m_momInit.
 * - Creates a genfit::RKTrackRep using the signed particle hypothesis.
 * - Instantiates a genfit::Track with the track representation,
 *   state vector, and m_covInit.
 * - Iterates over the hits selected for the fit (m_fitHits: all the hits, or the first arc of a looper).
 * - Converts each hit into a Genfit measurement:
 *     - edm4hep::TrackerHitPlane: planar measurement (detID = 0)
 *     - edm4hep::SenseWireHit: wire measurement   (detID = 1)
 * - Wraps each measurement into a genfit::TrackPoint
 *   and inserts it into the genfit::Track.
 * - Terminates execution if an unsupported hit type is found.
 *
 * @param particle_hypotesis PDG code hypothesis (sign adjusted internally).
 * @param debug_lvl Debug verbosity level passed to measurement constructors.
 *
 * @note InitializeTrack() must be called before this method.
 */
void GenfitTrack::CreateGenFitTrack(int particle_hypotesis, int debug_lvl) {

  m_genfitTrack.reset();

  m_signedParticleHypothesis = particle_hypotesis;
  if (particle_hypotesis == 11 || particle_hypotesis == 13) {
    m_signedParticleHypothesis = -m_chargeHypothesis * particle_hypotesis;

  } else {
    m_signedParticleHypothesis = m_chargeHypothesis * particle_hypotesis;
  }

  // Create stateVec
  TVectorD stateVec(6);

  stateVec[0] = m_posInit.X();
  stateVec[1] = m_posInit.Y();
  stateVec[2] = m_posInit.Z();

  stateVec[3] = m_momInit.X();
  stateVec[4] = m_momInit.Y();
  stateVec[5] = m_momInit.Z();

  auto trackRep = std::make_unique<genfit::RKTrackRep>(m_signedParticleHypothesis);
  m_genfitTrack = std::make_unique<genfit::Track>(trackRep.get(), stateVec, m_covInit);
  trackRep.release();

  auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();

  // Check whether the wire tracker has sectors.
  // m_dcDecoder is nullptr when the geometry has no wire tracker; in that case
  // there are no SenseWireHits either, so has_sectors is never actually consumed.
  const bool has_sectors = m_dcDecoder && (m_dcDecoder->fieldDescription().find("sector") != std::string::npos);

  int detID(-1);

  // The following code creates genfit::TrackPoints from EDM4hep hits and adds them to the genfit::Track.
  // Only the hits selected for the fit (m_fitHits: all of them, or the first arc of a looper) are used;
  // the hitId is the 1-based position of the hit in m_edm4hepTrack, so that Fit() can map it back.
  //
  // GENFIT manages ownership of the entire object chain:
  // - A genfit::Track owns its TrackPoints
  // - Deleting a Track automatically deletes all associated TrackPoints
  // - Deleting a TrackPoint also releases the underlying measurement object used to create it
  //
  // Measurement masters are cached per hit (per-track cache, shared only between the fit
  // attempts and hypotheses of this track) and cloned for every TrackPoint: the clone is
  // owned by the TrackPoint, the master by the cache. The clone carries the same raw
  // measurement and the same hitId as a directly constructed measurement.
  auto buildMeasurement = [&](std::size_t idx, int hitId) -> genfit::AbsMeasurement* {
    auto hit = hits_for_genfit[idx];

    if (hit.isA<edm4hep::TrackerHitPlane>()) {
      detID = 0;
      auto planar_hit = hit.as<edm4hep::TrackerHitPlane>();
      GenfitInterface::PlanarMeasurement measurement(planar_hit, detID, hitId, debug_lvl);
      return measurement.getGenFit();
    } else if (hit.isA<edm4hep::SenseWireHit>()) {
      detID = 1;
      auto wire_hit = hit.as<edm4hep::SenseWireHit>();
      GenfitInterface::WireMeasurement measurement(wire_hit, m_wireInfo, m_dcDecoder, has_sectors, detID, hitId,
                                                   debug_lvl);
      return measurement.getGenFit();
    }
    throw std::runtime_error("InitializeTrack: Unknown hit type encountered - Hit will be skipped.");
  };

  if (m_cache) {
    auto& masters = m_cache->measurements;
    if (masters.size() < hits_for_genfit.size()) {
      masters.resize(hits_for_genfit.size());
    }
    for (const std::size_t idx : m_fitHits) {
      const int hitId = static_cast<int>(idx) + 1;
      if (!masters[idx]) {
        masters[idx].reset(buildMeasurement(idx, hitId));
      }
      m_genfitTrack->insertPoint(new genfit::TrackPoint(masters[idx]->clone(), m_genfitTrack.get()));
    }
  } else {
    for (const std::size_t idx : m_fitHits) {
      const int hitId = static_cast<int>(idx) + 1;
      m_genfitTrack->insertPoint(new genfit::TrackPoint(buildMeasurement(idx, hitId), m_genfitTrack.get()));
    }
  }
}

/**
 * @brief Fit the Genfit track using a chosen track fitting algorithm.
 *
 * This method fits the `GenfitTrack` object using the specified fitting strategy
 * (Deterministic Annealing Filter, standard Kalman, or Kalman with reference track),
 * and updates the associated `edm4hep::Track` track states at key locations:
 * the Interaction Point (IP), the first hit, and the last hit.
 *
 * If DAF is used and `FilterHits` is enabled, track hits are filtered after fitting
 * based on their measurement weights, retaining only accepted hits in `m_fittedHits`.
 *
 * The method also computes helix parameters (d0, z0, phi, omega, tanLambda) for each
 * track state based on the fitted position and momentum, taking into account the
 * assumed charge hypothesis and magnetic field.
 *
 * @param fittedHits Output collection of fitted tracker hits after filtering (it is not updated if no filtering is
 * applied).
 * @param FitterType Fitting strategy to use (https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf):
 *        - "DAF"        : Deterministic Annealing Filter
 *        - "KALMAN"     : Standard Kalman filter
 *        - "KALMAN_REF" : Kalman filter with reference track
 * @param debug_lvl   Gaudi-compatible verbosity level for debug output
 *                       - 0 or >= 3 = silent
 *                       - 2 = fit results, track states, and fit failures
 *                       - 1 = level 2 output + detailed fitter information
 * @param Beta_init   Initial annealing parameter (only for DAF, default = 100)
 * @param Beta_final  Final annealing parameter (only for DAF, default = 0.1)
 * @param Beta_steps  Number of annealing steps (only for DAF, default = 10)
 * @param FilterHits  If true and DAF is used, applies hit filtering after fitting (default = true)
 *
 * @return true if the fit was successful, false otherwise.
 *
 * @note If any exception occurs during fitting or state extrapolation, the function
 *       returns false and does not update the track.
 */
bool GenfitTrack::Fit(std::string FitterType = "DAF", int debug_lvl = 0, std::optional<double> Beta_init = 100.,
                      std::optional<double> Beta_final = 0.1, std::optional<int> Beta_steps = 10,
                      std::optional<bool> FilterHits = true, bool NoMaterial = false) {

  const bool showFitDiagnostics = (debug_lvl == 1 || debug_lvl == 2);
  const bool showFitterDiagnostics = (debug_lvl == 1);

  // The navigator is process-global. GeoNavigatorStateGuard saves a backup of the navigator and
  // restores it when GeoNavigatorStateGuard is destroyed.
  // This step initalizes the navigator with the fit seed (position and momentum).
  GeoNavigatorStateGuard geoStateGuard;
  TGeoNavigator* navigator = gGeoManager->GetCurrentNavigator();
  navigator->ResetState();
  navigator->CdTop();
  navigator->FindNode(m_posInit.X(), m_posInit.Y(), m_posInit.Z()); // seed position, in cm

  // Initialize the genfit fitter
  std::unique_ptr<genfit::AbsKalmanFitter> genfitFitter = nullptr;

  if (FitterType == "DAF") {

    auto daf = std::make_unique<genfit::DAF>(true, 1e-3, 1e-3);
    daf->setAnnealingScheme(Beta_init.value(), Beta_final.value(), Beta_steps.value());
    daf->setProbCut(1e-5);
    daf->setConvergenceDeltaWeight(1e-2);

    genfitFitter = std::move(daf);

  } else if (FitterType == "KALMAN") {

    genfitFitter = std::make_unique<genfit::KalmanFitter>();

  } else if (FitterType == "KALMAN_REF") {

    genfitFitter = std::make_unique<genfit::KalmanFitterRefTrack>();

  } else {

    throw std::invalid_argument("Unknown fit method: " + FitterType);
  }

  genfitFitter->setDebugLvl(showFitterDiagnostics ? 1 : 0);

  // Recovery option: turn the material effects off for this fit. The step-size oscillations of
  // the energy-loss machinery at steep material boundaries can generate non positive-definite
  // transport noise and kill the Kalman recursion; for a first-arc fit the material budget is
  // tiny, so a fit without material effects is an acceptable last resort. Restored on exit.
  struct NoMaterialGuard {
    genfit::MaterialEffects* effects;
    bool active;
    explicit NoMaterialGuard(bool on) : effects(genfit::MaterialEffects::getInstance()), active(on) {
      if (active) {
        effects->setNoEffects(true);
      }
    }
    ~NoMaterialGuard() {
      if (active) {
        effects->setNoEffects(false);
      }
    }
  } noMaterialGuard(NoMaterial);

  // Process track
  genfit::Track genfitTrack = *m_genfitTrack;
  genfit::AbsTrackRep* trackRep = genfitTrack.getTrackRep(0);

  try {

    const auto tFitCoreStart = std::chrono::steady_clock::now();
    genfitFitter->processTrackWithRep(&genfitTrack, trackRep);
    m_fitCoreSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - tFitCoreStart).count();

  } catch (const std::exception& e) {

    if (showFitDiagnostics) {
      std::cerr << "Exception during track fitting: " << e.what() << std::endl;
    }
    m_edm4hepTrack.setChi2(-1);
    m_edm4hepTrack.setNdf(-1);

    m_trackWithFit.setChi2(-1);
    m_trackWithFit.setNdf(-1);

    return false;
  }

  const auto tStateStart = std::chrono::steady_clock::now();
  struct StateTimerGuard {
    double& acc;
    std::chrono::steady_clock::time_point t0;
    ~StateTimerGuard() { acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
  } stateTimerGuard{m_stateExtractionSeconds, tStateStart};

  if (genfitFitter->isTrackFitted(&genfitTrack, trackRep)) {

    // The fitted hit data is collected in m_pendingFittedHits (class member): the caller
    // decides the output collection via CreateFilteredHits().

    // Hit Filtering based on measurement weights (if FitterType == "DAF" and FilterHits == true)
    if (FilterHits.value() && FitterType == "DAF") {

      int numPoints = genfitTrack.getNumPoints();
      for (int idx = 0; idx < numPoints; idx++) {

        // Retrieve the genfit::TrackPoint for this index and get its associated fitter information.
        // Each TrackPoint can have one or more assigned measurements:
        // - Planar measurements: typically 1 measurement
        // - Drift chamber hits: typically 2 measurements (left and right positions)
        //
        // The loop iterates over all measurements to check their assigned weights and determine
        // whether the hit should be accepted or discarded.
        auto point_track = genfitTrack.getPoint(idx);
        auto fitterInfo = point_track->getFitterInfo(trackRep);
        genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(fitterInfo);

        unsigned int nMeas = kfi->getNumMeasurements();
        bool isAccepted = false;

        for (unsigned int j = 0; j < nMeas; j++) {
          genfit::MeasurementOnPlane* mop = kfi->getMeasurementOnPlane(j);
          double weight = mop->getWeight();

          // Weight threshold: 1e-5
          // - Genfit assigns 1e-10 to measurements it discards.
          // - Exception: when left and right hits are very close (within detector resolution),
          //   both may get weight 0.5. In this case, we resolve the ambiguity by selecting the left hit.
          if (weight > 1e-5) {
            isAccepted = true;
            break;
          }
        }

        if (isAccepted) {
          // Retrieving the fitted state averages the forward and backward Kalman
          // states. This can fail for an otherwise fitted track when one point has
          // an ill-conditioned covariance matrix. Do not let such a point abort the
          // full event; omit it from the filtered-hit output instead.
          genfit::MeasuredStateOnPlane measState;
          try {
            measState = kfi->getFittedState();
          } catch (const std::exception& e) {
            if (showFitDiagnostics) {
              std::cerr << "Exception retrieving fitted hit state: " << e.what() << std::endl;
            }
            continue;
          }

          // Extract the 3D position of the fitted state
          TVector3 pos = measState.getPos();

          // Extract the measurement plane's orientation vectors
          auto planeMeas = measState.getPlane();
          auto U = planeMeas->getU(); // Plane u-direction
          auto V = planeMeas->getV(); // Plane v-direction
          auto O = planeMeas->getO(); // Plane origin

          // Extract covariance matrix in Genfit's 6D (pos+momentum) space
          TMatrixDSym covMatrix_cm_gev = measState.get6DCov();

          // Scale the covariance matrix from cm/GeV to desired mm/GeV
          TMatrixDSym covMatrix = covMatrix_cm_gev;
          for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
              double scale = 1.0;
              bool pos_i = (i < 3);
              bool pos_j = (j < 3);
              if (pos_i && pos_j)
                scale = 100.0; // position-position
              else if (pos_i || pos_j)
                scale = 10.0; // position-momentum
              covMatrix(i, j) *= scale;
            }
          }

          // Extract only the 3x3 position covariance submatrix
          TMatrixDSym covXYZ(3);
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              covXYZ(i, j) = covMatrix(i, j);
            }
          }

          // Convert plane axes to TVector3 for error computation
          TVector3 u(U.X(), U.Y(), U.Z());
          TVector3 v(V.X(), V.Y(), V.Z());

          // Convert plane axes to TVectorD for similarity calculation with covariance
          TVectorD uVec(3);
          uVec[0] = U.X();
          uVec[1] = U.Y();
          uVec[2] = U.Z();
          TVectorD vVec(3);
          vVec[0] = V.X();
          vVec[1] = V.Y();
          vVec[2] = V.Z();

          // Compute variances along plane directions using covariance similarity
          double var_u = covXYZ.Similarity(uVec);
          double var_v = covXYZ.Similarity(vVec);

          // Convert variance to standard deviation (errors)
          double err_u = std::sqrt(var_u);
          double err_v = std::sqrt(var_v);

          // Identify the digi hit of this point: the measurement hitId is the 1-based
          // position of the hit in m_edm4hepTrack (see CreateGenFitTrack)
          const int hitId = point_track->getRawMeasurement(0)->getHitId();
          const uint64_t digiCellID = m_edm4hepTrack.getTrackerHits(hitId - 1).getCellID();

          m_pendingFittedHits.push_back({pos, covXYZ, err_u, err_v, u, v, digiCellID, hitId - 1});
        }
      }
    }

    // getFittedState() averages the forward/backward Kalman states and can throw
    // a std::exception (e.g. an ill-conditioned covariance)
    genfit::MeasuredStateOnPlane fittedState;

    // TrackState at the first hit. The reference point is the corresponding
    // digi-hit position; UpdateTrackState moves the fitted state to the
    // transverse PCA with respect to that pivot.
    edm4hep::TrackState trackStateFirstHit;
    try {
      fittedState = genfitTrack.getFittedState();
      trackStateFirstHit = UpdateTrackState(fittedState, m_firstHitReferencePoint, edm4hep::TrackState::AtFirstHit);
    } catch (const std::exception& e) {
      if (showFitDiagnostics) {
        std::cerr << "Exception creating TrackState at first hit: " << e.what() << std::endl;
      }
      m_edm4hepTrack.setChi2(-1);
      m_edm4hepTrack.setNdf(-1);
      m_trackWithFit.setChi2(-1);
      m_trackWithFit.setNdf(-1);

      return false;
    }

    // TrackState at the last hit, using the last digi-hit position as pivot.
    edm4hep::TrackState trackStateLastHit;
    try {
      fittedState = genfitTrack.getFittedState(genfitTrack.getNumPoints() - 1);
      trackStateLastHit = UpdateTrackState(fittedState, m_lastHitReferencePoint, edm4hep::TrackState::AtLastHit);
    } catch (const std::exception& e) {
      if (showFitDiagnostics) {
        std::cerr << "Exception creating TrackState at last hit: " << e.what() << std::endl;
      }
      m_edm4hepTrack.setChi2(-1);
      m_edm4hepTrack.setNdf(-1);
      m_trackWithFit.setChi2(-1);
      m_trackWithFit.setNdf(-1);

      return false;
    }

    // TrackState at the nominal IP. m_vpReferencePoint is used as pivot point
    genfit::TrackPoint* tp = genfitTrack.getPointWithFitterInfo(0);
    auto* fi = static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(trackRep));
    edm4hep::TrackState trackStateIP;

    try {
      fittedState = fi->getFittedState(true);
      trackStateIP = UpdateTrackState(
          fittedState, m_vpReferencePoint,
          edm4hep::TrackState::AtIP); // N.B. trackState does not carry meaningful information for displaced tracks
    } catch (const std::exception& e) {
      // The point-0 state can fail on ill-conditioned measurement averages; the globally
      // smoothed state (already used above) is a valid fallback instead of losing the track
      try {
        fittedState = genfitTrack.getFittedState();
        trackStateIP = UpdateTrackState(fittedState, m_vpReferencePoint, edm4hep::TrackState::AtIP);
      } catch (const std::exception&) {
        if (showFitDiagnostics) {
          std::cerr << "Exception creating TrackState at IP: " << e.what() << std::endl;
        }

        m_edm4hepTrack.setChi2(-1);
        m_edm4hepTrack.setNdf(-1);
        m_trackWithFit.setChi2(-1);
        m_trackWithFit.setNdf(-1);

        return false;
      }
    }

    // The filtered hits are committed by CreateFilteredHits(): they must be created in the
    // collection where they will live (an EDM4hep object cannot be shared between two
    // collections), and only the caller knows that collection.

    if (showFitDiagnostics) {

      std::cout << "GenfitTrack    DEBUG : TrackState at IP: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateIP.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateIP.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateIP.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateIP.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateIP.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateIP.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateIP.referencePoint.x << ", "
                << trackStateIP.referencePoint.y << ", " << trackStateIP.referencePoint.z << ") mm" << std::endl;

      std::cout << "\nGenfitTrack   DEBUG : TrackState at First Hit: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateFirstHit.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateFirstHit.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateFirstHit.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateFirstHit.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateFirstHit.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateFirstHit.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateFirstHit.referencePoint.x << ", "
                << trackStateFirstHit.referencePoint.y << ", " << trackStateFirstHit.referencePoint.z << ") mm"
                << std::endl;

      std::cout << "\nGenfitTrack   DEBUG : TrackState at Last Hit: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateLastHit.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateLastHit.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateLastHit.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateLastHit.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateLastHit.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateLastHit.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateLastHit.referencePoint.x << ", "
                << trackStateLastHit.referencePoint.y << ", " << trackStateLastHit.referencePoint.z << ") mm\n"
                << std::endl;
    }

    m_edm4hepTrack.addToTrackStates(trackStateIP);
    m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
    m_edm4hepTrack.addToTrackStates(trackStateLastHit);

    m_trackWithFit.addToTrackStates(trackStateIP);
    m_trackWithFit.addToTrackStates(trackStateFirstHit);
    m_trackWithFit.addToTrackStates(trackStateLastHit);

    m_edm4hepTrack.setChi2(genfitTrack.getFitStatus()->getChi2());
    m_edm4hepTrack.setNdf(genfitTrack.getFitStatus()->getNdf());

    m_trackWithFit.setChi2(genfitTrack.getFitStatus()->getChi2());
    m_trackWithFit.setNdf(genfitTrack.getFitStatus()->getNdf());

    return true;

  } else {
    if (showFitDiagnostics) {
      std::cerr << "Fit not converged: chi2/ndf = " << genfitTrack.getFitStatus(trackRep)->getChi2() << "/"
                << genfitTrack.getFitStatus(trackRep)->getNdf() << std::endl;
    }
    m_edm4hepTrack.setChi2(-1);
    m_edm4hepTrack.setNdf(-1);

    m_trackWithFit.setChi2(-1);
    m_trackWithFit.setNdf(-1);

    return false;
  }
}

/**
 * @brief Create the filtered hits collected by Fit() in the given output collection.
 *
 * The hits are created with the same fields Fit() computed (position, covariance, plane
 * axes, digi cellID and the digi hit index in getQuality()) and are related to the
 * filtered-hit track, exactly as when Fit() created them directly. The hits must be created
 * in the collection where they will live (an EDM4hep object cannot be shared between two
 * collections), which is why the caller passes it here. Call at most once per fitted track:
 * the pending hits are consumed.
 */
void GenfitTrack::CreateFilteredHits(edm4hep::TrackerHitPlaneCollection& out) {

  for (const auto& pendingHit : m_pendingFittedHits) {
    auto hit3D = out.create();
    hit3D.setPosition(edm4hep::Vector3d(pendingHit.position.X() / dd4hep::mm, pendingHit.position.Y() / dd4hep::mm,
                                        pendingHit.position.Z() / dd4hep::mm));
    hit3D.setCovMatrix({
        static_cast<float>(pendingHit.covariance(0, 0)), // xx
        static_cast<float>(pendingHit.covariance(1, 0)), // yx
        static_cast<float>(pendingHit.covariance(1, 1)), // yy
        static_cast<float>(pendingHit.covariance(2, 0)), // zx
        static_cast<float>(pendingHit.covariance(2, 1)), // zy
        static_cast<float>(pendingHit.covariance(2, 2))  // zz
    });
    hit3D.setDu(pendingHit.du);
    hit3D.setDv(pendingHit.dv);
    // plane axes as (theta, phi) of the unit vectors, following the edm4hep TrackerHitPlane convention
    hit3D.setU({static_cast<float>(pendingHit.u.Theta()), static_cast<float>(pendingHit.u.Phi())});
    hit3D.setV({static_cast<float>(pendingHit.v.Theta()), static_cast<float>(pendingHit.v.Phi())});
    hit3D.setCellID(pendingHit.cellID);
    hit3D.setQuality(pendingHit.hitIndex);

    m_trackWithFit.addToTrackerHits(hit3D);
  }
  m_pendingFittedHits.clear();
}

/**
 * @brief Propagation of the covariance matrix from helix-parameter representation to Cartesian coordinate
 * representation
 *
 * @param C_helix Covariance matrix in native EDM4hep helix units: (mm, rad, 1/mm, mm, dimensionless)
 * @param positionCm Nominal Cartesian position in GENFIT cm.
 * @param momentumGeV Nominal Cartesian momentum in GeV.
 * @param referencePointCm Reference point in GENFIT cm.
 * @param charge Charge hypothesis.
 * @param magneticFieldTesla Longitudinal magnetic field in T.
 *
 * @return Covariance matrix in cartesian-basis (cm and GeV units)
 */
TMatrixDSym GenfitTrack::InitialCovarianceMatrixHelixToCartesian(const TMatrixDSym& C_helix, const TVector3& positionCm,
                                                                 const TVector3& momentumGeV,
                                                                 const TVector3& referencePointCm, int charge,
                                                                 double magneticFieldTesla) {

  const double px = momentumGeV.X();
  const double py = momentumGeV.Y();
  const double pz = momentumGeV.Z();
  const double pt = momentumGeV.Perp();
  const double phi = momentumGeV.Phi();
  const double xMm = positionCm.X() / dd4hep::mm;
  const double yMm = positionCm.Y() / dd4hep::mm;
  const double referenceXMm = referencePointCm.X() / dd4hep::mm;
  const double referenceYMm = referencePointCm.Y() / dd4hep::mm;
  const double d0Mm = -(xMm - referenceXMm) * std::sin(phi) + (yMm - referenceYMm) * std::cos(phi);
  const double omegaPerMm =
      (magneticFieldTesla * charge < 0. ? -1. : 1.) * std::abs(ConversionUnits::a_lcio * magneticFieldTesla / pt);

  // --- Jacobian (6x5) ---
  TMatrixD J(6, 5);
  J.Zero();

  // These definitions are taken from
  // https://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf
  //
  // d0 = -(x_PCA - Xr)*sin(phi0) + (y_PCA - Yr)*cos(phi0)
  // => x_PCA = Xr - d0*sin(phi0), y_PCA = Yr + d0*cos(phi0)   [all in mm]

  J(0, 0) = -sin(phi);        // dx / dd0
  J(0, 1) = -d0Mm * cos(phi); // dx / dphi0
  J(0, 2) = 0.0;              // dx / domega
  J(0, 3) = 0.0;              // dx / dz0
  J(0, 4) = 0.0;              // dx / dtanLambda

  J(1, 0) = cos(phi);         // dy / dd0
  J(1, 1) = -d0Mm * sin(phi); // dy / dphi0
  J(1, 2) = 0.0;              // dy / domega
  J(1, 3) = 0.0;              // dy / dz0
  J(1, 4) = 0.0;              // dy / dtanLambda

  // z = z_PCA = z0 + P^r_z   [mm]
  J(2, 0) = 0.0; // dz / dd0
  J(2, 1) = 0.0; // dz / dphi0
  J(2, 2) = 0.0; // dz / domega
  J(2, 3) = 1.0; // dz / dz0
  J(2, 4) = 0.0; // dz / dtanLambda

  // px = pt * cos(phi)
  J(3, 0) = 0.0;              // dpx / dd0
  J(3, 1) = -py;              // dpx / dphi0
  J(3, 2) = -px / omegaPerMm; // dpx / domega
  J(3, 3) = 0.0;              // dpx / dz0
  J(3, 4) = 0.0;              // dpx / dtanLambda

  // py = pt * sin(phi)
  J(4, 0) = 0.0;              // dpy / dd0
  J(4, 1) = px;               // dpy / dphi0
  J(4, 2) = -py / omegaPerMm; // dpy / domega
  J(4, 3) = 0.0;              // dpy / dz0
  J(4, 4) = 0.0;              // dpy / dtanLambda

  // pz = pt * tanLambda
  J(5, 0) = 0.0;              // dpz / dd0
  J(5, 1) = 0.0;              // dpz / dphi0
  J(5, 2) = -pz / omegaPerMm; // dpz / domega
  J(5, 3) = 0.0;              // dpz / dz0
  J(5, 4) = pt;               // dpz / dtanLambda

  // --- Compute C_cart = J * C_helix * J^T   (positions in mm, momenta in GeV) ---
  TMatrixD Jt(TMatrixD::kTransposed, J);
  TMatrixD tmp = J * C_helix * Jt;

  // Convert position-related blocks back from mm to cm, to match Position_cm's units
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      tmp(i, j) *= dd4hep::mm * dd4hep::mm;
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 3; j < 6; ++j) {
      tmp(i, j) *= dd4hep::mm;
      tmp(j, i) *= dd4hep::mm;
    }
  }

  TMatrixDSym C_cart(6);
  C_cart.Zero();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j <= i; ++j) {
      C_cart(i, j) = tmp(i, j);
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      C_cart(i, j) = C_cart(j, i);
    }
  }

  return C_cart;
}

/**
 * @brief Propagation of the covariance matrix from cartesian-parameter representation to Helix coordinate
 * representation
 *
 * @param C_cartesian    Covariance Matrix in cartesian-basis (cm and GeV units)
 * @param Position_cm    Inizial Position in cm
 * @param Momentum_gev   Initial Momentum in gev/c
 * @param RefPoint_cm    Reference position (e.g. IP)
 * @param Charge         Charge hypothesis
 * @param Bz             Bz
 *
 * @return Covariance matrix in helix-basis (mm and GeV units)
 */
TMatrixDSym GenfitTrack::CovarianceMatrixCartesianToHelix(const TMatrixDSym& C_cartesian, // 6x6
                                                          TVector3 Position_cm, TVector3 Momentum_gev,
                                                          TVector3 RefPoint_cm, int Charge, double Bz) {

  double px = Momentum_gev.X();
  double py = Momentum_gev.Y();
  double pz = Momentum_gev.Z();

  double x_PCA_mm = Position_cm.X() / dd4hep::mm;
  double y_PCA_mm = Position_cm.Y() / dd4hep::mm;
  double RefX_mm = RefPoint_cm.X() / dd4hep::mm;
  double RefY_mm = RefPoint_cm.Y() / dd4hep::mm;

  double pt = Momentum_gev.Perp();
  double pt2 = std::pow(pt, 2);

  double phi0 = std::atan2(py, px);
  double tanLambda = pz / pt;
  double omega = (Bz * Charge < 0 ? -1.0 : 1.0) * std::abs(ConversionUnits::a_lcio * Bz / pt); // in 1/mm

  // Convert covariance matrix from cm/GeV to mm/GeV
  TMatrixDSym C_cart(C_cartesian);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C_cart(i, j) *= 1 / (dd4hep::mm * dd4hep::mm);
    }
  }
  for (int i = 0; i < 3; ++i) {
    for (int j = 3; j < 6; ++j) {
      C_cart(i, j) *= 1 / dd4hep::mm;
      C_cart(j, i) *= 1 / dd4hep::mm;
    }
  }

  // --- Jacobian (5x6) ---
  TMatrixD J(5, 6);
  J.Zero();

  // These definitions are taken from
  // https://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf

  // d0 = -(x - Xr)*sin(phi0) + (y - Yr)*cos(phi0)
  J(0, 0) = -sin(phi0); // dd0 / dx
  J(0, 1) = cos(phi0);  // dd0 / dy
  J(0, 2) = 0.0;        // dd0 / dz

  // chain rule via phi0 = atan2(py, px)
  double dd0_dphi = (RefX_mm - x_PCA_mm) * cos(phi0) + (RefY_mm - y_PCA_mm) * sin(phi0);

  J(0, 3) = dd0_dphi * (-py / pt2); // dd0 / dpx
  J(0, 4) = dd0_dphi * (px / pt2);  // dd0 / dpy
  J(0, 5) = 0.0;                    // dd0 / dpz

  // phi0 = atan2(py, px);
  J(1, 0) = 0.0;       // dphi0 / dx
  J(1, 1) = 0.0;       // dphi0 / dy
  J(1, 2) = 0.0;       // dphi0 / dz
  J(1, 3) = -py / pt2; // dphi0 / dpx
  J(1, 4) = px / pt2;  // dphi0 / dpy
  J(1, 5) = 0.0;       // dphi0 / dpz

  // omega = q * a * Bz / pt, with pt = sqrt(px^2 + py^2)
  J(2, 0) = 0.0;               // domega / dx_PCA
  J(2, 1) = 0.0;               // domega / dy_PCA
  J(2, 2) = 0.0;               // domega / dz_PCA
  J(2, 3) = -omega * px / pt2; // domega / dpx
  J(2, 4) = -omega * py / pt2; // domega / dpy
  J(2, 5) = 0.0;               // domega / dpz

  // z0 = P^0_z - P^r_z
  J(3, 0) = 0.0; // dz0 / dx_PCA
  J(3, 1) = 0.0; // dz0 / dy_PCA
  J(3, 2) = 1.0; // dz0 / dz_PCA
  J(3, 3) = 0.0; // dz0 / dpx
  J(3, 4) = 0.0; // dz0 / dpy
  J(3, 5) = 0.0; // dz0 / dpz

  // tanLambda = pz / pt
  J(4, 0) = 0.0;                   // dtanLambda / dx_PCA
  J(4, 1) = 0.0;                   // dtanLambda / dy_PCA
  J(4, 2) = 0.0;                   // dtanLambda / dz_PCA
  J(4, 3) = -tanLambda * px / pt2; // dtanLambda / dpx
  J(4, 4) = -tanLambda * py / pt2; // dtanLambda / dpy
  J(4, 5) = 1.0 / pt;              // dtanLambda / dpz

  // --- Compute C_helix = J * C_cart * J^T ---
  TMatrixD Jt(TMatrixD::kTransposed, J);
  TMatrixD tmp = J * C_cart * Jt;

  TMatrixDSym C_helix(5);
  C_helix.Zero();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      C_helix(i, j) = tmp(i, j);
    }
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      C_helix(i, j) = C_helix(j, i);
    }
  }

  return C_helix;
}

/**
 * @brief Update an edm4hep::TrackState using the fitted Genfit state parameters.
 *
 * This method extracts position, momentum, and covariance information from a
 * `genfit::MeasuredStateOnPlane` object and converts them into the corresponding
 * `edm4hep::TrackState` representation, including helix parameters.
 *
 * The procedure includes:
 *   - Extrapolation to the transverse point of closest approach to the line
 *     parallel to z through the requested reference point
 *   - Extraction of position, momentum, and covariance at that PCA
 *   - Retrieval of the local magnetic field (Bz) at the track position
 *   - Derivation of helix parameters:
 *       - d0          : transverse impact parameter (mm)
 *       - z0          : longitudinal impact parameter (mm)
 *       - phi         : azimuthal angle of momentum (rad)
 *       - omega       : curvature (1/mm), signed according to charge hypothesis
 *       - tanLambda   : dip angle (pz / pt)
 *   - Assignment of the reference point and track state location
 *
 * @param measuredState Genfit measured state containing fitted position and momentum
 * @param referencePointCm Reference point position [cm]
 * @param location Location identifier for the track state (e.g. AtIP, AtFirstHit, AtLastHit)
 *
 * @note The magnetic field is retrieved at the current position and converted to Tesla.
 * @note The curvature sign (omega) depends on the assumed particle charge hypothesis.
 */
edm4hep::TrackState GenfitTrack::UpdateTrackState(genfit::MeasuredStateOnPlane measuredState, TVector3 referencePointCm,
                                                  int location) {

  edm4hep::TrackState Edm4hepTrackState;

  // GENFIT extrapolates the fitted state to its transverse PCA relative to
  // the line parallel to z through referencePointCm. For AtIP, this is the beam
  // axis; for AtFirstHit and AtLastHit, it passes through the corresponding
  // digi-hit position.
  measuredState.extrapolateToLine(referencePointCm, TVector3(0., 0., 1.));

  TVector3 positionCm, momentumGeV;
  TMatrixDSym covarianceCmGeV(6);

  measuredState.getPosMomCov(positionCm, momentumGeV, covarianceCmGeV);

  double xCm = positionCm.X();
  double yCm = positionCm.Y();
  double zCm = positionCm.Z();
  double pzGeV = momentumGeV.Z();
  double ptGeV = momentumGeV.Perp();
  double phi = momentumGeV.Phi();

  double d0Mm =
      (-(xCm - referencePointCm.X()) * std::sin(phi) + (yCm - referencePointCm.Y()) * std::cos(phi)) / dd4hep::mm;
  double z0Mm = (zCm - referencePointCm.Z()) / dd4hep::mm;

  double bzTesla = m_fieldMap->getBz(positionCm) / (dd4hep::tesla / dd4hep::kilogauss);

  double tanLambda = pzGeV / ptGeV;
  double omegaPerMm = std::abs(ConversionUnits::a_lcio * bzTesla / ptGeV);
  if (m_chargeHypothesis < 0)
    omegaPerMm = -omegaPerMm;

  Edm4hepTrackState.D0 = d0Mm;
  Edm4hepTrackState.Z0 = z0Mm;
  Edm4hepTrackState.phi = phi;
  Edm4hepTrackState.omega = omegaPerMm;
  Edm4hepTrackState.tanLambda = tanLambda;
  Edm4hepTrackState.time = 0.;

  Edm4hepTrackState.referencePoint = edm4hep::Vector3f(
      referencePointCm.X() / dd4hep::mm, referencePointCm.Y() / dd4hep::mm, referencePointCm.Z() / dd4hep::mm);
  Edm4hepTrackState.location = location;

  TMatrixDSym covarianceHelix = CovarianceMatrixCartesianToHelix(covarianceCmGeV, positionCm, momentumGeV,
                                                                 referencePointCm, m_chargeHypothesis, bzTesla);

  // Conversion from TMatrixDSym(5x5) to lower-triangular packed format used in edm4hep::TrackState
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      Edm4hepTrackState.setCovMatrix(static_cast<float>(covarianceHelix(i, j)), static_cast<edm4hep::TrackParams>(i),
                                     static_cast<edm4hep::TrackParams>(j));
    }
  }

  return Edm4hepTrackState;
}

} // namespace GenfitInterface
