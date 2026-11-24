/*
 * Copyright (c) 2014-2024 Key4hep-Project.
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
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include <TGeoManager.h>
#include <TVector3.h>

#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>

#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/MaterialManager.h>
#include <DDRec/SurfaceManager.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>

#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "GenfitField.h"
#include "GenfitMaterialInterface.h"
#include "GenfitTrack.h"

#include "utils.h"

/** @struct GenfitTrackFitter
 *
 *  Gaudi MultiTransformer that refines the parameters of the reconstructed tracks using the GENFIT library.
 *  The fitting process minimizes the track chi2 while accounting for material effects, magnetic field, and detector
 * geometry.
 *
 *  If the fit fails for all hypotheses, the track is not saved to the output collections.
 *
 *  input:
 *    - initial track collection : edm4hep::TrackCollection
 *
 *  output:
 *    - fitted track collection : edm4hep::TrackCollection
 *    - fitted track collection with filtered hits (without L/R ambiguity if drift chamber is present):
 * edm4hep::TrackCollection
 *    - filtered hit collection : edm4hep::TrackerHitPlaneCollection
 *
 *  @author Andrea De Vita
 *  @date   2026-02
 *
 */

struct GenfitTrackFitter final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackCollection,
                                            edm4hep::TrackerHitPlaneCollection>(const edm4hep::TrackCollection&)> {
  GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {

                             KeyValues("InputTracks", {"InputTracks"})

                         },
                         {

                             KeyValues("OutputFittedTracks", {"Fitted_tracks"}),
                             KeyValues("OutputFittedTracksWithFilteredHits", {"Fitted_tracks_with_filtered_hits"}),
                             KeyValues("OutputFittedHits", {"Fitted_hits"})

                         }) {}

  StatusCode initialize() {

    // Initialize printout level for debug information
    m_printoutLevel = msgLevel();

    if (!gGeoManager) {
      error() << "Error: TGeoManager is not initialized!" << endmsg;
      return StatusCode::FAILURE;
    }

    // Initialize the Genfit: FieldManager and MaterialEffects
    const dd4hep::Detector* geo_detector = m_geoSvc->getDetector();
    auto fieldMap = geo_detector->field();
    m_genfitField = std::make_unique<GenfitInterface::GenfitField>(fieldMap);

    m_fieldManager = genfit::FieldManager::getInstance();
    m_fieldManager->init(m_genfitField.get());

    m_geoMaterial = GenfitInterface::GenfitMaterialInterface::getInstance(geo_detector);

    if (m_useBrems) {
      genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
      genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
    } else {
      genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
      genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
    }
    genfit::MaterialEffects::getInstance()->setMscModel(m_mscModel.value());
    genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(m_energyLossBetheBloch.value());
    genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(m_noiseBetheBloch.value());
    genfit::MaterialEffects::getInstance()->setNoiseCoulomb(m_noiseCoulomb.value());

    int debug_lvl_material = 0;
    if (m_printoutLevel == uint(MSG::VERBOSE)) {
      debug_lvl_material = 1;
    }
    genfit::MaterialEffects::getInstance()->setDebugLvl(debug_lvl_material);

    // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
    // These object are necessary to extract the drift chamber hits information, such as positions of the wire
    // extremities
    m_surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();
    // If the detector doesn't have a drift chamber, this part will be skipped
    try {

      std::string DCH_name(m_WireTracker_name.value());
      dd4hep::DetElement WireTracker_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
      m_wire_info = WireTracker_DE.extension<dd4hep::rec::WireTracker_info_struct>();
      dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
      dd4hep::Readout dch_readout = dch_sd.readout();
      m_dc_decoder = dch_readout.idSpec().decoder();

      std::string desc = m_dc_decoder->fieldDescription();

      std::set<std::string> fields;
      std::stringstream ss(desc);
      std::string token;

      while (std::getline(ss, token, ',')) {
        auto pos = token.find(':');
        if (pos != std::string::npos) {
          fields.insert(token.substr(0, pos));
        }
      }

      if (!(fields.count("superlayer") && fields.count("layer") && (fields.count("nphi") || fields.count("tube")))) {

        warning() << "WireTracker decoder missing required fields: " << desc << endmsg;
        throw std::runtime_error("Invalid WireTracker decoder");
      }

    } catch (const std::exception& e) {

      warning() << "WireTracker_info unexpected error: " << e.what()
                << ". This may indicate missing drift chamber or configuration issues." << endmsg;
    }

    // Calorimeter geometry
    // These parameters are used to propagate tracks to the
    // calorimeter surface (barrel and endcap geometry).
    if (m_runCalorimeterExtrapolation.value()) {

      dd4hep::Detector& detector = dd4hep::Detector::getInstance();

      // Check if there is a calorimeter. For simplicity, we assume that if a barrel calorimeter is present, then the
      // endcap is also present.
      const std::vector<dd4hep::DetElement> barrelECal = dd4hep::DetectorSelector(detector).detectors(
          dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL, dd4hep::DetType::AUXILIARY);

      if (barrelECal.empty()) {
        error()
            << "No barrel calorimeter found in detector description but runCalorimeterExtrapolation is set to true. "
               "Set it to false or define the detector types in the DD4hep geometry (dd4hep::DetType::CALORIMETER, ...)"
            << endmsg;
        return StatusCode::FAILURE;
      }

      // If IDEA_o1 we expect that the calorimeter element is not separated in barrel and endcap,
      // but we have a single calorimeter element with both barrel and endcap information.
      const std::vector<dd4hep::DetElement> dualReadoutElements = dd4hep::DetectorSelector(detector).detectors(
          dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY);

      if (!dualReadoutElements.empty()) {

        const dd4hep::rec::LayeredCalorimeterData* ecalExt = nullptr;
        try {

          ecalExt =
              getExtension(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL,
                           dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);

          if (!ecalExt) {
            error() << "ECal LayeredCalorimeterData is null: cannot retrieve calorimeter geometry." << endmsg;
            return StatusCode::FAILURE;
          }

        } catch (const std::exception& e) {

          error() << "Failed to retrieve ECal LayeredCalorimeterData: " << e.what()
                  << ". Calorimeter geometry for propagation is not available — please implement it or switch off "
                     "the extrapolation."
                  << endmsg;

          return StatusCode::FAILURE;
        }

        // Convert to mm (dd4hep::mm = 0.1)
        m_eCalBarrelInnerR = ecalExt->extent[0] / dd4hep::mm; // barrel rmin
        m_eCalBarrelMaxZ = ecalExt->extent[2] / dd4hep::mm;   // barrel zmax == endcap zmin

        m_eCalEndCapInnerR = ecalExt->extent[4] / dd4hep::mm; // endcap rmin
        m_eCalEndCapOuterR = ecalExt->extent[5] / dd4hep::mm; // endcap rmax
        m_eCalEndCapInnerZ = ecalExt->extent[2] / dd4hep::mm; // endcap zmin
        m_eCalEndCapOuterZ = ecalExt->extent[3] / dd4hep::mm; // endcap zmax

      } else {

        const dd4hep::rec::LayeredCalorimeterData* ecalBarrelExt = nullptr;
        try {
          ecalBarrelExt =
              getExtension(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL,
                           dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);
          if (!ecalBarrelExt) {
            error() << "ECal Barrel LayeredCalorimeterData is null: cannot retrieve calorimeter geometry." << endmsg;
            return StatusCode::FAILURE;
          }

        } catch (const std::exception& e) {

          error() << "Failed to retrieve ECal LayeredCalorimeterData: " << e.what()
                  << ". Calorimeter geometry for propagation is not available — please implement it or switch off "
                     "the extrapolation."
                  << endmsg;

          return StatusCode::FAILURE;
        }

        // Convert to mm (dd4hep::mm = 0.1)
        m_eCalBarrelInnerR = ecalBarrelExt->extent[0] / dd4hep::mm;
        m_eCalBarrelMaxZ = ecalBarrelExt->extent[3] / dd4hep::mm;

        const dd4hep::rec::LayeredCalorimeterData* ecalEndCapExt = nullptr;
        try {
          ecalEndCapExt =
              getExtension(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP,
                           dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD);

          if (!ecalEndCapExt) {
            error() << "ECal EndCap LayeredCalorimeterData is null: cannot retrieve calorimeter geometry." << endmsg;
            return StatusCode::FAILURE;
          }

        } catch (const std::exception& e) {

          error() << "Failed to retrieve ECal EndCap LayeredCalorimeterData: " << e.what()
                  << ". Calorimeter geometry for propagation is not available — please implement it or switch off "
                     "the extrapolation."
                  << endmsg;

          return StatusCode::FAILURE;
        }

        // Convert to mm (dd4hep::mm = 0.1)
        m_eCalEndCapInnerR = ecalEndCapExt->extent[0] / dd4hep::mm;
        m_eCalEndCapOuterR = ecalEndCapExt->extent[1] / dd4hep::mm;
        m_eCalEndCapInnerZ = ecalEndCapExt->extent[2] / dd4hep::mm;
        m_eCalEndCapOuterZ = ecalEndCapExt->extent[3] / dd4hep::mm;
      }
    }

    return StatusCode::SUCCESS;
  }

  // GenFit and ROOT geometry navigation use process-global mutable state.
  // Do not allow the scheduler to execute this algorithm concurrently for
  // different events.
  // If changing to multi-threaded be careful with side effects when
  // resetting the TGeoManager, that will apply globally. In addition
  // dereferencing gGeoManager->GetListOfNavigators() may crash
  // for any thread that is not the one that initialized the geometry
  bool isReEntrant() const override { return false; }

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection>
  operator()(const edm4hep::TrackCollection& tracks_input) const override {

    debug() << "Event number: " << event_counter++ << endmsg;
    num_tracks_event = 0;

    // These collections store the output of the fit
    edm4hep::TrackCollection FittedTracks;
    edm4hep::TrackCollection FittedTracksWithFilteredHits;
    edm4hep::TrackerHitPlaneCollection FittedHits;

    // Loop over the tracks created by the pattern recognition step
    for (const auto& track : tracks_input) {

      num_tracks++;
      num_tracks_event++;

      // Skip unmatched tracks if the option is enabled
      if (m_ListOfTypesToSkip.size() > 0 && std::find(m_ListOfTypesToSkip.begin(), m_ListOfTypesToSkip.end(),
                                                      track.getType()) != m_ListOfTypesToSkip.end()) {
        num_skip += 1;
        debug() << "Skipping track " << num_tracks_event << " with type " << track.getType() << "\n" << endmsg;
        continue;
      }

      // skip tracks with less then 3 hits (seed initialization needs 3 hits)
      if (track.getTrackerHits().size() < 3) {
        num_skip += 1;
        debug() << "Track " << num_tracks_event << ": less than 3 hits, skipping fit.\n" << endmsg;
        continue;
      }

      num_processed_tracks += 1;

      // Per-track cache and attempt ledger: both live only for this track's processing and
      // are never shared with other tracks or events, so the results cannot depend on the
      // order of tracks in the event or on the order of events.
      auto cache = std::make_shared<GenfitInterface::GenfitTrack::TrackFitCache>();
      std::vector<AttemptKey> attempted;

      // Attempt chain: all hits (with the automatic first-arc selection for loopers), the forced
      // first-arc selection, the silicon-only fallback, a thinned hit set, a scan of alternative
      // momentum seeds on the thinned set, and finally fits without material effects
      auto tryFit = [&](bool forceFirstArc, bool siOnly, int thinTarget, double seedOverride, bool noMaterial,
                        bool ransacInliers) -> int {
        ransacInliers = ransacInliers && m_ransacFallback;
        std::optional<FittedResult> result;
        if (m_singleEvaluation) {
          result = FitTrack(track, forceFirstArc, siOnly, ransacInliers, thinTarget, seedOverride, noMaterial,
                            m_particleHypothesis[0], cache, attempted);
        } else {
          result = FitBestHypothesis(track, forceFirstArc, siOnly, ransacInliers, thinTarget, seedOverride,
                                     noMaterial, cache, attempted);
        }
        if (!result) {
          return 0;
        }
        FillOutput(std::move(*result), track, FittedTracks, FittedTracksWithFilteredHits, FittedHits,
                   m_runCalorimeterExtrapolation.value());
        return 1;
      };

      int isSuccess = tryFit(false, false, 0, 0., false, false);

      if (!isSuccess) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying with the first-arc hit selection." << endmsg;
        isSuccess = tryFit(true, false, 0, 0., false, false);
      }

      if (!isSuccess && m_siOnlyFallback) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying with the silicon hits only." << endmsg;
        isSuccess = tryFit(false, true, 0, 0., false, false);
      }

      if (!isSuccess && m_ransacFallback) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying with the dominant-circle inliers." << endmsg;
        isSuccess = tryFit(true, false, 0, 0., false, true);
      }

      if (!isSuccess) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying with a thinned hit set." << endmsg;
        isSuccess = tryFit(true, false, 50, 0., false, true);
      }

      for (double seedP : {0.2, 0.7, 2.0}) {
        if (isSuccess) {
          break;
        }
        debug() << "Track " << num_tracks_event << ": fit failed, trying with a " << seedP << " GeV/c momentum seed."
                << endmsg;
        isSuccess = tryFit(true, false, 50, seedP, false, true);
      }

      if (!isSuccess) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying without material effects." << endmsg;
        isSuccess = tryFit(true, false, 0, 0., true, true);
      }

      if (!isSuccess) {
        debug() << "Track " << num_tracks_event << ": fit failed, trying a thinned set without material effects."
                << endmsg;
        isSuccess = tryFit(true, false, 50, 0., true, true);
      }

      if (!isSuccess) {
        number_failures += 1;
        debug() << "Track " << num_tracks_event << ": fit failed for all attempts, skipping track." << endmsg;
        continue;
      }
    }

    return std::make_tuple(std::move(FittedTracks), std::move(FittedTracksWithFilteredHits), std::move(FittedHits));
  }

  StatusCode finalize() {

    info() << "Run report:" << endmsg;
    info() << "Number of tracks: " << num_tracks << endmsg;
    info() << "Number of successes: " << (num_processed_tracks - number_failures) << "/" << num_processed_tracks
           << endmsg;
    info() << "Number of skipped tracks: " << num_skip << "/" << num_tracks << endmsg;
    info() << "----------------\n" << endmsg;

    return StatusCode::SUCCESS;
  }

public:
  // Debug / monitoring counters (not part of core algorithm logic).
  // These are mainly used to keep track of what happens during processing,
  // e.g. for logging, sanity checks, and performance evaluation.

  mutable int event_counter = 0; // Counts how many events have been processed

  // Num_tracks = num_processed_tracks + num_skip
  mutable int num_tracks = 0;           // Total number of tracks seen (including skipped ones)
  mutable int num_tracks_event = 0;     // Total number of tracks seen in one event (including skipped ones)
  mutable int num_skip = 0;             // Number of tracks skipped (e.g. failing pre-selection)
  mutable int num_processed_tracks = 0; // Number of tracks actually processed (i.e. not skipped)
  mutable int number_failures = 0;      // Number of track fits that failed

private:
  //////////////////////
  // Debug & Printout //
  //////////////////////

  // Debug level for track fitting and initialization printouts
  // 0       : no printouts
  // INFO    : prints the final fit summary
  // DEBUG   : INFO + initial track parameters, measurements, track states, and fit failures
  // VERBOSE : DEBUG + detailed internal fitter information
  uint m_printoutLevel;

  /////////////////////////
  // Geometry & Detector //
  /////////////////////////

  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};

  std::unique_ptr<GenfitInterface::GenfitField> m_genfitField;
  genfit::FieldManager* m_fieldManager;

  GenfitInterface::GenfitMaterialInterface* m_geoMaterial;

  dd4hep::rec::SurfaceManager* m_surfMan;
  dd4hep::rec::WireTracker_info_struct* m_wire_info{nullptr};
  dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder{nullptr};

  Gaudi::Property<std::string> m_WireTracker_name{
      this, "WireTrackerName", "DCH_v2",
      "WireTrackerName in the detector description (used to retrieve Wire Tracker geometry and material information)"};

  // ECAL Geometry Parameters
  double m_eCalBarrelInnerR;
  double m_eCalBarrelMaxZ;
  double m_eCalEndCapInnerR;
  double m_eCalEndCapOuterR;
  double m_eCalEndCapInnerZ;
  double m_eCalEndCapOuterZ;

  /////////////////////
  // Fitter Settings //
  /////////////////////

  Gaudi::Property<bool> m_useBrems{this, "UseBrems", false, "Include Bremsstrahlung energy loss and noise in the fit"};

  Gaudi::Property<std::string> m_mscModel{
      this, "MscModel", "Highland",
      "Multiple-scattering model used for the material noise (genfit MSC formulas)"};
  Gaudi::Property<bool> m_energyLossBetheBloch{
      this, "EnergyLossBetheBloch", true, "Include Bethe-Bloch energy loss in the material effects"};
  Gaudi::Property<bool> m_noiseBetheBloch{
      this, "NoiseBetheBloch", true, "Include Bethe-Bloch energy-loss noise in the material effects"};
  Gaudi::Property<bool> m_noiseCoulomb{
      this, "NoiseCoulomb", true, "Include multiple-scattering (Coulomb) noise in the material effects"};

  Gaudi::Property<std::string> m_Fitter_type{
      this, "FitterType", "DAF",
      "Fitter type to use [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]: 'DAF', 'KALMAN', "
      "or 'KALMAN_REF'"};
  Gaudi::Property<double> m_Beta_init{this, "BetaInit", 100.,
                                      "Initial annealing parameter (if m_Fitter_type == 'DAF') "
                                      "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};
  Gaudi::Property<double> m_Beta_final{this, "BetaFinal", 0.05,
                                       "Final annealing parameter (if m_Fitter_type == 'DAF') "
                                       "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};
  Gaudi::Property<int> m_Beta_steps{this, "BetaSteps", 15,
                                    "Number of steps in the annealing schedule (if m_Fitter_type == 'DAF') "
                                    "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};

  //////////////////////////
  // Track Initialization //
  //////////////////////////

  Gaudi::Property<std::vector<int>> m_particleHypothesis{
      this,
      "ParticleHypothesisList",
      {211},
      "List of particle hypotheses to consider: {11,13,211,321,2212} -> e, mu, pi, K, p"};

  Gaudi::Property<bool> m_useFirstHitAsReference{
      this, "UseFirstHitAsReference", false,
      "If true, the first hit is used as the reference point for track initialization"};

  Gaudi::Property<int> m_initializationType{
      this, "InitializationType", 1,
      "Method for initializing track parameters before fitting: "
      "0: first two hits (position from first hit, momentum from direction between first and second hit); "
      "1: refined initialization using ComputeInitialParameters(Bz); "
      "2: initialize from a track state (position from reference point, momentum from helix parameters); "
      "3: use user-provided InitPosition and InitMomentum"};

  Gaudi::Property<double> m_sigma_d0{
      this, "Sigma_d0", 0.05,
      "Initial uncertainty on d0 (mm) for the initial covariance matrix when InitializationType = 0,1,3."};
  Gaudi::Property<double> m_sigma_phi{
      this, "Sigma_phi", 0.1,
      "Initial uncertainty on phi (rad) for the initial covariance matrix when InitializationType = 0,1,3."};
  Gaudi::Property<double> m_omega_factor{
      this, "OmegaFactor", 0.5,
      "Scaling factor for omega uncertainty for the initial covariance matrix when InitializationType = 0,1,3."
      "The actual sigma_omega is computed as OmegaFactor * |omega|"};
  Gaudi::Property<double> m_z0_factor{
      this, "Z0Factor", 0.1,
      "Scaling factor for z0 uncertainty for the initial covariance matrix when InitializationType = 0,1,3."
      "The actual sigma_z0 is computed as Z0Factor * |z0 (in mm)|"};
  Gaudi::Property<double> m_sigma_tanLambda{
      this, "Sigma_tanLambda", 0.1,
      "Initial uncertainty on tanLambda for the initial covariance matrix when InitializationType = 0,1,3."};

  Gaudi::Property<int> m_trackStateLocation{this, "TrackStateLocation", edm4hep::TrackState::AtFirstHit,
                                            "TrackState location used for initialization when InitializationType = 2. "
                                            "Defines where the reference point and helix parameters are taken from"};

  Gaudi::Property<std::vector<double>> m_init_position{
      this, "InitPosition", {0., 0., 0.}, "User-defined initial position [in mm] for InitializationType = 3"};

  Gaudi::Property<std::vector<double>> m_init_momentum{
      this, "InitMomentum", {0., 0., 0.}, "User-defined initial momentum [in GeV/c] for InitializationType = 3"};

  Gaudi::Property<double> m_looperPTThreshold{
      this, "LooperPTThreshold", 0.7,
      "Transverse momentum threshold [GeV/c] below which a track is treated as a looper: the circle is fitted in the "
      "xy plane, phi_local (the angle around the fitted circle center, zeroed at the hit closest to the z-origin) is "
      "computed, and only the hits before the first peak of the sawtooth in the (z, phi_local) plane are seeded and "
      "fitted, i.e. the first arc of the loop. Tracks curling inside the tracker scale with the magnetic field "
      "(pT ~ 0.3 * B[T] * R[m]): adjust the threshold when the field is changed (0.7 GeV/c corresponds to B = 2 T)."};

  ////////////////////////////////////
  // Track Filtering and Evaluation //
  ////////////////////////////////////

  Gaudi::Property<bool> m_skipTrackOrdering{this, "SkipTrackOrdering", false, "Skip hit ordering before fitting"};

  Gaudi::Property<std::vector<int>> m_ListOfTypesToSkip{
      this, "ListOfTypesToSkip", std::vector<int>(),
      "List of track types to skip during fitting (e.g. unmatched tracks with type = 0)"};

  Gaudi::Property<bool> m_singleEvaluation{
      this, "RunSingleEvaluation", false,
      "If true, only the first particle hypothesis is evaluated. "
      "If false, all hypotheses are scanned and the best fit is chosen based on chi2/ndf"};

  Gaudi::Property<bool> m_filterTrackHits{
      this, "FilterTrackHits", true,
      "Filter track hits after fitting based on their quality (weights assigned by the fitter). "
      "This also resolves the left/right ambiguity in drift chamber hits: only the hit with the higher weight is "
      "considered."};

  Gaudi::Property<bool> m_siOnlyFallback{
      this, "SiliconOnlyFallback", true,
      "If true, tracks whose fit fails both with all the hits and with the first-arc hits of a looper are retried "
      "fitting only the silicon (planar) hits of the selected hit set: those are much more precise than the drift "
      "chamber ones (whose z resolution is of the order of cm) and the fit often converges on them alone."};

  Gaudi::Property<bool> m_ransacFallback{
      this, "CircleInlierFallback", true,
      "If true, tracks whose fit fails with the other hit sets are retried on the inliers of a robust (RANSAC) "
      "circle fit, ordered along the dominant circle: this isolates one turn of one particle in impure (merged) "
      "tracks and gives the fitter a proper hit order."};

  Gaudi::Property<bool> m_runCalorimeterExtrapolation{
      this, "RunCalorimeterExtrapolation", true,
      "If true, the track is extrapolated to the calorimeter surfaces (barrel and endcap)"
      "and the corresponding track states are stored in the output track."
      "This is applied only if a non-zero magnetic field (Bz > 0) is found at the last hit position."};

  //////////////////////////
  // Track Classification //
  //////////////////////////

  Gaudi::Property<double> m_RadialThresholdPromptTrack{
      this, "RadialThresholdPromptTrack", 100.,
      "Radius [mm] defining a spherical region centered at {0,0,0}. "
      "Tracks produced within this radius are classified as prompt, while those outside are considered displaced. "
      "This parameter is used to select different track initialization strategies: prompt tracks are initialized using "
      "the PCA of the hits with respect to a reference point (IP or VP), whereas displaced tracks are initialized "
      "using the PCA "
      "of the first hit."};

  /////////////////
  // Fit Helpers //
  /////////////////

  // Result of a successful fit: the interface owns the fitted EDM4hep tracks, and the
  // filtered hits were produced by this fit only. Both are local to one track's processing:
  // nothing is shared between tracks or events, so results cannot depend on the processing
  // order of tracks or events.
  struct FittedResult {
    std::unique_ptr<GenfitInterface::GenfitTrack> interface;
  };

  // Identity of a fit configuration: the effective hit selection, the seed, the material
  // setting and the particle hypothesis. The fit is deterministic, so a configuration that
  // has already been attempted in this track's fallback chain (and failed, since the chain
  // stops at the first success) cannot succeed on a second execution and is skipped.
  struct AttemptKey {
    std::vector<std::size_t> fitHits;
    double seed = 0.;
    bool noMaterial = false;
    int pdg = 0;
    bool operator==(const AttemptKey& other) const {
      return seed == other.seed && noMaterial == other.noMaterial && pdg == other.pdg && fitHits == other.fitHits;
    }
  };

  /**
   * @brief Initialize and fit one configuration of one track, returning the fitted result.
   *
   * This is the per-configuration unit of work of ProcessTrack: initialization (seeded with
   * the cached geometry products of the track), hit restrictions, GENFIT track creation
   * (cloned cached measurements) and the fit itself. Identical configurations are executed
   * at most once per track (deterministic fit). Everything is local to the call apart from
   * the per-track cache and the per-track attempt ledger passed in by the caller.
   *
   * @return the fitted result, or std::nullopt if the configuration is restricted away,
   *         duplicated, or the fit fails or throws.
   */
  std::optional<FittedResult>
  FitTrack(const edm4hep::Track& track, bool ForceFirstArc, bool SiOnly, bool RansacInliers, int ThinTarget,
           double SeedMomentumOverride, bool NoMaterial, int particleHypothesis,
           const std::shared_ptr<GenfitInterface::GenfitTrack::TrackFitCache>& cache,
           std::vector<AttemptKey>& attempted) const {

    auto track_interface = std::make_unique<GenfitInterface::GenfitTrack>(track, m_skipTrackOrdering, m_wire_info,
                                                                         m_dc_decoder, m_genfitField.get());
    track_interface->SetTrackCache(cache);

    TVector3 Init_position(m_init_position.value()[0] * dd4hep::mm, m_init_position.value()[1] * dd4hep::mm,
                           m_init_position.value()[2] * dd4hep::mm);

    TVector3 Init_momentum(m_init_momentum.value()[0], m_init_momentum.value()[1], m_init_momentum.value()[2]);

    double RadialThresholdPromptTrack_cm = m_RadialThresholdPromptTrack.value() * dd4hep::mm;

    track_interface->InitializeTrack(RadialThresholdPromptTrack_cm, m_useFirstHitAsReference, ForceFirstArc,
                                     m_initializationType, m_trackStateLocation.value(), Init_position, Init_momentum,
                                     m_looperPTThreshold.value(), SeedMomentumOverride, m_sigma_d0.value(),
                                     m_sigma_phi.value(), m_omega_factor.value(), m_z0_factor.value(),
                                     m_sigma_tanLambda.value());

    if (SiOnly && !track_interface->RestrictToSiliconHits()) {
      return std::nullopt;
    }

    if (RansacInliers && !track_interface->RestrictToCircleInliers(3.0)) {
      return std::nullopt;
    }

    if (ThinTarget > 0) {
      track_interface->RestrictToEvenlySpacedHits(static_cast<std::size_t>(ThinTarget));
    }

    const AttemptKey key{track_interface->GetFitHits(), SeedMomentumOverride, NoMaterial, particleHypothesis};
    if (std::find(attempted.begin(), attempted.end(), key) != attempted.end()) {
      return std::nullopt; // identical configuration already attempted and failed: skip the duplicate fit
    }
    attempted.push_back(key);

    auto track_init = track_interface->GetInitialization();
    const bool showDebugOutput = m_printoutLevel == uint(MSG::DEBUG) || m_printoutLevel == uint(MSG::VERBOSE);

    debug() << "Track " << num_tracks_event << " with " << track.getTrackerHits().size()
            << " hits: initial seed for track fit:" << endmsg;

    debug() << "  Initial position [mm]: (" << track_init.Position.X() / dd4hep::mm << ", "
            << track_init.Position.Y() / dd4hep::mm << ", " << track_init.Position.Z() / dd4hep::mm << ")" << endmsg;

    debug() << "  Initial momentum [GeV/c]: (" << track_init.Momentum.X() << ", " << track_init.Momentum.Y() << ", "
            << track_init.Momentum.Z() << ")" << endmsg;

    debug() << "  Charge hypothesis: " << track_init.Charge << endmsg;
    debug() << "  Hits used for the fit: " << track_init.NumHits << " / " << track.getTrackerHits().size() << endmsg;

    if (showDebugOutput) {
      debug() << "  Initial covariance matrix:" << endmsg;
      track_init.CovMatrix.Print();
    }

    debug() << endmsg;

    int debug_track = showDebugOutput ? 1 : 0;

    track_interface->CreateGenFitTrack(particleHypothesis, debug_track);

    FittedResult result;
    result.interface = std::move(track_interface);

    bool isFit = false;
    try {
      isFit = result.interface->Fit(m_Fitter_type.value(), m_printoutLevel, m_Beta_init, m_Beta_final, m_Beta_steps,
                                    m_filterTrackHits, NoMaterial);
    } catch (const std::exception& e) {
      debug() << "Track fit threw for track " << num_tracks_event << ": " << e.what() << endmsg;
      return std::nullopt;
    }

    if (!isFit) {
      debug() << "Track fit FAILED for track " << num_tracks_event << endmsg;
      return std::nullopt;
    }

    // Acceptance criterion shared by all fit modes: a converged fit with a non-positive chi2
    // or ndf carries no fit-quality information and is not reported (it falls through to the
    // fallback attempts of the chain instead).
    const auto& fitted_track = result.interface->GetTrack_edm4hep();
    if (fitted_track.getChi2() <= 0 || fitted_track.getNdf() <= 0) {
      debug() << "Track fit produced invalid chi2/ndf for track " << num_tracks_event << endmsg;
      return std::nullopt;
    }

    return result;
  }

  /**
   * @brief Store the result of a successful fit in the output collections.
   *
   * The post-fit stage of the old ProcessTrack: field query at the last hit, optional
   * calorimeter extrapolation, and the output copies (fitted track, filtered-hit track and
   * filtered hits of this fit). The values written are the same ones the old code produced
   * from the equivalent re-fit of the winning hypothesis.
   */
  void FillOutput(FittedResult&& result, const edm4hep::Track& track, edm4hep::TrackCollection& FittedTracks,
                  edm4hep::TrackCollection& FittedTracksWithFilteredHits,
                  edm4hep::TrackerHitPlaneCollection& FittedHits, bool runCalorimeterExtrapolation) const {

    auto& track_interface = *result.interface;
    auto& edm4hep_track = track_interface.GetTrack_edm4hep();
    const bool showDebugOutput = m_printoutLevel == uint(MSG::DEBUG) || m_printoutLevel == uint(MSG::VERBOSE);

    double Bz = 0.;
    for (auto ts : edm4hep_track.getTrackStates()) {
      if (ts.location == edm4hep::TrackState::AtLastHit) {

        auto ref = ts.referencePoint;
        Bz = m_genfitField->getBz(TVector3(ref.x * dd4hep::mm, ref.y * dd4hep::mm, ref.z * dd4hep::mm)) /
             (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla;
      }
    }

    if (runCalorimeterExtrapolation && Bz > 0) {

      FillTrackWithCalorimeterExtrapolation(edm4hep_track, Bz, track_interface.GetCharge(), m_eCalBarrelInnerR,
                                            m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR,
                                            m_eCalEndCapInnerZ);

      if (showDebugOutput) {
        auto trackStates = edm4hep_track.getTrackStates();
        edm4hep::TrackState trackStateCalo;
        for (const auto& ts : trackStates) {
          if (ts.location == edm4hep::TrackState::AtCalorimeter) {
            trackStateCalo = ts;
            break;
          }
        }

        debug() << ": TrackState at Calo: " << endmsg;
        debug() << ":  D0: " << trackStateCalo.D0 << " mm" << endmsg;
        debug() << ":  Z0: " << trackStateCalo.Z0 << " mm" << endmsg;
        debug() << ":  phi: " << trackStateCalo.phi << " rad" << endmsg;
        debug() << ":  omega: " << trackStateCalo.omega << " 1/mm" << endmsg;
        debug() << ":  tanLambda: " << trackStateCalo.tanLambda << endmsg;
        debug() << ":  location: " << trackStateCalo.location << endmsg;
        debug() << ":  reference point: (" << trackStateCalo.referencePoint.x << ", " << trackStateCalo.referencePoint.y
                << ", " << trackStateCalo.referencePoint.z << ") mm\n"
                << endmsg;
      }
    }

    FittedTracks.push_back(edm4hep_track);

    if (m_filterTrackHits) {

      // The filtered hits of this fit are created directly in the output collection and
      // related to the filtered-hit track (an EDM4hep object cannot be shared between two
      // collections). Only the winning configuration reaches this point, so the output
      // receives exactly the hits of the fit being reported.
      track_interface.CreateFilteredHits(FittedHits);

      auto edm4hep_track_with_fit = track_interface.GetTrackWithFit_edm4hep();

      if (runCalorimeterExtrapolation && Bz > 0) {

        FillTrackWithCalorimeterExtrapolation(edm4hep_track_with_fit, Bz, track_interface.GetCharge(),
                                              m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR,
                                              m_eCalEndCapOuterR, m_eCalEndCapInnerZ);

        if (showDebugOutput) {
          auto trackStates = edm4hep_track.getTrackStates();
          edm4hep::TrackState trackStateCalo;
          for (const auto& ts : trackStates) {
            if (ts.location == edm4hep::TrackState::AtCalorimeter) {
              trackStateCalo = ts;
              break;
            }
          }
        }
      }

      FittedTracksWithFilteredHits.push_back(edm4hep_track_with_fit);
    }
  }

  /**
   * @brief Fit all particle hypotheses of one configuration and keep the best result.
   *
   * Replaces the old FindBestHypothesis + re-fit pair: the hypotheses are fitted through
   * FitTrack (which shares the per-track initialization products) and the fitted result of
   * the best chi2/ndf is returned instead of being discarded and recomputed by a second fit.
   * The chi2/ndf comparison is unchanged and the fit itself is identical for every
   * hypothesis (hit filtering is a post-fit step and does not enter the comparison).
   *
   * @return the fitted result of the winning hypothesis, or std::nullopt if no hypothesis
   *         produced a valid fit.
   */
  std::optional<FittedResult>
  FitBestHypothesis(const edm4hep::Track& track, bool ForceFirstArc, bool SiOnly, bool RansacInliers, int ThinTarget,
                    double SeedMomentumOverride, bool NoMaterial,
                    const std::shared_ptr<GenfitInterface::GenfitTrack::TrackFitCache>& cache,
                    std::vector<AttemptKey>& attempted) const {

    std::optional<FittedResult> best;
    double winning_chi2_ndf = std::numeric_limits<double>::max();

    for (int pdgCode : m_particleHypothesis) {

      auto result = FitTrack(track, ForceFirstArc, SiOnly, RansacInliers, ThinTarget, SeedMomentumOverride,
                             NoMaterial, pdgCode, cache, attempted);

      if (!result) {
        continue;
      }

      auto& edm4hep_track = result->interface->GetTrack_edm4hep();

      float chi2 = edm4hep_track.getChi2();
      int ndf = edm4hep_track.getNdf();

      double chi2_ndf = chi2 / ndf;

      if (chi2_ndf < winning_chi2_ndf) {

        winning_chi2_ndf = chi2_ndf;
        best = std::move(result);
      }
    }

    return best;
  }
};

DECLARE_COMPONENT(GenfitTrackFitter)