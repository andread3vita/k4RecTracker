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

//=== Standard Library ===
#include <algorithm>
#include <cmath>
#include <cstdlib>  
#include <fstream>  
#include <filesystem>  
#include <iostream>
#include <iterator> 
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

//=== ROOT / C++ ABI ===
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <cxxabi.h>
#include <Eigen/Dense>

//=== GenFit ===
#include <ConstField.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <PlanarMeasurement.h>
#include <RKTrackRep.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <TGeoMaterialInterface.h>

//=== Gaudi / k4FWCore ===
#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

//=== DD4hep / DDRec / DDSegmentation ===
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>
#include "DD4hep/Fields.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/MaterialManager.h"

//=== podio / edm4hep ===
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"

//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "GenfitWireMeasurement.hpp"
#include "GenfitPlanarMeasurement.hpp"
#include "GenfitField.hpp"
#include "GenfitTrack.hpp"
#include "GenfitMaterialInterface.hpp"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>

// Define collection types
#include "podio/UserDataCollection.h"

#include "utils.hpp"
#include "FastCircleSeed.hpp"

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TFile.h"


/** @struct GenfitTrackFitter
 *
 *  Gaudi MultiTransformer that refines the parameters of the reconstructed tracks using the GENFIT library.  
 *  For each input track, the module performs a fit under five different particle hypotheses (electron, muon, pion, kaon, proton), 
 *  producing five distinct output collections of fitted tracks.  
 *  The fitting process minimizes the track χ² while accounting for material effects, magnetic field, and detector geometry.
 *
 *  If the fit fails for a given hypothesis, a fallback track is still generated with χ² = -1 and ndf = -1 to preserve collection integrity.
 *  Each fitted track also includes the extrapolated state at the calorimeter (barrel or endcap).
 *
 *  input:
 *    - initial track collection : edm4hep::TrackCollection
 *
 *  output:
 *    - fitted track collection under electron hypothesis : edm4hep::TrackCollection
 *    - fitted track collection under muon hypothesis : edm4hep::TrackCollection
 *    - fitted track collection under pion hypothesis : edm4hep::TrackCollection
 *    - fitted track collection under kaon hypothesis : edm4hep::TrackCollection
 *    - fitted track collection under proton hypothesis : edm4hep::TrackCollection
 *
 *
 *  @author Andrea De Vita
 *  @date   2025-11
 *
*/

struct Point3D {
    double x;
    double y;
    double z;
};

struct TrackInitialState
{

    TVector3 init_pos;
    TVector3 init_mom;
    int charge;
   
};

std::vector<Point3D> preparePoints(const edm4hep::Track& track, int dir = 1)
{
    const auto hits = track.getTrackerHits();

    // Track endpoints in z
    double zMin =  std::numeric_limits<double>::max();
    double zMax = -std::numeric_limits<double>::max();

    double xAtZMin = 0., yAtZMin = 0.;
    double xAtZMax = 0., yAtZMax = 0.;

    // Find hits with minimum and maximum z
    for (const auto& hit : hits) {
        const auto& p = hit.getPosition();

        if (p.z < zMin) {
            zMin = p.z;
            xAtZMin = p.x;
            yAtZMin = p.y;
        }
        if (p.z > zMax) {
            zMax = p.z;
            xAtZMax = p.x;
            yAtZMax = p.y;
        }
    }

    // Choose the starting point:
    // - by default, the hit closer to z = 0
    // - reversed if dir < 0
    const bool minCloserToZero = std::abs(zMin) < std::abs(zMax);
    const bool takeMin = (dir > 0) ? minCloserToZero : !minCloserToZero;

    double firstX = takeMin ? xAtZMin : xAtZMax;
    double firstY = takeMin ? yAtZMin : yAtZMax;
    double firstZ = takeMin ? zMin    : zMax;

    // Estimate track inclination in the yz-plane
    const double dy = std::abs(yAtZMax - yAtZMin);
    const double dz = std::abs(zMax    - zMin);
    const double cosTheta = dz / std::hypot(dy, dz);

    // If the track is almost perpendicular to z,
    // fall back to the hit with smaller radius
    if (std::abs(cosTheta) < 0.01) {
        const double rMin = std::hypot(xAtZMin, yAtZMin, zMin);
        const double rMax = std::hypot(xAtZMax, yAtZMax, zMax);

        const bool minHasSmallerR = (dir > 0) ? (rMin < rMax) : (rMin > rMax);

        firstX = minHasSmallerR ? xAtZMin : xAtZMax;
        firstY = minHasSmallerR ? yAtZMin : yAtZMax;
        firstZ = minHasSmallerR ? zMin    : zMax;
    }

    std::cout<< "GenfitTrackFitter    DEBUG First point chosen at (x, y, z) = (" << firstX << ", " << firstY << ", " << firstZ << ")" << std::endl;

    // Compute distances from the chosen starting point
    std::vector<std::pair<float, std::size_t>> distIndex;
    distIndex.reserve(hits.size());

    for (std::size_t i = 0; i < hits.size(); ++i) {
        const auto& p = hits[i].getPosition();
        const float d = std::hypot(p.x - firstX, p.y - firstY, p.z - firstZ);
        distIndex.emplace_back(d, i);
    }

    // Sort hits along the track
    std::ranges::sort(distIndex, {}, &std::pair<float, std::size_t>::first);

    // Build (z, y) points in track order
    std::vector<Point3D> xyzPoints;
    xyzPoints.reserve(hits.size());

    for (const auto& [_, idx] : distIndex) {
        const auto& p = hits[idx].getPosition();
        xyzPoints.push_back({p.x, p.y, p.z});
    }

    return xyzPoints;
}

int findFirstExtremum(const std::vector<Point3D>& points,
                      double epsilon = 0.05,
                      int smoothWindow = 5)
{
    int n = points.size();
    if (n < smoothWindow || n < 3) return -1; // not enough points

    int W = smoothWindow / 2;

    // Compute smoothed y-values
    std::vector<double> ys(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        int count = 0;
        for (int k = -W; k <= W; ++k) {
            int j = i + k;
            if (j >= 0 && j < n) {
                sum += points[j].y;
                count++;
            }
        }
        ys[i] = sum / count;
    }

    // Compute approximate derivative along x
    std::vector<double> d(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        double dz = points[i+1].z - points[i-1].z;
        if (dz != 0.0) {
            d[i] = (ys[i+1] - ys[i-1]) / dz;
        }
    }

    // Find first zero-crossing of derivative above threshold
    for (int i = 1; i < n - 1; ++i) {
        if ((d[i-1] >  epsilon && d[i] < -epsilon) ||  // local max
            (d[i-1] < -epsilon && d[i] >  epsilon))    // local min
        {
            return i;
        }
    }

    return -1; // no extremum found
}


TrackInitialState inizialize_seed(std::vector<Point3D> points_raw, int maxHit = 50, double B = 2.0)
{

    // FIT CIRCLE TO XY PROJECTION
    std::vector<Point3D> points;

    int idx_fill = 0;
    for (const auto& p : points_raw) {
        
        if (idx_fill > maxHit) continue;
        idx_fill +=1;

        points.push_back(Point3D(p.x, p.y, p.z));
    }

    std::vector<Point2D_xy> points_xy;
    for (const auto& p : points) {
        points_xy.push_back(Point2D_xy(p.x, p.y));
    }

    // FIT CIRCLE TO XY PROJECTION
    FastCircleFit circle(points_xy);
    Point2D_xy closestPoint = circle.closestPointTo(points_xy[0]);
    Point2D_xy tangent_xy = circle.tangentAtPCA(closestPoint, points_xy[1]);

    double rho = circle.rho();
    double init_pT = std::abs(rho * 0.3 * B) / 1000; 

    TVector3 init_mom = TVector3(tangent_xy.x * init_pT, tangent_xy.y * init_pT, 0);

    // FIT LINE TO SZ PROJECTION: s_j = r | phi_j - phi_(j-1)|
    double phi1 = std::atan2(closestPoint.y - circle.y0(), closestPoint.x - circle.x0());
    double phi2 = std::atan2(points_xy[1].y - circle.y0(), points_xy[1].x - circle.x0());

    if (phi2 - phi1 > M_PI) phi2 -= 2*M_PI;
    if (phi2 - phi1 < - M_PI) phi2 += 2*M_PI;

    double s1 = 0.0;
    double s2 = rho * std::abs(phi2 - phi1);
    double z1 = points[0].z;
    double z2 = points[1].z;

    double sums = s1 + s2;
    double sumz = z1 + z2;
    double sumsz = s1*z1 + s2*z2;
    double sums2 = s1*s1 + s2*s2;
    const int n = 2;
    
    // Calculate slope and intercept
    double denominator = n * sums2 - sums * sums;
    double slope = 0.0;
    // double intercept = 0.0;
    if (std::abs(denominator) < 1e-6) {
        
        slope = 1e6;
        // intercept = sumz / n;
    } else {
        slope = (n * sumsz - sums * sumz) / denominator;
        // intercept = (sumz - slope * sums) / n;
    }

    double pZ = slope * init_pT;
    init_mom.SetZ(pZ);


    // charge
    double x0_ = circle.x0();
    double y0_ = circle.y0();
    TVector3 B_field = TVector3(0., 0., 2.);

    TVector3 pos(closestPoint.x, closestPoint.y, 0);
    TVector3 center(x0_, y0_, 0);

    TVector3 curvDir = center - pos;
    TVector3 lorentzDir = init_mom.Cross(B_field);
    lorentzDir.SetZ(0);
    curvDir.SetZ(0);

    int charge = (lorentzDir.Dot(curvDir) > 0) ? +1 : -1;


    return TrackInitialState{TVector3(closestPoint.x, closestPoint.y, points[0].z), init_mom, charge};

};

struct GenfitTrackFitter final : 
        k4FWCore::MultiTransformer< std::tuple< edm4hep::TrackCollection>(const edm4hep::TrackCollection&)>                                                                         
{
    GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("InputTracks", {"InputTracks"})
                
            },
            {   

                KeyValues("OutputFittedTracks", {"Fitted_tracks"})

            }) {}
    
            
    StatusCode initialize() {     
        
        if (!gGeoManager) {
            std::cerr << "Error: TGeoManager is not initialized!" << std::endl;
            return StatusCode::FAILURE;
        }
        
        // Initialize the Genfit: FieldManager and MaterialEffects
        m_detector = m_geoSvc->getDetector();
        m_field = m_detector->field();
        m_genfitField=new GenfitInterface::GenfitField(m_field);

        m_fieldManager = genfit::FieldManager::getInstance();
        m_fieldManager->init(m_genfitField);

        m_geoMaterial=GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);      
        // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
        // genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
        // genfit::MaterialEffects::getInstance()->setMscModel("Highland");

        genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
        genfit::MaterialEffects::getInstance()->setNoiseBrems(false);   

        int debug_lvl_material = 0;
        if (m_debug_lvl > 2)
        {
            debug_lvl_material = 2;
        }
        genfit::MaterialEffects::getInstance()->setDebugLvl(debug_lvl_material);         

        // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
        // These object are necessary to extract the drift chamber hits information, such as positions of the wire extremities
        m_surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();
        // If the detector doesn't have a drift chamber, this part will be skipped
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            m_dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            m_dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}

        // Retrive calorimeter information
        // These parameters are necessary to propagate the track to the calorimeter surface
        dd4hep::Detector& detector = dd4hep::Detector::getInstance();
        const std::vector< dd4hep::DetElement>& isDuaReadout = dd4hep::DetectorSelector(detector).detectors(  dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY );
        if( isDuaReadout.size() > 0 )
        {

            const dd4hep::rec::LayeredCalorimeterData * DualReadoutExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

            // mm (dd4hep::mm = 0.1)
            m_eCalBarrelInnerR  = DualReadoutExtension->extent[0] / dd4hep::mm;         // barrel rmin
            m_eCalBarrelMaxZ    = DualReadoutExtension->extent[2] / dd4hep::mm;         // barrel zmax == endcap zmin

            m_eCalEndCapInnerR  = DualReadoutExtension->extent[4] / dd4hep::mm;         // endcap rmin
            m_eCalEndCapOuterR  = DualReadoutExtension->extent[5] / dd4hep::mm;         // endcap rmax
            m_eCalEndCapInnerZ  = DualReadoutExtension->extent[2] / dd4hep::mm;         // endcap zmin
            m_eCalEndCapOuterZ  = DualReadoutExtension->extent[3] / dd4hep::mm;         // endcap zmax
        }
        else
        {
            const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

            // mm (dd4hep::mm = 0.1)                                                                                  
            m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
            m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
                    
            const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));


            // mm (dd4hep::mm = 0.1)                                                                                  
            m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
            m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
            m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
            m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
        }

        // N.B. we are assuming that the magnetic field is uniform and along the z direction
        // Z-component of the magnetic field at the center of the detector
        // This component is used to propagate the track to the calorimeter surface, but it is also necessary to
        // compute the pt component given omega
        dd4hep::Position center(0, 0, 0);
        dd4hep::Direction bfield = m_field.magneticField(center);
        m_Bz = bfield.z() / dd4hep::tesla;

        return StatusCode::SUCCESS;

    }
    
    std::tuple<edm4hep::TrackCollection> operator()( const edm4hep::TrackCollection& tracks_input) const override                                                 
    {
        
        // These 5 collections store the output of the fit for 5 particle hypotesis: e, mu, pi, K, p
        edm4hep::TrackCollection FittedTracks;

        info() << "Event number: " << event_counter++ << endmsg;
            
        // Loop over the tracks created by the pattern recognition step
        for (const auto& track : tracks_input)
        {

            num_tracks  +=1;

            // Skip background tracks if the option is enabled
            // Consider background tracks those with type = 0
            if (m_skip_background && track.getType() == 0) 
            {
                num_skip += 1;
                warning() << "Track " << num_tracks - 1<< ": background track (type = 0), skipping fit." << endmsg;
                warning() << "" << endmsg;
                continue;        // skip background        
            } 
            
            if (track.getTrackerHits().size() <= 3) 
            {   
                num_skip += 1;
                warning() << "Track " << num_tracks - 1 << ": less than 3 hits, skipping fit." << endmsg;
                warning() << "" << endmsg;
                continue;        // skip empty tracks and tracks with less then 3 hits (seed initialization needs 3 hits)
            }

            if (num_tracks > 3)
            {
                continue;
            }

            // Prepare the points for the extremum finding
            debug() << "Forward track " << num_tracks - 1 << " with " << track.getTrackerHits().size() << " hits summary:" << endmsg;
            auto xyzPoints = preparePoints(track, 1);

            int maxHitForLoopers = findFirstExtremum(xyzPoints, 0.01);
            if (maxHitForLoopers < 0) maxHitForLoopers = track.getTrackerHits().size();

            TrackInitialState init_seed_temp = inizialize_seed(xyzPoints, maxHitForLoopers);
            TVector3 init_mom_temp = init_seed_temp.init_mom;   // GeV/c

            if (init_mom_temp.Perp() > 1)
            {
                maxHitForLoopers = track.getTrackerHits().size();
            }

            // Compute seed for initial position and inizial momentum
            TrackInitialState init_seed = inizialize_seed(xyzPoints, maxHitForLoopers);
            TVector3 init_pos = init_seed.init_pos;             // mm
            TVector3 init_mom = init_seed.init_mom;             // GeV/c

            TVector3 init_pos_cm(
                init_seed.init_pos.X() * 0.1,
                init_seed.init_pos.Y() * 0.1,
                init_seed.init_pos.Z() * 0.1
            ); // cm
            TVector3 init_mom_gev = init_seed.init_mom;         // GeV/c

            // Summary
            debug() << "  Initial position: (" << init_pos.x() << ", " << init_pos.y() << ", " << init_pos.z() << ")" << endmsg;
            debug() << "  Initial momentum: (" << init_mom.x() << ", " << init_mom.y() << ", " << init_mom.z() << ")" << endmsg;
            debug() << "  Charge hypothesis: " << init_seed.charge << endmsg;
            debug() << "  Max hit for loopers: " << maxHitForLoopers << endmsg;
            debug() << "\n" << endmsg;

            num_processed_tracks += 1;

            if (m_singleEvaluation)
            {
                int pdgCode = m_particleHypotesis[0];

                int pdg_with_charge = pdgCode;
                int charge = 0;
                if (pdgCode == 11 || pdgCode == 13)
                {
                    pdg_with_charge = - init_seed.charge * pdgCode;
                    charge = - init_seed.charge;
                }
                else
                {   
                    pdg_with_charge = init_seed.charge * pdgCode;
                    charge = init_seed.charge;
                }

                GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(    track, 1,
                                                                                                pdg_with_charge, 
                                                                                                m_dch_info, m_dc_decoder, 
                                                                                                maxHitForLoopers, 
                                                                                                init_pos_cm, init_mom_gev); 
                
                track_interface.createGenFitTrack(true, m_debug_lvl);
                bool isFit = track_interface.fit(charge, m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl);     
                
                if (isFit)
                {
                    auto edm4hep_track = track_interface.getTrack_edm4hep();

                    // Propagate to calorimeter and store the state at calorimeter
                    FillTrackWithCalorimeterExtrapolation(edm4hep_track, m_Bz, charge, m_eCalBarrelInnerR, m_eCalBarrelMaxZ,
                                                                                m_eCalEndCapInnerR, m_eCalEndCapOuterR,
                                                                                m_eCalEndCapInnerZ, m_eCalEndCapOuterZ);
                   

                    // Add the fitted track to the output collection
                    FittedTracks.push_back(edm4hep_track);

                }
                else
                {

                    debug() << "Track " << num_tracks - 1 << ": fit failed for single evaluation hypothesis, skipping track." << endmsg;
                    number_failures += 1;

                    auto failedTrack = FittedTracks.create();
                    failedTrack.setChi2(-1);
                    failedTrack.setNdf(-1);
                    continue;
                }

            }

            int winning_hypothesis = -1; double winning_chi2_ndf = std::numeric_limits<double>::max();
            for (int pdgCode : m_particleHypotesis)
            {

                //Create trackInterface, initialize genfit track and fit it
                int pdg_with_charge = pdgCode;
                int charge = 0;
                if (pdgCode == 11 || pdgCode == 13)
                {
                    pdg_with_charge = - init_seed.charge * pdgCode;
                    charge = - init_seed.charge;
                }
                else
                {   
                    pdg_with_charge = init_seed.charge * pdgCode;
                    charge = init_seed.charge;
                }

                GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(    track, 1,
                                                                                                pdg_with_charge, 
                                                                                                m_dch_info, m_dc_decoder, 
                                                                                                maxHitForLoopers, 
                                                                                                init_pos_cm, init_mom_gev); 

                track_interface.createGenFitTrack(true, 0);
                bool isFit = track_interface.fit(charge, m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, 0); 
                
                
                if (isFit)
                {

                    auto edm4hep_track = track_interface.getTrack_edm4hep();
                    float genfit_chi2_val = edm4hep_track.getChi2();
                    int genfit_ndf_val = edm4hep_track.getNdf();

                    if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0)  continue; // skip invalid fits
                    
                    double chi2_ndf = genfit_chi2_val / genfit_ndf_val;
                    if (chi2_ndf < winning_chi2_ndf)
                    {
                        winning_chi2_ndf = chi2_ndf;
                        winning_hypothesis = pdgCode;
                    }
                
                }

            }   
            
            if (winning_hypothesis == -1)
            {
                debug() << "Track " << num_tracks - 1 << ": fit failed for all hypotheses, skipping track." << endmsg;
                number_failures += 1;

                auto failedTrack = FittedTracks.create();
                failedTrack.setChi2(-1);
                failedTrack.setNdf(-1);
            }
            else
            {

                debug() << "Track " << num_tracks - 1 << ": winning hypothesis is " << winning_hypothesis << " with chi2/ndf = " << winning_chi2_ndf << endmsg;

                // Create trackInterface, initialize genfit track and fit it
                int pdg_with_charge = winning_hypothesis;
                int charge = 0;
                if (winning_hypothesis == 11 || winning_hypothesis == 13)
                {
                    pdg_with_charge = - init_seed.charge * winning_hypothesis;
                    charge = - init_seed.charge;
                }
                else
                {   
                    pdg_with_charge = init_seed.charge * winning_hypothesis;
                    charge = init_seed.charge;
                }

                GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(    track, 1,
                                                                                                pdg_with_charge, 
                                                                                                m_dch_info, m_dc_decoder, 
                                                                                                maxHitForLoopers, 
                                                                                                init_pos_cm, init_mom_gev); 

                track_interface.createGenFitTrack(true, m_debug_lvl);
                bool _ = track_interface.fit(charge, m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl); 
                
                auto edm4hep_track = track_interface.getTrack_edm4hep();

                // Propagate to calorimeter and store the state at calorimeter
                FillTrackWithCalorimeterExtrapolation(edm4hep_track, m_Bz, charge, m_eCalBarrelInnerR, m_eCalBarrelMaxZ,
                                                                            m_eCalEndCapInnerR, m_eCalEndCapOuterR,
                                                                            m_eCalEndCapInnerZ, m_eCalEndCapOuterZ);
               

                // Add the fitted track to the output collection
                FittedTracks.push_back(edm4hep_track);

            }
        }

        
        return std::make_tuple( std::move(FittedTracks));
        
    } 
    
    StatusCode finalize() {     
        
        info() << "Run report:" << endmsg;
        info() << "Number of tracks: " << num_tracks << endmsg;
        info() << "Number of successes: " << (num_processed_tracks - number_failures) <<  "/" << num_processed_tracks << endmsg;
        info() << "Number of skipped tracks: " << num_skip << "/" << num_tracks << endmsg;
        info() << "----------------\n" << endmsg;

        // std::string filename = "/afs/cern.ch/work/a/adevita/public/workDir/testFolder/output_preparePoints_track_" + std::to_string(num_tracks) + "_forward.txt";
        // std::ofstream out(filename);
        // if (!out) {
        //     std::cerr << "Error opening file: " << filename << std::endl;
        // }

        // if (!out) {
        //     throw std::runtime_error("Cannot open output file");
        // }

        // // One point per line: z y
        // for (const auto& p : zy_pos) {
        //     out << p.z << " " << p.y << '\n';
        // }
        // out.close();



        return StatusCode::SUCCESS;

    }

    public:
        
        mutable int event_counter = 0;

        // Num_tracks = num_processed_tracks + num_skip
        mutable int num_tracks = 0;             // Total number of tracks
        mutable int num_skip = 0;               // Number of skipped tracks
        mutable int num_processed_tracks = 0;   // Number of tracks that have been processed (not skipped)
        mutable int number_failures = 0;        // Number of failed fits

        // mutable std::vector<Point2D_zy> zy_pos;


    private:

        

        ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
        dd4hep::Detector* m_detector{nullptr};  // Detector instance
        dd4hep::OverlayedField m_field;         // Magnetic field

        GenfitInterface::GenfitField* m_genfitField;
        genfit::FieldManager* m_fieldManager;

        GenfitInterface::GenfitMaterialInterface* m_geoMaterial;
        genfit::MaterialEffects* m_materialEffects;

        dd4hep::rec::SurfaceManager* m_surfMan;
        dd4hep::rec::DCH_info* m_dch_info;
        dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder;

        double m_Bz;

        double m_eCalBarrelInnerR;
        double m_eCalBarrelMaxZ;
        double m_eCalEndCapInnerR;
        double m_eCalEndCapOuterR;
        double m_eCalEndCapInnerZ;
        double m_eCalEndCapOuterZ;

        std::vector<int> m_particleHypotesis = {211};   // {11,13,211,321,2212} -> e, mu, pi, K, p

        Gaudi::Property<double> m_Beta_init{this, "Beta_init", 100, "Beta Initial value"};
        Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.05, "Beta Final value"};
        Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 15, "Beta number of Steps"};


        Gaudi::Property<bool> m_skip_background{this, "skip_background", true, "Skip background track"};
        Gaudi::Property<bool> m_singleEvaluation{this, "single_evaluation", false, "Single evaluation mode"};
        


        Gaudi::Property<int> m_debug_lvl{this, "debug_lvl", 0, "Debug level: Genfit - 0: none, 1: hit info, 2: + fitter info, 3: + material effects info"}; 

};

DECLARE_COMPONENT(GenfitTrackFitter)
