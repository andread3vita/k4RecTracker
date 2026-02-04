// /*
//  * Copyright (c) 2014-2024 Key4hep-Project.
//  *
//  * This file is part of Key4hep.
//  * See https://key4hep.github.io/key4hep-doc/ for further info.
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  *     http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
// */

// //=== Standard Library ===
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>  
// #include <fstream>  
// #include <filesystem>  
// #include <iostream>
// #include <iterator> 
// #include <map>
// #include <memory>
// #include <numeric>
// #include <queue>
// #include <sstream>
// #include <string>
// #include <typeinfo>
// #include <vector>

// //=== ROOT / C++ ABI ===
// #include <TEveManager.h>
// #include <TGeoManager.h>
// #include <TVector3.h>
// #include <cxxabi.h>
// #include <Eigen/Dense>

// //=== GenFit ===
// #include <ConstField.h>
// #include <DAF.h>
// #include <EventDisplay.h>
// #include <Exception.h>
// #include <FieldManager.h>
// #include <KalmanFitterInfo.h>
// #include <KalmanFitterRefTrack.h>
// #include <MaterialEffects.h>
// #include <PlanarMeasurement.h>
// #include <RKTrackRep.h>
// #include <StateOnPlane.h>
// #include <Track.h>
// #include <TrackPoint.h>
// #include <TGeoMaterialInterface.h>

// //=== Gaudi / k4FWCore ===
// #include "Gaudi/Property.h"
// #include "k4FWCore/DataHandle.h"
// #include "k4FWCore/Transformer.h"

// //=== DD4hep / DDRec / DDSegmentation ===
// #include "DD4hep/DD4hepUnits.h"
// #include "DD4hep/DetType.h"
// #include "DD4hep/DetElement.h"
// #include "DD4hep/Detector.h"
// #include "DD4hep/DetectorSelector.h"
// #include "DD4hep/Readout.h"
// #include <DDRec/DetectorData.h>
// #include <DDRec/Vector3D.h>
// #include <DDSegmentation/BitFieldCoder.h>
// #include "DD4hep/Fields.h"
// #include "DDRec/SurfaceManager.h"
// #include "DDRec/MaterialManager.h"

// //=== podio / edm4hep ===
// #include "podio/Frame.h"
// #include "podio/ROOTReader.h"
// #include "edm4hep/TrackerHitPlaneCollection.h"
// #include "edm4hep/SenseWireHitCollection.h"
// #include "edm4hep/TrackCollection.h"
// #include "edm4hep/MCParticleCollection.h"
// #include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
// #include "edm4hep/TrackState.h"
// #include "edm4hep/TrackMCParticleLinkCollection.h"

// //=== k4Interface ===
// #include "k4Interface/IGeoSvc.h"
// #include "k4Interface/IUniqueIDGenSvc.h"

// //=== GenfitInterface ===
// #include "GenfitWireMeasurement.hpp"
// #include "GenfitPlanarMeasurement.hpp"
// #include "GenfitField.hpp"
// #include "GenfitTrack.hpp"
// #include "GenfitMaterialInterface.hpp"

// //=== Others ===
// #include <marlinutil/HelixClass_double.h>
// #include <Objects/Helix.h>

// // Define collection types
// #include "podio/UserDataCollection.h"

// #include "utils.hpp"
// #include "FastCircleSeed.hpp"

// /** @struct GenfitTrackFitter
//  *
//  *  Gaudi MultiTransformer that refines the parameters of the reconstructed tracks using the GENFIT library.  
//  *  For each input track, the module performs a fit under five different particle hypotheses (electron, muon, pion, kaon, proton), 
//  *  producing five distinct output collections of fitted tracks.  
//  *  The fitting process minimizes the track χ² while accounting for material effects, magnetic field, and detector geometry.
//  *
//  *  If the fit fails for a given hypothesis, a fallback track is still generated with χ² = -1 and ndf = -1 to preserve collection integrity.
//  *  Each fitted track also includes the extrapolated state at the calorimeter (barrel or endcap).
//  *
//  *  input:
//  *    - initial track collection : edm4hep::TrackCollection
//  *
//  *  output:
//  *    - fitted track collection under electron hypothesis : edm4hep::TrackCollection
//  *    - fitted track collection under muon hypothesis : edm4hep::TrackCollection
//  *    - fitted track collection under pion hypothesis : edm4hep::TrackCollection
//  *    - fitted track collection under kaon hypothesis : edm4hep::TrackCollection
//  *    - fitted track collection under proton hypothesis : edm4hep::TrackCollection
//  *
//  *
//  *  @author Andrea De Vita
//  *  @date   2025-11
//  *
// */

// struct Point3D
// {
//     double x;
//     double y;
//     double z;

// };

// struct TrackInitialState
// {

//     TVector3 init_pos;
//     TVector3 init_mom;
//     int charge;
   
// };

// std::vector<double> smooth(const std::vector<double>& v, int window = 5) {
    
//     std::vector<double> out(v.size()); 
//     int W = window / 2;
//     for (size_t i = 0; i < v.size(); ++i) {
//         double sum = 0;
//         int count = 0;
//         for (int k = -W; k <= W; ++k) {
//             int j = i + k;
//             if (j >= 0 && j < (int)v.size()) {
//                 sum += v[j];
//                 count++;
//             }
//         }
//         out[i] = sum / count;
//     }
//     return out;
// }

// int findFirstExtremum(const std::vector<Point2D>& p,
//                       double epsilon = 0.05,
//                       int smoothWindow = 5)
// {
//     int n = p.size();
//     if (n < smoothWindow) return -1;

//     std::vector<double> y(n);
//     for (int i = 0; i < n; ++i) y[i] = p[i].y;

//     std::vector<double> ys = smooth(y, smoothWindow);

//     // std::cout << "Smoothed x y values:\n";
//     // for (auto k = 0; k < n; k++)
//     // {
//     //     std::cout << "Points: " << p[k].x << " " << ys[k] << std::endl;
//     // }


//     std::vector<double> d(n, 0.0);
//     for (int i = 1; i < n - 1; ++i) {
//         double dz = p[i+1].x - p[i-1].x;
//         if (dz != 0)
//         {
//             d[i] = (ys[i+1] - ys[i-1]) / dz;
//             // std::cout << "d[" << i << "] = " << d[i] << std::endl;
//         }
//     }

//     for (int i = 1; i < n - 1; ++i) {

//         if (d[i-1] >  epsilon && d[i] < -epsilon)
//             return i;

//         if (d[i-1] < -epsilon && d[i] >  epsilon)
//             return i;
//     }

//     return -1;
// }

// int findLooperLimit(const std::vector<Point2D>& p)
// {
//     if (p.size() < 3)
//         return -1;

//     constexpr int W = 1;

//     for (size_t i = W; i + W < p.size(); ++i) {

//         // direzione smooth
//         double vx = p[i+W].x - p[i-W].x;
//         double vy = p[i+W].y - p[i-W].y;

//         // vettore radiale
//         double rx = p[i].x;
//         double ry = p[i].y;

//         double dot = vx * rx + vy * ry;

//         if (dot < 0.0) {
//             return static_cast<int>(i);
//         }
//     }

//     return -1;
// }

// const std::vector<Point2D> preparePoints(const edm4hep::Track& track, int dir = 1)
// {

//     // // Sort the hits by distance from the origin
//     // std::vector<std::pair<float, int>> hitDistIndices{};
//     // int index = 0;
//     // auto hits_in_track = track.getTrackerHits();
//     // for (auto hit : hits_in_track) {

//     //     const auto pos_siHit = hit.getPosition();
//     //     const auto distance = std::sqrt(pos_siHit.x * pos_siHit.x + pos_siHit.y * pos_siHit.y + pos_siHit.z * pos_siHit.z);
//     //     hitDistIndices.emplace_back(distance, index++);

//     // }

//     // std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);

//     auto hits_in_track = track.getTrackerHits();
//     double z_max = -1e9; double x_withZmax = 0.0; double y_withZmax = 0.0;
//     double z_min = 1e9; double x_withZmin = 0.0; double y_withZmin = 0.0;

//     auto index = 0;
//     for (auto hit : hits_in_track) {

//         const auto pos_Hit = hit.getPosition();
//         if (pos_Hit.z > z_max) {
            
//             z_max = pos_Hit.z;
//             x_withZmax = pos_Hit.x;
//             y_withZmax = pos_Hit.y;
//         }
//         if (pos_Hit.z < z_min) {
                
//             z_min = pos_Hit.z;
//             x_withZmin = pos_Hit.x;
//             y_withZmin = pos_Hit.y;

//         }
//         index++;

//     }


//     double first_x = 0;
//     double first_y = 0;
//     double first_z = 0;

//     // Take x y z correspoding to the Z which is closer to zero
//     if (std::abs(z_min) < std::abs(z_max)) {
        
//         if (dir > 0) {
//             first_x = x_withZmin;
//             first_y = y_withZmin;
//             first_z = z_min;
//         }
//         else {
//             first_x = x_withZmax;
//             first_y = y_withZmax;
//             first_z = z_max;
//         }
        
//     } else {
        
//         if (dir > 0) {
//             first_x = x_withZmax;
//             first_y = y_withZmax;
//             first_z = z_max;
//         }
//         else {
//             first_x = x_withZmin;
//             first_y = y_withZmin;
//             first_z = z_min;
//         }
        
//     }

//     double delta_y = std::abs(y_withZmax - y_withZmin);
//     double delta_z = std::abs(z_max - z_min);
//     double cos_theta = delta_z / std::sqrt(std::pow(delta_y,2) + std::pow(delta_z,2));

//     if ( std::abs(cos_theta) < 0.01 ) {

//         double z_min_R = std::sqrt( x_withZmin * x_withZmin + y_withZmin * y_withZmin + z_min * z_min);
//         double z_max_R = std::sqrt( x_withZmax * x_withZmax + y_withZmax * y_withZmax + z_max * z_max);

//         if (z_min_R < z_max_R)
//         {   
//             if (dir > 0) {
//                 first_x = x_withZmin;
//                 first_y = y_withZmin;
//                 first_z = z_min;
//             }
//             else
//             {
//                 first_x = x_withZmax;
//                 first_y = y_withZmax;
//                 first_z = z_max;
//             }
//         }
//         else
//         {   
//             if (dir > 0) {
//                 first_x = x_withZmax;
//                 first_y = y_withZmax;
//                 first_z = z_max;
//             }
//             else
//             {
//                 first_x = x_withZmin;
//                 first_y = y_withZmin;
//                 first_z = z_min;
//             }
//         }
            
//     }


//     std::vector<std::pair<float, int>> hitDistIndices{};
//     index = 0;
//     for (auto hit : hits_in_track) {

//         const auto pos_siHit = hit.getPosition();
//         const auto distance = std::sqrt(std::pow(pos_siHit.x - first_x, 2) + std::pow(pos_siHit.y - first_y, 2) + std::pow(pos_siHit.z - first_z, 2));
//         hitDistIndices.emplace_back(distance, index++);

//     }
        
//     // Fill the internal track with sorted hits
//     std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);
//     // std::ranges::reverse(hitDistIndices);

//     // limit track to one round (for loopers)
//     std::vector<Point2D> zy_pos;
//     for (const auto& [_, idx] : hitDistIndices) {

//         auto hit = hits_in_track[idx];
//         const auto pos_siHit = hit.getPosition();

//         Point2D posForMax;
//         posForMax.x = pos_siHit.z;
//         posForMax.y = pos_siHit.y;

//         // std::cout << pos_siHit.z << " " << pos_siHit.y << std::endl;

//         zy_pos.push_back(posForMax);

//     }

//     return zy_pos;

// }

// TrackInitialState inizialize_seed(edm4hep::Track track, int dir, int maxHit = 50, double B = 2.0)
// {

//     // // Sort the hits by distance from the origin
//     // std::vector<std::pair<float, int>> hitDistIndices{};
//     // int index = 0;

//     // auto hits_in_track = track.getTrackerHits();
//     // for (auto hit : hits_in_track) {

//     //     const auto pos_siHit = hit.getPosition();
//     //     const auto distance = std::sqrt(pos_siHit.x * pos_siHit.x + pos_siHit.y * pos_siHit.y + pos_siHit.z * pos_siHit.z);
//     //     hitDistIndices.emplace_back(distance, index++);

//     // }

//     // std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);

//     auto hits_in_track = track.getTrackerHits();
//     double z_max = -1e9; double x_withZmax = 0.0; double y_withZmax = 0.0;
//     double z_min = 1e9; double x_withZmin = 0.0; double y_withZmin = 0.0;

//     auto index = 0;
//     for (auto hit : hits_in_track) {

//         const auto pos_Hit = hit.getPosition();
//         if (pos_Hit.z > z_max) {
            
//             z_max = pos_Hit.z;
//             x_withZmax = pos_Hit.x;
//             y_withZmax = pos_Hit.y;
//         }
//         if (pos_Hit.z < z_min) {
                
//             z_min = pos_Hit.z;
//             x_withZmin = pos_Hit.x;
//             y_withZmin = pos_Hit.y;

//         }
//         index++;

//     }


//     double first_x = 0;
//     double first_y = 0;
//     double first_z = 0;

//     // Take x y z correspoding to the Z which is closer to zero
//     if (std::abs(z_min) < std::abs(z_max)) {
        
//         if (dir > 0) {
//             first_x = x_withZmin;
//             first_y = y_withZmin;
//             first_z = z_min;
//         }
//         else {
//             first_x = x_withZmax;
//             first_y = y_withZmax;
//             first_z = z_max;
//         }
        
//     } else {
        
//         if (dir > 0) {
//             first_x = x_withZmax;
//             first_y = y_withZmax;
//             first_z = z_max;
//         }
//         else {
//             first_x = x_withZmin;
//             first_y = y_withZmin;
//             first_z = z_min;
//         }
        
//     }

//     double delta_y = std::abs(y_withZmax - y_withZmin);
//     double delta_z = std::abs(z_max - z_min);
//     double cos_theta = delta_z / std::sqrt(std::pow(delta_y,2) + std::pow(delta_z,2));

//     if ( std::abs(cos_theta) < 0.01 ) {

//         double z_min_R = std::sqrt( x_withZmin * x_withZmin + y_withZmin * y_withZmin + z_min * z_min);
//         double z_max_R = std::sqrt( x_withZmax * x_withZmax + y_withZmax * y_withZmax + z_max * z_max);

//         if (z_min_R < z_max_R)
//         {   
//             if (dir > 0) {
//                 first_x = x_withZmin;
//                 first_y = y_withZmin;
//                 first_z = z_min;
//             }
//             else
//             {
//                 first_x = x_withZmax;
//                 first_y = y_withZmax;
//                 first_z = z_max;
//             }
//         }
//         else
//         {   
//             if (dir > 0) {
//                 first_x = x_withZmax;
//                 first_y = y_withZmax;
//                 first_z = z_max;
//             }
//             else
//             {
//                 first_x = x_withZmin;
//                 first_y = y_withZmin;
//                 first_z = z_min;
//             }
//         }
            
//     }


//     std::vector<std::pair<float, int>> hitDistIndices{};
//     index = 0;
//     for (auto hit : hits_in_track) {

//         const auto pos_siHit = hit.getPosition();
//         const auto distance = std::sqrt(std::pow(pos_siHit.x - first_x, 2) + std::pow(pos_siHit.y - first_y, 2) + std::pow(pos_siHit.z - first_z, 2));
//         hitDistIndices.emplace_back(distance, index++);

//     }
        
//     // Fill the internal track with sorted hits
//     std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);

//     // limit track to one round (for loopers)
//     std::vector<Point3D> points;
//     std::vector<HitError> errors;
//     int idx_fill = 0;
//     for (const auto& [_, idx] : hitDistIndices) {

//         idx_fill +=1;
      

//         auto hit = hits_in_track[idx];
//         auto pos = hit.getPosition();

//         // if (idx_fill > 50) continue;
//         // std::cout << "Hit position: " << pos.x << " " << pos.y << " " << pos.z << std::endl;

//         if (idx_fill > maxHit) continue;

       
//         points.push_back(Point3D(pos.x, pos.y, pos.z));
        
//     }

//     std::vector<Point2D> points_xy;
//     for (const auto& p : points) {
//         points_xy.push_back(Point2D(p.x, p.y));
//     }

//     // FIT CIRCLE TO XY PROJECTION
//     FastCircleFit circle(points_xy);
//     Point2D closestPoint = circle.closestPointTo(points_xy[0]);
//     Point2D tangent_xy = circle.tangentAtPCA(closestPoint, points_xy[1]);

//     double rho = circle.rho();
//     double init_pT = std::abs(rho * 0.3 * B) / 1000; 

//     TVector3 init_mom = TVector3(tangent_xy.x * init_pT, tangent_xy.y * init_pT, 0);

//     // FIT LINE TO SZ PROJECTION: s_j = r | phi_j - phi_(j-1)|
//     double phi1 = std::atan2(closestPoint.y - circle.y0(), closestPoint.x - circle.x0());
//     double phi2 = std::atan2(points_xy[1].y - circle.y0(), points_xy[1].x - circle.x0());

//     if (phi2 - phi1 > M_PI) phi2 -= 2*M_PI;
//     if (phi2 - phi1 < - M_PI) phi2 += 2*M_PI;

//     double s1 = 0.0;
//     double s2 = rho * std::abs(phi2 - phi1);
//     double z1 = points[0].z;
//     double z2 = points[1].z;

//     double sums = s1 + s2;
//     double sumz = z1 + z2;
//     double sumsz = s1*z1 + s2*z2;
//     double sums2 = s1*s1 + s2*s2;
//     const int n = 2;
    
//     // Calculate slope and intercept
//     double denominator = n * sums2 - sums * sums;
//     double slope = 0.0;
//     // double intercept = 0.0;
//     if (std::abs(denominator) < 1e-6) {
        
//         slope = 1e6;
//         // intercept = sumz / n;
//     } else {
//         slope = (n * sumsz - sums * sumz) / denominator;
//         // intercept = (sumz - slope * sums) / n;
//     }

//     double pZ = slope * init_pT;
//     init_mom.SetZ(pZ);


//     // charge
//     double x0_ = circle.x0();
//     double y0_ = circle.y0();
//     TVector3 B_field = TVector3(0., 0., 2.);

//     TVector3 pos(closestPoint.x, closestPoint.y, 0);
//     TVector3 center(x0_, y0_, 0);

//     TVector3 curvDir = center - pos;
//     TVector3 lorentzDir = init_mom.Cross(B_field);
//     lorentzDir.SetZ(0);
//     curvDir.SetZ(0);

//     int charge = (lorentzDir.Dot(curvDir) > 0) ? +1 : -1;


//     return TrackInitialState{TVector3(closestPoint.x / 10, closestPoint.y / 10, points[0].z / 10), init_mom, charge};

// };

// struct GenfitTrackFitter final : 
//         k4FWCore::MultiTransformer< std::tuple< edm4hep::TrackCollection,
//                                                 edm4hep::TrackCollection,
//                                                 edm4hep::TrackCollection,
//                                                 edm4hep::TrackCollection,
//                                                 edm4hep::TrackCollection>(const edm4hep::TrackCollection&)>                                                                         
// {
//     GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : 
//         MultiTransformer ( name, svcLoc,
//             {   
                
//                 KeyValues("InputTracks", {"InputTracks"})
                
//             },
//             {   
//                 KeyValues("OutputFittedTracksElectronHypotesis", {"Fitted_tracks"}),
//                 KeyValues("OutputFittedTracksMuonHypotesis", {"Fitted_tracks"}),
//                 KeyValues("OutputFittedTracksPionHypotesis", {"Fitted_tracks"}),
//                 KeyValues("OutputFittedTracksKaonHypotesis", {"Fitted_tracks"}),
//                 KeyValues("OutputFittedTracksProtonHypotesis", {"Fitted_tracks"})
            
//             }) {}
    
            
//     StatusCode initialize() {     
        
//         if (!gGeoManager) {
//             std::cerr << "Error: TGeoManager is not initialized!" << std::endl;
//             return StatusCode::FAILURE;
//         }
        
//         // Initialize the Genfit: FieldManager and MaterialEffects
//         m_detector = m_geoSvc->getDetector();
//         m_field = m_detector->field();
//         m_genfitField=new GenfitInterface::GenfitField(m_field);

//         fieldManager = genfit::FieldManager::getInstance();
//         fieldManager->init(m_genfitField);

//         m_geoMaterial=GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);      
//         // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
//         // genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
//         // genfit::MaterialEffects::getInstance()->setMscModel("Highland");

//         genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
//         genfit::MaterialEffects::getInstance()->setNoiseBrems(false);   

//         genfit::MaterialEffects::getInstance()->setDebugLvl(m_debug_lvl);           

//         // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
//         // These object are necessary to extract the drift chamber hits information, such as positions of the wire extremities
//         surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();

//         // If the detector doesn't have a drift chamber, this part will be skipped
//         try {
//             std::string DCH_name("DCH_v2");
//             dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
//             dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
//             dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
//             dd4hep::Readout dch_readout = dch_sd.readout();
//             dc_decoder = dch_readout.idSpec().decoder();
//         } catch (const std::out_of_range& e) {}

//         // Z-component of the magnetic field at the center of the detector
//         // This component is used to propagate the track to the calorimeter surface, but it is also necessary to
//         // compute the pt component given omega
//         dd4hep::Position center(0, 0, 0);
//         dd4hep::Direction bfield = m_field.magneticField(center);
//         m_Bz = bfield.z() / dd4hep::tesla;

//         // Retrive calorimeter information
//         // These parameters are necessary to propagate the track to the calorimeter surface
//         dd4hep::Detector& detector = dd4hep::Detector::getInstance();
//         const std::vector< dd4hep::DetElement>& isDuaReadout = dd4hep::DetectorSelector(detector).detectors(  dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY );
//         if( isDuaReadout.size() > 0 )
//         {

//             const dd4hep::rec::LayeredCalorimeterData * DualReadoutExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

//             // mm
//             m_eCalBarrelInnerR = DualReadoutExtension->extent[0] / dd4hep::mm;     // barrel rmin
//             m_eCalBarrelMaxZ = DualReadoutExtension->extent[2] / dd4hep::mm;       // barrel zmax == endcap zmin

//             m_eCalEndCapInnerR = DualReadoutExtension->extent[4] / dd4hep::mm;     // endcap rmin
//             m_eCalEndCapOuterR = DualReadoutExtension->extent[5] / dd4hep::mm;     // endcap rmax
//             m_eCalEndCapInnerZ = DualReadoutExtension->extent[2] / dd4hep::mm;     // endcap zmin
//             m_eCalEndCapOuterZ = DualReadoutExtension->extent[3] / dd4hep::mm;     // endcap zmax
//         }
//         else
//         {
//             const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
//                                                                                               ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));
//             m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
//             m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
                    
//             const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
//                                                                                               ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));
//             m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
//             m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
//             m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
//             m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
//         }

//         // std::cout << m_eCalBarrelInnerR << " " << m_eCalBarrelMaxZ << " " << m_eCalEndCapInnerR << " " << m_eCalEndCapOuterR << " " << m_eCalEndCapInnerZ << " " << m_eCalEndCapOuterZ << std::endl;
//         return StatusCode::SUCCESS;

//     }
    
//     std::tuple< edm4hep::TrackCollection,
//                 edm4hep::TrackCollection,
//                 edm4hep::TrackCollection,
//                 edm4hep::TrackCollection,
//                 edm4hep::TrackCollection> operator()( const edm4hep::TrackCollection& tracks_input) const override                                                 
//     {
        
//         // These 5 collections store the output of the fit for 5 particle hypotesis: e, mu, pi, K, p
//         edm4hep::TrackCollection FittedTracks_electron;
//         edm4hep::TrackCollection FittedTracks_muon;
//         edm4hep::TrackCollection FittedTracks_pion;
//         edm4hep::TrackCollection FittedTracks_kaon;
//         edm4hep::TrackCollection FittedTracks_proton;

//         info() << "Event number: " << index_counter++ << endmsg;
            
//         // Loop over the tracks created by the pattern recognition step
//         for (const auto& track : tracks_input)
//         {
            
//             if (track.getType() == 0) 
//             {
//                 num_skip += 1;
//                 warning() << "Track " << num_tracks << ": background track (type = 0), skipping fit." << endmsg;
//                 warning() << "" << endmsg;
//                 continue;        // skip background        
//             }    

//             num_tracks  +=1;
            
//             if (track.getTrackerHits().size() <= 3) 
//             {   
//                 num_skip += 1;
//                 warning() << "Track " << num_tracks << ": less than 3 hits, skipping fit." << endmsg;
//                 warning() << "" << endmsg;
//                 continue;        // skip empty tracks and tracks with less then 3 hits (seed initialization needs 3 hits)
//             }

//             if (num_tracks != 1)
//             {
//                 continue;
//             }

//             /////////////////
//             //// FORWARD ////
//             /////////////////

//             // Prepare the points for the extremum finding
//             auto zy_pos = preparePoints(track, 1);
//             int maxHitForLoopers = findFirstExtremum(zy_pos, 0.01); // - 12;
//             if (maxHitForLoopers < 0) maxHitForLoopers = track.getTrackerHits().size();

//             TrackInitialState init_seed_temp = inizialize_seed(track, 1, maxHitForLoopers);
//             TVector3 init_pos_temp = init_seed_temp.init_pos;
//             TVector3 init_mom_temp = init_seed_temp.init_mom;

//             if (init_mom_temp.Perp() > 1)
//             {
//                 maxHitForLoopers = track.getTrackerHits().size();
//             }
//             info() << "Fitting forward track " << num_tracks << " with " << track.getTrackerHits().size() << " hits, maxHitForLoopers = " << maxHitForLoopers << endmsg;

//             // Compute seed for initial position and inizial momentum
//             TrackInitialState init_seed = inizialize_seed(track, 1, maxHitForLoopers);
//             TVector3 init_pos = init_seed.init_pos; // cm
//             TVector3 init_mom = init_seed.init_mom; // GeV/c

//             //////////////////
//             //// BACKWARD ////
//             //////////////////

//             // Prepare the points for the extremum finding
//             auto zy_pos_back = preparePoints(track, -1);
//             int maxHitForLoopers_back = findFirstExtremum(zy_pos_back, 0.01); // - 12;
//             if (maxHitForLoopers_back < 0) maxHitForLoopers_back = track.getTrackerHits().size();

//             TrackInitialState init_seed_back_temp = inizialize_seed(track, -1, maxHitForLoopers_back);
//             TVector3 init_pos_back_temp = init_seed_back_temp.init_pos;
//             TVector3 init_mom_back_temp = init_seed_back_temp.init_mom;

//             if (init_mom_back_temp.Perp() > 1)
//             {
//                 maxHitForLoopers_back = track.getTrackerHits().size();
//             }
//             info() << "Fitting backward track " << num_tracks << " with " << track.getTrackerHits().size() << " hits, maxHitForLoopers_back = " << maxHitForLoopers_back << endmsg;

//             // Compute seed for initial position and inizial momentum
//             TrackInitialState init_seed_back = inizialize_seed(track, -1, maxHitForLoopers_back);
//             TVector3 init_pos_back = init_seed_back.init_pos;
//             TVector3 init_mom_back = init_seed_back.init_mom;

//             // Loop over the particle hypotesis
//             num_processed_tracks += 1;
//             for (int pdgCode : m_particleHypotesis)
//             {

//                 // Create trackInterface, initialize genfit track and fit it
//                 int pdg_with_charge = pdgCode;
//                 if (pdgCode == 11 || pdgCode == 13)
//                 {
//                     pdg_with_charge = - init_seed.charge * pdgCode;
//                 }
//                 else
//                 {   
//                     pdg_with_charge = init_seed.charge * pdgCode;
//                 }
//                 int charge = getHypotesisCharge(pdg_with_charge);  

//                 GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(    track, 
//                                                                                                 pdg_with_charge, 
//                                                                                                 dch_info, dc_decoder, 
//                                                                                                 maxHitForLoopers, maxHitForLoopers_back, 
//                                                                                                 init_pos, init_mom, 
//                                                                                                 init_pos_back, init_mom_back); 

//                 track_interface.createGenFitTrack(1, m_debug_lvl, init_mom.Z(), init_mom.Perp());
//                 track_interface.createGenFitTrack(-1, m_debug_lvl, init_mom_back.Z(), init_mom_back.Perp());

//                 bool isFit = track_interface.fit(1, m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl); 
          
//                 if (isFit)
//                 {

//                     auto edm4hep_track = track_interface.getTrack_edm4hep();
//                     auto IP_trackstate = edm4hep_track.getTrackStates()[0];
                    
//                     double genfit_omega = IP_trackstate.omega;
//                     double genfit_pt_val = charge * m_Bz / abs(genfit_omega);
//                     double genfit_phi_val = IP_trackstate.phi;
//                     double genfit_pz = IP_trackstate.tanLambda * genfit_pt_val;
//                     double genfit_px = genfit_pt_val * std::cos(genfit_phi_val);
//                     double genfit_py = genfit_pt_val * std::sin(genfit_phi_val);
//                     double genfit_p = std::sqrt(genfit_px * genfit_px + genfit_py * genfit_py + genfit_pz * genfit_pz);
//                     TVector3 genfit_mom = TVector3(genfit_px, genfit_py, genfit_pz);

//                     float genfit_chi2_val = edm4hep_track.getChi2();
//                     int genfit_ndf_val = edm4hep_track.getNdf();

//                     // If the fitter suppresses the hits (ndf or chi2 <= 0) the fit is considered as failed
//                     if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0) {
                        
//                         number_failures++;
//                         if (pdgCode == 11)
//                         {
//                             auto failedTrack = FittedTracks_electron.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 13)
//                         {
//                             auto failedTrack = FittedTracks_muon.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);;
//                         }
//                         else if (pdgCode == 211)
//                         {
//                             auto failedTrack = FittedTracks_pion.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 321)
//                         {
//                             auto failedTrack = FittedTracks_kaon.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 2212)
//                         {
//                             auto failedTrack = FittedTracks_proton.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }

//                         debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                         debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                         debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;
//                         debug() << "" << endmsg;
//                         continue; 

//                     }

//                     double genfit_chi2_ndf_val = genfit_chi2_val / genfit_ndf_val;
//                     auto genfit_hits_in_track = edm4hep_track.getTrackerHits();

//                     debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                     debug() << "GENFIT Initialization: CMS" << endmsg;
                  
//                     track_interface.getTrack_init();

//                     debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                     debug() << "GENFIT Omega: " << genfit_omega << endmsg;
//                     debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
//                     debug() << "GENFIT tanLambda: " << IP_trackstate.tanLambda << endmsg;

//                     debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
//                     debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
//                     debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

//                     debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
//                     debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
//                     debug() << "" << endmsg;
                    
//                     debug() << "Track State at IP: " << endmsg;
//                     debug() << "  D0 [mm]: " << IP_trackstate.D0 << endmsg;
//                     debug() << "  Z0 [mm]: " << IP_trackstate.Z0 << endmsg;
//                     debug() << "  Phi [rad]: " << IP_trackstate.phi << endmsg;
//                     debug() << "  Omega [q/pt]: " << IP_trackstate.omega << endmsg;
//                     debug() << "  TanLambda: " << IP_trackstate.tanLambda << endmsg;
//                     debug() << "" << endmsg;

//                     // Propagation of the track to the calorimeter surface: add atCalorimeter trackState to the track         
//                     // FillTrackWithCalorimeterExtrapolation(
//                     //     edm4hep_track,
//                     //     m_Bz,
//                     //     charge,
//                     //     m_eCalBarrelInnerR,
//                     //     m_eCalBarrelMaxZ,
//                     //     m_eCalEndCapInnerR,
//                     //     m_eCalEndCapOuterR,
//                     //     m_eCalEndCapInnerZ
//                     // );
                    
//                     // Fill track collections
//                     if (std::abs(pdgCode) == 11) {
//                         FittedTracks_electron.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 13) {
//                         FittedTracks_muon.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 211) {
//                         FittedTracks_pion.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 321) {
//                         FittedTracks_kaon.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 2212) {
//                         FittedTracks_proton.push_back(edm4hep_track);
//                     }     

//                 } else {  
                    

//                     // If the fit fails, it returns a track with chi2=ndf=-1
//                     if (std::abs(pdgCode) == 11) {
//                         auto failedTrack = FittedTracks_electron.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 13) {
//                         auto failedTrack = FittedTracks_muon.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 211) {
//                         auto failedTrack = FittedTracks_pion.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 321) {
//                         auto failedTrack = FittedTracks_kaon.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 2212) {
//                         auto failedTrack = FittedTracks_proton.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     }

                        
//                     debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                     debug() << "GENFIT Initialization: CMS" << endmsg; 
//                     track_interface.getTrack_init();

//                     debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                     debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;
//                     debug() << "" << endmsg;

//                     number_failures++;

//                 }

//                 bool isFit_back = track_interface.fit(-1, m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl); 
//                 if (isFit_back)
//                 {

                    
//                     auto edm4hep_track = track_interface.getTrack_edm4hep(-1);
//                     auto genfit_trackstate = edm4hep_track.getTrackStates()[0];
                    
//                     double genfit_omega = genfit_trackstate.omega;
//                     double genfit_pt_val = charge * m_Bz / abs(genfit_omega);
//                     double genfit_phi_val = genfit_trackstate.phi;
//                     double genfit_pz = genfit_trackstate.tanLambda * genfit_pt_val;
//                     double genfit_px = genfit_pt_val * std::cos(genfit_phi_val);
//                     double genfit_py = genfit_pt_val * std::sin(genfit_phi_val);
//                     double genfit_p = std::sqrt(genfit_px * genfit_px + genfit_py * genfit_py + genfit_pz * genfit_pz);
//                     TVector3 genfit_mom = TVector3(genfit_px, genfit_py, genfit_pz);

//                     float genfit_chi2_val = edm4hep_track.getChi2();
//                     int genfit_ndf_val = edm4hep_track.getNdf();

//                     // If the fitter suppresses the hits (ndf or chi2 <= 0) the fit is considered as failed
//                     if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0) {
                        
//                         number_failures++;
//                         if (pdgCode == 11)
//                         {
//                             auto failedTrack = FittedTracks_electron.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 13)
//                         {
//                             auto failedTrack = FittedTracks_muon.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);;
//                         }
//                         else if (pdgCode == 211)
//                         {
//                             auto failedTrack = FittedTracks_pion.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 321)
//                         {
//                             auto failedTrack = FittedTracks_kaon.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }
//                         else if (pdgCode == 2212)
//                         {
//                             auto failedTrack = FittedTracks_proton.create();
//                             failedTrack.setChi2(-1);
//                             failedTrack.setNdf(-1);
//                         }

//                         debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                         debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                         debug() << "GENFIT Initialization: CMS" << endmsg;
//                         track_interface.getTrack_init(-1);
//                         debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;
//                         debug() << "" << endmsg;
//                         continue; 
//                     }

//                     double genfit_chi2_ndf_val = genfit_chi2_val / genfit_ndf_val;
//                     auto genfit_hits_in_track = edm4hep_track.getTrackerHits();

//                     debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                     debug() << "GENFIT Initialization: CMS" << endmsg;
                  
//                     track_interface.getTrack_init(-1);

//                     debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                     debug() << "GENFIT Omega: " << genfit_omega << endmsg;
//                     debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
//                     debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;

//                     debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
//                     debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
//                     debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

//                     debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
//                     debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
//                     debug() << "" << endmsg;

//                     // // Propagation of the track to the calorimeter surface: add atCalorimeter trackState to the track
// 		            // int charge = getHypotesisCharge(pdg_with_charge);                   
//                     // FillTrackWithCalorimeterExtrapolation(
//                     //     edm4hep_track,
//                     //     m_Bz,
//                     //     charge,
//                     //     a,
//                     //     m_eCalBarrelInnerR,
//                     //     m_eCalBarrelMaxZ,
//                     //     m_eCalEndCapInnerR,
//                     //     m_eCalEndCapOuterR,
//                     //     m_eCalEndCapInnerZ
//                     // );
                    
//                     // Fill track collections
//                     if (std::abs(pdgCode) == 11) {
//                         FittedTracks_electron.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 13) {
//                         FittedTracks_muon.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 211) {
//                         FittedTracks_pion.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 321) {
//                         FittedTracks_kaon.push_back(edm4hep_track);
//                     } else if (std::abs(pdgCode) == 2212) {
//                         FittedTracks_proton.push_back(edm4hep_track);
//                     }     

//                 } else {  

//                     // If the fit fails, it returns a track with chi2=ndf=-1
//                     if (std::abs(pdgCode) == 11) {

//                         auto failedTrack = FittedTracks_electron.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);

//                     } else if (std::abs(pdgCode) == 13) {

//                         auto failedTrack = FittedTracks_muon.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 211) {

//                         auto failedTrack = FittedTracks_pion.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 321) {

//                         auto failedTrack = FittedTracks_kaon.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     } else if (std::abs(pdgCode) == 2212) {

//                         auto failedTrack = FittedTracks_proton.create();
//                         failedTrack.setChi2(-1);
//                         failedTrack.setNdf(-1);
//                     }

                        
//                     debug() << "Number of hits in the track: " << track.getTrackerHits().size() << endmsg;
//                     debug() << "GENFIT Initialization: CMS" << endmsg; 
//                     track_interface.getTrack_init();

//                     debug() << "GENFIT PDG: " << pdg_with_charge << endmsg;
//                     debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;
//                     debug() << "" << endmsg;

//                     number_failures++;

//                 }
            
//             }
        
//         }

        
//         return std::make_tuple( std::move(FittedTracks_electron),
//                                 std::move(FittedTracks_muon),
//                                 std::move(FittedTracks_pion),
//                                 std::move(FittedTracks_kaon),
//                                 std::move(FittedTracks_proton));
        
//     } 
    
//     StatusCode finalize() {     
        
//         info() << "Run report:" << endmsg;
//         info() << "Number of tracks: " << num_tracks << endmsg;
//         info() << "Number of successes: " << (num_processed_tracks - number_failures) <<  "/" << num_processed_tracks << endmsg;
//         info() << "Number of skipped tracks: " << num_skip - 1 << "/" << num_tracks<< endmsg;
//         info() << "----------------\n" << endmsg;



//         return StatusCode::SUCCESS;

//     }

//     public:
        
//         mutable int index_counter = 0;
//         mutable int number_failures = 0;
//         mutable int num_tracks = 0;             // Total number of tracks processed
//         mutable int num_skip = 0;               // Number of skipped tracks
//         mutable int num_processed_tracks = 0;   // Number of tracks that have been processed (not skipped)


//     private:

//         ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
//         dd4hep::Detector* m_detector{nullptr};  // Detector instance
//         dd4hep::OverlayedField m_field;         // Magnetic field

//         GenfitInterface::GenfitField* m_genfitField;
//         genfit::FieldManager* fieldManager;

//         GenfitInterface::GenfitMaterialInterface* m_geoMaterial;
//         genfit::MaterialEffects* materialEffects;

//         dd4hep::rec::SurfaceManager* surfMan;
//         dd4hep::rec::DCH_info* dch_info;
//         dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;

//         double m_Bz;

//         double m_eCalBarrelInnerR;
//         double m_eCalBarrelMaxZ;
//         double m_eCalEndCapInnerR;
//         double m_eCalEndCapOuterR;
//         double m_eCalEndCapInnerZ;
//         double m_eCalEndCapOuterZ;

//         std::vector<int> m_particleHypotesis = {211}; // {11,13,211,321,2212} -> e, mu, pi, K, p

//         Gaudi::Property<double> m_Beta_init{this, "Beta_init", 100, "Beta Initial value"};
//         Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.05, "Beta Final value"};
//         Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 15, "Beta number of Steps"};
//         Gaudi::Property<int> m_debug_lvl{this, "debug_lvl", 0, "Debug level: Genfit"};

// };

// DECLARE_COMPONENT(GenfitTrackFitter)
