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
//  */

// #include <algorithm>
// #include <cmath>
// #include <iostream>
// #include <map>
// #include <memory> 
// #include <numeric>
// #include <queue>
// #include <sstream>
// #include <string>
// #include <typeinfo>
// #include <vector>
// #include <fstream>  // Per std::ifstream
// #include <vector>   // Per std::vector
// #include <iterator> // Per std::istreambuf_iterator
// #include <typeinfo>

// #include <iostream>
// #include <typeinfo>
// #include <cxxabi.h>
// #include <memory>

// #include <ATen/ATen.h>
// #include <torch/torch.h>
// #include "onnxruntime_cxx_api.h"
// #include "onnxruntime_run_options_config_keys.h"

// #include "dbscan.hpp"
// #include "utils.hpp"
// #include "SI_measurement.hpp"
// #include "CLDtrack.hpp"

// #include "Gaudi/Property.h"
// #include "k4FWCore/Transformer.h"

// #include "DD4hep/Detector.h"
// #include "DDRec/Vector3D.h"
// #include "DDSegmentation/BitFieldCoder.h"
// #include "k4FWCore/DataHandle.h"
// #include "k4Interface/IGeoSvc.h"
// #include <iostream>
// #include <string>
// #include <vector>
// #include <cstdlib>  // For getenv
// #include "k4Interface/IGeoSvc.h"
// #include "k4Interface/IUniqueIDGenSvc.h"
// #include <filesystem>  // For std::filesystem::path

// //genfit
// #include <ConstField.h>
// #include <Exception.h>
// #include <FieldManager.h>
// #include <KalmanFitterRefTrack.h>
// #include <DAF.h>
// #include <StateOnPlane.h>
// #include <Track.h>
// #include <TrackPoint.h>

// #include <MaterialEffects.h>
// #include <RKTrackRep.h>
// #include <TGeoMaterialInterface.h>

// #include <EventDisplay.h>

// #include <PlanarMeasurement.h>

// #include <TEveManager.h>
// #include <TGeoManager.h>
// #include <TVector3.h>
// #include <vector>

// #include "podio/ROOTReader.h"
// #include "podio/Frame.h"
// #include <KalmanFitterInfo.h>

// // Define collection types
// #include "podio/UserDataCollection.h"
// using DoubleColl = podio::UserDataCollection<double>;
// using IntColl = podio::UserDataCollection<int>;
// using FloatColl = podio::UserDataCollection<float>;


// #include "edm4hep/TrackCollection.h"
// #include "edm4hep/TrackerHitPlaneCollection.h"
// #include "edm4hep/ReconstructedParticleCollection.h"
// #include "edm4hep/TrackMCParticleLinkCollection.h"

// #include <marlinutil/HelixClass_double.h>
// #include <Objects/Helix.h>

// #include <DDRec/DetectorData.h>
// #include "DD4hep/Detector.h"
// #include "DD4hep/DD4hepUnits.h"
// #include "DD4hep/DetType.h"
// #include "DD4hep/DetectorSelector.h"
// #include "DD4hep/Readout.h"

// #include "utils.hpp"
// /** @struct CLD_fitter_validation
// *
// *
// */

// struct CLD_fitter_validation final : 
//         k4FWCore::MultiTransformer< std::tuple< edm4hep::TrackCollection,IntColl, DoubleColl, DoubleColl, DoubleColl, DoubleColl,IntColl,DoubleColl,DoubleColl,DoubleColl, DoubleColl,IntColl,DoubleColl,DoubleColl,DoubleColl>(const edm4hep::TrackMCParticleLinkCollection&)> 
            
                                                                                            
// {
//     CLD_fitter_validation(const std::string& name, ISvcLocator* svcLoc) : 
//         MultiTransformer ( name, svcLoc,
//             {   
                
//                 KeyValues("TrackMCLinks", {"TrackMCLinks"})
//             },
//             {   
//                 KeyValues("Tracks_Genfit", {"Tracks_Genfit"}),  
//                 KeyValues("t_pdg", {"t_pdg"}),
//                 KeyValues("t_pt", {"t_pt"}),
//                 KeyValues("t_costheta", {"t_costheta"}),
//                 KeyValues("t_phi", {"t_phi"}),

//                 KeyValues("cld_chi2", {"cld_chi2"}),
//                 KeyValues("cld_ndf", {"cld_ndf"}),
//                 KeyValues("cld_pt", {"cld_pt"}),
//                 KeyValues("cld_costheta", {"cld_costheta"}),
//                 KeyValues("cld_phi", {"cld_phi"}),

//                 KeyValues("genfit_chi2", {"genfit_chi2"}),
//                 KeyValues("genfit_ndf", {"genfit_ndf"}),
//                 KeyValues("genfit_pt", {"genfit_pt"}),
//                 KeyValues("genfit_costheta", {"genfit_costheta"}),
//                 KeyValues("genfit_phi", {"genfit_phi"})
            
//             }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
//     StatusCode initialize() {

//         // geoManager = new TGeoManager("Geometry", "CLD geometry");
//         // geoManager->Import(geoPath_m.value().c_str());

//         materialEffects = genfit::MaterialEffects::getInstance();
//         materialEffects->init(new genfit::TGeoMaterialInterface());

//         fieldManager = genfit::FieldManager::getInstance();
//         fieldManager->init(new genfit::ConstField(0., 0., m_Bz.value())); // kGauss

//         dd4hep::Detector* detector = (m_geoSvc->getDetector());
//         surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
//         // surfaceMap_vertex = surfMan->map("Vertex");
//         // surfaceMap_InnerTrackers = surfMan->map("InnerTrackers");
//         // surfaceMap_OuterTrackers = surfMan->map("OuterTrackers");

//         const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
//                                                                                           ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) , detector);
//         m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
//         m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
                  
//         const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
//                                                                                           ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) , detector);
//         m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
//         m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
//         m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
//         m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
                    
//         return StatusCode::SUCCESS;

//     }

    
//     std::tuple<edm4hep::TrackCollection,IntColl, DoubleColl, DoubleColl, DoubleColl, DoubleColl,IntColl,DoubleColl,DoubleColl,DoubleColl, DoubleColl,IntColl,DoubleColl,DoubleColl,DoubleColl> operator()( const edm4hep::TrackMCParticleLinkCollection& TrackMCLinks) const override 
//     {
        
//         new TGeoManager("Geometry", "CLD geometry");
//         TGeoManager::Import("/afs/cern.ch/user/a/adevita/public/workDir/test/fitting_validation/CLD_validation/TGeo_CLD.root");
        
//         std::cout << "m_eCalBarrelInnerR: " << m_eCalBarrelInnerR << std::endl;
//         std::cout << "m_eCalBarrelMaxZ: " << m_eCalBarrelMaxZ << std::endl;
//         std::cout << "m_eCalEndCapInnerR: " << m_eCalEndCapInnerR << std::endl;
//         std::cout << "m_eCalEndCapOuterR: " << m_eCalEndCapOuterR << std::endl;
//         std::cout << "m_eCalEndCapInnerZ: " << m_eCalEndCapInnerZ << std::endl;
//         std::cout << "m_eCalEndCapOuterZ: " << m_eCalEndCapOuterZ << std::endl;
       
//         IntColl t_pdg;
//         DoubleColl t_pt;
//         DoubleColl t_costheta;
//         DoubleColl t_phi;

//         DoubleColl cld_chi2;
//         IntColl cld_ndf;
//         DoubleColl cld_pt;
//         DoubleColl cld_costheta;
//         DoubleColl cld_phi;

//         DoubleColl genfit_chi2;
//         IntColl genfit_ndf;
//         DoubleColl genfit_pt;
//         DoubleColl genfit_costheta;
//         DoubleColl genfit_phi;
        
//         edm4hep::TrackCollection trackColl;
//         for (auto link : TrackMCLinks)
//         {
//             auto track = link.getFrom();
//             auto particle = link.getTo();

//             auto genStatus = particle.getGeneratorStatus();
//             auto vertex = particle.getVertex();
            
//             if (genStatus == 1 && vertex.x == 0. && vertex.y == 0. && vertex.z == 0.) // make sure that you are taking only pions with vertex at the IP
//             {
//                 auto p = particle.getMomentum();
//                 TVector3 mom(p.x, p.y, p.z);
//                 auto pt = mom.Pt();
//                 auto costheta = mom.CosTheta();
//                 auto phi = mom.Phi();
                
//                 debug() << "True pdg: " << particle.getPDG() << endmsg;
//                 debug() << "True Momentum: [" << mom.X() << " " << mom.Y() << " " << mom.Z() << " ]-> " << mom.Mag() << endmsg;
//                 debug() << "True pt: " << pt << endmsg;
//                 debug() << "True cos(theta): " << costheta << endmsg;
//                 debug() << "True phi: " << phi << endmsg;
//                 debug() << "" << endmsg;

//                 t_pdg.push_back(particle.getPDG());
//                 t_pt.push_back(pt);
//                 t_costheta.push_back(costheta);
//                 t_phi.push_back(phi);

//                 auto hits_in_track = track.getTrackerHits();

//                 ////////////////////////////////////////
//                 ////////// CLD RECONSTRUCTION //////////
//                 ////////////////////////////////////////

//                 double c_light = 2.99792458e8;
//                 double a = c_light * 1e3 * 1e-15;
//                 double Bz = 2.0;

//                 auto trackstate = track.getTrackStates()[0];
          
//                 double omega = trackstate.omega;
//                 double cld_pt_val = a * Bz / abs(omega);
//                 double cld_phi_val = trackstate.phi;
//                 double pz = trackstate.tanLambda * cld_pt_val;
//                 double px = cld_pt_val * std::cos(cld_phi_val);
//                 double py = cld_pt_val * std::sin(cld_phi_val);
//                 double cld_p = std::sqrt(px * px + py * py + pz * pz);
//                 TVector3 mom_CLD = TVector3(px, py, pz);

//                 float cld_chi2_val = track.getChi2();
//                 int cld_ndf_val = track.getNdf();
//                 double cld_chi2_ndf_val = (cld_ndf_val > 0 ? cld_chi2_val / cld_ndf_val : -1);

//                 debug() << "CLDreco Omega: " << omega << endmsg;
//                 debug() << "CLDreco Phi: " << cld_phi_val << endmsg;
//                 debug() << "CLDreco tanLambda: " << trackstate.tanLambda << endmsg;
                
//                 debug() << "CLDreco Momentum: " << "[ " << mom_CLD.X() << " " << mom_CLD.Y() << " " << mom_CLD.Z() << " ] -> " << cld_p << endmsg;
//                 debug() << "CLDreco pt: " << mom_CLD.Pt() << endmsg;
//                 debug() << "CLDreco cosTheta: " << mom_CLD.CosTheta() << endmsg;

//                 debug() << "CLD Chi2: " << cld_chi2_val << ", CLD NDF: " << cld_ndf_val << ", CLD Chi2/NDF: " << cld_chi2_ndf_val << endmsg;
//                 debug() << "Number of hits in CLD track: " << hits_in_track.size() << endmsg;
//                 debug() << "" << endmsg;
                

//                 cld_chi2.push_back(cld_chi2_val);
//                 cld_ndf.push_back(cld_ndf_val);
//                 cld_costheta.push_back(std::abs(pz / cld_p));
//                 cld_pt.push_back(cld_pt_val);
//                 cld_phi.push_back(cld_phi_val);

//                 //////////////////////////
//                 //////// GENFIT //////////
//                 //////////////////////////

//                 GENFIT::CLDtrack track_interface = GENFIT::CLDtrack(track,surfMan);
//                 track_interface.createGenFitTrack();

//                 bool isFit = track_interface.fit();

//                 if (isFit)
//                 {

//                     auto edm4hep_track = track_interface.getTrack_edm4hep();
//                     auto genfit_trackstate = edm4hep_track.getTrackStates()[0];
          
//                     double genfit_omega = genfit_trackstate.omega;
//                     double genfit_pt_val = a * Bz / abs(genfit_omega);
//                     double genfit_phi_val = genfit_trackstate.phi;
//                     double genfit_pz = genfit_trackstate.tanLambda * genfit_pt_val;
//                     double genfit_px = genfit_pt_val * std::cos(genfit_phi_val);
//                     double genfit_py = genfit_pt_val * std::sin(genfit_phi_val);
//                     double genfit_p = std::sqrt(genfit_px * genfit_px + genfit_py * genfit_py + genfit_pz * genfit_pz);
//                     TVector3 genfit_mom = TVector3(genfit_px, genfit_py, genfit_pz);

//                     float genfit_chi2_val = edm4hep_track.getChi2();
//                     int genfit_ndf_val = edm4hep_track.getNdf();
//                     double genfit_chi2_ndf_val = (genfit_ndf_val > 0 ? genfit_chi2_val / genfit_ndf_val : -1);

//                     auto genfit_hits_in_track = edm4hep_track.getTrackerHits();
//                     debug() << "GENFIT Omega: " << genfit_omega << endmsg;
//                     debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
//                     debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;
                    
//                     debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
//                     debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
//                     debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

//                     debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
//                     debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
//                     debug() << "" << endmsg;
                    
//                     trackColl.push_back(edm4hep_track);
//                     genfit_chi2.push_back(genfit_chi2_val);
//                     genfit_ndf.push_back(genfit_ndf_val);
//                     genfit_costheta.push_back(std::abs(genfit_pz / genfit_p));
//                     genfit_pt.push_back(genfit_mom.Pt());
//                     genfit_phi.push_back(genfit_phi_val);


//                     // trackState at Calo
//                     auto trackStateLastHit = edm4hep_track.getTrackStates()[2];

//                     double omega_lastHit = trackStateLastHit.omega;
//                     double pt_lasthit = a * Bz / abs(omega_lastHit);
//                     double phi_lasthit = trackStateLastHit.phi;
//                     double pz_lasthit = trackStateLastHit.tanLambda * pt_lasthit;
//                     double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
//                     double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
//                     auto ref_lastHit = trackStateLastHit.referencePoint;

//                     // produce new helix at last hit position
//                     double posAtLastHit[] = {ref_lastHit[0], ref_lastHit[1], ref_lastHit[2]};
//                     double momAtLastHit[] = {px_lasthit, py_lasthit, pz_lasthit};
//                     auto helixAtLastHit = HelixClass_double();
//                     helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, 1, 2.);

//                     // TrackState at Calorimeter
//                     if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {

//                         pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
//                         pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
//                         float minGenericTime(std::numeric_limits<float>::max());
            
//                         // create helix to project
//                         // rather than using parameters at production, better to use those from
//                         // last hit
//                         pandora::CartesianVector pos_lasthit(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
//                         pandora::CartesianVector mom_lasthit(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);
//                         const pandora::Helix helix(pos_lasthit, mom_lasthit, 1,2.);
//                         const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
//                         const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
            
//                         // First project to endcap
//                         pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
//                         bool hasEndCapProjection(false);
//                         if (m_eCalEndCapInnerR>0) {
//                             float genericTime(std::numeric_limits<float>::max());
//                             const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint,
//                                                                 endCapProjection, genericTime));
//                             float x = endCapProjection.GetX();
//                             float y = endCapProjection.GetY();
//                             float r = std::sqrt(x*x+y*y);
//                             if (
//                                 (pandora::STATUS_CODE_SUCCESS == statusCode) &&
//                                 (genericTime < minGenericTime) &&
//                                 (r >= m_eCalEndCapInnerR) &&
//                                 (r <= m_eCalEndCapOuterR)
//                             ) {
//                                 minGenericTime = genericTime;
//                                 bestECalProjection = endCapProjection;
//                                 hasEndCapProjection = true;
//                             }
//                         }
                        
                        
//                         // Then project to barrel surface(s), and keep projection
//                         // if extrapolation is within the z acceptance of the detector
//                         pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
//                         bool hasBarrelProjection = false;
//                         if (m_eCalBarrelInnerR>0) {
//                             float genericTime(std::numeric_limits<float>::max());
//                             const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint,
//                                                                 barrelProjection, genericTime));
//                             if (
//                                 (pandora::STATUS_CODE_SUCCESS == statusCode) &&
//                                 (std::fabs(barrelProjection.GetZ())<= m_eCalBarrelMaxZ)
//                             ) {
//                                 hasBarrelProjection = true;
//                                 if (genericTime < minGenericTime) {
//                                 minGenericTime = genericTime;
//                                 secondBestECalProjection = bestECalProjection;
//                                 bestECalProjection = barrelProjection;
//                                 }
//                                 else {
//                                 secondBestECalProjection = barrelProjection;
//                                 }
//                             }
//                         }
 
//                         // store extrapolation to calo
//                         // by default, store extrapolation with lower arrival time
//                         // get extrapolated position
//                         edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit);

//                         omega_lastHit = trackState_AtCalorimeter.omega;
//                         pt_lasthit = a * Bz / abs(omega_lastHit);
//                         phi_lasthit = trackState_AtCalorimeter.phi;
//                         pz_lasthit = trackState_AtCalorimeter.tanLambda * pt_lasthit;
//                         px_lasthit = pt_lasthit * std::cos(phi_lasthit);
//                         py_lasthit = pt_lasthit * std::sin(phi_lasthit);
//                         ref_lastHit = trackState_AtCalorimeter.referencePoint;

//                         // attach the TrackState to the track
//                         edm4hep_track.addToTrackStates(trackState_AtCalorimeter);

//                     }
                        
//                 }
//                 else
//                 {
//                     auto failedTrack = trackColl.create();
//                     failedTrack.setChi2(-1);
//                     failedTrack.setNdf(-1);

//                     genfit_chi2.push_back(-1);
//                     genfit_ndf.push_back(-1);
//                     genfit_costheta.push_back(-1);
//                     genfit_pt.push_back(-1);
//                     genfit_phi.push_back(-1);

//                     debug() << "Chi2: " << failedTrack.getChi2() << " NDF: " << failedTrack.getNdf() << endmsg;
//                 }
                

//                 debug() << "----------------\n" << endmsg;
//             }
            
//         }

//         return std::make_tuple(     std::move(trackColl),
//                                     std::move(t_pdg),
//                                     std::move(t_pt),
//                                     std::move(t_costheta),
//                                     std::move(t_phi),
//                                     std::move(cld_chi2),
//                                     std::move(cld_ndf),
//                                     std::move(cld_pt),
//                                     std::move(cld_costheta),
//                                     std::move(cld_phi),
//                                     std::move(genfit_chi2),
//                                     std::move(genfit_ndf),
//                                     std::move(genfit_pt),
//                                     std::move(genfit_costheta),
//                                     std::move(genfit_phi));
        
//     } 

//     private:

//         // TGeoManager* geoManager;
//         genfit::MaterialEffects* materialEffects;
//         genfit::FieldManager* fieldManager;

//         // Gaudi::Property<std::string> geoPath_m{this, "geoPath_m", "/afs/cern.ch/user/a/adevita/public/workDir/test/fitting_validation/CLD_validation/TGeo_CLD.root", "geoPath_m"};
//         Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
//         Gaudi::Property<double> m_Bz{this, "Bz", 20., "Bz value (kilogauss)"};

//         SmartIF<IGeoSvc> m_geoSvc;
//         const dd4hep::rec::SurfaceManager* surfMan;

//         double m_eCalBarrelInnerR = 0;
//         double m_eCalBarrelMaxZ = 0;

//         double m_eCalEndCapInnerR = 0;
//         double m_eCalEndCapOuterR = 0;
//         double m_eCalEndCapInnerZ = 0;
//         double m_eCalEndCapOuterZ = 0;

//         // const dd4hep::rec::SurfaceMap* surfaceMap_vertex;
//         // const dd4hep::rec::SurfaceMap* surfaceMap_InnerTrackers;
//         // const dd4hep::rec::SurfaceMap* surfaceMap_OuterTrackers;  

// };

// DECLARE_COMPONENT(CLD_fitter_validation)