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

//=== podio / edm4hep ===
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "Wire_measurement.hpp"
#include "Planar_measurement.hpp"
#include "GenfitField.hpp"
#include "GenfitTrack.hpp"
#include "GenfitMaterialInterface.hpp"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>

// Define collection types
#include "podio/UserDataCollection.h"

#include "utils.hpp"


/** @struct TrackFitter_Genfit

*
*
*/

struct TrackFitter_Genfit final : 
        k4FWCore::MultiTransformer< std::tuple<extension::TrackCollection>(  const extension::TrackCollection&,
                                                                             const edm4hep::MCParticleCollection&, 
                                                                             const podio::UserDataCollection<int>&)>                                                                         
{
    TrackFitter_Genfit(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("tracks_input", {"tracks_input"}),
                KeyValues("MCParticles", {"MCParticles"}),
                KeyValues("MCparticleIndex", {"MCparticleIndex"}),
            },
            {   
                KeyValues("Fitted_tracks", {"Fitted_tracks"})
            
            }) {}
         
    StatusCode initialize() {       
        
        // Initialize the Genfit
        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());
        
        m_detector = m_geoSvc->getDetector();
        m_field = m_detector->field();
        m_genfitField=new GenfitField(m_field);

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(m_genfitField);
    
        dd4hep::Position center(0, 0, 0);
        dd4hep::Direction bfield = m_field.magneticField(center);
        m_Bz = bfield.z() / dd4hep::tesla;

        
        // Initialize the surface manager
        surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();

        // Retrieve the drift chamber information and decoder
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}
        

        // TODO: Retrive calorimeter information
        dd4hep::Detector& detector = dd4hep::Detector::getInstance();
        const std::vector< dd4hep::DetElement>& isDuaReadout = dd4hep::DetectorSelector(detector).detectors(  dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY );
        if( isDuaReadout.size() > 0 )
        {

            const dd4hep::rec::LayeredCalorimeterData * DualReadoutExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

            m_eCalBarrelInnerR = DualReadoutExtension->extent[0] / dd4hep::mm;     // barrel rmin
            m_eCalBarrelMaxZ = DualReadoutExtension->extent[2] / dd4hep::mm;       // barrel zmax == endcap zmin

            m_eCalEndCapInnerR = DualReadoutExtension->extent[4] / dd4hep::mm;     // endcap rmin
            m_eCalEndCapOuterR = DualReadoutExtension->extent[5] / dd4hep::mm;     // endcap rmax
            m_eCalEndCapInnerZ = DualReadoutExtension->extent[2] / dd4hep::mm;     // endcap zmin
            m_eCalEndCapOuterZ = DualReadoutExtension->extent[3] / dd4hep::mm;     // endcap zmax
        }
        else
        {
            const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));
            m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
            m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
                    
            const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));
            m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
            m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
            m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
            m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
        }


        return StatusCode::SUCCESS;

    }

    
    std::tuple<extension::TrackCollection> operator()(  const extension::TrackCollection& tracks_input,
                                                        const edm4hep::MCParticleCollection& mcParticles,
                                                        const podio::UserDataCollection<int>& MCparticleIndex) const override                                                 
    {
        

        info() << "Event number: " << index_counter++ << endmsg;
        extension::TrackCollection FittedTracks;


        // Loop over the extension::tracks
        int track_idx = 0;
        for (const auto& track : tracks_input)
        {
            
            auto particle = mcParticles[MCparticleIndex[track_idx]];

            auto p = particle.getMomentum();
            TVector3 mom(p.x, p.y, p.z);
            auto pt = mom.Pt();
            auto costheta = mom.CosTheta();
            auto phi = mom.Phi();
                
            debug() << "True pdg: " << particle.getPDG() << endmsg;
            debug() << "True Momentum: [" << mom.X() << " " << mom.Y() << " " << mom.Z() << " ]-> " << mom.Mag() << endmsg;
            debug() << "True pt: " << pt << endmsg;
            debug() << "True cos(theta): " << costheta << endmsg;
            debug() << "True phi: " << phi << endmsg;
            debug() << "" << endmsg;

            track_idx += 1;

            for (auto partPDG : m_particleHypotesis)
            {
                GENFIT::GenfitTrack track_interface = GENFIT::GenfitTrack(track, dch_info, dc_decoder, partPDG);
                track_interface.createGenFitTrack();
                bool isFit = track_interface.fit(m_Beta_init,m_Beta_final,m_Beta_steps);  

                if (isFit)
                {

                    auto edm4hep_track = track_interface.getTrack_edm4hep();
                    auto genfit_trackstate = edm4hep_track.getTrackStates()[0];
                
                    double genfit_omega = genfit_trackstate.omega;
                    double genfit_pt_val = a * m_Bz / abs(genfit_omega);
                    double genfit_phi_val = genfit_trackstate.phi;
                    double genfit_pz = genfit_trackstate.tanLambda * genfit_pt_val;
                    double genfit_px = genfit_pt_val * std::cos(genfit_phi_val);
                    double genfit_py = genfit_pt_val * std::sin(genfit_phi_val);
                    double genfit_p = std::sqrt(genfit_px * genfit_px + genfit_py * genfit_py + genfit_pz * genfit_pz);
                    TVector3 genfit_mom = TVector3(genfit_px, genfit_py, genfit_pz);

                    float genfit_chi2_val = edm4hep_track.getChi2();
                    int genfit_ndf_val = edm4hep_track.getNdf();
                    double genfit_chi2_ndf_val = (genfit_ndf_val > 0 ? genfit_chi2_val / genfit_ndf_val : -1);

                    auto genfit_hits_in_track = edm4hep_track.getTrackerHits();
                    debug() << "GENFIT PDG: " << partPDG << endmsg;
                    debug() << "GENFIT Omega: " << genfit_omega << endmsg;
                    debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
                    debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;
                            
                    debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
                    debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
                    debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

                    debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
                    debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
                    debug() << "" << endmsg;
                            
                    // // TrackState at Calorimeter
                    // auto trackStateLastHit = edm4hep_track.getTrackStates()[2];

                    // double omega_lastHit = trackStateLastHit.omega;
                    // double pt_lasthit = a * m_Bz / abs(omega_lastHit);
                    // double phi_lasthit = trackStateLastHit.phi;
                    // double pz_lasthit = trackStateLastHit.tanLambda * pt_lasthit;
                    // double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
                    // double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
                    // auto ref_lastHit = trackStateLastHit.referencePoint;

                    // // produce new helix at last hit position
                    // double posAtLastHit[] = {ref_lastHit[0], ref_lastHit[1], ref_lastHit[2]};
                    // double momAtLastHit[] = {px_lasthit, py_lasthit, pz_lasthit};
                    // auto helixAtLastHit = HelixClass_double();
                    // helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, 1, m_Bz);

                    // // Propagation to Endcap
                    // if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {

                    //     pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
                    //     pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
                    //     float minGenericTime(std::numeric_limits<float>::max());
                    
                    //     // create helix to project
                    //     // rather than using parameters at production, better to use those from
                    //     // last hit
                    //     pandora::CartesianVector pos_lasthit(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
                    //     pandora::CartesianVector mom_lasthit(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);

                    //     const pandora::Helix helix(pos_lasthit, mom_lasthit, 1,m_Bz);
                    //     const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
                    //     const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
                    
                    //     // First project to endcap
                    //     pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
                    //     bool hasEndCapProjection(false);
                    //     if (m_eCalEndCapInnerR>0) {
                    //         float genericTime(std::numeric_limits<float>::max());
                    //         const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint, endCapProjection, genericTime));
                    //         float x = endCapProjection.GetX();
                    //         float y = endCapProjection.GetY();
                    //         float r = std::sqrt(x*x+y*y);
                    //         if (
                    //             (pandora::STATUS_CODE_SUCCESS == statusCode) &&
                    //             (genericTime < minGenericTime) &&
                    //             (r >= m_eCalEndCapInnerR) &&
                    //             (r <= m_eCalEndCapOuterR)
                    //         ) {
                    //             minGenericTime = genericTime;
                    //             bestECalProjection = endCapProjection;
                    //             hasEndCapProjection = true;
                    //         }
                    //     }
                                
                                
                    //     // Then project to barrel surface(s), and keep projection
                    //     // if extrapolation is within the z acceptance of the detector
                    //     pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
                    //     bool hasBarrelProjection = false;
                    //     if (m_eCalBarrelInnerR>0) {

                    //         float genericTime(std::numeric_limits<float>::max());
                    //         const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint,
                    //                                                 barrelProjection, genericTime));
                    //         if (
                    //             (pandora::STATUS_CODE_SUCCESS == statusCode) &&
                    //             (std::fabs(barrelProjection.GetZ())<= m_eCalBarrelMaxZ)
                    //         ) {
                    //             hasBarrelProjection = true;
                    //             if (genericTime < minGenericTime) {
                    //             minGenericTime = genericTime;
                    //             secondBestECalProjection = bestECalProjection;
                    //             bestECalProjection = barrelProjection;
                    //             }
                    //             else {
                    //             secondBestECalProjection = barrelProjection;
                    //             }
                    //         }
                    //     }
        
                    //     // store extrapolation to calo
                    //     // by default, store extrapolation with lower arrival time
                    //     // get extrapolated position
                    //     edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit,2.0);
                    //     omega_lastHit = trackState_AtCalorimeter.omega;
                    //     pt_lasthit = a * m_Bz / abs(omega_lastHit);
                    //     phi_lasthit = trackState_AtCalorimeter.phi;
                    //     pz_lasthit = trackState_AtCalorimeter.tanLambda * pt_lasthit;
                    //     px_lasthit = pt_lasthit * std::cos(phi_lasthit);
                    //     py_lasthit = pt_lasthit * std::sin(phi_lasthit);
                    //     ref_lastHit = trackState_AtCalorimeter.referencePoint;
                    //     // attach the TrackState to the track
                    //     edm4hep_track.addToTrackStates(trackState_AtCalorimeter);

                    // }
                
                    // FittedTracks.push_back(edm4hep_track);

                }
                else
                {
                        
                    // auto failedTrack = FittedTracks.create();
                    // failedTrack.setChi2(-1);
                    // failedTrack.setNdf(-1);

                    // debug() << "Chi2: " << failedTrack.getChi2() << " NDF: " << failedTrack.getNdf() << endmsg;
                    debug() << "Fit Failed!" << endmsg;
                }
                debug() << "----------------\n" << endmsg;

            }
            
 
        }

        return std::make_tuple(std::move(FittedTracks));
        
    } 
    
    public:
        mutable int index_counter = 0;

    private:

        mutable genfit::MaterialEffects* materialEffects;
        mutable genfit::FieldManager* fieldManager;

        ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
        mutable dd4hep::Detector* m_detector{nullptr};  // Detector instance
        mutable dd4hep::OverlayedField m_field;         // Magnetic field
        mutable GenfitField* m_genfitField;
        
        mutable dd4hep::rec::SurfaceManager* surfMan;
        mutable dd4hep::rec::DCH_info* dch_info;

        mutable dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;

        Gaudi::Property<double> m_Beta_init{this, "Beta_init", 10., "Beta Initial value"};
        Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.1, "Beta Final value"};
        Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 10, "Beta number of Steps"};

        mutable double m_eCalBarrelInnerR = 0;
        mutable double m_eCalBarrelMaxZ = 0;

        mutable double m_eCalEndCapInnerR = 0;
        mutable double m_eCalEndCapOuterR = 0;
        mutable double m_eCalEndCapInnerZ = 0;
        mutable double m_eCalEndCapOuterZ = 0;

        double c_light = 2.99792458e8;
        double a = c_light * 1e3 * 1e-15;
        mutable double m_Bz = 2.0; // T

        // std::vector<int> m_particleHypotesis = {11,-11,13,-13,211,-211, 321,-321, 2212, -2212}; // e, mu, pi, K, p
        std::vector<int> m_particleHypotesis = {221}; // pi+


};

DECLARE_COMPONENT(TrackFitter_Genfit)
