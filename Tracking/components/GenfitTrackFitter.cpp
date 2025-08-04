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
#include "DDRec/SurfaceManager.h"
#include "DDRec/MaterialManager.h"

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

#include "utils.h"


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
 *    - initial track collection : extension::TrackCollection
 *
 *  output:
 *    - fitted track collection under electron hypothesis : extension::TrackCollection
 *    - fitted track collection under muon hypothesis : extension::TrackCollection
 *    - fitted track collection under pion hypothesis : extension::TrackCollection
 *    - fitted track collection under kaon hypothesis : extension::TrackCollection
 *    - fitted track collection under proton hypothesis : extension::TrackCollection
 *
 *
 *
 *  @author Andrea De Vita
 *  @date   2025-06
 *
 */


struct GenfitTrackFitter final : 
        k4FWCore::MultiTransformer< std::tuple< extension::TrackCollection,
                                                extension::TrackCollection,
                                                extension::TrackCollection,
                                                extension::TrackCollection,
                                                extension::TrackCollection>(const extension::TrackCollection&)>                                                                         
{
    GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("tracks_input", {"tracks_input"})
                
            },
            {   
                KeyValues("Fitted_tracks_electron", {"Fitted_tracks"}),
                KeyValues("Fitted_tracks_muon", {"Fitted_tracks"}),
                KeyValues("Fitted_tracks_pion", {"Fitted_tracks"}),
                KeyValues("Fitted_tracks_kaon", {"Fitted_tracks"}),
                KeyValues("Fitted_tracks_proton", {"Fitted_tracks"})
            
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

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(m_genfitField);

        m_geoMaterial=GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);      
        genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
        genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
        genfit::MaterialEffects::getInstance()->setMscModel("Highland");

        // genfit::MaterialEffects::getInstance()->setDebugLvl(2); // Set the debug level for material effects


        // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
        // These object are necessary to extract the drift chamber hits information, such as positions of the wire extremities
        surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();

        // If the detector doesn't have a drift chamber, this part will be skipped
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}

        // Z-component of the magnetic field at the center of the detector
        // This component is used to propagate the track to the calorimeter surface, but it is also necessary to
        // compute the pt component given omega
        dd4hep::Position center(0, 0, 0);
        dd4hep::Direction bfield = m_field.magneticField(center);
        m_Bz = bfield.z() / dd4hep::tesla;

        // Retrive calorimeter information
        // These parameters are necessary to propagate the track to the calorimeter surface
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

    
    std::tuple< extension::TrackCollection,
                extension::TrackCollection,
                extension::TrackCollection,
                extension::TrackCollection,
                extension::TrackCollection> operator()( const extension::TrackCollection& tracks_input) const override                                                 
    {
        
        // These 5 collections store the output of the fit for 5 particle hypotesis: e, mu, pi, K, p
        extension::TrackCollection FittedTracks_electron;
        extension::TrackCollection FittedTracks_muon;
        extension::TrackCollection FittedTracks_pion;
        extension::TrackCollection FittedTracks_kaon;
        extension::TrackCollection FittedTracks_proton;

        info() << "Event number: " << index_counter++ << endmsg;
            
        // Loop over the tracks created by the pattern recognition step
        for (const auto& track : tracks_input)
        {
            

            if (track.getTrackerHits().size() == 0) continue; // skip empty tracks
            num_tracks  +=1;
            
            // Loop over the particle hypotesis
            for (int pdgCode : m_particleHypotesis)
            {

                // Create trackInterface, initialize genfit trakc and fit it
                GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, dch_info, dc_decoder,pdgCode,TVector3(m_init_pos_x, m_init_pos_y, m_init_pos_z),TVector3(m_init_mom_x, m_init_mom_y, m_init_mom_z));
                track_interface.createGenFitTrack(m_debug_lvl);
                bool isFit = track_interface.fit(m_Beta_init,m_Beta_final,m_Beta_steps, m_Bz); 
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

                    // If the fitter suppresses the hits (ndf or chi2 <= 0) the fit is considered as failed
                    if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0) {
                        
                        number_failures++;
                        if (pdgCode == 11)
                        {
                            auto failedTrack = FittedTracks_electron.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == -11)
                        {
                            auto failedTrack = FittedTracks_electron.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);

                        }
                        else if (pdgCode == 13)
                        {
                            auto failedTrack = FittedTracks_muon.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);;
                        }
                        else if (pdgCode == -13)
                        {
                            auto failedTrack = FittedTracks_muon.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == 211)
                        {
                            auto failedTrack = FittedTracks_pion.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == -211)
                        {
                            auto failedTrack = FittedTracks_pion.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == 321)
                        {
                            auto failedTrack = FittedTracks_kaon.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == -321)
                        {
                            auto failedTrack = FittedTracks_kaon.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == 2212)
                        {
                            auto failedTrack = FittedTracks_proton.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }
                        else if (pdgCode == -2212)
                        {
                            auto failedTrack = FittedTracks_proton.create();
                            failedTrack.setChi2(-1);
                            failedTrack.setNdf(-1);
                        }

                        debug() << "Number of hits in the track: " << track.getTrackerHits().size() << std::endl;
                        debug() << "GENFIT PDG: " << pdgCode << endmsg;
                        debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;
                        continue; 
                    }


                    double genfit_chi2_ndf_val = genfit_chi2_val / genfit_ndf_val;
                    auto genfit_hits_in_track = edm4hep_track.getTrackerHits();

                    debug() << "Number of hits in the track: " << track.getTrackerHits().size() << std::endl;
                    debug() << "GENFIT PDG: " << pdgCode << endmsg;
                    debug() << "GENFIT Omega: " << genfit_omega << endmsg;
                    debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
                    debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;

                    debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
                    debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
                    debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

                    debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
                    debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
                    debug() << "" << endmsg;

                    // Propagation of the track to the calorimeter surface: add atCalorimeter trackState to the track
		            int charge = getHypotesisCharge(pdgCode);                   
                    FillTrackWithCalorimeterExtrapolation(
                        edm4hep_track,
                        m_Bz,
                        charge,
                        a,
                        m_eCalBarrelInnerR,
                        m_eCalBarrelMaxZ,
                        m_eCalEndCapInnerR,
                        m_eCalEndCapOuterR,
                        m_eCalEndCapInnerZ
                    );
                    
                    // Fill track collections
                    if (pdgCode == 11)
                    {
                        FittedTracks_electron.push_back(edm4hep_track);
                    }
                    else if (pdgCode == -11)
                    {
                        FittedTracks_electron.push_back(edm4hep_track);
                    }
                    else if (pdgCode == 13)
                    {
                        FittedTracks_muon.push_back(edm4hep_track);
                    }
                    else if (pdgCode == -13)
                    {
                        FittedTracks_muon.push_back(edm4hep_track);
                    }
                    else if (pdgCode == 211)
                    {
                        FittedTracks_pion.push_back(edm4hep_track);
                    }
                    else if (pdgCode == -211)
                    {
                        FittedTracks_pion.push_back(edm4hep_track);
                    }
                    else if (pdgCode == 321)
                    {
                        FittedTracks_kaon.push_back(edm4hep_track);
                    }
                    else if (pdgCode == -321)
                    {
                        FittedTracks_kaon.push_back(edm4hep_track);
                    }
                    else if (pdgCode == 2212)
                    {
                        FittedTracks_proton.push_back(edm4hep_track);
                    }
                    else if (pdgCode == -2212)
                    {
                        FittedTracks_proton.push_back(edm4hep_track);
                    }
                       
                }
                else
                {   
                    // If the fit fails, it returns a track with chi2=ndf=-1
                    if (pdgCode == 11)
                    {
                        auto failedTrack = FittedTracks_electron.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == -11)
                    {
                        auto failedTrack = FittedTracks_electron.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);

                    }
                    else if (pdgCode == 13)
                    {
                        auto failedTrack = FittedTracks_muon.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);;
                    }
                    else if (pdgCode == -13)
                    {
                        auto failedTrack = FittedTracks_muon.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == 211)
                    {
                        auto failedTrack = FittedTracks_pion.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == -211)
                    {
                        auto failedTrack = FittedTracks_pion.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == 321)
                    {
                        auto failedTrack = FittedTracks_kaon.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == -321)
                    {
                         auto failedTrack = FittedTracks_kaon.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == 2212)
                    {
                        auto failedTrack = FittedTracks_proton.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }
                    else if (pdgCode == -2212)
                    {
                        auto failedTrack = FittedTracks_proton.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                    }

                    
                    debug() << "Number of hits in the track: " << track.getTrackerHits().size() << std::endl;
                    debug() << "GENFIT PDG: " << pdgCode << endmsg;
                    debug() << "GENFIT Chi2: " << -1 << ", GENFIT NDF: " << -1 << ", GENFIT Chi2/NDF: " << -1 << endmsg;


                    number_failures++;
                }
                
            }
            debug() << "----------------\n" << endmsg;
        }

        
        return std::make_tuple( std::move(FittedTracks_electron),
                                std::move(FittedTracks_muon),
                                std::move(FittedTracks_pion),
                                std::move(FittedTracks_kaon),
                                std::move(FittedTracks_proton));
        
    } 
    
    StatusCode finalize() {     
        
        info() << "Run report:" << endmsg;
        info() << "Number of failed tracks: " << number_failures <<  "/" << num_tracks << endmsg;
        info() << "----------------\n" << endmsg;



        return StatusCode::SUCCESS;

    }

    public:
        mutable int index_counter = 0;
        mutable int number_failures = 0;
        mutable int num_tracks = 0; // Total number of tracks processed

    private:

        ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
        genfit::MaterialEffects* materialEffects;
        genfit::FieldManager* fieldManager;

       
        dd4hep::Detector* m_detector{nullptr};  // Detector instance
        dd4hep::OverlayedField m_field;         // Magnetic field
        GenfitInterface::GenfitField* m_genfitField;
        GenfitInterface::GenfitMaterialInterface* m_geoMaterial;

        dd4hep::rec::SurfaceManager* surfMan;
        dd4hep::rec::DCH_info* dch_info;

        dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;

        Gaudi::Property<double> m_Beta_init{this, "Beta_init", 100, "Beta Initial value"};
        Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.1, "Beta Final value"};
        Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 10, "Beta number of Steps"};

        Gaudi::Property<int> m_debug_lvl{this, "debug_lvl", 0, "Debug level for GenfitTrack"};

        Gaudi::Property<int> m_init_pos_x{this, "init_pos_x", -1, "m_init_pos_x"};
        Gaudi::Property<int> m_init_pos_y{this, "init_pos_y", -1, "m_init_pos_y"};
        Gaudi::Property<int> m_init_pos_z{this, "init_pos_z", -1, "m_init_pos_z"};
        Gaudi::Property<int> m_init_mom_x{this, "init_mom_x", 0, "init_mom_x"};
        Gaudi::Property<int> m_init_mom_y{this, "init_mom_y", 0, "init_mom_y"};
        Gaudi::Property<int> m_init_mom_z{this, "init_mom_z", 0, "init_mom_z"};

        double c_light = 2.99792458e8;
        double a = c_light * 1e3 * 1e-15;
        double m_Bz = 2.0; // T

        double m_eCalBarrelInnerR = 0;
        double m_eCalBarrelMaxZ = 0;

        double m_eCalEndCapInnerR = 0;
        double m_eCalEndCapOuterR = 0;
        double m_eCalEndCapInnerZ = 0;
        double m_eCalEndCapOuterZ = 0;

        // std::vector<int> m_particleHypotesis = {11,13,211,321,2212}; // e, mu, pi, K, p
        std::vector<int> m_particleHypotesis = {211}; //pi


};

DECLARE_COMPONENT(GenfitTrackFitter)