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

//=== podio / edm4hep ===
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackState.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "Wire_measurement.hpp"
#include "GenfitTrack.hpp"
#include "Planar_measurement.hpp"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>
#include "utils.hpp"

// Define collection types
#include "podio/UserDataCollection.h"
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;


/** @struct FitterWithXML

*
*
*/

struct FitterWithXML final : 
        k4FWCore::MultiTransformer< std::tuple< extension::TrackCollection,IntColl,IntColl >(  const edm4hep::MCParticleCollection&,
                                                                                                                                        const extension::SenseWireHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> 
            
                                                                                            
{
    FitterWithXML(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("mcParticles", {"mcParticles"}),
                KeyValues("DC_associations", {"DC_associations"}),
                KeyValues("VTXD_links", {"VTXD_links"}),
                KeyValues("VTXB_links", {"VTXIB_links"}),
                KeyValues("wrapperB_links", {"wrapperB_links"}),
                KeyValues("wrapperD_links", {"wrapperD_links"})
            },
            {   
                KeyValues("fittedTracks", {"fittedTracks"}),
                KeyValues("indexTrack", {"indexTrack"}),
                KeyValues("indexMCparticle", {"indexMCparticle"})
            
            }) {}
         
    StatusCode initialize() {

        

        return StatusCode::SUCCESS;

    }

    
    std::tuple<extension::TrackCollection,IntColl,IntColl> operator()(  const edm4hep::MCParticleCollection& mcParticles,
                                                        const extension::SenseWireHitSimTrackerHitLinkCollection& DC_asso,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXD_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperD_links) const override 
    {
        
        IntColl indexTrack;
        IntColl indexMCparticle;
        extension::TrackCollection trackColl;
        
        std::cout << "Event number: " << index_counter++ << std::endl;
        if (index_counter != 15)
        {
            
            return std::make_tuple(std::move(trackColl),std::move(indexTrack),std::move(indexMCparticle));
        }
       
        // Initialization
        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(new genfit::ConstField(0., 0., m_Bz.value())); // kGauss
        

        dd4hep::Detector* detector = (m_geoSvc->getDetector());
        surfMan = detector->extension<dd4hep::rec::SurfaceManager>();

        //-----------------
        // Retrieve the subdetector
        std::string DCH_name("DCH_v2");
        dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);

        dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();

        // Retrieve the readout associated with the detector element (subdetector)
        dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
        dd4hep::Readout dch_readout = dch_sd.readout();
        // set the cellID decoder
        dc_decoder = dch_readout.idSpec().decoder();

       
        for (auto particle : mcParticles)
        {
            auto genStatus = particle.getGeneratorStatus();
            auto vertex = particle.getVertex();

            if (genStatus == 1 && vertex.x == 0. && vertex.y == 0. && vertex.z == 0.) // make sure that you are taking only pions with vertex at the IP
            {

                auto particle_id = particle.getObjectID().index;

                auto p = particle.getMomentum();
                TVector3 mom(p.x, p.y, p.z);
                auto pt = mom.Pt();
                auto costheta = mom.CosTheta();
                auto phi = mom.Phi();
                
                
                
                indexMCparticle.push_back(particle.getObjectID().index);
                extension::MutableTrack track;
                // fill track with digi hits
                for (const auto link : VTXD_links)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();
                    if (part_id == particle_id && !flag)
                    {
                        auto digi_hit = link.getFrom().as<edm4hep::TrackerHitPlane>();
                        track.addToTrackerHits(digi_hit);
                    }
                }

                for (const auto link : VTXB_links)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();

                    if (part_id == particle_id && !flag)
                    {   
                        auto digi_hit = link.getFrom().as<edm4hep::TrackerHitPlane>();
                        track.addToTrackerHits(digi_hit);
                    }
                }

                for (const auto link : wrapperB_links)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();
                    
                    if (part_id == particle_id && !flag)
                    {
                        auto digi_hit = link.getFrom().as<edm4hep::TrackerHitPlane>();
                        track.addToTrackerHits(digi_hit);
                    }
                }

                for (const auto link : wrapperD_links)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();

                    if (part_id == particle_id && !flag)
                    {
                        auto digi_hit = link.getFrom().as<edm4hep::TrackerHitPlane>();
                        track.addToTrackerHits(digi_hit);
                    }
                }

                for (const auto& link : DC_asso)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();

                    if (part_id == particle_id && !flag)
                    {
                        auto digi_hit = link.getFrom();
                        track.addToTrackerHits(digi_hit);

                    }
                        

                }
            
                //////////////////////////
                //////// GENFIT //////////
                //////////////////////////

                GENFIT::GenfitTrack track_interface = GENFIT::GenfitTrack(track, dch_info, dc_decoder,211);
                track_interface.createGenFitTrack();

                auto test_track = track_interface.getTrack_edm4hep();
                std::cout << "Number of hits: " << test_track.getTrackerHits().size() << std::endl;

                bool isFit = track_interface.fit(m_Beta_init,m_Beta_final,m_Beta_steps);

                if (isFit)
                {

                    double c_light = 2.99792458e8;
                    double a = c_light * 1e3 * 1e-15;
                    double Bz = 2.0;

                    auto edm4hep_track = track_interface.getTrack_edm4hep();
                    auto genfit_trackstate = edm4hep_track.getTrackStates()[0];
            
                    double genfit_omega = genfit_trackstate.omega;
                    double genfit_pt_val = a * Bz / abs(genfit_omega);
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

                    debug() << "True pdg: " << particle.getPDG() << endmsg;
                    debug() << "True Momentum: [" << mom.X() << " " << mom.Y() << " " << mom.Z() << " ]-> " << mom.Mag() << endmsg;
                    debug() << "True pt: " << pt << endmsg;
                    debug() << "True cos(theta): " << costheta << endmsg;
                    debug() << "True phi: " << phi << endmsg;
                    debug() << "" << endmsg;

                    debug() << "GENFIT Omega: " << genfit_omega << endmsg;
                    debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
                    debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;
                        
                    debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
                    debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
                    debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

                    debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
                    debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
                    debug() << "" << endmsg;
                        
                    
                    indexTrack.push_back(0);
                    trackColl.push_back(edm4hep_track);

                }
                else
                {
                    
                    auto failedTrack = trackColl.create();
                    failedTrack.setChi2(-1);
                    failedTrack.setNdf(-1);
                    indexTrack.push_back(0);

                    debug() << "Chi2: " << failedTrack.getChi2() << " NDF: " << failedTrack.getNdf() << endmsg;

                    failures++;
                }
                debug() << "----------------\n" << endmsg;
            }
        }
        
        debug() << "Number of failures: " << failures << endmsg;
        return std::make_tuple(std::move(trackColl),std::move(indexTrack),std::move(indexMCparticle));
        
    } 

    public:
        mutable int index_counter = 0;
        mutable int failures = 0;
       
        // TGeoManager* geoManager;
        mutable genfit::MaterialEffects* materialEffects;
        mutable genfit::FieldManager* fieldManager;

        mutable Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        mutable Gaudi::Property<double> m_Bz{this, "Bz", 20., "Bz value (kilogauss)"};

        mutable Gaudi::Property<double> m_Beta_init{this, "Beta_init", 10., "Beta Initial value"};
        mutable Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.1, "Beta Final value"};
        mutable Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 10, "Beta number of Steps"};

        mutable ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
        mutable const dd4hep::rec::SurfaceManager* surfMan;
        mutable const dd4hep::rec::DCH_info* dch_info;

        mutable dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;


};

DECLARE_COMPONENT(FitterWithXML)
