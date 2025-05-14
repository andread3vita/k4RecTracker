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
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "DC_measurement.hpp"
#include "IDEAtrack.hpp"
#include "SI_measurement.hpp"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>
#include "utils.hpp"

// Define collection types
#include "podio/UserDataCollection.h"
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;


/** @struct IDEA_fitter_validation

*
*
*/

struct IDEA_fitter_validation final : 
        k4FWCore::MultiTransformer< std::tuple< extension::TrackCollection>(const edm4hep::MCParticleCollection&,
                                                                            const extension::SenseWireHitSimTrackerHitLinkCollection&,
                                                                            const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                            const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                            const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                            const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> 
            
                                                                                            
{
    IDEA_fitter_validation(const std::string& name, ISvcLocator* svcLoc) : 
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
                KeyValues("fittedTracks", {"fittedTracks"})
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
         
    StatusCode initialize() {

        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(new genfit::ConstField(0., 0., m_Bz.value())); // kGauss
        // fieldManager->init(m_geoSvc->lcdd()->field()); // kGauss
        

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

        // put an if statement to choose the detector
        dd4hep::Detector& mainDetector = *detector;
        const dd4hep::DetElement& detDRC = mainDetector.detector("DRcalo"); //put a gaudi property to choose the subetector
        const dd4hep::rec::LayeredCalorimeterData* DRCextension = detDRC.extension<dd4hep::rec::LayeredCalorimeterData>();

        m_eCalBarrelInnerR = DRCextension->extent[0] / dd4hep::mm;     // barrel rmin
        m_eCalBarrelMaxZ = DRCextension->extent[2] / dd4hep::mm;       // barrel zmax == endcap zmin

        m_eCalEndCapInnerR = DRCextension->extent[4] / dd4hep::mm;     // endcap rmin
        m_eCalEndCapOuterR = DRCextension->extent[5] / dd4hep::mm;     // endcap rmax
        m_eCalEndCapInnerZ = DRCextension->extent[2] / dd4hep::mm;     // endcap zmin
        m_eCalEndCapOuterZ = DRCextension->extent[3] / dd4hep::mm;     // endcap zmax

        return StatusCode::SUCCESS;

    }

    
    std::tuple<extension::TrackCollection> operator()(  const edm4hep::MCParticleCollection& mcParticles,
                                                        const extension::SenseWireHitSimTrackerHitLinkCollection& DC_asso,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXD_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperD_links) const override 
    {
        
        index_counter++;
        // if (index_counter != 11)
        // {
            
        //     return std::make_tuple(extension::TrackCollection());
        // }

        std::cout << "Event number: " << index_counter << std::endl;

        new TGeoManager("Geometry", "IDEA geometry");
        TGeoManager::Import("/eos/user/a/adevita/saveSpace/k4RecTracker/Tracking/TGeo_IDEA.root");
       
        
        extension::TrackCollection trackColl;
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
                
                debug() << "True pdg: " << particle.getPDG() << endmsg;
                debug() << "True Momentum: [" << mom.X() << " " << mom.Y() << " " << mom.Z() << " ]-> " << mom.Mag() << endmsg;
                debug() << "True pt: " << pt << endmsg;
                debug() << "True cos(theta): " << costheta << endmsg;
                debug() << "True phi: " << phi << endmsg;
                debug() << "" << endmsg;

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

                GENFIT::IDEAtrack track_interface = GENFIT::IDEAtrack(track, surfMan, dch_info, dc_decoder);
                track_interface.createGenFitTrack();

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
                    debug() << "GENFIT Omega: " << genfit_omega << endmsg;
                    debug() << "GENFIT Phi: " << genfit_phi_val << endmsg;
                    debug() << "GENFIT tanLambda: " << genfit_trackstate.tanLambda << endmsg;
                        
                    debug() << "GENFIT Momentum: " << "[ " << genfit_mom.X() << " " << genfit_mom.Y() << " " << genfit_mom.Z() << " ] -> " << genfit_p << endmsg;
                    debug() << "GENFIT pt: " << genfit_mom.Pt() << endmsg;
                    debug() << "GENFIT cosTheta: " << genfit_mom.CosTheta() << endmsg;

                    debug() << "GENFIT Chi2: " << genfit_chi2_val << ", GENFIT NDF: " << genfit_ndf_val << ", GENFIT Chi2/NDF: " << genfit_chi2_ndf_val << endmsg;
                    debug() << "Number of hits in GENFIT track: " << genfit_hits_in_track.size() << endmsg;
                    debug() << "" << endmsg;
                        
                    // // trackState at Calo
                    // auto trackStateLastHit = edm4hep_track.getTrackStates()[2];

                    // double omega_lastHit = trackStateLastHit.omega;
                    // double pt_lasthit = a * Bz / abs(omega_lastHit);
                    // double phi_lasthit = trackStateLastHit.phi;
                    // double pz_lasthit = trackStateLastHit.tanLambda * pt_lasthit;
                    // double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
                    // double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
                    // auto ref_lastHit = trackStateLastHit.referencePoint;

                    // // produce new helix at last hit position
                    // double posAtLastHit[] = {ref_lastHit[0], ref_lastHit[1], ref_lastHit[2]};
                    // double momAtLastHit[] = {px_lasthit, py_lasthit, pz_lasthit};
                    // auto helixAtLastHit = HelixClass_double();
                    // helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, 1, 2.);

                    // // TrackState at Calorimeter
                    // if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {

                    //     pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
                    //     pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
                    //     float minGenericTime(std::numeric_limits<float>::max());
                
                    //     // create helix to project
                    //     // rather than using parameters at production, better to use those from
                    //     // last hit
                    //     pandora::CartesianVector pos_lasthit(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
                    //     pandora::CartesianVector mom_lasthit(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);
                    //     const pandora::Helix helix(pos_lasthit, mom_lasthit, 1,2.);
                    //     const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
                    //     const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
                
                    //     // First project to endcap
                    //     pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
                    //     bool hasEndCapProjection(false);
                    //     if (m_eCalEndCapInnerR>0) {
                    //         float genericTime(std::numeric_limits<float>::max());
                    //         const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint,
                    //                                             endCapProjection, genericTime));
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
                    //                                             barrelProjection, genericTime));
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
                    //     edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit);
                    //     omega_lastHit = trackState_AtCalorimeter.omega;
                    //     pt_lasthit = a * Bz / abs(omega_lastHit);
                    //     phi_lasthit = trackState_AtCalorimeter.phi;
                    //     pz_lasthit = trackState_AtCalorimeter.tanLambda * pt_lasthit;
                    //     px_lasthit = pt_lasthit * std::cos(phi_lasthit);
                    //     py_lasthit = pt_lasthit * std::sin(phi_lasthit);
                    //     ref_lastHit = trackState_AtCalorimeter.referencePoint;
                    //     // attach the TrackState to the track
                    //     edm4hep_track.addToTrackStates(trackState_AtCalorimeter);

                    // }
                    
                    trackColl.push_back(edm4hep_track);
                }
                else
                {
                    auto failedTrack = trackColl.create();
                    failedTrack.setChi2(-1);
                    failedTrack.setNdf(-1);

                    debug() << "Chi2: " << failedTrack.getChi2() << " NDF: " << failedTrack.getNdf() << endmsg;
                }
                debug() << "----------------\n" << endmsg;
            }
        }
        
        return std::make_tuple(std::move(trackColl));
        
    } 

    public:
        mutable int index_counter = 0;

    private:

       
        // TGeoManager* geoManager;
        genfit::MaterialEffects* materialEffects;
        genfit::FieldManager* fieldManager;

        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        Gaudi::Property<double> m_Bz{this, "Bz", 20., "Bz value (kilogauss)"};

        Gaudi::Property<double> m_Beta_init{this, "Beta_init", 10., "Beta Initial value"};
        Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.1, "Beta Final value"};
        Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 10, "Beta number of Steps"};

        SmartIF<IGeoSvc> m_geoSvc;
        const dd4hep::rec::SurfaceManager* surfMan;
        const dd4hep::rec::DCH_info* dch_info;

        dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;

        double m_eCalBarrelInnerR = 0;
        double m_eCalBarrelMaxZ = 0;

        double m_eCalEndCapInnerR = 0;
        double m_eCalEndCapOuterR = 0;
        double m_eCalEndCapInnerZ = 0;
        double m_eCalEndCapOuterZ = 0;


};

DECLARE_COMPONENT(IDEA_fitter_validation)
