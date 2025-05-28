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
#include <iostream>
#include <vector>
#include <tuple>

//=== Gaudi / k4FWCore ===
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

//=== podio / edm4hep ===
#include "podio/UserDataCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "extension/TrackCollection.h"

//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"

/** @struct IDEATracksFromGenParticles
 *
 *  Gaudi MultiTransformer that generates a collection of MC truth tracks by associating tracker and drift chamber hits 
 *  to MCParticles at the interaction point (IP). The algorithm loops over the MCParticle collection and builds track-like 
 *  objects by collecting all hits linked to primary particles through simulation truth associations.
 *
 *  Only particles with generator status equal to 1 and originating at the origin (0,0,0) are selected for track building.
 *  For each selected MCParticle, the corresponding hits are gathered from different subdetectors:
 *    - Vertex detector (VTXD and VTXB)
 *    - Barrel and disk wrappers (wrapperB, wrapperD)
 *    - Drift Chamber (DC)
 *
 *  Each hit is associated via a SimTrackerHitLink, and only hits produced by the primary particle (not secondaries) are kept.
 *
 *  input:
 *    - MCParticles collection : mcParticles
 *    - DC hit-to-MC links : DC_links
 *    - VTXD hit-to-MC links : VTXD_links
 *    - VTXB hit-to-MC links : VTXB_links
 *    - Wrapper barrel hit-to-MC links : wrapperB_links
 *    - Wrapper disk hit-to-MC links : wrapperD_links
 *
 *  output:
 *    - Track collection built from truth information : MCTracks
 *    - Index mapping from tracks to MCParticles : MCParticleIndex
 *
 *
 *  @author Andrea De Vita
 *  @date   2025-05
 *
 */


struct IDEATracksFromGenParticles final : 
        k4FWCore::MultiTransformer< std::tuple<extension::TrackCollection,podio::UserDataCollection<int>>(  const edm4hep::MCParticleCollection&,
                                                                                                                                        const extension::SenseWireHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> 
            
                                                                                            
{
    IDEATracksFromGenParticles(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("mcParticles", {"mcParticles"}),
                KeyValues("DC_links", {"DC_links"}),
                KeyValues("VTXD_links", {"VTXD_links"}),
                KeyValues("VTXB_links", {"VTXIB_links"}),
                KeyValues("wrapperB_links", {"wrapperB_links"}),
                KeyValues("wrapperD_links", {"wrapperD_links"})
            },
            {   
                KeyValues("MCTracks", {"MCTracks"}),
                KeyValues("MCParticleIndex", {"MCParticleIndex"})
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
         
    StatusCode initialize() {

        return StatusCode::SUCCESS;

    }

    
    std::tuple<extension::TrackCollection,podio::UserDataCollection<int>> operator()(   const edm4hep::MCParticleCollection& mcParticles,
                                                                                        const extension::SenseWireHitSimTrackerHitLinkCollection& DC_links,
                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXD_links,
                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXB_links,
                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperB_links,
                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperD_links) const override 
    {

        extension::TrackCollection trackColl;
        std::vector<int> MCParticleIndex;
        for (auto particle : mcParticles)
        {
            auto genStatus = particle.getGeneratorStatus();
            auto vertex = particle.getVertex();


            if (genStatus == 1 && vertex.x == 0. && vertex.y == 0. && vertex.z == 0.) // make sure that you are taking only particles with vertex at the IP
            {

                auto particle_id = particle.getObjectID().index;
                MCParticleIndex.push_back(particle_id);

                auto track = trackColl.create();
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

                for (const auto& link : DC_links)
                {

                    auto part_id = link.getTo().getParticle().getObjectID().index;
                    auto flag = link.getTo().isProducedBySecondary();

                    if (part_id == particle_id && !flag)
                    {
                        auto digi_hit = link.getFrom();
                        track.addToTrackerHits(digi_hit);

                    }
                        

                }
            
            }
            
            debug() << "OK" << endmsg;
        }
        
        return std::make_tuple(std::move(trackColl),std::move(MCParticleIndex));
        
    }

    private:

        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        SmartIF<IGeoSvc> m_geoSvc;
        


};

DECLARE_COMPONENT(IDEATracksFromGenParticles)
