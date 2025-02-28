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
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <thread> 
#include <chrono> 
#include <typeinfo>

#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "utils.hpp"

// Define collection types
#include "podio/UserDataCollection.h"
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;

#include "edm4hep/MCParticleCollection.h"
using ParticleColl = edm4hep::MCParticleCollection;

#include "edm4hep/SimTrackerHitCollection.h"
using SimHits = edm4hep::SimTrackerHitCollection;

#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHit3D.h"
#include "edm4hep/MutableTrackerHit3D.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
using VertexHitsColl = edm4hep::TrackerHit3DCollection;
using VTX_links = edm4hep::TrackerHitSimTrackerHitLinkCollection;

#include "extension/DriftChamberDigiCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociation.h"
using DCHitsColl = extension::DriftChamberDigiCollection;
using DC_associations = extension::MCRecoDriftChamberDigiAssociationCollection;

#include "edm4hep/TrackCollection.h"
#include "extension/TrackCollection.h"
#include "extension/TrackerHit.h"
using TrackColl = extension::TrackCollection;
using TrackHit = extension::TrackerHit;


/** @struct GGTF_performance_IDEAv2
 *
 *  Gaudi MultiTransformer...
 *
 *  input: 
*
 *
 *  output:
 *
 *
 *
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2024-08
 *
 */

struct GGTF_performance_IDEAv2 final : 
        k4FWCore::MultiTransformer<std::tuple<DoubleColl,DoubleColl,DoubleColl,DoubleColl,FloatColl,FloatColl,IntColl,IntColl,IntColl,IntColl,IntColl,IntColl>(          
                                                                                                                                                            const TrackColl&, 
                                                                                                                                                            const ParticleColl&,
                                                                                                                                                            const DC_associations&,
                                                                                                                                                            const VTX_links&,
                                                                                                                                                            const VTX_links&,
                                                                                                                                                            const VTX_links&)> 
{
    GGTF_performance_IDEAv2(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("InputCollectionTracks", {"inputTracks"}),
                KeyValues("InputCollectionParticles", {"inputMCparticles"}),
                KeyValues("DC_associations", {"DC_associations"}),
                KeyValues("VTXD_links", {"VTXD_links"}),
                KeyValues("VTXIB_links", {"VTXIB_links"}),
                KeyValues("VTXOB_links", {"VTXOB_links"})
            
            },
            {   
                KeyValues("costheta_mc", {"costheta_mc"}),
                KeyValues("pt_mc", {"pt_mc"}),
                KeyValues("vertex_mc", {"vertex_mc"}),
                KeyValues("deltaMC_mc", {"deltaMC_mc"}),

                KeyValues("purity_mc", {"purity_mc"}),
                KeyValues("efficiency_mc", {"efficiency_mc"}),
                KeyValues("isReco_mc", {"isReco_mc"}),
                KeyValues("isMatched_mc", {"isMatched_mc"}),
                KeyValues("matched_track_index", {"matched_track_index"}),

                KeyValues("matched_butNoReco_particles_index", {"matched_butNoReco_particles_index"}),
                KeyValues("fakeTracks_index", {"fakeTracks_index"}),
                KeyValues("tracking_eff", {"tracking_eff"})

            }) {}
            

    std::tuple<DoubleColl,DoubleColl,DoubleColl,DoubleColl,FloatColl,FloatColl,IntColl,IntColl,IntColl,IntColl,IntColl,IntColl> operator()(    const TrackColl& inputTracks, 
                                                                                                                                          const ParticleColl& inputMCparticles,
                                                                                                                                          const DC_associations& DC_asso,
                                                                                                                                          const VTX_links& VTXD_links,
                                                                                                                                          const VTX_links& VTXIB_links,
                                                                                                                                          const VTX_links& VTXOB_links) const override 
    {
        IntColl tracking_eff;

        DoubleColl costheta_mc = {};
        DoubleColl pt_mc = {};
        DoubleColl phi_mc = {};
        IntColl  pdg_mc = {};
        DoubleColl vertex_mc = {};
        DoubleColl deltaMC_mc = {};
        
        // MC particle physics parameters
        for (const auto MC_par : inputMCparticles) {

            double         px     = MC_par.getMomentum().x;
            double         py     = MC_par.getMomentum().y;
            double         pz     = MC_par.getMomentum().z;
            std::vector<double> info = computeCosThetaPtAndPhi(px, py, pz);

            auto vertex_vec = MC_par.getVertex();
            double R = sqrt(pow(vertex_vec[0],2)+pow(vertex_vec[1],2)+pow(vertex_vec[2],2));
            
            double costheta = info[0];
            double pt       = info[1];
            double phi = info[2];
            int pdg = MC_par.getPDG();

            pt_mc.push_back(pt);
            costheta_mc.push_back(costheta);
            pdg_mc.push_back(pdg);
            phi_mc.push_back(phi);
            vertex_mc.push_back(R);

            //TO-DO
            deltaMC_mc.push_back(-1.);

        }
        
        // Number of Hits per particle in each track
        // MC particle total number of Hits
        std::vector<std::map<int, int>> particles_in_tracks;
        std::vector<int> particles_numHits(inputMCparticles.size(),0);
        for (size_t i = 0; i< inputTracks.size(); i++)
        {
            
            auto hits_in_track = inputTracks[i].getTrackerHits();
            std::map<int, int> particles_track;
            for (const auto &hit : hits_in_track)
            {
                auto hit_id = hit.getObjectID();
                if (hit.isA<edm4hep::TrackerHit3D>()) {
                
                    for (const auto link : VTXD_links)
                    {
                        auto digi_hit_id = link.getFrom().getObjectID();

                        if (digi_hit_id == hit_id)
                        {
                            auto sim_hit = link.getTo();
                            auto particle_id = sim_hit.getParticle().getObjectID().index;
                            particles_numHits[particle_id] +=1;

                            // Increment the value in the map for the particle_id
                            if (particles_track.find(particle_id) == particles_track.end())
                            {
                                // If the particle_id is not in the map, initialize it with 1
                                particles_track[particle_id] = 1;
                            }
                            else
                            {
                                // If the particle_id is already in the map, increment its value
                                particles_track[particle_id]++;
                            }

                                continue;
                        }

                    }

                    for (const auto link : VTXIB_links)
                    {
                        auto digi_hit_id = link.getFrom().getObjectID();

                        if (digi_hit_id == hit_id)
                        {
                            auto sim_hit = link.getTo();
                            auto particle_id = sim_hit.getParticle().getObjectID().index;
                            particles_numHits[particle_id] +=1;

                            // Increment the value in the map for the particle_id
                            if (particles_track.find(particle_id) == particles_track.end())
                            {
                                // If the particle_id is not in the map, initialize it with 1
                                particles_track[particle_id] = 1;
                            }
                            else
                            {
                                // If the particle_id is already in the map, increment its value
                                particles_track[particle_id]++;
                            }

                                continue;
                        }
                    }

                    for (const auto link : VTXOB_links)
                    {
                        auto digi_hit_id = link.getFrom().getObjectID();

                        if (digi_hit_id == hit_id)
                        {
                            auto sim_hit = link.getTo();
                            auto particle_id = sim_hit.getParticle().getObjectID().index;
                            particles_numHits[particle_id] +=1;

                            // Increment the value in the map for the particle_id
                            if (particles_track.find(particle_id) == particles_track.end())
                            {
                                // If the particle_id is not in the map, initialize it with 1
                                particles_track[particle_id] = 1;
                            }
                            else
                            {
                                // If the particle_id is already in the map, increment its value
                                particles_track[particle_id]++;
                            }

                                continue;
                        }
                    }


                }
               
                if (hit.isA<extension::DriftChamberDigi>()) {
                    
                    for (const auto& association : DC_asso)
                    {
                        auto digi_hit_id = association.getDigi().getObjectID();

                        if (digi_hit_id == hit_id)
                        {
                            auto sim_hit = association.getSim();
                            auto particle_id = sim_hit.getParticle().getObjectID().index;
                            particles_numHits[particle_id] +=1;

                            // Increment the value in the map for the particle_id
                            if (particles_track.find(particle_id) == particles_track.end())
                            {
                                // If the particle_id is not in the map, initialize it with 1
                                particles_track[particle_id] = 1;
                            }
                            else
                            {
                                // If the particle_id is already in the map, increment its value
                                particles_track[particle_id]++;
                            }

                            continue;
                        }

                    }
                }

            }

            particles_in_tracks.push_back(particles_track);

        }

        IntColl isReco_mc = {};
        IntColl isMatched_mc;
        IntColl matched_butNoReco_mc = {};
        IntColl matched_track_mc = {};
        FloatColl purity_mc;
        FloatColl efficiency_mc;
       
        
        for (size_t i = 0; i < inputMCparticles.size(); i++)
        {
            int RECO = int((pt_mc[i]> pt_thr) && (costheta_mc[i] < cos_thr) && (particles_numHits[i] > hits_thr));
            isReco_mc.push_back(RECO);

            std::vector<int> isMATCHED;
            std::vector<float> track_purity;
            std::vector<float> track_efficiency;
            for (size_t j = 1; j < inputTracks.size(); j++)
            {
                    
                std::map<int, int> particles_track = particles_in_tracks[j];
                if (particles_track.find(i) != particles_track.end()) {
                
                    float purity = (particles_track[i]*1.0) / (float(inputTracks[j].getTrackerHits().size())*1.0);
                    float efficiency = (particles_track[i]*1.0) / (particles_numHits[i]*1.0);
                    int MATCHED = int((purity > pur_thr) && (efficiency > eff_thr));

                    isMATCHED.push_back(MATCHED);
                    track_purity.push_back(purity);
                    track_efficiency.push_back(efficiency);

                }
                else
                {
                    isMATCHED.push_back(0);
                    track_purity.push_back(-1.);
                    track_efficiency.push_back(-1.);

               }

            }

            int best_index = -1;
            float max_purity = std::numeric_limits<float>::lowest(); 
            float max_efficiency = std::numeric_limits<float>::lowest();

            for (size_t k = 0; k < isMATCHED.size(); ++k) {
                if (isMATCHED[k] > 0) {
                    if (track_purity[k] > max_purity || (track_purity[k] == max_purity && track_efficiency[k] > max_efficiency)) {
                        max_purity = track_purity[k];
                        max_efficiency = track_efficiency[k];
                        best_index = static_cast<int>(k+1);
                    }   
                }
            }

            if (best_index > 0)
            {
                isMatched_mc.push_back(1);
                purity_mc.push_back(max_purity);
                efficiency_mc.push_back(max_efficiency);

                if (!RECO)
                {
                    matched_butNoReco_mc.push_back(i);
                }
            }
            else
            {
                isMatched_mc.push_back(0);
                purity_mc.push_back(-1.);
                efficiency_mc.push_back(-1.);
            }
            
            matched_track_mc.push_back(best_index);

        }


        int number_of_matched_particles = 0;
        int number_of_reconstructable_particles = 0;
        for (size_t i = 0; i < inputMCparticles.size(); i++)
        {
            int matched_track_index = matched_track_mc[i];

            if (isReco_mc[i] > 0)
            {
                number_of_reconstructable_particles += 1;
                if (matched_track_index > -1)
                {
                    number_of_matched_particles += 1;
                }
            }
        }
        tracking_eff.push_back(number_of_matched_particles);
        tracking_eff.push_back(number_of_reconstructable_particles);

        std::map<int, int> countMap;
        for (int value : matched_track_mc) {
            countMap[value]++;
        }

        IntColl fakeTracks_index;
        for (size_t i = 1; i < inputTracks.size(); ++i) {
            if (!countMap.contains(i)) {
                fakeTracks_index.push_back(i);
            }
        }

        float denom = tracking_eff[1]*1;
        if (denom < 1.)
        {
            denom = 1.;
        }

        int num_track = inputTracks.size();
        if (num_track > 0)
        {
            num_track -= 1;
        }



        std::cout << "Number of Reconstructable MC particles: " << tracking_eff[1] << std::endl;
        std::cout << "Number of Tracks: " << num_track << std::endl;
        std::cout << "Tracking efficiency (overall): " << (tracking_eff[0]*1.)/(denom) << std::endl;
        std::cout << "Number of fake tracks (overall): " << fakeTracks_index.size() << std::endl;
        std::cout << "-----------------------" <<  std::endl;
        

        return std::make_tuple( std::move(costheta_mc),
                                std::move(pt_mc),
                                std::move(vertex_mc),
                                std::move(deltaMC_mc),

                                std::move(purity_mc),
                                std::move(efficiency_mc),
                                std::move(isReco_mc),
                                std::move(isMatched_mc),
                                std::move(matched_track_mc),
                                
                                std::move(matched_butNoReco_mc),
                                std::move(fakeTracks_index),
                                std::move(tracking_eff));
    }
    
    private:

        Gaudi::Property< double > pt_thr{this, "pt_thr", 0.1, "pt_thr"};
        Gaudi::Property< double > cos_thr{this, "cos_thr", 0.99, "cos_thr"};
        Gaudi::Property< int > hits_thr{this, "hits_thr", 4 , "hits_thr"};

        Gaudi::Property< double > pur_thr{this, "pur_thr",  0.5, "pur_thr"};
        Gaudi::Property< double > eff_thr{this, "eff_thr", 0.5, "eff_thr"};

};

DECLARE_COMPONENT(GGTF_performance_IDEAv2)
