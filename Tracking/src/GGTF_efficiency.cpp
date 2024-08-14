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

#include "Gaudi/Property.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHit3D.h"
#include "podio/UserDataCollection.h"
#include "k4FWCore/Transformer.h"

#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/TrackCollection.h"
#include "extension/MutableTrackerHit3D.h"

#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep { using TrackerHit3DCollection = edm4hep::TrackerHitCollection;}  // namespace edm4hep
#endif

#include "utils.hpp"

// Define collection types
using DoubleColl = podio::UserDataCollection<double>;
using FloatColl = podio::UserDataCollection<float>;
using IntColl = podio::UserDataCollection<int>;
using ParticleColl = edm4hep::MCParticleCollection;
using DCTrackerHitColl = extension::DriftChamberDigiCollection;
using VertexColl = edm4hep::TrackerHit3DCollection;
using DCTrackerHitColl_sim = edm4hep::SimTrackerHitCollection;
using VertexColl_sim = edm4hep::SimTrackerHitCollection;
using TrackColl = extension::TrackCollection;

struct GGTF_efficiency final : 
        k4FWCore::MultiTransformer<std::tuple<DoubleColl, DoubleColl, DoubleColl, IntColl, IntColl,DoubleColl,DoubleColl,IntColl,IntColl,IntColl,IntColl>(          
                                                                                                                            const TrackColl&, 
                                                                                                                            const ParticleColl&,
                                                                                                                            const DCTrackerHitColl_sim&,
                                                                                                                            const VertexColl_sim&,
                                                                                                                            const VertexColl_sim&,
                                                                                                                            const VertexColl_sim&)> 
{
    GGTF_efficiency(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("InputCollectionTracks", {"inputTracks"}),
                KeyValues("InputCollectionParticles", {"inputMCparticles"}),
                KeyValues("inputHits_CDC_sim", {"inputHits_CDC_sim"}),
                KeyValues("inputHits_VTXIB_sim", {"inputHits_VTXIB_sim"}),
                KeyValues("inputHits_VTXD_sim", {"inputHits_VTXD_sim"}),
                KeyValues("inputHits_VTXOB_sim", {"inputHits_VTXOB_sim"})
            
            },
            {   
                KeyValues("out_costheta", {"out_costheta"}),
                KeyValues("out_pt", {"out_pt"}),
                KeyValues("out_phi", {"out_phi"}),
                KeyValues("out_pdg", {"out_pdg"}),
                KeyValues("out_num_hits", {"out_num_hits"}),
                KeyValues("out_pur", {"out_pur"}),
                KeyValues("out_eff", {"out_eff"}),
                KeyValues("assigned_track_mc", {"assigned_track_mc"}),
                KeyValues("numberFakes", {"numberFakes"}),
                KeyValues("isReco", {"isReco"}),
                KeyValues("isTrack", {"isTrack"}),

            }) {}
            

    std::tuple<DoubleColl, DoubleColl, DoubleColl, IntColl, IntColl,DoubleColl,DoubleColl,IntColl,IntColl,IntColl,IntColl> operator()( const TrackColl& inputTracks, 
                                                                                             const ParticleColl& inputMCparticles,
                                                                                             const DCTrackerHitColl_sim& inputHits_CDC_sim,
                                                                                             const VertexColl_sim& inputHits_VTXIB_sim,
                                                                                             const VertexColl_sim& inputHits_VTXD_sim,
                                                                                             const VertexColl_sim& inputHits_VTXOB_sim) const override 
    {

        DoubleColl costheta_mc = {};
        DoubleColl pt_mc = {};
        DoubleColl phi_mc = {};
        IntColl  pdg_mc = {};
        
        for (const auto MC_par : inputMCparticles) {

            // define a lorentz vector and compute theta and pt
            double         px     = MC_par.getMomentum().x;
            double         py     = MC_par.getMomentum().y;
            double         pz     = MC_par.getMomentum().z;
            std::vector<double> info = computeCosThetaPtAndPhi(px, py, pz);
            
            double costheta = info[0];
            double pt       = info[1];
            double phi = info[2];
            int pdg = MC_par.getPDG();

            pt_mc.push_back(pt);
            costheta_mc.push_back(costheta);
            pdg_mc.push_back(pdg);
            phi_mc.push_back(phi);

        }

        /// compute the number of hits for each particle
        std::vector<int> particle_hits(phi_mc.size(), 0);
        for (const auto hit : inputHits_VTXD_sim) {

            int part_index = hit.getParticle().getObjectID().index;
            particle_hits[part_index] +=1;
        }

        for (const auto hit : inputHits_VTXIB_sim) {
        int part_index = hit.getParticle().getObjectID().index;
        particle_hits[part_index] += 1;
        }

        for (const auto hit : inputHits_VTXOB_sim) {
        int part_index = hit.getParticle().getObjectID().index;
        particle_hits[part_index] += 1;
        }

        for (const auto hit : inputHits_CDC_sim) {
        int part_index = hit.getParticle().getObjectID().index;
        particle_hits[part_index] += 1;
        }
        size_t numPart = particle_hits.size();

        /// Efficiency
        int total_hits = 0;
        std::vector<std::map<int, double>> eff_mc;
        std::vector<std::map<int, int>>    hits_Tracks;

        std::vector<int> hits_num;
        std::vector<int> track_labels;
        for (const auto track : inputTracks) {
            
            int label = track.getDEdx();
            track_labels.push_back(label);

            auto hits_in_track = track.getTrackerHits();
            int numHits_track = hits_in_track.size();
            hits_num.push_back(numHits_track);

            std::vector<int> hit_idx;
            hit_idx.reserve(numHits_track);

            for (auto hit : hits_in_track) {

                int idx = hit.getEDep();
                hit_idx.push_back(idx);
                total_hits += 1;
            }

            std::map<int, int> counter;
            for (int idx : hit_idx) {
                counter[idx]++;
            }

            std::map<int, double> EFF;
            std::map<int, int>    hitsPart;
            for (size_t value = 0; value < numPart; ++value) {
                
                EFF[value]      = counter[value] / static_cast<double>(numHits_track);
                hitsPart[value] = counter[value];
            
            }

            eff_mc.push_back(EFF);
            hits_Tracks.push_back(hitsPart);

        }

        /// Purity
        std::vector<std::map<int, double>> pur_tracks;
        int temp_count = 0;
        for (const auto part : inputMCparticles) {
            int numHitsTracker = 0;

            for (const auto& hits_map : hits_Tracks) {
            auto it = hits_map.find(temp_count);
            if (it != hits_map.end()) {
                numHitsTracker += it->second;
            }
            }

            numHitsTracker = (numHitsTracker > 0) ? numHitsTracker : 1;

            std::map<int, double> PUR;
            for (size_t value = 0; value < hits_Tracks.size(); ++value) {
                
                auto it       = hits_Tracks[value].find(temp_count);
                int  hitCount = (it != hits_Tracks[value].end()) ? it->second : 0;
                PUR[value]    = static_cast<double>(hitCount) / numHitsTracker;
            }

            temp_count += 1;
            pur_tracks.push_back(PUR);
        }

        DoubleColl efficiency_mc;
        DoubleColl purity_mc;
        IntColl assigned_track_mc;
        IntColl isReco;
        IntColl isTrack;
        for (size_t i = 0; i < pur_tracks.size(); ++i) {

            double pur;
            double eff;
            double assigned_eff;
            double assigned_pur;
            int assigned_track;

            int RECO;
            int TRACK;
            for (size_t s = 0; s < eff_mc.size(); ++s) {
                
                auto pur_it = pur_tracks[i].find(s);
                if (pur_it != pur_tracks[i].end()) {
                    pur = pur_it->second;
                }

                auto eff_it = eff_mc[s].find(i);
                if (eff_it != eff_mc[s].end()) {
                    eff = eff_it->second;

                }
                
                if (pur > pur_thr__ && eff > eff_thr__)
                {   

                    assigned_track = track_labels[s];
                    assigned_eff = eff;
                    assigned_pur = pur;
                    TRACK = 1;

                    break;
                    
                }
                else
                {   
                    assigned_track = 0;
                    assigned_eff = -1.0;
                    assigned_pur = -1.0;
                    TRACK = 0;
                }
                
            
            }

            RECO = int((pt_mc[i]> pt_thr) && (costheta_mc[i] < cos_thr) && (particle_hits[i] > hits_thr));

            isReco.push_back(RECO);
            isTrack.push_back(TRACK && RECO);

            assigned_track_mc.push_back(assigned_track);
            purity_mc.push_back(assigned_pur);
            efficiency_mc.push_back(assigned_eff);

        }

        std::set<int> unique_tracks;    
        for (int num : assigned_track_mc) 
        {        
            if (num != 0) 
            {            
                unique_tracks.insert(num);        
            }   
        }

        int numberOfFakes = inputTracks.size()-unique_tracks.size();

        IntColl numberFakes;
        numberFakes.push_back(numberOfFakes);

        return std::make_tuple( std::move(costheta_mc), 
                                std::move(pt_mc), 
                                std::move(phi_mc), 
                                std::move(pdg_mc),
                                std::move(particle_hits),
                                std::move(purity_mc),
                                std::move(efficiency_mc),
                                std::move(assigned_track_mc),
                                std::move(numberFakes),
                                std::move(isReco),
                                std::move(isTrack));
    }
    
    private:

        Gaudi::Property< double > pt_thr{this, "pt_thr", 0.1, "pt_thr"};
        Gaudi::Property< double > cos_thr{this, "cos_thr", 0.99, "cos_thr"};
        Gaudi::Property< int > hits_thr{this, "hits_thr", 4 , "hits_thr"};

        Gaudi::Property< double > pur_thr__{this, "modelPath",  0.66, "pur_thr__"};
        Gaudi::Property< double > eff_thr__{this, "eff_thr__", 0.05, "eff_thr__"};

};

DECLARE_COMPONENT(GGTF_efficiency)