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

// Define threshold constants
const double pt_thr_ = 0.100;
const double cos_thr_ = 0.99;
const double hits_thr_ = 4.0;

const double pur_thr_ = 0.66;
const double eff_thr_ = 0.05;

struct Efficiency_calc_gun final : 
        k4FWCore::MultiTransformer<std::tuple<IntColl, DoubleColl, IntColl, IntColl, IntColl, IntColl, IntColl>( const TrackColl&, 
                                                                 const ParticleColl&,
                                                                 const DCTrackerHitColl_sim&,
                                                                 const VertexColl_sim&,
                                                                 const VertexColl_sim&,
                                                                 const VertexColl_sim&)> 
{
    Efficiency_calc_gun(const std::string& name, ISvcLocator* svcLoc) : 
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
                KeyValues("pdg_MCParticles", {"pdg_MCParticles"}),
                KeyValues("pt_MCParticles", {"pt_MCParticles"}),
                KeyValues("index_MCParticles", {"index_MCParticles"}),
                KeyValues("isReco", {"isReco"}),
                KeyValues("isPrimary", {"isPrimary"}),
                KeyValues("outputTrackingEff_primary", {"outputTrackingEff_primary"}),
                KeyValues("outputTrackingEff_secondary", {"outputTrackingEff_secondary"})
            
            }) {}
            

    std::tuple<IntColl, DoubleColl, IntColl, IntColl, IntColl, IntColl, IntColl> operator()( const TrackColl& inputTracks, 
                                                                                             const ParticleColl& inputMCparticles,
                                                                                             const DCTrackerHitColl_sim& inputHits_CDC_sim,
                                                                                             const VertexColl_sim& inputHits_VTXIB_sim,
                                                                                             const VertexColl_sim& inputHits_VTXD_sim,
                                                                                             const VertexColl_sim& inputHits_VTXOB_sim) const override 
    {

        std::vector<double> costheta_mc = {};
        std::vector<double> pt_mc = {};
        std::vector<int>  pdg_mc = {};
        std::vector<int>  index_mc = {};
        std::vector<int>  isPrimary= {};

        for (const auto MC_par : inputMCparticles) {

            // define a lorentz vector and compute theta and pt
            double         px     = MC_par.getMomentum().x;
            double         py     = MC_par.getMomentum().y;
            double         pz     = MC_par.getMomentum().z;
            std::vector<double> info = computeCosThetaAndPt(px, py, pz);

            double pt       = info[1];
            double costheta = info[0];

            int pdg = MC_par.getPDG();
            auto index_MC = static_cast<int>(MC_par.getObjectID().index);

            pt_mc.push_back(pt);
            costheta_mc.push_back(costheta);
            pdg_mc.push_back(pdg);
            index_mc.push_back(index_MC);

            int isPrim = (MC_par.getGeneratorStatus() == 1);
            isPrimary.push_back(isPrim);

        }

        /// compute the number of hits for each particle
        std::vector<int> particle_hits(isPrimary.size(), 0);
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

        /// check which particle is reconstrutable
        std::vector<int> isReco;
        for (size_t i = 0; i < numPart; i++) {

            int isRec = (pt_mc[i] > pt_thr_) && (costheta_mc[i] < cos_thr_) && (particle_hits[i] > hits_thr_);
            isReco.push_back(isRec);
        
        }

        // std::cout << "MC particles properties:" << std::endl;
        // for (size_t i = 0; i < numPart; i++) {
        //             std::cout << "Particle:" << i << 
        //             "\tindex:" << index_mc[i] <<
        //             "\tpdg_id:" << pdg_mc[i] <<
        //             "\tIsPrim:" << isPrimary[i] <<
        //             "\tcostheta:" << costheta_mc[i] << 
        //             "\tpt:" << pt_mc[i] << 
        //             "\tnumHits:" << particle_hits[i] <<
        //             "\tIsReco:" << isReco[i] << std::endl;
        // }

        /// Efficiency
        int total_hits = 0;
        std::vector<std::map<int, double>> eff_mc;
        std::vector<std::map<int, int>>    hits_Tracks;

        std::vector<int> track_isCluster;
        for (const auto track : inputTracks) {

            auto hits_in_track = track.getTrackerHits();
            int  numHits       = hits_in_track.size() > 0 ? hits_in_track.size() : 1;
            
            int label = hits_in_track[0].getType();
            if (label >= 0)
            {
                track_isCluster.push_back(1);
            }
            else
            {
                track_isCluster.push_back(0);
            }

            std::vector<int> hit_idx;
            hit_idx.reserve(hits_in_track.size());
            for (int k = 0; k < numHits; ++k) {
            int idx = static_cast<int>(hits_in_track[k].getEDep());
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
                
                EFF[value]      = counter[value] / static_cast<double>(numHits);
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

    
        int numberOfMatchedPrimary = 0;
        int numberOfMatchedSecondary = 0;

        for (size_t i = 0; i < pur_tracks.size(); ++i) {
            for (size_t s = 0; s < eff_mc.size(); ++s) {
                double pur = 0.0;
                double eff = 0.0;

                // Find the purity
                auto pur_it = pur_tracks[i].find(s);
                if (pur_it != pur_tracks[i].end()) {
                    pur = pur_it->second;
                }

                // Find the efficiency
                auto eff_it = eff_mc[s].find(i);
                if (eff_it != eff_mc[s].end()) {
                    eff = eff_it->second;
                }

                // Check if purity and efficiency are above the thresholds
                if (pur > pur_thr_ && eff > eff_thr_ && isReco[i] == 1 && isPrimary[i] == 1 && track_isCluster[s] == 1) {
                    numberOfMatchedPrimary++;
                }
                else if (pur > pur_thr_ && eff > eff_thr_  && isReco[i] == 1 && isPrimary[i] == 0 && track_isCluster[s] == 1) {
                numberOfMatchedSecondary++;
                }
            }
        }
        

        int tracking_eff_denom_primary = 0;
        int tracking_eff_denom_secondary = 0;
        for(size_t res = 0; res < isPrimary.size(); res++)
        {
            if(isPrimary[res] == 1 && isReco[res] == 1)
            {
                tracking_eff_denom_primary += 1;
            }
            else if(isPrimary[res] == 0 && isReco[res] == 1)
            {
                tracking_eff_denom_secondary += 1;
            }

        }

        int tracking_eff_num_primary                       = numberOfMatchedPrimary;
        int tracking_eff_num_secondary                     = numberOfMatchedSecondary;

        IntColl primary_eff;
        IntColl secondary_eff;

        primary_eff.push_back(tracking_eff_num_primary);
        primary_eff.push_back(tracking_eff_denom_primary);

        secondary_eff.push_back(tracking_eff_num_secondary);
        secondary_eff.push_back(tracking_eff_denom_secondary);

        // std::cout << "\n\nResults:" << std::endl;
        // std::cout << "Efficiency for primary particles:" << tracking_eff_num_primary << "/" << tracking_eff_denom_primary << std::endl;
        // std::cout << "Efficiency for secondary particles:" << tracking_eff_num_secondary << "/" << tracking_eff_denom_secondary << std::endl;

        return std::make_tuple(std::move(pdg_mc), std::move(pt_mc), std::move(index_mc), std::move(isReco),std::move(isPrimary),std::move(primary_eff), std::move(secondary_eff));
    }
    

};

DECLARE_COMPONENT(Efficiency_calc_gun)