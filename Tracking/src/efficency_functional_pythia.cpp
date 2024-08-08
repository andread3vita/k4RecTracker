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
const double pt_thr__ = 0.100;
const double cos_thr__ = 0.99;
const double hits_thr__ = 4.0;

const double pur_thr__ = 0.66;
const double eff_thr__ = 0.05;

struct Efficiency_calc_pythia final : 
        k4FWCore::MultiTransformer<std::tuple<DoubleColl, 
                                              DoubleColl, 
                                              DoubleColl, 
                                              IntColl, 
                                              DoubleColl,
                                              DoubleColl,
                                              IntColl, 
                                              IntColl,
                                              IntColl,
                                              IntColl,
                                              DoubleColl>(          const TrackColl&, 
                                                                 const ParticleColl&,
                                                                 const DCTrackerHitColl_sim&,
                                                                 const VertexColl_sim&,
                                                                 const VertexColl_sim&,
                                                                 const VertexColl_sim&)> 
{
    Efficiency_calc_pythia(const std::string& name, ISvcLocator* svcLoc) : 
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
                KeyValues("out_phi", {"out_phi"}),
                KeyValues("out_pt", {"out_pt"}),
                KeyValues("out_pdg", {"out_pdg"}),
                KeyValues("out_purity", {"out_purity"}),
                KeyValues("out_eff", {"out_eff"}),
                KeyValues("out_index_part", {"out_index_part"}),
                KeyValues("out_nHits", {"out_nHits"}),
                KeyValues("out_isReco_tracks", {"out_isReco_tracks"}),
                KeyValues("out_isReco_particles", {"out_isReco_particles"}),
                KeyValues("out_pt_particles", {"out_pt_particles"}),
            
            }) {}
            

    std::tuple<DoubleColl, DoubleColl, DoubleColl, IntColl, DoubleColl,DoubleColl, IntColl, IntColl,IntColl,IntColl,DoubleColl> operator()( const TrackColl& inputTracks, 
                                                                                             const ParticleColl& inputMCparticles,
                                                                                             const DCTrackerHitColl_sim& inputHits_CDC_sim,
                                                                                             const VertexColl_sim& inputHits_VTXIB_sim,
                                                                                             const VertexColl_sim& inputHits_VTXD_sim,
                                                                                             const VertexColl_sim& inputHits_VTXOB_sim) const override 
    {

        std::vector<double> costheta_mc = {};
        std::vector<double> phi_mc = {};
        std::vector<int>  pdg_mc = {};
        std::vector<double> pt_mc = {};

        std::vector<int>  index_mc = {};
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
            auto index_MC = static_cast<int>(MC_par.getObjectID().index);

            pt_mc.push_back(pt);
            costheta_mc.push_back(costheta);
            pdg_mc.push_back(pdg);
            phi_mc.push_back(phi);
            index_mc.push_back(index_MC);

        }

        /// compute the number of hits for each particle
        std::vector<int> particle_hits(index_mc.size(), 0);
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

            int isRec = (pt_mc[i] > pt_thr__) && (costheta_mc[i] < cos_thr__) && (particle_hits[i] > hits_thr__);
            isReco.push_back(isRec);
        
        }

        /// Efficiency
        int total_hits = 0;
        std::vector<std::map<int, double>> eff_mc;
        std::vector<std::map<int, int>>    hits_Tracks;

        std::vector<int> track_isCluster;

        std::vector<int> hits_num;
        for (const auto track : inputTracks) {

            auto hits_in_track = track.getTrackerHits();
            int  numHits       = hits_in_track.size() > 0 ? hits_in_track.size() : 1;

            hits_num.push_back(hits_in_track.size());

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

        DoubleColl out_costheta;
        DoubleColl out_phi;
        DoubleColl out_pt;
        IntColl out_pdg;
      
        DoubleColl out_purity;
        DoubleColl out_eff;
        IntColl out_index_part;
        IntColl out_nHits;
        IntColl out_isReco_tracks;

        for (size_t i = 0; i < eff_mc.size(); ++i) {
            
            double pur = -1.0;
            int idx_assigned_part = -1;

            for (size_t s = 0; s < pur_tracks.size(); ++s) {
                
                double pur_temp;
                auto pur_it = pur_tracks[s].find(i);
                if (pur_it != pur_tracks[s].end()) {
                    pur_temp = pur_it->second;
                }

                if(pur_temp>pur)
                {
                    pur = pur_temp;
                    idx_assigned_part = s;
                }

            }

            // Find the efficiency
            double eff = 0.0;
            auto eff_it = eff_mc[i].find(idx_assigned_part);
            if (eff_it != eff_mc[i].end()) {
                eff = eff_it->second;

            }

            double cos = -1e3;
            double phi = -1e3;
            double pt = -1e3;
            int pdg = -1e3;
            double pur_ass = -1e3;
            double eff_ass = -1e3;
            int idx =-1e3;
            int number_of_hits = hits_num[i];
            int isReco_val = -1e3;

            // Check if purity and efficiency are above the thresholds
            if (pur > pur_thr__ && eff > eff_thr__)
            {   
                cos = costheta_mc[idx_assigned_part];
                phi = phi_mc[idx_assigned_part];
                pt = pt_mc[idx_assigned_part];
                pdg = pdg_mc[idx_assigned_part];
                pur_ass = pur;
                eff_ass = eff;
                idx = index_mc[idx_assigned_part];

                if(isReco[idx_assigned_part])
                {
                    isReco_val = 1;
                }
            }

            out_costheta.push_back(cos);
            out_phi.push_back(phi);
            out_pt.push_back(pt);
            out_pdg.push_back(pdg);
            out_purity.push_back(pur_ass);
            out_eff.push_back(eff_ass);
            out_index_part.push_back(idx);
            out_nHits.push_back(number_of_hits);
            out_isReco_tracks.push_back(isReco_val);

        }

        // IntColl out_isReco_particles = isReco;

        return std::make_tuple(std::move(out_costheta), std::move(out_phi), std::move(out_pt), std::move(out_pdg),std::move(out_purity),std::move(out_eff),std::move(out_index_part),std::move(out_nHits),std::move(out_isReco_tracks),std::move(isReco),std::move(pt_mc));
    }
    

};

DECLARE_COMPONENT(Efficiency_calc_pythia)