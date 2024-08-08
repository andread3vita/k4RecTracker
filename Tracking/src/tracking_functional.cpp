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
#include "dbscan.hpp"
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include <ATen/ATen.h>
#include <queue>
#include <typeinfo>

#include "Gaudi/Property.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHit3D.h"
#include "podio/UserDataCollection.h"
#include "k4FWCore/Transformer.h"

#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/TrackCollection.h"
#include "extension/MutableTrackerHit3D.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "extension/TrackerHit3DCollection.h"

#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep { using TrackerHit3DCollection = edm4hep::TrackerHitCollection; }  // namespace edm4hep
#endif


// Define collection types
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using ParticleColl = edm4hep::MCParticleCollection;
using DCTrackerHitColl = extension::DriftChamberDigiCollection;
using HitsColl = extension::TrackerHit3DCollection;
using DCTrackerHitColl_sim = edm4hep::SimTrackerHitCollection;
using VertexColl_sim = edm4hep::SimTrackerHitCollection;
using TrackColl = extension::TrackCollection;
using SimHits = edm4hep::SimTrackerHitCollection;

struct tracking_func final : 
        k4FWCore::MultiTransformer<std::tuple<TrackColl, HitsColl, DoubleColl, DoubleColl>(
                                                                                            const DCTrackerHitColl&, 
                                                                                            const HitsColl&,
                                                                                            const HitsColl&,
                                                                                            const HitsColl&,

                                                                                            const SimHits&,
                                                                                            const SimHits&,
                                                                                            const SimHits&,
                                                                                            const SimHits&)> 
            
                                                                                            
{
    tracking_func(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputHits_CDC", {"inputHits_CDC"}),
                KeyValues("inputHits_VTXIB", {"inputHits_VTXIB"}),
                KeyValues("inputHits_VTXD", {"inputHits_VTXD"}),
                KeyValues("inputHits_VTXOB", {"inputHits_VTXOB"}),

                KeyValues("inputHits_CDC_sim", {"inputHits_CDC_sim"}),
                KeyValues("inputHits_VTXIB_sim", {"inputHits_VTXIB_sim"}),
                KeyValues("inputHits_VTXD_sim", {"inputHits_VTXD_sim"}),
                KeyValues("inputHits_VTXOB_sim", {"inputHits_VTXOB_sim"})
            
            },
            {   
                
                KeyValues("outputTracks", {"outputTracks"}),
                KeyValues("outputHits", {"outputHits"}),
                KeyValues("clustering_space", {"clustering_space"}),
                KeyValues("clustering_space_tracks", {"clustering_space_tracks"})
            
            }
            ) {}
            

    std::tuple<TrackColl, HitsColl, DoubleColl, DoubleColl> operator()(
                                                                        const DCTrackerHitColl& inputHits_CDC, 
                                                                        const HitsColl& inputHits_VTXIB,
                                                                        const HitsColl& inputHits_VTXD,
                                                                        const HitsColl& inputHits_VTXOB,
                                                                        const SimHits& inputHits_CDC_sim,
                                                                        const SimHits& inputHits_VTXIB_sim,
                                                                        const SimHits& inputHits_VTXD_sim,
                                                                        const SimHits& inputHits_VTXOB_sim) const override 
    {

        /// Pointer to the ONNX enviroment
        std::unique_ptr<Ort::Env> fEnv;
        /// Pointer to the ONNX inference session
        std::unique_ptr<Ort::Session> fSession;
        /// ONNX settings
        Ort::SessionOptions fSessionOptions;
        /// ONNX memory info
        const OrtMemoryInfo* fInfo;
        struct MemoryInfo;

        /// the input names represent the names given to the model
        /// when defining  the model's architecture (if applicable)
        /// they can also be retrieved from model.summary()
        std::vector<const char*> fInames;
        std::vector<const char*> fOnames;

        std::string modelPath="/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx";

        fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

        auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
        fEnv          = std::move(envLocal);
        
        fSessionOptions.SetIntraOpNumThreads(1);
        fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

        // printf("Using Onnxruntime C++ API\n");
        auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.c_str(), fSessionOptions);	
        fSession          = std::move(sessionLocal);
        // printf("Starting to run inference\n");
        Ort::AllocatorWithDefaultOptions allocator;
        std::size_t i = 0;
        const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();
        const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

        fInames.push_back(input_name);
        fOnames.push_back(output_names);

        // For now there is a 1to1 correspondence between sim and digi so we can get the list of mc for each hit
        std::vector <float> ListHitMC_VTXD; 
        for (const auto input_sim_hit : inputHits_VTXD_sim) {
          auto MC_particle = input_sim_hit.getParticle();
          auto object_id_MC = MC_particle.getObjectID();
          auto index_MC = object_id_MC.index;
          ListHitMC_VTXD.push_back(index_MC);                 
        }
        std::vector <float> ListHitMC_VTXIB; 
        for (const auto input_sim_hit : inputHits_VTXIB_sim) {
          auto MC_particle = input_sim_hit.getParticle();
          auto object_id_MC = MC_particle.getObjectID();
          auto index_MC = object_id_MC.index;
          ListHitMC_VTXIB.push_back(index_MC);                 
        }
        std::vector <float> ListHitMC_VTXOB; 
        for (const auto input_sim_hit : inputHits_VTXOB_sim) {
          auto MC_particle = input_sim_hit.getParticle();
          auto object_id_MC = MC_particle.getObjectID();
          auto index_MC = object_id_MC.index;
          ListHitMC_VTXOB.push_back(index_MC);                 
        }
        std::vector <float> ListHitMC_CDC; 
        for (const auto input_sim_hit : inputHits_CDC_sim) {
          auto MC_particle = input_sim_hit.getParticle();
          auto object_id_MC = MC_particle.getObjectID();
          auto index_MC = object_id_MC.index;
          ListHitMC_CDC.push_back(index_MC);                 
        }

  

        // size_t size_total = size_CDC+size_VTXD+size_VTXIB+size_VTXOB;
        std::vector <float> ListGlobalInputs; 
        int it = 0;

        std::vector <float> ListHitType_VTXD;
        int it_0 = 0;
        for (const auto input_sim_hit : inputHits_VTXD) {
          ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
          ListGlobalInputs.push_back(1.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0); 
          ListHitType_VTXD.push_back(it);
          it += 1;  
          it_0 += 1;                        
        }
        torch::Tensor ListHitType_VTXD_tensor = torch::from_blob(ListHitType_VTXD.data(), {it_0}, torch::kFloat32);


        std::vector <float> ListHitType_VTXIB;
        int it_1 = 0;
        for (const auto input_sim_hit : inputHits_VTXIB) {
          ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
          ListGlobalInputs.push_back(1.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0); 
          ListHitType_VTXIB.push_back(it);
          it += 1; 
          it_1 += 1;                      
        }
        torch::Tensor ListHitType_VTXIB_tensor = torch::from_blob(ListHitType_VTXIB.data(), {it_1}, torch::kFloat32);


        std::vector <float> ListHitType_VTXOB;
        int it_2 = 0;
        for (const auto input_sim_hit : inputHits_VTXOB) {
          ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
          ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
          ListGlobalInputs.push_back(1.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(0.0);
          ListHitType_VTXOB.push_back(it);
          it += 1;  
          it_2 += 1;                   
        }
        torch::Tensor ListHitType_VTXOB_tensor = torch::from_blob(ListHitType_VTXOB.data(), {it_2}, torch::kFloat32);


        std::vector <float> ListHitType_CDC;
        int it_3 = 0;
        for (const auto input_sim_hit : inputHits_CDC) {
          ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().x);
          ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().y);
          ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().z);
          ListGlobalInputs.push_back(0.0);
          ListGlobalInputs.push_back(input_sim_hit.getRightPosition().x - input_sim_hit.getLeftPosition().x);
          ListGlobalInputs.push_back(input_sim_hit.getRightPosition().y - input_sim_hit.getLeftPosition().y);
          ListGlobalInputs.push_back(input_sim_hit.getRightPosition().z - input_sim_hit.getLeftPosition().z); 
          ListHitType_CDC.push_back(it);
          it += 1;    
          it_3 += 1;                     
        }
        torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_3}, torch::kFloat32);

        size_t total_size = it*7;
        std::vector<int64_t> tensor_shape = {it, 7};
        std::vector<Ort::Value> input_tensors;
        input_tensors.emplace_back( Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));

        auto output_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());

        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        std::vector<float> output_vector(floatarr, floatarr + it*4);

        DoubleColl output_cluster;
        for (int64_t t = 0; t < it; ++t) {
            for (int64_t j = 0; j < 4; ++j) {
              output_cluster.push_back(output_vector[t * 4 + j]);  
            }
        }

        torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {it, 4}, torch::kFloat32).clone(); 

        std::vector<point3w> points;

        for (int t = 0; t < output_model_tensor.size(0); ++t) {
            auto x = output_model_tensor[t][0].item<float>();
            auto y = output_model_tensor[t][1].item<float>();
            auto z = output_model_tensor[t][2].item<float>();
            auto weight = output_model_tensor[t][3].item<float>();
            points.push_back({x, y, z, weight});
        }

        float eps =  0.50; 
        int min_pts = 4;
        auto clusters = dbscan(points, eps, min_pts);

        std::vector<int> labels(points.size(), -1);
        for (size_t cluster_idx = 0; cluster_idx < clusters.size(); ++cluster_idx) {
            for (auto point_idx : clusters[cluster_idx]) {
                labels[point_idx] = static_cast<int>(cluster_idx);
            }
        }

        auto clustering = torch::from_blob(labels.data(), {static_cast<long>(labels.size())}, torch::kInt32).clone();
        DoubleColl output_cluster_track;
        for (int64_t t = 0; t < it; ++t) {
            
            auto el = clustering.index({t}).item<double>();
            output_cluster_track.push_back(el);  
  
        }

        torch::Tensor unique_tensor;
        torch::Tensor inverse_indices;
        std::tie(unique_tensor, inverse_indices)  = at::_unique(clustering, true, true);
        

        TrackColl output_tracks;
        HitsColl output_hits;
        int64_t number_of_tracks = unique_tensor.numel(); 
        for (int t = 0; t < number_of_tracks; ++t) {

            auto id_of_track = unique_tensor.index({t});
            auto output_track  = output_tracks.create();
            output_track.setChi2(1.);
            output_track.setNdf(1);
            output_track.setDEdx(1.);
            torch::Tensor mask = (clustering == id_of_track);
            torch::Tensor indices = torch::nonzero(mask);
            int64_t number_of_hits = indices.numel();

            for (int j = 0; j < number_of_hits; ++j) {
                
                auto index_id = indices.index({j});
                torch::Tensor mask_VTXD = (ListHitType_VTXD_tensor == index_id);
                torch::Tensor mask_VTXIB = (ListHitType_VTXIB_tensor == index_id);
                torch::Tensor mask_VTOB = (ListHitType_VTXOB_tensor == index_id);
                torch::Tensor mask_CDC = (ListHitType_CDC_tensor == index_id);

                if ((torch::sum(mask_VTXD)>0).item<bool>()){

                    auto hit = inputHits_VTXD.at(index_id.item<int>());
                    auto hit_extension  = output_hits.create();
                    hit_extension.setCellID(hit.getCellID());
                    hit_extension.setType(id_of_track.item<int>());
                    hit_extension.setEDep(ListHitMC_VTXD[index_id.item<int>()]);
                    hit_extension.setPosition(hit.getPosition());
                    output_track.addToTrackerHits(hit_extension);
                
                } else if ((torch::sum(mask_VTXIB)>0).item<bool>()){
                    index_id = index_id-it_0;
                    auto hit = inputHits_VTXIB.at(index_id.item<int>());
                    auto hit_extension  = output_hits.create();
                    hit_extension.setCellID(hit.getCellID());
                    hit_extension.setType(id_of_track.item<int>());
                    hit_extension.setEDep(ListHitMC_VTXIB[index_id.item<int>()]);
                    hit_extension.setPosition(hit.getPosition());
                    output_track.addToTrackerHits(hit_extension);

                } else if ((torch::sum(mask_VTOB)>0).item<bool>()){
                    index_id = index_id-(it_1+it_0);
                    auto hit = inputHits_VTXOB.at(index_id.item<int>());
                    auto hit_extension  = output_hits.create();
                    hit_extension.setCellID(hit.getCellID());
                    hit_extension.setType(id_of_track.item<int>());
                    hit_extension.setEDep(ListHitMC_VTXOB[index_id.item<int>()]);
                    hit_extension.setPosition(hit.getPosition());
                    output_track.addToTrackerHits(hit_extension);

                } else if ((torch::sum(mask_CDC)>0).item<bool>()){
                    index_id = index_id-(it_1+it_2 +it_0);
                    auto hit = inputHits_CDC.at(index_id.item<int>());
                    auto hit_extension  = output_hits.create();

                    hit_extension.setCellID(hit.getCellID());
                    hit_extension.setType(id_of_track.item<int>());
                    hit_extension.setEDep(ListHitMC_CDC[index_id.item<int>()]);
                    hit_extension.setPosition(hit.getLeftPosition());
                    output_track.addToTrackerHits(hit_extension);
                }
            }
        }

        return std::make_tuple(std::move(output_tracks), std::move(output_hits), std::move(output_cluster), std::move(output_cluster_track));

    }

};

DECLARE_COMPONENT(tracking_func)