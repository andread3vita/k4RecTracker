// /*
//  * Copyright (c) 2014-2024 Key4hep-Project.
//  *
//  * This file is part of Key4hep.
//  * See https://key4hep.github.io/key4hep-doc/ for further info.
//  *
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  *
//  *     http://www.apache.org/licenses/LICENSE-2.0
//  *
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  */

// #include <algorithm>
// #include <cmath>
// #include <iostream>
// #include <map>
// #include <numeric>
// #include <queue>
// #include <sstream>
// #include <string>
// #include <typeinfo>
// #include <vector>

// #include <ATen/ATen.h>
// #include <torch/torch.h>
// #include "onnxruntime_cxx_api.h"

// #include "dbscan.hpp"
// #include "utils.hpp"

// #include "Gaudi/Property.h"
// #include "edm4hep/MCParticleCollection.h"
// #include "edm4hep/SimTrackerHitCollection.h"
// #include "edm4hep/TrackCollection.h"
// #include "edm4hep/TrackerHit3D.h"
// #include "podio/UserDataCollection.h"
// #include "k4FWCore/Transformer.h"

// #include "extension/DriftChamberDigiCollection.h"
// #include "extension/DriftChamberDigiLocalCollection.h"
// #include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
// #include "extension/TrackCollection.h"
// #include "extension/MutableTrackerHit3D.h"
// #include "extension/TrackerHit3DCollection.h"
// #include "edm4hep/SimTrackerHitCollection.h"

// #if __has_include("edm4hep/TrackerHit3DCollection.h")
// #include "edm4hep/TrackerHit3DCollection.h"
// #else
// #include "edm4hep/TrackerHitCollection.h"
// namespace edm4hep { using TrackerHit3DCollection = edm4hep::TrackerHitCollection; }  // namespace edm4hep
// #endif

// // Define collection types
// using DoubleColl = podio::UserDataCollection<double>;
// using IntColl = podio::UserDataCollection<int>;
// using ParticleColl = edm4hep::MCParticleCollection;
// using DCTrackerHitColl = extension::DriftChamberDigiCollection;
// using HitsColl = extension::TrackerHit3DCollection;
// using DCTrackerHitColl_sim = edm4hep::SimTrackerHitCollection;
// using VertexColl_sim = edm4hep::SimTrackerHitCollection;
// using TrackColl = extension::TrackCollection;
// using SimHits = edm4hep::SimTrackerHitCollection;


// /** @struct GGTF_tracking
//  *
//  *  Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the GGTF_tracking. 
//  *  The first step takes the raw hits and it returns a collection of 4-dimensional points inside an embedding space.
//  *  Eeach 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described intuitively by a potential, 
//  *  which attracts hits belonging to the same cluster and drives away those that do not.
//  *  This collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the same track.
//  *
//  *  input: 
//  *    - MC hits from DC and vertex : SimHits
//  *    - digitalized hits from DC (global coordinates) : DCTrackerHitColl 
//  *    - digitalized hits from vertex (global coordinates) : HitsColl
//  *
//  *  output:
//  *    - Clustering labels : IntColl -> VTXD , VTXIB , VTXOB , CDC
//  *    - Track collection : TrackColl
//  *
//  *
//  *
//  *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
//  *  @date   2024-08
//  *
//  */

// struct GGTF_tracking final : 
//         k4FWCore::MultiTransformer<std::tuple<IntColl, TrackColl>(
//                                                                     const DCTrackerHitColl&, 
//                                                                     const HitsColl&,
//                                                                     const HitsColl&,
//                                                                     const HitsColl&,

//                                                                     const SimHits&,
//                                                                     const SimHits&,
//                                                                     const SimHits&,
//                                                                     const SimHits&)> 
            
                                                                                            
// {
//     GGTF_tracking(const std::string& name, ISvcLocator* svcLoc) : 
//         MultiTransformer ( name, svcLoc,
//             {   
                 
//                 KeyValues("inputHits_CDC", {"inputHits_CDC"}),
//                 KeyValues("inputHits_VTXIB", {"inputHits_VTXIB"}),
//                 KeyValues("inputHits_VTXD", {"inputHits_VTXD"}),
//                 KeyValues("inputHits_VTXOB", {"inputHits_VTXOB"}),

//                 KeyValues("inputHits_CDC_sim", {"inputHits_CDC_sim"}),
//                 KeyValues("inputHits_VTXIB_sim", {"inputHits_VTXIB_sim"}),
//                 KeyValues("inputHits_VTXD_sim", {"inputHits_VTXD_sim"}),
//                 KeyValues("inputHits_VTXOB_sim", {"inputHits_VTXOB_sim"})
            
//             },
//             {   

//                 KeyValues("clustering_space_tracks", {"clustering_space_tracks"}),
//                 KeyValues("outputTracks", {"outputTracks"})      
            
//             }
//             ) {}
            
//     StatusCode initialize() {
      
//         // Initialize the ONNX memory info object for CPU memory allocation.
//         // This specifies that the memory will be allocated using the Arena Allocator on the CPU.
//         fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

//         // Create and initialize the ONNX environment with a logging level set to WARNING.
//         // This environment handles logging and runtime configuration for the ONNX session.
//         auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
//         fEnv          = std::move(envLocal);

//         // Set the session options to configure the ONNX inference session.
//         // Set the number of threads used for intra-op parallelism to 1 (single-threaded execution).
//         fSessionOptions.SetIntraOpNumThreads(1);

//         // Disable all graph optimizations to keep the model execution as close to the original as possible.
//         fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

//         // Create an ONNX inference session using the configured environment and session options.
//         // The session is used to load the model specified by the `modelPath`.
//         auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);
//         fSession          = std::move(sessionLocal);

//         // Create an ONNX allocator with default options to manage memory allocations during runtime.
//         Ort::AllocatorWithDefaultOptions allocator;

//         // Get the name of the first input node (index i) from the ONNX model and store it.
//         // This retrieves the name from the model and releases the memory after storing it in fInames.
//         std::size_t i = 0;
//         const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();

//         // Get the name of the first output node (index i) from the ONNX model and store it.
//         // This retrieves the name from the model and releases the memory after storing it in fOnames.
//         const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

//         // Store the retrieved input and output names in the respective vectors.
//         fInames.push_back(input_name);
//         fOnames.push_back(output_names);
      
      
//         return StatusCode::SUCCESS;

//     }


//     std::tuple<IntColl,TrackColl> operator()(
//                                                 const DCTrackerHitColl& inputHits_CDC, 
//                                                 const HitsColl& inputHits_VTXIB,
//                                                 const HitsColl& inputHits_VTXD,
//                                                 const HitsColl& inputHits_VTXOB,
//                                                 const SimHits& inputHits_CDC_sim,
//                                                 const SimHits& inputHits_VTXIB_sim,
//                                                 const SimHits& inputHits_VTXD_sim,
//                                                 const SimHits& inputHits_VTXOB_sim) const override 
//     {

//         ////////////////////////////////////////
//         ////////// DATA PREPROCESSING //////////
//         ////////////////////////////////////////
        
//         // Although there is a link between digi hits and MC hits, we cannot use it because the output tracks 
//         // do not contain the original hits, but rather copies. This step is necessary because it is not possible 
//         // to use hits from both the drift chamber and the vertex in the same collection. Therefore, for each hit, 
//         // a MutableTrackerHit3D is created, and its properties are defined based on the original hit.
//         //
//         // These for loops are used to store the index of the MC particle for each hit. For now, this information 
//         // will be saved in the EDep of the hit stored in each track.

//         // Vector to store MC particle indices for the VTXD detector hits
//         std::vector <float> ListHitMC_VTXD; 
//         for (const auto input_sim_hit : inputHits_VTXD_sim) {
//             auto MC_particle = input_sim_hit.getParticle(); // Retrieve the associated MC particle
//             auto object_id_MC = MC_particle.getObjectID();  // Get the object ID of the MC particle
//             auto index_MC = object_id_MC.index;             // Extract the index from the object ID
//             ListHitMC_VTXD.push_back(index_MC);             // Store the index in the vector
//         }

//         // Vector to store MC particle indices for the VTXIB detector hits
//         std::vector <float> ListHitMC_VTXIB; 
//         for (const auto input_sim_hit : inputHits_VTXIB_sim) {
//             auto MC_particle = input_sim_hit.getParticle(); // Retrieve the associated MC particle
//             auto object_id_MC = MC_particle.getObjectID();  // Get the object ID of the MC particle
//             auto index_MC = object_id_MC.index;             // Extract the index from the object ID
//             ListHitMC_VTXIB.push_back(index_MC);            // Store the index in the vector
//         }

//         // Vector to store MC particle indices for the VTXOB detector hits
//         std::vector <float> ListHitMC_VTXOB; 
//         for (const auto input_sim_hit : inputHits_VTXOB_sim) {
//             auto MC_particle = input_sim_hit.getParticle(); // Retrieve the associated MC particle
//             auto object_id_MC = MC_particle.getObjectID();  // Get the object ID of the MC particle
//             auto index_MC = object_id_MC.index;             // Extract the index from the object ID
//             ListHitMC_VTXOB.push_back(index_MC);            // Store the index in the vector
//         }

//         // Vector to store MC particle indices for the CDC detector hits
//         std::vector <float> ListHitMC_CDC; 
//         for (const auto input_sim_hit : inputHits_CDC_sim) {
//             auto MC_particle = input_sim_hit.getParticle(); // Retrieve the associated MC particle
//             auto object_id_MC = MC_particle.getObjectID();  // Get the object ID of the MC particle
//             auto index_MC = object_id_MC.index;             // Extract the index from the object ID
//             ListHitMC_CDC.push_back(index_MC);              // Store the index in the vector
//         }

        
//         // Vector to store the global input values for all hits.
//         // This will contain position and other hit-specific data to be used as input for the model.
//         std::vector<float> ListGlobalInputs; 
//         int it = 0;  // Iterator to keep track of the index in ListGlobalInputs.

//         /// Processing hits from the VTXD (Vertex Detector).
//         std::vector<float> ListHitType_VTXD; // Vector to store hit indices for VTXD.
//         int it_0 = 0;  // Iterator to keep track of the number of VTXD hits.
//         for (const auto input_hit : inputHits_VTXD) {
//             // Add the 3D position of the hit to the global input list.
//             ListGlobalInputs.push_back(input_hit.getPosition().x);
//             ListGlobalInputs.push_back(input_hit.getPosition().y);
//             ListGlobalInputs.push_back(input_hit.getPosition().z);
            
//             // Add placeholder values for additional input dimensions.
//             ListGlobalInputs.push_back(1.0); // Placeholder for detector type label (1 for vertex detector, 0 for drift chamber)
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0); 
            
//             // Store the current index in ListHitType_VTXD and increment the global iterator.
//             ListHitType_VTXD.push_back(it);
//             it += 1;  
//             it_0 += 1;                        
//         }

//         // Convert ListHitType_VTXD to a Torch tensor for use in PyTorch models.
//         torch::Tensor ListHitType_VTXD_tensor = torch::from_blob(ListHitType_VTXD.data(), {it_0}, torch::kFloat32);

//         /// Processing hits from the VTXIB (Vertex Inner Barrel).
//         std::vector<float> ListHitType_VTXIB; // Vector to store hit indices for VTXIB.
//         int it_1 = 0;  // Iterator to keep track of the number of VTXIB hits.
//         for (const auto input_hit : inputHits_VTXIB) {
//             // Add the 3D position of the hit to the global input list.
//             ListGlobalInputs.push_back(input_hit.getPosition().x);
//             ListGlobalInputs.push_back(input_hit.getPosition().y);
//             ListGlobalInputs.push_back(input_hit.getPosition().z);
            
//             // Add placeholder values for additional input dimensions.
//             ListGlobalInputs.push_back(1.0); // Placeholder for detector type label (1 for vertex detector, 0 for drift chamber)
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0); 
            
//             // Store the current index in ListHitType_VTXIB and increment the global iterator.
//             ListHitType_VTXIB.push_back(it);
//             it += 1; 
//             it_1 += 1;                      
//         }

//         // Convert ListHitType_VTXIB to a Torch tensor for use in PyTorch models.
//         torch::Tensor ListHitType_VTXIB_tensor = torch::from_blob(ListHitType_VTXIB.data(), {it_1}, torch::kFloat32);

//         /// Processing hits from the VTXOB (Vertex Outer Barrel).
//         std::vector<float> ListHitType_VTXOB; // Vector to store hit indices for VTXOB.
//         int it_2 = 0;  // Iterator to keep track of the number of VTXOB hits.
//         for (const auto input_hit : inputHits_VTXOB) {
//             // Add the 3D position of the hit to the global input list.
//             ListGlobalInputs.push_back(input_hit.getPosition().x);
//             ListGlobalInputs.push_back(input_hit.getPosition().y);
//             ListGlobalInputs.push_back(input_hit.getPosition().z);
            
//             // Add placeholder values for additional input dimensions.
//             ListGlobalInputs.push_back(1.0); // Placeholder for detector type label (1 for vertex detector, 0 for drift chamber)
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0);
//             ListGlobalInputs.push_back(0.0); 
            
//             // Store the current index in ListHitType_VTXOB and increment the global iterator.
//             ListHitType_VTXOB.push_back(it);
//             it += 1;  
//             it_2 += 1;                   
//         }

//         // Convert ListHitType_VTXOB to a Torch tensor for use in PyTorch models.
//         torch::Tensor ListHitType_VTXOB_tensor = torch::from_blob(ListHitType_VTXOB.data(), {it_2}, torch::kFloat32);

//         /// Processing hits from the CDC (Central Drift Chamber).
//         std::vector<float> ListHitType_CDC; // Vector to store hit indices for CDC.
//         int it_3 = 0;  // Iterator to keep track of the number of CDC hits.
//         for (const auto input_hit : inputHits_CDC) {
//             // Add the 3D position of the left hit to the global input list.
//             ListGlobalInputs.push_back(input_hit.getLeftPosition().x);
//             ListGlobalInputs.push_back(input_hit.getLeftPosition().y);
//             ListGlobalInputs.push_back(input_hit.getLeftPosition().z);
            
//             // Add the difference between the right and left hit positions to the global input list.
//             ListGlobalInputs.push_back(0.0); // Placeholder for detector type label (1 for vertex detector, 0 for drift chamber)
//             ListGlobalInputs.push_back(input_hit.getRightPosition().x - input_hit.getLeftPosition().x);
//             ListGlobalInputs.push_back(input_hit.getRightPosition().y - input_hit.getLeftPosition().y);
//             ListGlobalInputs.push_back(input_hit.getRightPosition().z - input_hit.getLeftPosition().z); 
            
//             // Store the current index in ListHitType_CDC and increment the global iterator.
//             ListHitType_CDC.push_back(it);
//             it += 1;    
//             it_3 += 1;                     
//         }

//         // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
//         torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_3}, torch::kFloat32);

//         /////////////////////////////
//         ////////// ML STEP //////////
//         /////////////////////////////

//         // Calculate the total size of the input tensor, based on the number of hits (it) and the 
//         // number of features per hit (7: x, y, z, and four placeholders).
//         size_t total_size = it * 7;

//         // Define the shape of the input tensor. The tensor will have `it` rows (one for each hit)
//         // and 7 columns (representing the features of each hit).
//         std::vector<int64_t> tensor_shape = {it, 7};

//         // Create a vector to store the input tensors that will be fed into the ONNX model.
//         std::vector<Ort::Value> input_tensors;

//         // Create an ONNX tensor from the input data (ListGlobalInputs) and add it to the input_tensors vector.
//         // The tensor is created using the memory information (fInfo), the data pointer, the total size,
//         // and the tensor shape defined earlier.
//         input_tensors.emplace_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));

//         // Run the ONNX inference session with the provided input tensor.
//         // The model will process the input data and produce output tensors.
//         // The input tensor names (fInames), input tensors (input_tensors), output tensor names (fOnames), 
//         // and the sizes of these vectors are passed to the Run method.
//         auto output_model_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());

//         // Extract the raw output data from the first output tensor returned by the model.
//         // This data is stored as a float array.
//         float* floatarr = output_model_tensors.front().GetTensorMutableData<float>();

//         // Convert the raw output data into a std::vector of floats for easier handling.
//         // The size of the vector is determined by the number of hits (it) and the number of 
//         // output features per hit (3 coordinates and charge).
//         std::vector<float> output_model_vector(floatarr, floatarr + it * 4);

//         /////////////////////////////////////
//         ////////// CLUSTERING STEP //////////
//         /////////////////////////////////////

//         // Perform clustering on the data using the get_clustering function. 
//         // The result is a tensor where each element represents the cluster assignment for each point.
//         auto clustering = get_clustering(output_model_vector, it, 0.6, 0.05);

//         // Create a container to store the cluster labels and Iterate over each point in the clustering result.
//         IntColl labels;
//         for (int64_t t = 0; t < it; ++t) {
            
//             // Retrieve the cluster label for the current point and store it in 'el'.
//             auto el = clustering.index({t}).item<int>();
            
//             // Add the label to the 'labels' collection.
//             labels.push_back(el);  
//         }

//         // Declare tensors to hold the unique cluster labels and their corresponding inverse indices.
//         torch::Tensor unique_tensor;
//         torch::Tensor inverse_indices;
//         // Compute the unique cluster labels and the inverse indices, which map each original index to its unique label.
//         std::tie(unique_tensor, inverse_indices) = at::_unique(clustering, true, true);
        

//         /////////////////////////////////
//         ////////// OUTPUT STEP //////////
//         /////////////////////////////////

//         // Initialize the collections to store output tracks and hits
//         TrackColl output_tracks;
//         HitsColl output_hits;

//         // Determine the total number of unique tracks to process
//         int64_t number_of_tracks = unique_tensor.numel(); 

//         // Loop over each unique track
//         for (int i = 0; i < number_of_tracks; ++i) {
//             // Retrieve the track ID from the unique tensor
//             auto id_of_track = unique_tensor.index({i});
//             // Create a new track in the output tracks collection
//             auto output_track  = output_tracks->create();

//             // Set global properties for the track
//             output_track.setChi2(1.);
//             output_track.setNdf(1);
//             output_track.setDEdx(id_of_track.item<int>());

//             // Create a mask to identify hits belonging to the current track ID
//             torch::Tensor mask = (clustering == id_of_track);
//             torch::Tensor indices = torch::nonzero(mask);
//             int64_t number_of_hits = indices.numel();  // Determine the number of hits associated with this track

//             // Loop over each hit associated with the current track
//             for (int j = 0; j < number_of_hits; ++j) {
//                 // Retrieve the hit index
//                 auto index_id = indices.index({j});

//                 // Create masks to check which detector (VTXD, VTXIB, VTXOB, or CDC) the hit belongs to
//                 torch::Tensor mask_VTXD = (ListHitType_VTXD_tensor == index_id);
//                 torch::Tensor mask_VTXIB = (ListHitType_VTXIB_tensor == index_id);
//                 torch::Tensor mask_VTOB = (ListHitType_VTXOB_tensor == index_id);
//                 torch::Tensor mask_CDC = (ListHitType_CDC_tensor == index_id);

//                 // If the hit belongs to VTXD
//                 if ((torch::sum(mask_VTXD) > 0).item<bool>()){
//                     // Retrieve the hit from VTXD
//                     auto hit = inputHits_VTXD.at(index_id.item<int>());

//                     // Create an hit in the output hits collection
//                     auto hit_extension  = output_hits.create();
//                     // Set properties for the hit
//                     hit_extension.setCellID(hit.getCellID());
//                     hit_extension.setType(1);
//                     hit_extension.setEDep(ListHitMC_VTXD[index_id.item<int>()]); // EDep is the index of the MC particle
//                     hit_extension.setPosition(hit.getPosition());
//                     // Add the hit to the current track
//                     output_track.addToTrackerHits(hit_extension);

//                 // If the hit belongs to VTXIB
//                 } else if ((torch::sum(mask_VTXIB) > 0).item<bool>()){
//                     // Adjust the index and retrieve the hit from VTXIB
//                     index_id = index_id - it_0;
//                     auto hit = inputHits_VTXIB.at(index_id.item<int>());
//                     auto hit_extension  = output_hits.create();
//                     // Set properties for the hit
//                     hit_extension.setCellID(hit.getCellID());
//                     hit_extension.setType(1);
//                     hit_extension.setEDep(ListHitMC_VTXIB[index_id.item<int>()]); // EDep is the index of the MC particle
//                     hit_extension.setPosition(hit.getPosition());
//                     // Add the hit to the current track
//                     output_track.addToTrackerHits(hit_extension);

//                 // If the hit belongs to VTXOB
//                 } else if ((torch::sum(mask_VTOB) > 0).item<bool>()){
//                     // Adjust the index and retrieve the hit from VTXOB
//                     index_id = index_id - (it_1 + it_0);
//                     auto hit = inputHits_VTXOB.at(index_id.item<int>());
//                     auto hit_extension  = output_hits.create();
//                     // Set properties for the extended hit
//                     hit_extension.setCellID(hit.getCellID());
//                     hit_extension.setType(1);
//                     hit_extension.setEDep(ListHitMC_VTXOB[index_id.item<int>()]); // EDep is the index of the MC particle
//                     hit_extension.setPosition(hit.getPosition());
//                     // Add the extended hit to the current track
//                     output_track.addToTrackerHits(hit_extension);

//                 // If the hit belongs to the CDC
//                 } else if ((torch::sum(mask_CDC) > 0).item<bool>()){
//                     // Adjust the index and retrieve the hit from CDC
//                     index_id = index_id - (it_1 + it_2 + it_0);
//                     auto hit = inputHits_CDC.at(index_id.item<int>());
//                     auto hit_extension  = output_hits.create();
//                     // Set properties for the extended hit
//                     hit_extension.setCellID(hit.getCellID());
//                     hit_extension.setType(0);
//                     hit_extension.setEDep(ListHitMC_CDC[index_id.item<int>()]); // EDep is the index of the MC particle
//                     hit_extension.setPosition(hit.getLeftPosition());
//                     // Add the extended hit to the current track
//                     output_track.addToTrackerHits(hit_extension);
//                 }
//             }
//         }


//         // return std::make_tuple(std::move(output_tracks), std::move(output_hits), std::move(output_cluster), std::move(output_cluster_track));
//         return std::make_tuple(std::move(labels),std::move(output_tracks));
//     } 

//     private:

//         /// Pointer to the ONNX environment.
//         /// This object manages the global state of the ONNX runtime, such as logging and threading.
//         std::unique_ptr<Ort::Env> fEnv;

//         /// Pointer to the ONNX inference session.
//         /// This session is used to execute the model for inference.
//         std::unique_ptr<Ort::Session> fSession;

//         /// ONNX session options.
//         /// These settings control the behavior of the inference session, such as optimization level, 
//         /// execution providers, and other configuration parameters.
//         Ort::SessionOptions fSessionOptions;

//         /// ONNX memory info.
//         /// This object provides information about memory allocation and is used during the creation of 
//         /// ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
//         const OrtMemoryInfo* fInfo;

//         /// Stores the input and output names for the ONNX model.
//         /// These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
//         std::vector<const char*> fInames;
//         std::vector<const char*> fOnames;

//         /// Property to specify the path to the ONNX model file.
//         /// This is a configurable property that defines the location of the ONNX model file on the filesystem.
//         Gaudi::Property< std::string > modelPath{this, "modelPath", "/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx", "modelPath"};


// };

// DECLARE_COMPONENT(GGTF_tracking)