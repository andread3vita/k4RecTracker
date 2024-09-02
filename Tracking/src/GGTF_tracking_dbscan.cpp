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
#include <memory> 
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include <fstream>  // Per std::ifstream
#include <vector>   // Per std::vector
#include <iterator> // Per std::istreambuf_iterator

#include <iostream>
#include <typeinfo>
#include <cxxabi.h>
#include <memory>

#include <ATen/ATen.h>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"

#include "dbscan.hpp"
#include "utils.hpp"

#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"


// Define collection types
#include "podio/UserDataCollection.h"
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;

#include "edm4hep/MCParticleCollection.h"
using ParticleColl = edm4hep::MCParticleCollection;

#include "extension/TrackerHit3DCollection.h"
using VertexHitsColl = extension::TrackerHit3DCollection;

#include "extension/DriftChamberDigiV2Collection.h"
using DCHitColl = edm4hep::DriftChamberDigiV2Collection;

#include "extension/TrackCollection.h"
using TrackColl = extension::TrackCollection;

struct GGTF_tracking_dbscan final : 
        k4FWCore::MultiTransformer< std::tuple<TrackColl>( 
                                                                    const DCHitColl&, 
                                                                    const VertexHitsColl&,
                                                                    const VertexHitsColl&,
                                                                    const VertexHitsColl&)> 
            
                                                                                            
{
    GGTF_tracking_dbscan(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputHits_CDC", {"inputHits_CDC"}),
                KeyValues("inputHits_VTXIB", {"inputHits_VTXIB"}),
                KeyValues("inputHits_VTXD", {"inputHits_VTXD"}),
                KeyValues("inputHits_VTXOB", {"inputHits_VTXOB"})

            },
            {   
                KeyValues("outputTracks", {"outputTracks"})        
            
            }) {}
    
    StatusCode initialize() {

        // Leggi il file .onnx in modalità binaria
        std::ifstream model_file(modelPath.value().c_str(), std::ios::binary | std::ios::ate);
        if (!model_file) {
            std::cerr << "Failed to open the ONNX model file: " << modelPath << std::endl;
        }

        // Ottieni la dimensione del file
        std::streamsize model_size = model_file.tellg();
        model_file.seekg(0, std::ios::beg);

        // Crea un buffer per contenere il file modello
        std::vector<char> model_data_temp(model_size);
        if (!model_file.read(model_data_temp.data(), model_size)) {
            std::cerr << "Failed to read the ONNX model file: " << modelPath << std::endl;
        }

        model_data = model_data_temp;

        // Initialize the ONNX memory info object for CPU memory allocation.
        // This specifies that the memory will be allocated using the Arena Allocator on the CPU.
        fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

        // Create and initialize the ONNX environment with a logging level set to WARNING.
        // This environment handles logging and runtime configuration for the ONNX session.
        auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
        fEnv          = std::move(envLocal);

        // Set the session options to configure the ONNX inference session.
        // Set the number of threads used for intra-op parallelism to 1 (single-threaded execution).
        fSessionOptions.SetIntraOpNumThreads(1);

        // Disable all graph optimizations to keep the model execution as close to the original as possible.
        fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

        fSessionOptions.DisableMemPattern();

        // Create an ONNX inference session using the configured environment and session options.
        // The session is used to load the model specified by the `modelPath`.
        auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);
        fSession          = std::move(sessionLocal);

        // Create an ONNX allocator with default options to manage memory allocations during runtime.
        Ort::AllocatorWithDefaultOptions allocator;

        // Get the name of the first input node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fInames.
        std::size_t i = 0;
        const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();

        // Get the name of the first output node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fOnames.
        const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

        // Store the retrieved input and output names in the respective vectors.
        fInames.push_back(input_name);
        fOnames.push_back(output_names);


        return StatusCode::SUCCESS;

    }

    
    std::tuple<TrackColl> operator()(
                                                const DCHitColl& inputHits_CDC, 
                                                const VertexHitsColl& inputHits_VTXIB,
                                                const VertexHitsColl& inputHits_VTXD,
                                                const VertexHitsColl& inputHits_VTXOB) const override 
    {

        torch::NoGradGuard no_grad;

        ////////////////////////////////////////
        ////////// DATA PREPROCESSING //////////
        ////////////////////////////////////////
        
        // Although there is a link between digi hits and MC hits, we cannot use it because the output tracks 
        // do not contain the original hits, but rather copies. This step is necessary because it is not possible 
        // to use hits from both the drift chamber and the vertex in the same collection. Therefore, for each hit, 
        // a MutableTrackerHit3D is created, and its properties are defined based on the original hit.
        //
        // These for loops are used to store the index of the MC particle for each hit. For now, this information 
        // will be saved in the EDep of the hit stored in each track.

        std::cout << "Input Hit collection size VTXD: " << inputHits_VTXD.size() << std::endl;
        std::cout << "Input Hit collection size VTXIB: " << inputHits_VTXIB.size() << std::endl;
        std::cout << "Input Hit collection size VTXOB: " << inputHits_VTXOB.size() << std::endl;
        std::cout << "Input Hit collection size CDC: " << inputHits_CDC.size() << std::endl;
        std::cout << "______________________________________________________" << std::endl;
        std::cout << "Tot hits: " << inputHits_CDC.size() + inputHits_VTXOB.size() + inputHits_VTXIB.size() + inputHits_VTXD.size() << std::endl;
        std::cout << "  " << std::endl;

        // Vector to store the global input values for all hits.
        // This will contain position and other hit-specific data to be used as input for the model.
        std::vector<float> ListGlobalInputs; 
        int it = 0;  // Iterator to keep track of the index in ListGlobalInputs.

        /// Processing hits from the VTXD (Vertex Detector).
        std::vector<float> ListHitType_VTXD; // Vector to store hit indices for VTXD.
        int it_0 = 0;  // Iterator to keep track of the number of VTXD hits.
        for (const auto input_hit : inputHits_VTXD) {
            // Add the 3D position of the hit to the global input list.
            ListGlobalInputs.push_back(input_hit.getPosition().x);
            ListGlobalInputs.push_back(input_hit.getPosition().y);
            ListGlobalInputs.push_back(input_hit.getPosition().z);
            
            // Add placeholder values for additional input dimensions.
            ListGlobalInputs.push_back(1.0); 
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0); 
            
            // Store the current index in ListHitType_VTXD and increment the global iterator.
            ListHitType_VTXD.push_back(it);
            it += 1;  
            it_0 += 1;                        
        }

        // Convert ListHitType_VTXD to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_VTXD_tensor = torch::from_blob(ListHitType_VTXD.data(), {it_0}, torch::kFloat32);

        /// Processing hits from the VTXIB (Vertex Inner Barrel).
        std::vector<float> ListHitType_VTXIB; // Vector to store hit indices for VTXIB.
        int it_1 = 0;  // Iterator to keep track of the number of VTXIB hits.
        for (const auto input_hit : inputHits_VTXIB) {
            // Add the 3D position of the hit to the global input list.
            ListGlobalInputs.push_back(input_hit.getPosition().x);
            ListGlobalInputs.push_back(input_hit.getPosition().y);
            ListGlobalInputs.push_back(input_hit.getPosition().z);
            
            // Add placeholder values for additional input dimensions.
            ListGlobalInputs.push_back(1.0);
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0); 
            
            // Store the current index in ListHitType_VTXIB and increment the global iterator.
            ListHitType_VTXIB.push_back(it);
            it += 1; 
            it_1 += 1;                      
        }

        // Convert ListHitType_VTXIB to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_VTXIB_tensor = torch::from_blob(ListHitType_VTXIB.data(), {it_1}, torch::kFloat32);

        /// Processing hits from the VTXOB (Vertex Outer Barrel).
        std::vector<float> ListHitType_VTXOB; // Vector to store hit indices for VTXOB.
        int it_2 = 0;  // Iterator to keep track of the number of VTXOB hits.
        for (const auto input_hit : inputHits_VTXOB) {
            // Add the 3D position of the hit to the global input list.
            ListGlobalInputs.push_back(input_hit.getPosition().x);
            ListGlobalInputs.push_back(input_hit.getPosition().y);
            ListGlobalInputs.push_back(input_hit.getPosition().z);
            
            // Add placeholder values for additional input dimensions.
            ListGlobalInputs.push_back(1.0); 
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0);
            ListGlobalInputs.push_back(0.0); 
            
            // Store the current index in ListHitType_VTXOB and increment the global iterator.
            ListHitType_VTXOB.push_back(it);
            it += 1;  
            it_2 += 1;                   
        }

        // Convert ListHitType_VTXOB to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_VTXOB_tensor = torch::from_blob(ListHitType_VTXOB.data(), {it_2}, torch::kFloat32);

        // /// Processing hits from the CDC (Central Drift Chamber).
        std::vector<float> ListHitType_CDC; // Vector to store hit indices for CDC.
        int it_3 = 0;  // Iterator to keep track of the number of CDC hits.
        // for (const auto input_hit : inputHits_CDC) {
        //     // Add the 3D position of the left hit to the global input list.
        //     ListGlobalInputs.push_back(input_hit.getLeftPosition().x);
        //     ListGlobalInputs.push_back(input_hit.getLeftPosition().y);
        //     ListGlobalInputs.push_back(input_hit.getLeftPosition().z);
            
        //     // Add the difference between the right and left hit positions to the global input list.
        //     ListGlobalInputs.push_back(0.0); 
        //     ListGlobalInputs.push_back(input_hit.getRightPosition().x - input_hit.getLeftPosition().x);
        //     ListGlobalInputs.push_back(input_hit.getRightPosition().y - input_hit.getLeftPosition().y);
        //     ListGlobalInputs.push_back(input_hit.getRightPosition().z - input_hit.getLeftPosition().z); 
            
        //     // Store the current index in ListHitType_CDC and increment the global iterator.
        //     ListHitType_CDC.push_back(it);
        //     it += 1;    
        //     it_3 += 1;                     
        // }

        // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_3}, torch::kFloat32);

        /////////////////////////////
        ////////// ML STEP //////////
        /////////////////////////////

         // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
        TrackColl* output_tracks = new TrackColl;
        std::vector<int> labels(it, -1);
        if (it > 0 && it < 20000)
        {
            // Calculate the total size of the input tensor, based on the number of hits (it) and the 
            // number of features per hit (7: x, y, z, and four placeholders).
            size_t total_size = it * 7;

            // Define the shape of the input tensor. The tensor will have `it` rows (one for each hit)
            // and 7 columns (representing the features of each hit).
            std::vector<int64_t> tensor_shape = {it, 7};

            // // Create a vector to store the input tensors that will be fed into the ONNX model.
            std::vector<Ort::Value> input_tensors;

            // // Create an ONNX tensor from the input data (ListGlobalInputs) and add it to the input_tensors vector.
            // // The tensor is created using the memory information (fInfo), the data pointer, the total size,
            // // and the tensor shape defined earlier.
            // input_tensors.emplace_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));
            
            input_tensors.push_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));

            // // Run the ONNX inference session with the provided input tensor.
            // // The model will process the input data and produce output tensors.
            // // The input tensor names (fInames), input tensors (input_tensors), output tensor names (fOnames), 
            // // and the sizes of these vectors are passed to the Run method.
            auto output_model_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());
            
            // Extract the raw output data from the first output tensor returned by the model.
            // This data is stored as a float array.
            double* floatarr = output_model_tensors.front().GetTensorMutableData<double>();

            // Convert the raw output data into a std::vector of floats for easier handling.
            // The size of the vector is determined by the number of hits (it) and the number of 
            // output features per hit (3 coordinates and charge).
            std::vector<double> output_model_vector(floatarr, floatarr + it * 4);

            /////////////////////////////////////
            ////////// CLUSTERING STEP //////////
            /////////////////////////////////////

            // Convert the output vector from the model to a Torch tensor, specifying the shape and data type
            torch::Tensor output_model_tensor = torch::from_blob(output_model_vector.data(), {it, 4}, torch::kFloat32).clone();    

            // Initialize a vector to store the 3D points with an associated weight
            std::vector<point3w> points;
            for (int i = 0; i < output_model_tensor.size(0); ++i) {
                auto x = output_model_tensor[i][0].item<double>();
                auto y = output_model_tensor[i][1].item<double>();
                auto z = output_model_tensor[i][2].item<double>();
                auto weight = output_model_tensor[i][3].item<double>();
                
                // Store the point in the points vector
                points.push_back({x, y, z, weight});
            }


            // Apply the DBSCAN clustering algorithm to the points with the given step size and minimum points per cluster
            // NOTE: dbscan doesn't take into account the weights, so the fourth entry of each point will be ignored.
            auto clusters = dbscan(points, step_size, min_points);

            // Initialize a vector to store the cluster labels, with default label -1 (indicating noise or unclustered points)
            // Assign cluster labels to the corresponding points
            for (size_t cluster_idx = 0; cluster_idx < clusters.size(); ++cluster_idx) {
                for (auto point_idx : clusters[cluster_idx]) {
                    labels[point_idx] = static_cast<int>(cluster_idx)+1;
                }
            }

            // Convert the vector of cluster labels back to a Torch tensor
            auto clustering = torch::from_blob(labels.data(), {static_cast<long>(labels.size())}, torch::kInt32).clone();

            // Find unique cluster labels and the corresponding indices that map the original points to these unique labels
            torch::Tensor unique_tensor;
            torch::Tensor inverse_indices;
            std::tie(unique_tensor, inverse_indices) = at::_unique(clustering, true, true);

            /////////////////////////////////
            ////////// OUTPUT STEP //////////
            /////////////////////////////////

            // Get the total number of unique tracks based on the unique_tensor size
            int64_t number_of_tracks = unique_tensor.numel(); 

            // Loop through each unique track ID
            for (int i = 0; i < number_of_tracks; ++i) {

                // Retrieve the current track ID
                auto id_of_track = unique_tensor.index({i});
                
                // Create a new track in the output collection and set its type to the current track ID
                auto output_track = output_tracks->create();
                output_track.setType(id_of_track.item<int>());

                // std::cout << typeid(output_track).name() << std::endl;

                // Create a mask to select all hits belonging to the current track
                torch::Tensor mask = (clustering == id_of_track);
                
                // Find the indices of the hits that belong to the current track
                torch::Tensor indices = torch::nonzero(mask);
                int64_t number_of_hits = indices.numel();

                // Loop through each hit index for the current track
                for (int j = 0; j < number_of_hits; ++j) {

                    // Get the current hit index
                    auto index_id = indices.index({j});
                    
                    // Check which detector the hit belongs to (VTXD, VTXIB, VTOB, CDC)
                    torch::Tensor mask_VTXD = (ListHitType_VTXD_tensor == index_id);
                    torch::Tensor mask_VTXIB = (ListHitType_VTXIB_tensor == index_id);
                    torch::Tensor mask_VTOB = (ListHitType_VTXOB_tensor == index_id);
                    torch::Tensor mask_CDC = (ListHitType_CDC_tensor == index_id);

                    // If the hit belongs to the VTXD detector
                    if ((torch::sum(mask_VTXD) > 0).item<bool>()) {

                        auto hit = inputHits_VTXD.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);

                    } 

                    // If the hit belongs to the VTXIB detector
                    else if ((torch::sum(mask_VTXIB) > 0).item<bool>()) {
                        index_id = index_id - it_0;

                        auto hit = inputHits_VTXIB.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);

                    } 
                    // If the hit belongs to the VTOB detector
                    else if ((torch::sum(mask_VTOB) > 0).item<bool>()) {
                        index_id = index_id - (it_1 + it_0);

                        auto hit = inputHits_VTXOB.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);

                    }
                    // If the hit belongs to the CDC detector
                    else if ((torch::sum(mask_CDC) > 0).item<bool>()) {
                        index_id = index_id - (it_1 + it_2 + it_0);

                        auto hit = inputHits_CDC.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);
                    }
                }
            }

            inverse_indices.reset();
            unique_tensor.reset();
            clustering.reset();
            
            input_tensors.clear();
            output_model_tensors.clear();
            
            
        }

        ListHitType_VTXIB_tensor.reset();
        ListHitType_VTXOB_tensor.reset();
        ListHitType_VTXD_tensor.reset();
        ListHitType_CDC_tensor.reset();

        std::vector<float>().swap(ListHitType_VTXD);
        std::vector<float>().swap(ListHitType_VTXIB);
        std::vector<float>().swap(ListHitType_VTXOB);
        std::vector<float>().swap(ListHitType_CDC);

        // Return the output collections as a tuple
        return std::make_tuple(std::move(*output_tracks));

    } 

    private:
        
        std::vector<char> model_data;

        /// Pointer to the ONNX environment.
        /// This object manages the global state of the ONNX runtime, such as logging and threading.
        std::unique_ptr<Ort::Env> fEnv;

        /// Pointer to the ONNX inference session.
        /// This session is used to execute the model for inference.
        std::unique_ptr<Ort::Session> fSession;

        /// ONNX session options.
        /// These settings control the behavior of the inference session, such as optimization level, 
        /// execution providers, and other configuration parameters.
        Ort::SessionOptions fSessionOptions;

        /// ONNX memory info.
        /// This object provides information about memory allocation and is used during the creation of 
        /// ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
        const OrtMemoryInfo* fInfo;
        struct MemoryInfo;

        /// Stores the input and output names for the ONNX model.
        /// These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
        std::vector<const char*> fInames;
        std::vector<const char*> fOnames;

        /// Property to specify the path to the ONNX model file.
        /// This is a configurable property that defines the location of the ONNX model file on the filesystem.
        Gaudi::Property<std::string> modelPath{this, "modelPath", "/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx", "modelPath"};

        /// Property to configure the step size for the DBSCAN clustering algorithm.
        /// This parameter controls the maximum distance between points in a cluster.
        Gaudi::Property<double> step_size{this, "step_size", 0.5, "step_size"};

        /// Property to configure the minimum number of points required to form a cluster in the DBSCAN algorithm.
        /// This parameter defines the density threshold for identifying clusters.
        Gaudi::Property<int> min_points{this, "min_points", 10, "min_points"};
        

};

DECLARE_COMPONENT(GGTF_tracking_dbscan)