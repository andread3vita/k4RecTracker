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
// #include <memory> 
// #include <numeric>
// #include <queue>
// #include <sstream>
// #include <string>
// #include <typeinfo>
// #include <vector>
// #include <fstream>  // Per std::ifstream
// #include <vector>   // Per std::vector
// #include <iterator> // Per std::istreambuf_iterator

// #include <ATen/ATen.h>
// #include <torch/torch.h>
// #include "onnxruntime_cxx_api.h"
// #include "onnxruntime_run_options_config_keys.h"

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
// using FloatColl = podio::UserDataCollection<float>;

// using ParticleColl = edm4hep::MCParticleCollection;
// using DCTrackerHitColl = extension::DriftChamberDigiCollection;
// using DCTrackerHitColl_sim = edm4hep::SimTrackerHitCollection;
// using VertexColl_sim = edm4hep::SimTrackerHitCollection;
// using SimHits = edm4hep::SimTrackerHitCollection;

// using HitsColl = extension::TrackerHit3DCollection;
// using TrackColl = extension::TrackCollection;

// #include <iostream>
// #include <typeinfo>
// #include <cxxabi.h>
// #include <memory>


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
//  *    - Hits collection : HitsColl (VTXD , VTXIB , VTXOB , CDC)
//  *    - Track collection : TrackColl
//  *
//  *
//  *
//  *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
//  *  @date   2024-08
//  *
//  */

// struct GGTF_tracking_dbscan_check final : 
//         k4FWCore::MultiTransformer<std::tuple<FloatColl,DoubleColl,IntColl>(
//                                                                     const DCTrackerHitColl&, 
//                                                                     const HitsColl&,
//                                                                     const HitsColl&,
//                                                                     const HitsColl&,

//                                                                     const SimHits&,
//                                                                     const SimHits&,
//                                                                     const SimHits&,
//                                                                     const SimHits&)> 
            
                                                                                            
// {
//     GGTF_tracking_dbscan_check(const std::string& name, ISvcLocator* svcLoc) : 
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
//                 KeyValues("inputModel", {"inputModel"}),
//                 KeyValues("outputModel", {"outputModel"}),
//                 KeyValues("outputClustering", {"outputClustering"})           
            
//             }
//             ) {}
    

//     // StatusCode initialize() {

//     //     // Leggi il file .onnx in modalità binaria
//     //     std::ifstream model_file(modelPath.value().c_str(), std::ios::binary | std::ios::ate);
//     //     if (!model_file) {
//     //         std::cerr << "Failed to open the ONNX model file: " << modelPath << std::endl;
//     //     }

//     //     // Ottieni la dimensione del file
//     //     std::streamsize model_size = model_file.tellg();
//     //     model_file.seekg(0, std::ios::beg);

//     //     // Crea un buffer per contenere il file modello
//     //     std::vector<char> model_data_temp(model_size);
//     //     if (!model_file.read(model_data_temp.data(), model_size)) {
//     //         std::cerr << "Failed to read the ONNX model file: " << modelPath << std::endl;
//     //     }

//     //     model_data = model_data_temp;

//     //     // Initialize the ONNX memory info object for CPU memory allocation.
//     //     // This specifies that the memory will be allocated using the Arena Allocator on the CPU.
//     //     fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

//     //     // Create and initialize the ONNX environment with a logging level set to WARNING.
//     //     // This environment handles logging and runtime configuration for the ONNX session.
//     //     auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
//     //     fEnv          = std::move(envLocal);

//     //     // Set the session options to configure the ONNX inference session.
//     //     // Set the number of threads used for intra-op parallelism to 1 (single-threaded execution).
//     //     fSessionOptions.SetIntraOpNumThreads(1);

//     //     // Disable all graph optimizations to keep the model execution as close to the original as possible.
//     //     fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

//     //     fSessionOptions.DisableMemPattern();

//     //     // Create an ONNX inference session using the configured environment and session options.
//     //     // The session is used to load the model specified by the `modelPath`.
//     //     auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);
//     //     fSession          = std::move(sessionLocal);

//     //     // Create an ONNX allocator with default options to manage memory allocations during runtime.
//     //     Ort::AllocatorWithDefaultOptions allocator;

//     //     // Get the name of the first input node (index i) from the ONNX model and store it.
//     //     // This retrieves the name from the model and releases the memory after storing it in fInames.
//     //     std::size_t i = 0;
//     //     const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();

//     //     // Get the name of the first output node (index i) from the ONNX model and store it.
//     //     // This retrieves the name from the model and releases the memory after storing it in fOnames.
//     //     const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

//     //     // Store the retrieved input and output names in the respective vectors.
//     //     fInames.push_back(input_name);
//     //     fOnames.push_back(output_names);


//     //     return StatusCode::SUCCESS;

//     // }


//     std::tuple<FloatColl,DoubleColl,IntColl> operator()(
//                                                 const DCTrackerHitColl& inputHits_CDC, 
//                                                 const HitsColl& inputHits_VTXIB,
//                                                 const HitsColl& inputHits_VTXD,
//                                                 const HitsColl& inputHits_VTXOB,
//                                                 const SimHits& inputHits_CDC_sim,
//                                                 const SimHits& inputHits_VTXIB_sim,
//                                                 const SimHits& inputHits_VTXD_sim,
//                                                 const SimHits& inputHits_VTXOB_sim) const override 
//     {

//         // torch::NoGradGuard no_grad;

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

//         fSessionOptions.DisableMemPattern();

//         /// ONNX memory info.
//         /// This object provides information about memory allocation and is used during the creation of 
//         /// ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
//         const OrtMemoryInfo* fInfo;

//         /// Stores the input and output names for the ONNX model.
//         /// These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
//         std::vector<const char*> fInames;
//         std::vector<const char*> fOnames;
        
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
//         const auto input_name = fSession->GetInputNameAllocated(0, allocator).release();

//         // Get the name of the first output node (index i) from the ONNX model and store it.
//         // This retrieves the name from the model and releases the memory after storing it in fOnames.
//         const auto output_names = fSession->GetOutputNameAllocated(0, allocator).release();

//         // Store the retrieved input and output names in the respective vectors.
//         fInames.push_back(input_name);
//         fOnames.push_back(output_names);

//         //////////////////////
//         ///////////////////////
//         //////////////////////
//         std::ifstream file("/afs/cern.ch/user/a/adevita/public/check_file/matrix_python.csv");
//         std::vector<std::vector<double>> matrix;
//         std::string line;
//         while (std::getline(file, line)) {
//             std::istringstream iss(line);
//             std::vector<double> row;
//             double value;
//             while (iss >> value) {
//                 row.push_back(value);
//                 if (iss.peek() == ',') {
//                     iss.ignore();
//                 }
//             }
//             matrix.push_back(row);
//         }
        
//         // Obtain the shape 
//         int64_t num_rows = matrix.size();
//         int64_t num_cols = matrix[0].size();

//         std::cout << num_rows << " " << num_cols << std::endl;

//         std::vector<int64_t> t_shape = {num_rows, num_cols};
//         // Convert vector of vectors to a single vector (because the CreateTensor only takes vectors)
//         std::vector<float> flatVector;
//         for (const auto& row : matrix) {
//             flatVector.insert(flatVector.end(), row.begin(), row.end());
//         }

//         Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

//         // Create an Ort::Env object
//         Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
//         std::cout << "Environment created" << std::endl;
//         // Calculate total size of the tensor
//         size_t t_size = 1;
//         for (int64_t dim : t_shape) {
//             t_size *= dim;
//         }
//         // Create a tensor from the data
//         std::vector<Ort::Value> in_tensors;
//         in_tensors.emplace_back( Ort::Value::CreateTensor<float>(memory_info, flatVector.data(), t_size, t_shape.data(), t_shape.size()));
        
//         // Run the model
//         Ort::SessionOptions session_options;
//         session_options.SetIntraOpNumThreads(1);
//         session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
//         printf("Using Onnxruntime C++ API\n");
//         Ort::Session session(env, modelPath.value().c_str(), session_options);	
//         printf("Starting to run inference\n");
//         Ort::AllocatorWithDefaultOptions allo;

//         // get inputs and outputs
//         std::vector<std::string> input_names;
//         std::vector<std::int64_t> input_shapes;
//         for (std::size_t i = 0; i < session.GetInputCount(); i++) {
//             input_names.emplace_back(session.GetInputNameAllocated(i, allo).get());
//             input_shapes = session.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
//         }
//         for (auto& s : input_shapes) {
//             if (s < 0) {
//             s = 1;
//             }
//         }

//         // Extract the data from the tensor
//         float* tensor_data = in_tensors.front().GetTensorMutableData<float>();

//         // Copy tensor data into a std::vector<float>
//         std::vector<float> tensor_vector(tensor_data, tensor_data + t_size);
//         std::ofstream fil_IN("/afs/cern.ch/user/a/adevita/public/check_file/input.csv");   
//         // Writing the tensor to the CSV file
//         for (int64_t i = 0; i < num_rows; ++i) {
//             for (int64_t j = 0; j < num_cols; ++j) {
//             fil_IN << tensor_vector[i * num_cols + j];
//             if (j < num_cols - 1) {
//                 fil_IN << ",";  // Add comma between values
//             }
//         }
//         fil_IN << "\n";  // New line after each row
//         }

//         fil_IN.close();

//         std::vector<std::string> output_name_collection;
//         for (std::size_t i = 0; i < session.GetOutputCount(); i++) {
//             output_name_collection.emplace_back(session.GetOutputNameAllocated(i, allo).get());
//             auto output_shapes = session.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
//         }
        
//         std::vector<const char*> input_names_char(input_names.size(), nullptr);
//         std::transform(std::begin(input_names), std::end(input_names), std::begin(input_names_char),
//                     [&](const std::string& str) { return str.c_str(); });
                    
//         std::vector<const char*> output_names_char(output_name_collection.size(), nullptr);
//         std::transform(std::begin(output_name_collection), std::end(output_name_collection), std::begin(output_names_char),
//                     [&](const std::string& str) { return str.c_str(); });

//         auto output_tensors = session.Run(Ort::RunOptions{nullptr}, input_names_char.data(), in_tensors.data(),
//                                         input_names_char.size(), output_names_char.data(), output_names_char.size());
//         std::cout << "Done!" << std::endl;

//         // Get the outputs and covert to vector 
//         float* floatar = output_tensors.front().GetTensorMutableData<float>();
//         std::vector<float> output_vector(floatar, floatar + num_rows*4);

//         std::ofstream fil("/afs/cern.ch/user/a/adevita/public/check_file/out.csv");   
         
//         // Scrittura del vettore nel file come righe con 4 colonne    
//         for (int i = 0; i < num_rows; ++i) {        
//             for (int j = 0; j < 4; ++j) {            
//                 fil << output_vector[i * 4 + j];            
                
//                 if (j < 3) {                
//                     fil << ",";         
//                 }   
//                 else{
//                     fil << "\n";     
//                 }     
//             }        
                
//         }

//         fil.close();

//         // ////////////////////////////////////////
//         // ////////// DATA PREPROCESSING //////////
//         // ////////////////////////////////////////
        
//         // Although there is a link between digi hits and MC hits, we cannot use it because the output tracks 
//         // do not contain the original hits, but rather copies. This step is necessary because it is not possible 
//         // to use hits from both the drift chamber and the vertex in the same collection. Therefore, for each hit, 
//         // a MutableTrackerHit3D is created, and its properties are defined based on the original hit.
//         //
//         // These for loops are used to store the index of the MC particle for each hit. For now, this information 
//         // will be saved in the EDep of the hit stored in each track.

//         std::cout << "Input Hit collection size VTXD: " << inputHits_VTXD.size() << std::endl;
//         std::cout << "Input Hit collection size VTXIB: " << inputHits_VTXIB.size() << std::endl;
//         std::cout << "Input Hit collection size VTXOB: " << inputHits_VTXOB.size() << std::endl;
//         std::cout << "Input Hit collection size CDC: " << inputHits_CDC.size() << std::endl;
//         std::cout << "______________________________________________________" << std::endl;
//         std::cout << "Tot hits: " << inputHits_CDC.size() + inputHits_VTXOB.size() + inputHits_VTXIB.size() + inputHits_VTXD.size() << std::endl;
//         std::cout << "  " << std::endl;

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
//             ListGlobalInputs.push_back(1.0); 
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
//             ListGlobalInputs.push_back(1.0);
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
//             ListGlobalInputs.push_back(1.0); 
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
//             ListGlobalInputs.push_back(0.0); 
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

//         // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
//         DoubleColl outputVec_model;
//         std::vector<int> labels(it, -1);
//         if (it > 0 && it < 20000)
//         {
//             // Calculate the total size of the input tensor, based on the number of hits (it) and the 
//             // number of features per hit (7: x, y, z, and four placeholders).
//             size_t total_size = it * 7;

//             // Define the shape of the input tensor. The tensor will have `it` rows (one for each hit)
//             // and 7 columns (representing the features of each hit).
//             std::vector<int64_t> tensor_shape = {it, 7};

//             // // Create a vector to store the input tensors that will be fed into the ONNX model.
//             std::vector<Ort::Value> input_tensors;

//             // // Create an ONNX tensor from the input data (ListGlobalInputs) and add it to the input_tensors vector.
//             // // The tensor is created using the memory information (fInfo), the data pointer, the total size,
//             // // and the tensor shape defined earlier.
//             // input_tensors.emplace_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));
            
//             input_tensors.push_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));

//             // // Run the ONNX inference session with the provided input tensor.
//             // // The model will process the input data and produce output tensors.
//             // // The input tensor names (fInames), input tensors (input_tensors), output tensor names (fOnames), 
//             // // and the sizes of these vectors are passed to the Run method.
//             auto output_model_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());
            
//             // Extract the raw output data from the first output tensor returned by the model.
//             // This data is stored as a float array.
//             double* floatarr = output_model_tensors.front().GetTensorMutableData<double>();

//             // Convert the raw output data into a std::vector of floats for easier handling.
//             // The size of the vector is determined by the number of hits (it) and the number of 
//             // output features per hit (3 coordinates and charge).
//             std::vector<double> output_model_vector(floatarr, floatarr + it * 4);

//             /////////////////////////////////////
//             ////////// CLUSTERING STEP //////////
//             /////////////////////////////////////

//             // Convert the output vector from the model to a Torch tensor, specifying the shape and data type
//             torch::Tensor output_model_tensor = torch::from_blob(output_model_vector.data(), {it, 4}, torch::kFloat32).clone();    

//             // Initialize a vector to store the 3D points with an associated weight
//             std::vector<point3w> points;
//             for (int i = 0; i < output_model_tensor.size(0); ++i) {
//                 auto x = output_model_tensor[i][0].item<double>();
//                 auto y = output_model_tensor[i][1].item<double>();
//                 auto z = output_model_tensor[i][2].item<double>();
//                 auto weight = output_model_tensor[i][3].item<double>();
                
//                 // Store the point in the points vector
//                 points.push_back({x, y, z, weight});
//             }

//             for (int i = 0; i < it; ++i) {

//                 auto x = output_model_tensor[i][0].item<double>();
//                 auto y = output_model_tensor[i][1].item<double>();
//                 auto z = output_model_tensor[i][2].item<double>();
//                 auto weight = output_model_tensor[i][3].item<double>();
                
//                 outputVec_model.push_back(x);
//                 outputVec_model.push_back(y);
//                 outputVec_model.push_back(z);
//                 outputVec_model.push_back(weight);
//             }


//             // Apply the DBSCAN clustering algorithm to the points with the given step size and minimum points per cluster
//             // NOTE: dbscan doesn't take into account the weights, so the fourth entry of each point will be ignored.
//             auto clusters = dbscan(points, step_size, min_points);

//             // Initialize a vector to store the cluster labels, with default label -1 (indicating noise or unclustered points)
//             // Assign cluster labels to the corresponding points
//             for (size_t cluster_idx = 0; cluster_idx < clusters.size(); ++cluster_idx) {
//                 for (auto point_idx : clusters[cluster_idx]) {
//                     labels[point_idx] = static_cast<int>(cluster_idx)+1;
//                 }
//             }

//             // Convert the vector of cluster labels back to a Torch tensor
//             auto clustering = torch::from_blob(labels.data(), {static_cast<long>(labels.size())}, torch::kInt32).clone();

//             clustering.reset();
            
//             input_tensors.clear();
//             output_model_tensors.clear();
            
            
//         }

//         ListHitType_VTXIB_tensor.reset();
//         ListHitType_VTXOB_tensor.reset();
//         ListHitType_VTXD_tensor.reset();
//         ListHitType_CDC_tensor.reset();

//         std::vector<float>().swap(ListHitType_VTXD);
//         std::vector<float>().swap(ListHitType_VTXIB);
//         std::vector<float>().swap(ListHitType_VTXOB);
//         std::vector<float>().swap(ListHitType_CDC);

//         // Return the output collections as a tuple
//         return std::make_tuple(std::move(ListGlobalInputs),std::move(outputVec_model), std::move(labels));


//     } 

//     private:
        
//         // std::vector<char> model_data;

//         // /// Pointer to the ONNX environment.
//         // /// This object manages the global state of the ONNX runtime, such as logging and threading.
//         // std::unique_ptr<Ort::Env> fEnv;

//         // /// Pointer to the ONNX inference session.
//         // /// This session is used to execute the model for inference.
//         // std::unique_ptr<Ort::Session> fSession;

//         // /// ONNX session options.
//         // /// These settings control the behavior of the inference session, such as optimization level, 
//         // /// execution providers, and other configuration parameters.
//         // Ort::SessionOptions fSessionOptions;

//         // /// ONNX memory info.
//         // /// This object provides information about memory allocation and is used during the creation of 
//         // /// ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
//         // const OrtMemoryInfo* fInfo;
//         // struct MemoryInfo;

//         // /// Stores the input and output names for the ONNX model.
//         // /// These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
//         // std::vector<const char*> fInames;
//         // std::vector<const char*> fOnames;

//         /// Property to specify the path to the ONNX model file.
//         /// This is a configurable property that defines the location of the ONNX model file on the filesystem.
//         Gaudi::Property<std::string> modelPath{this, "modelPath", "/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx", "modelPath"};

//         /// Property to configure the step size for the DBSCAN clustering algorithm.
//         /// This parameter controls the maximum distance between points in a cluster.
//         Gaudi::Property<double> step_size{this, "step_size", 0.5, "step_size"};

//         /// Property to configure the minimum number of points required to form a cluster in the DBSCAN algorithm.
//         /// This parameter defines the density threshold for identifying clusters.
//         Gaudi::Property<int> min_points{this, "min_points", 10, "min_points"};



// };

// DECLARE_COMPONENT(GGTF_tracking_dbscan_check)