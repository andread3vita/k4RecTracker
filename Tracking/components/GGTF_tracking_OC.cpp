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

// Standard Library
#include <algorithm>
#include <cmath>
#include <cstdlib>         // For getenv
#include <filesystem>      // For std::filesystem::path
#include <fstream>         // For std::ifstream
#include <iostream>
#include <iterator>        // For std::istreambuf_iterator
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// ONNX & Torch
#include <ATen/ATen.h>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"

// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TGeoMatrix.h"
#include "TRandom3.h"
#include "TVector3.h"

// === Gaudi Framework ===
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// === k4FWCore / k4Interface ===
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// === EDM4HEP & PODIO ===
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "podio/UserDataCollection.h"

// === EDM4HEP Extensions ===
#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "extension/TrackCollection.h"
#include "extension/TrackerHit.h"

// === DD4hep ===
#include "DD4hep/Detector.h"       
#include "DDRec/DCH_info.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

// === Project-specific ===
#include "utils.hpp"

/** @struct GGTF_tracking_OC_IDEAv3
 *
 *  Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the GGTF_tracking. 
 *  The first step takes the raw hits and it returns a collection of 4-dimensional points inside an embedding space.
 *  Eeach 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described intuitively by a potential, 
 *  which attracts hits belonging to the same cluster and drives away those that do not.
 *  This collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the same track.
 *
 *  input: 
 *    - digitalized hits from DC (global coordinates) : extension::SenseWireHitCollection
 *    - digitalized hits from vertex (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *      - digitalized hits from silicon wrapper (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *
 *  output:
 *    - Track collection : extension::TrackCollection
 *
 *
 *
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2025-06
 *
 */

struct GGTF_tracking_OC_IDEAv3 final : 
        k4FWCore::MultiTransformer< std::tuple<extension::TrackCollection>(   
                                                                    const extension::SenseWireHitCollection&, 
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&)>                                                                                            
{
    GGTF_tracking_OC_IDEAv3(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputHits_CDC", {"inputHits_CDC"}),
                KeyValues("inputHits_VTXB", {"inputHits_VTXB"}),
                KeyValues("inputHits_VTXD", {"inputHits_VTXD"}),
                KeyValues("inputHits_SWB", {"inputHits_SWB"}),
                KeyValues("inputHits_SWD", {"inputHits_SWD"})

            },
            {   
               
               KeyValues("outputTracks", {"outputTracks"})      
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
    StatusCode initialize() {

        ///////////////////////////////
        ///// ONNX Initialization /////
        ///////////////////////////////

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
        // fSessionOptions.DisableMemPattern();

        // Create an ONNX inference session using the configured environment and session options.
        // The session is used to load the model specified by the `modelPath`.
        auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);
        fSession          = std::move(sessionLocal);

        // Create an ONNX allocator with default options to manage memory allocations during runtime.
        Ort::AllocatorWithDefaultOptions allocator;

        // Get the name of the first input node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fInames.
        // Get the name of the first output node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fOnames.
        std::size_t i = 0;
        const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();
        const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

        // Store the retrieved input and output names in the respective vectors.
        fInames.push_back(input_name);
        fOnames.push_back(output_names);

        //////////////////////////////
        ///// Drift Chamber Info /////
        //////////////////////////////

        // Retrieve the drift chamber information and decoder
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}



        return StatusCode::SUCCESS;

   }

    
    std::tuple<extension::TrackCollection> operator()(      const extension::SenseWireHitCollection& inputHits_CDC, 
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_VTXB,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_VTXD,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_SWB,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_SWD) const override 
    {

        ////////////////////////////////////////
        ////////// DATA PREPROCESSING //////////
        ////////////////////////////////////////


        // Vector to store the global input values for all hits.
        // This will contain position and other hit-specific data to be used as input for the model.
        std::vector<float> ListGlobalInputs; 
        int it = 0;

        /// Processing hits from the VTXD (Vertex Disk).
        std::vector<float> ListHitType_VTXD;
        int it_0 = 0; 
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



        /// Processing hits from the VTXB (Vertex Barrel).
        std::vector<float> ListHitType_VTXB; 
        int it_1 = 0; 
        for (const auto input_hit : inputHits_VTXB) {

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
            ListHitType_VTXB.push_back(it);
            it += 1; 
            it_1 += 1;                      
        }
        // Convert ListHitType_VTXD to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_VTXB_tensor = torch::from_blob(ListHitType_VTXB.data(), {it_1}, torch::kFloat32);

        /// Processing hits from the VTXB (Vertex Barrel).
        std::vector<float> ListHitType_SWB; 
        int it_2 = 0; 
        for (const auto input_hit : inputHits_SWB) {

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
            ListHitType_SWB.push_back(it);
            it += 1; 
            it_2 += 1;                      
        }
        // Convert ListHitType_VTXD to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_SWB_tensor = torch::from_blob(ListHitType_SWB.data(), {it_2}, torch::kFloat32);

        /// Processing hits from the VTXB (Vertex Barrel).
        std::vector<float> ListHitType_SWD; 
        int it_3 = 0; 
        for (const auto input_hit : inputHits_SWD) {

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
            ListHitType_VTXB.push_back(it);
            it += 1; 
            it_3 += 1;                      
        }
        // Convert ListHitType_VTXD to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_SWD_tensor = torch::from_blob(ListHitType_SWD.data(), {it_3}, torch::kFloat32);
    
    

        std::vector<float> ListHitType_CDC;
        int it_4 = 0;
        for (const auto input_hit : inputHits_CDC) {
            
            // position along the wire, wire direction and drift distance
            edm4hep::Vector3d wirePos = input_hit.getPosition();   
            TVector3 wire_pos(static_cast<float>(wirePos.x),static_cast<float>(wirePos.y),static_cast<float>(wirePos.z));

            double distanceToWire = input_hit.getDistanceToWire(); 
            double wireAzimuthalAngle = input_hit.getWireAzimuthalAngle();  
            double wireStereoAngle = input_hit.getWireStereoAngle();  
            
            // Direction of the wire: z' axis
            double d_x = std::sin(wireStereoAngle) * std::sin(wireAzimuthalAngle);
            double d_y = -std::sin(wireStereoAngle) * std::cos(wireAzimuthalAngle);
            double d_z = std::cos(wireStereoAngle);
            TVector3 z_prime(d_x, d_y, d_z);
            z_prime = z_prime.Unit();

            // x' axis, orthogonal to z'
            TVector3 x_prime(1.0, 0.0, -d_x / d_z);
            x_prime = x_prime.Unit();

            // y' axis = z' × x'
            TVector3 y_prime = z_prime.Cross(x_prime);
            y_prime = y_prime.Unit();

            // Define the local left/right positions
            TVector3 left_local_pos(-distanceToWire, 0.0, 0.0);
            TVector3 right_local_pos(distanceToWire, 0.0, 0.0);

            // Transform to global

            TVector3 left_global_pos = x_prime * left_local_pos.X() + y_prime * left_local_pos.Y() + z_prime * left_local_pos.Z() + wire_pos;
            TVector3 right_global_pos = x_prime * right_local_pos.X() + y_prime * right_local_pos.Y() + z_prime * right_local_pos.Z() + wire_pos;
    
            // Add the 3D position of the left hit to the global input list.
            ListGlobalInputs.push_back(left_global_pos.X());
            ListGlobalInputs.push_back(left_global_pos.Y());
            ListGlobalInputs.push_back(left_global_pos.Z());
            
            // Add the difference between the right and left hit positions to the global input list.
            ListGlobalInputs.push_back(0.0); 
            ListGlobalInputs.push_back(right_global_pos.X()-left_global_pos.X());
            ListGlobalInputs.push_back(right_global_pos.y()-left_global_pos.Y());
            ListGlobalInputs.push_back(right_global_pos.Z()-left_global_pos.Z());
            
            // Store the current index in ListHitType_CDC and increment the global iterator.
            ListHitType_CDC.push_back(it);
            it += 1; 
            it_4 += 1;                      
        }
        // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_4}, torch::kFloat32);


        /////////////////////////////
        ////////// ML STEP //////////
        /////////////////////////////

        // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
        extension::TrackCollection output_tracks;
        if (it > 0 && it < 20000)
        {
            // Calculate the total size of the input tensor, based on the number of hits (it) and the 
            // number of features per hit (7: x, y, z, and four placeholders).
            size_t total_size = it * 7;
            std::vector<int64_t> tensor_shape = {it, 7};

            // Create a vector to store the input tensors that will be fed into the ONNX model.
            std::vector<Ort::Value> input_tensors;
            input_tensors.emplace_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));
            
            // Run the ONNX inference session with the provided input tensor.
            auto output_model_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());
            float* floatarr = output_model_tensors.front().GetTensorMutableData<float>();
            std::vector<float> output_model_vector(floatarr, floatarr + it * 4);

            /////////////////////////////////////
            ////////// CLUSTERING STEP //////////
            /////////////////////////////////////

            auto clustering = get_clustering(output_model_vector, it,0.6,0.3);
            torch::Tensor unique_tensor;
            torch::Tensor inverse_indices;
            std::tie(unique_tensor, inverse_indices) = at::_unique(clustering, true, true);
          
            /////////////////////////////////
            ////////// OUTPUT STEP //////////
            /////////////////////////////////

            // Get the total number of unique tracks based on the unique_tensor size
            int64_t number_of_tracks = unique_tensor.numel(); 
            
            bool has_zero = (unique_tensor == 0).any().item<bool>();
            if (!has_zero)
            {
                auto output_track = output_tracks.create();
                output_track.setType(0);
            }

            // Loop through each unique track ID
            for (int i = 0; i < number_of_tracks; ++i) {

                // Retrieve the current track ID
                auto id_of_track = unique_tensor.index({i});


                // Create a new track in the output collection and set its type to the current track ID
                auto output_track = output_tracks.create();
                output_track.setType(id_of_track.item<int>());

                // Create a mask to select all hits belonging to the current track
                torch::Tensor mask = (clustering == id_of_track);
                
                // Find the indices of the hits that belong to the current track
                torch::Tensor indices = torch::nonzero(mask);
                int64_t number_of_hits = indices.numel();

                // Loop through each hit index for the current track
                for (int j = 0; j < number_of_hits; ++j) {

                    // Get the current hit index
                    auto index_id = indices.index({j});
                    
                    // Check which detector the hit belongs to (VTXD, VTXB, CDC)
                    torch::Tensor mask_VTXD = (ListHitType_VTXD_tensor == index_id);
                    torch::Tensor mask_VTXB = (ListHitType_VTXB_tensor == index_id);
                    torch::Tensor mask_CDC = (ListHitType_CDC_tensor == index_id);

                    // If the hit belongs to the VTXD detector
                    if ((torch::sum(mask_VTXD) > 0).item<bool>()) {

                        auto hit = inputHits_VTXD.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);

                    } 
                    // If the hit belongs to the VTOB detector
                    else if ((torch::sum(mask_VTXB) > 0).item<bool>()) {
                        index_id = index_id - (it_0);

                        auto hit = inputHits_VTXB.at(index_id.item<int>());
                        output_track.addToTrackerHits(hit);

                    } 
                    // If the hit belongs to the CDC detector
                    else if ((torch::sum(mask_CDC) > 0).item<bool>()) {
                        index_id = index_id - (it_1 + it_0);

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

        ListHitType_VTXB_tensor.reset();
        ListHitType_VTXD_tensor.reset();
        ListHitType_CDC_tensor.reset();

        std::vector<float>().swap(ListHitType_VTXD);
        std::vector<float>().swap(ListHitType_VTXB);
        std::vector<float>().swap(ListHitType_CDC);

        // Return the output collections as a tuple
        return std::make_tuple(std::move(output_tracks));


    } 

    private:

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


        //------------------------------------------------------------------
        //          machinery for geometry

        /// Geometry service name
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

        /// Detector name
        Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

        /// Pointer to the geometry service
        SmartIF<IGeoSvc> m_geoSvc;

        dd4hep::rec::DCH_info* dch_info;
        dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;
        

};

DECLARE_COMPONENT(GGTF_tracking_OC_IDEAv3)