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
#include <typeinfo>

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

#include "edm4hep/SimTrackerHitCollection.h"
using SimHits = edm4hep::SimTrackerHitCollection;

#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
using VertexHitsColl = edm4hep::TrackerHit3DCollection;
using VTX_links = edm4hep::TrackerHitSimTrackerHitLinkCollection;

#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
using DCHitsColl = extension::SenseWireHitCollection;
using DC_links = extension::SenseWireHitSimTrackerHitLinkCollection;

#include "extension/TrackCollection.h"
using TrackColl = extension::TrackCollection;

#include "extension/TrackerHit.h"
using TrackHit = extension::TrackerHit;

#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>  // For getenv
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"
#include <filesystem>  // For std::filesystem::path

#include <Track.h>

/** @struct GGTF_tracking_dbscan_IDEAv3
 *
 *  Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the GGTF_tracking. 
 *  The first step takes the raw hits and it returns a collection of 4-dimensional points inside an embedding space.
 *  Eeach 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described intuitively by a potential, 
 *  which attracts hits belonging to the same cluster and drives away those that do not.
 *  This collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the same track.
 *
 *  input: 
 *    - digitalized hits from DC (global coordinates) : DCHitsColl 
 *    - digitalized hits from vertex (global coordinates) : VertexHitsColl
 *
 *  output:
 *    - Track collection : TrackColl
 *
 *
 *
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2025-02
 *
 */

struct GGTF_fitter final : 
        k4FWCore::MultiTransformer< std::tuple<TrackColl>( 
                                                                    
                                                                    const DCHitsColl&, 
                                                                    const VertexHitsColl&,
                                                                    const VertexHitsColl&)> 
            
                                                                                            
{
    GGTF_fitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputHits_CDC", {"inputHits_CDC"}),
                KeyValues("inputHits_VTXB", {"inputHits_VTXB"}),
                KeyValues("inputHits_VTXD", {"inputHits_VTXD"})
            },
            {   
                KeyValues("outputTracks", {"outputTracks"})      
            
            }) {}
    
    StatusCode initialize() {


        return StatusCode::SUCCESS;

    }

    
    std::tuple<TrackColl> operator()(   const DCHitsColl& inputHits_CDC, 
                                                const VertexHitsColl& inputHits_VTXB,
                                                const VertexHitsColl& inputHits_VTXD) const override 
    {

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
        std::cout << "Input Hit collection size VTXB: " << inputHits_VTXB.size() << std::endl;
        std::cout << "Input Hit collection size CDC: " << inputHits_CDC.size() << std::endl;
        std::cout << "______________________________________________________" << std::endl;
        std::cout << "Tot hits: " << inputHits_CDC.size() + inputHits_VTXB.size() + inputHits_VTXD.size() << std::endl;
        std::cout << "  " << std::endl;

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
        std::vector<float> ListHitType_CDC;
        int it_2 = 0;
        for (const auto input_hit : inputHits_CDC) {

            
            edm4hep::Vector3d wirePos = input_hit.getPosition();    // position along the wire
            std::vector<float> wire_pos = {static_cast<float>(wirePos.x),static_cast<float>(wirePos.y),static_cast<float>(wirePos.z)};

            double distanceToWire = input_hit.getDistanceToWire();   // drift distance
            double wire_azimuthal_angle = input_hit.getWireAzimuthalAngle();  
            double wire_stereo_angle = input_hit.getWireStereoAngle();  

            // Wire direction
            float d_x = std::sin(wire_stereo_angle) * std::cos(wire_azimuthal_angle);
            float d_y = std::sin(wire_stereo_angle) * std::sin(wire_azimuthal_angle);
            float d_z = std::cos(wire_stereo_angle);

            // z_prime
            std::vector<float> z_prime = {d_x, d_y, d_z};
            
            // x_prime
            std::vector<float> x_prime = {d_z, 0.0f, -d_z};
            float norm_x_prime = std::sqrt(x_prime[0] * x_prime[0] + x_prime[1] * x_prime[1] + x_prime[2] * x_prime[2]);
            x_prime[0] /= norm_x_prime;
            x_prime[1] /= norm_x_prime;
            x_prime[2] /= norm_x_prime;
            
            // y_prime
            std::vector<float> y_prime(3);
            y_prime[0] = z_prime[1] * x_prime[2] - z_prime[2] * x_prime[1]; 
            y_prime[1] = z_prime[2] * x_prime[0] - z_prime[0] * x_prime[2]; 
            y_prime[2] = z_prime[0] * x_prime[1] - z_prime[1] * x_prime[0];
            float norm_y_prime = std::sqrt(y_prime[0] * y_prime[0] + y_prime[1] * y_prime[1] + y_prime[2] * y_prime[2]);
            y_prime[0] /= norm_y_prime;
            y_prime[1] /= norm_y_prime;
            y_prime[2] /= norm_y_prime;
            
            // Conversion from local to global
            std::vector<float> leftHitLocalPosition = {-float(distanceToWire), 0.0f, 0.0f};
            std::vector<float> rightHitLocalPosition = {float(distanceToWire), 0.0f, 0.0f};
            
            std::vector<float> leftHitGlobalPosition = local_to_global(leftHitLocalPosition,x_prime,y_prime,z_prime,wire_pos);
            std::vector<float> rightHitGlobalPosition = local_to_global(rightHitLocalPosition,x_prime,y_prime,z_prime,wire_pos); 

            // Add the 3D position of the left hit to the global input list.
            ListGlobalInputs.push_back(leftHitGlobalPosition[0]);
            ListGlobalInputs.push_back(leftHitGlobalPosition[1]);
            ListGlobalInputs.push_back(leftHitGlobalPosition[2]);
            
            // Add the difference between the right and left hit positions to the global input list.
            ListGlobalInputs.push_back(0.0); 
            ListGlobalInputs.push_back(rightHitGlobalPosition[0]-leftHitGlobalPosition[0]);
            ListGlobalInputs.push_back(rightHitGlobalPosition[1]-leftHitGlobalPosition[1]);
            ListGlobalInputs.push_back(rightHitGlobalPosition[2]-leftHitGlobalPosition[2]);
            
            // Store the current index in ListHitType_CDC and increment the global iterator.
            ListHitType_CDC.push_back(it);
            it += 1; 
            it_2 += 1;                      
        }
        // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_2}, torch::kFloat32);


        /////////////////////////////
        ////////// ML STEP //////////
        /////////////////////////////

        // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
        extension::TrackCollection* output_tracks = new extension::TrackCollection();
        genfit::Track track;

        // Return the output collections as a tuple
        return std::make_tuple(std::move(*output_tracks));

    } 

    private:

        
        

};

DECLARE_COMPONENT(GGTF_fitter)