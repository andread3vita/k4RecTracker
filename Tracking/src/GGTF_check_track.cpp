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

#include "extension/TrackCollection.h"
using TrackColl = extension::TrackCollection;

#include "extension/DriftChamberDigiCollection.h"
using DCHitColl = extension::DriftChamberDigiCollection;

struct GGTF_check final : 
        k4FWCore::MultiTransformer< std::tuple<IntColl>( 
                                                                    const TrackColl&)>
            
                                                                                            
{
    GGTF_check(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("outputTracks", {"outputTracks"})

            },
            {   
                KeyValues("test", {"test"})        
            
            }) {}

    
    std::tuple<IntColl> operator()(
                                                const TrackColl& outputTracks) const override 
    {
        
        for (auto track : outputTracks)
        {
            auto hitColl = track.getTrackerHits();
            
            for (auto& hit : hitColl)
            {
                    std::cout << hit.getType() << std::endl;
            }
        
        }
       
        IntColl test;
        return std::make_tuple(std::move(test));

    } 

    private:
        
        
        

};

DECLARE_COMPONENT(GGTF_check)