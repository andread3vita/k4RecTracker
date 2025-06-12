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
#include <cstdlib>  // getenv
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TVector3.h"

// ONNX / Torch
#include <ATen/ATen.h>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"

// Gaudi Framework
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "k4FWCore/Transformer.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4hep & PODIO
#include "podio/UserDataCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/ParticleIDData.h"

// EDM4hep extensions
#include "extension/TrackCollection.h"
#include "extension/TrackerHit.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DDRec/DCH_info.h"
#include "DDSegmentation/BitFieldCoder.h"

// Custom modules
#include "utils.hpp"
#include "WeaverInterface.hpp"

/** @struct GGTF_tracking_OC
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

struct GGTF_tracking_OC_weiver final : 
        k4FWCore::MultiTransformer< std::tuple<extension::TrackCollection>(   
                                                                    const extension::SenseWireHitCollection&, 
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&,
                                                                    const edm4hep::TrackerHitPlaneCollection&)>                                                                                            
{
    GGTF_tracking_OC_weiver(const std::string& name, ISvcLocator* svcLoc) : 
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

        
        // Create the WeaverInterface object
        // weaver = std::make_unique<WeaverInterface>(model_path, json_path, vars);

        return StatusCode::SUCCESS;

   }

    
    std::tuple<extension::TrackCollection> operator()(      const extension::SenseWireHitCollection& inputHits_CDC, 
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_VTXB,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_VTXD,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_SWB,
                                                            const edm4hep::TrackerHitPlaneCollection& inputHits_SWD) const override 
    {

        



        extension::TrackCollection outputTracks;
        return std::make_tuple(std::move(outputTracks));

    } 

    private:

        SmartIF<IGeoSvc> m_geoSvc;
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

        Gaudi::Property<std::string> model_path{this, "model_path", "/eos/experiment/fcc/ee/GGTF_tracking/model.onnx","Path to the ONNX model"};
        mutable std::unique_ptr<WeaverInterface> weaver;
        

};

DECLARE_COMPONENT(GGTF_tracking_OC_weiver)