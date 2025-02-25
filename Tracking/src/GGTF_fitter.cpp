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
#include "DC_measurement.hpp"
#include "VTX_measurement.hpp"

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


//genfit
#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <PlanarMeasurement.h>

#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <vector>


#include "IDEAtrackFitter.hpp"

/** @struct GGTF_fitter
 *
 *
 */

struct GGTF_fitter final : 
        k4FWCore::MultiTransformer< std::tuple<IntColl>( 
                                                                    
                                                                    const TrackColl&)> 
            
                                                                                            
{
    GGTF_fitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("input_tracks", {"input_tracks"})
            },
            {   
                KeyValues("test", {"test"})      
            
            }) {}
    
    StatusCode initialize() {


        return StatusCode::SUCCESS;

    }

    
    std::tuple<IntColl> operator()(   const TrackColl& GGTF_tracks) const override 
    {

        // init geometry and mag. field
        new TGeoManager("Geometry", "IDEA geometry");
        TGeoManager::Import("/eos/user/a/adevita/workDir/k4RecTracker/Tracking/TGeo_IDEA.root");
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,0., 20.)); // 2 T

        
        for (auto track : GGTF_tracks)
        {
            // particle pdg code; pion hypothesis
            const int pdg = 211;

            // start values for the fit, e.g. from pattern recognition
            TVector3 pos(0, 0, 0);
            TVector3 mom(0, 0, 3);

            IDEAtracking::IDEAtrackFitter fitter_framework(pos,mom,pdg);
            
            auto hits_in_track = track.getTrackerHits();

            int vtx_idx = 0;
            int dc_idx = 0;
            for (auto hit : hits_in_track)
            {   
                if (hit.isA<edm4hep::TrackerHit3D>())
                {
                    int detID = 0;
                    int hitID = vtx_idx;
                    const auto vtx_hit =  hit.as<edm4hep::TrackerHit3D>();
                    IDEAtracking::VTX_measurement vtx_measure = IDEAtracking::VTX_measurement(vtx_hit,detID,hitID);
                    vtx_idx += 1;

                    auto measurement = vtx_measure.getGenFit();
                    fitter_framework.insertPoint(measurement);

           
                }
                else if (hit.isA<extension::SenseWireHit>())
                {
                    int detID = 1;
                    int hitID = dc_idx;
                    const auto dc_hit =  hit.as<extension::SenseWireHit>();
                    IDEAtracking::DC_measurement dc_measure = IDEAtracking::DC_measurement(dc_hit,detID,hitID);
                    dc_idx += 1;

                    auto measurement = dc_measure.getGenFit();
                    fitter_framework.insertPoint(measurement);
                }


            }

        }

       
        IntColl test;
        return std::make_tuple(std::move(test));

    } 

    private:

        
        

};

DECLARE_COMPONENT(GGTF_fitter)