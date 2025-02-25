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
        k4FWCore::MultiTransformer< std::tuple<IntColl>( 
                                                                    
                                                                    const DCHitsColl&)> 
            
                                                                                            
{
    GGTF_fitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputHits_CDC", {"inputHits_CDC"})
            },
            {   
                KeyValues("test", {"test"})      
            
            }) {}
    
    StatusCode initialize() {


        return StatusCode::SUCCESS;

    }

    
    std::tuple<IntColl> operator()(   const DCHitsColl& inputHits_CDC) const override 
    {

        // init geometry and mag. field
        new TGeoManager("Geometry", "IDEA geometry");
        TGeoManager::Import("/eos/user/a/adevita/workDir/k4RecTracker/Tracking/TGeo_IDEA.root");
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,0., 20.)); // 2 T

        // init fitter
        genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();


        // particle pdg code; pion hypothesis
        const int pdg = 211;

        // start values for the fit, e.g. from pattern recognition
        TVector3 pos(0, 0, 0);
        TVector3 mom(0, 0, 3);

        // trackrep
        genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

        // create track
        genfit::Track fitTrack(rep, pos, mom);
        

        int dc_idx = 0;
        for (const auto hit : inputHits_CDC)
        {   
            
            int detID = 1;
            int hitID = dc_idx;

            IDEAtracking::DC_measurement dc_hit = IDEAtracking::DC_measurement(hit,detID,hitID);
            dc_idx += 1;

            auto measurement = dc_hit.getGenFit();
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

        }

        //check
        fitTrack.checkConsistency();

        // // do the fit
        // fitter->processTrack(&fitTrack);

        // fitTrack.checkConsistency();

        delete fitter;

        // Return the output collections as a tuple
        IntColl test;
        return std::make_tuple(std::move(test));

    } 

    private:

        
        

};

DECLARE_COMPONENT(GGTF_fitter)