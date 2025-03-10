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

#include "edm4hep/TrackerHitPlaneCollection.h"
using VertexHitsColl = edm4hep::TrackerHitPlaneCollection;

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
#include <DAF.h>
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




#include "IDEAtrack.hpp"

/** @struct GGTF_fitter_IDEAv3
 *
 *
 */
 

struct GGTF_fitter_IDEAv3 final : 
        k4FWCore::MultiTransformer< std::tuple<IntColl>(const VertexHitsColl&,const VertexHitsColl&)> 
            
                                                                                            
{
    GGTF_fitter_IDEAv3(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("vtx_barrel", {"vtx_barrel"}),
                KeyValues("vtx_endcap", {"vtx_endcap"})
            },
            {   
                KeyValues("test", {"test"})      
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
    StatusCode initialize() {

        geoManager = new TGeoManager("Geometry", "IDEA geometry");

        std::string geoPath = geoPath_m.value();
        geoManager->Import(geoPath.c_str());

        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(new genfit::ConstField(0., 0., 20.)); // 2 T

        const auto detector = m_geoSvc->getDetector();
        const auto surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
        surfaceMap = surfMan->map(m_subDetName.value());


        return StatusCode::SUCCESS;

    }

    
    std::tuple<IntColl> operator()( const VertexHitsColl& vtx_barrel, const VertexHitsColl& vtx_endcap) const override 
    {
        
   
        genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
        // genfit::AbsKalmanFitter* fitter = new genfit::DAF();

        // particle pdg code; pion hypothesis
        const int pdg = 211;

        // start values for the fit, e.g. from pattern recognition
        TVector3 pos(0, 0, 0);
        TVector3 mom(0.,0.,80);

        // trackrep and create track
        genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
        genfit::Track fitTrack(rep, pos, mom);

        int vtx_idx(0);
        for (auto hit : vtx_barrel)
        {

            int detID = 1;
            auto vtx_struct = IDEAtracking::VTX_measurement(hit,surfaceMap,detID,++vtx_idx);

            auto measurement = vtx_struct.getGenFit();
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

        }
        vtx_idx = 0;
        for (auto hit : vtx_endcap)
        {
            int detID = 2;
            auto vtx_struct = IDEAtracking::VTX_measurement(hit,surfaceMap,detID,++vtx_idx);

            auto measurement = vtx_struct.getGenFit();
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

        }

        // int track_idx = 0;
        // for (auto track : GGTF_tracks)
        // {   
        //     if (track_idx > 0)
        //     {
        //         auto hits_in_track = track.getTrackerHits();
        //         std::vector<float> mom_estimation;
        //         for (auto hit : hits_in_track)
        //         {
                    
        //             int mom_idx = 0;      
        //             if (hit.isA<edm4hep::TrackerHitPlane>())
        //             {
        //                 auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
        //                 auto pos = vtx_hit.getPosition();

        //                 if (mom_idx == 0 || mom_idx == 1)
        //                 {
        //                     mom_estimation.push_back(pos[0]);
        //                     mom_estimation.push_back(pos[1]);
        //                     mom_estimation.push_back(pos[2]);
        //                 }

        //                 mom_idx+=1;
        //             }    
        //         }

        //         std::vector<float> mom_dir;
        //         mom_dir.push_back(mom_estimation[3]-mom_estimation[0]);
        //         mom_dir.push_back(mom_estimation[4]-mom_estimation[1]);
        //         mom_dir.push_back(mom_estimation[5]-mom_estimation[2]);

        //         // init fitter
        //         genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
        //         // genfit::AbsKalmanFitter* fitter = new genfit::DAF();

        //         // particle pdg code; pion hypothesis
        //         const int pdg = 211;

        //         // start values for the fit, e.g. from pattern recognition
        //         TVector3 pos(0, 0, 0);
        //         TVector3 mom(mom_dir[0]*10, mom_dir[1]*10, mom_dir[2]*10);

        //         // trackrep and create track
        //         genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
        //         genfit::Track fitTrack(rep, pos, mom);


        //         int dc_idx = 0;
        //         for (auto hit : hits_in_track)
        //         {
                    
        //             int vtx_idx(0);
        //             if (hit.isA<edm4hep::TrackerHitPlane>())
        //             {   
        //                 std::cout << vtx_idx << std::endl;
        //                 auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
        //                 int detID = 0; //TODO: set 0 for barrel and 1 for edncap and 2 for dc
        //                 auto vtx_struct = IDEAtracking::VTX_measurement(vtx_hit,surfaceMap,detID,++vtx_idx);
                        
                        
        //                 auto measurement = vtx_struct.getGenFit();
        //                 fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
                        
        //             }

                    
        //             if (hit.isA<extension::SenseWireHit>())
        //             {   
        //                 auto dc_hit =  hit.as<extension::SenseWireHit>();
        //                 int detID = 2; //TODO: set 0 for barrel and 1 for edncap and 2 for dc
        //                 auto dc_struct = IDEAtracking::DC_measurement(dc_hit,detID,vtx_idx);

        //                 auto measurement = dc_struct.getGenFit();
        //                 // fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

        //                 dc_idx+=1;
        //             }


        //         }

        //         for (int k=0; k< dc_idx; k++)
        //         {
        //             auto point = fitTrack.getPoint(k);
        //             point->Print();
        //         }

        //         fitTrack.checkConsistency();
        //         fitter->processTrack(&fitTrack);
        //         fitTrack.checkConsistency();

        //     }
        // }


        IntColl test;
        return std::make_tuple(std::move(test));

    } 

    private:

        TGeoManager* geoManager;
        genfit::MaterialEffects* materialEffects;
        genfit::FieldManager* fieldManager;

        Gaudi::Property<std::string> geoPath_m{this, "geoPath_m", "/eos/user/a/adevita/workDir/k4RecTracker/Tracking/TGeo_IDEA.root", "geoPath_m"};
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

        Gaudi::Property<std::string>  m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
        SmartIF<IGeoSvc> m_geoSvc;

        const dd4hep::rec::SurfaceMap* surfaceMap;

        
        

};

DECLARE_COMPONENT(GGTF_fitter_IDEAv3)