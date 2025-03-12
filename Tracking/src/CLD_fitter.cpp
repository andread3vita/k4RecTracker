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

#include "edm4hep/TrackCollection.h"
using CLD_trackColl = edm4hep::TrackCollection;

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



/** @struct GGTF_fitter_IDEAv3
 *
 *
 */
 

struct GGTF_fitter_IDEAv3 final : 
        k4FWCore::MultiTransformer< std::tuple<IntColl>(const VertexHitsColl&,const VertexHitsColl&,const VertexHitsColl&,const VertexHitsColl&,const VertexHitsColl&,const VertexHitsColl&,const CLD_trackColl&)> 
            
                                                                                            
{
    GGTF_fitter_IDEAv3(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("InnerTrackerBarrelCollection", {"Tracks from CLD"}),
                KeyValues("InnerTrackerEndcapCollection", {"Tracks from CLD"}),
                KeyValues("OuterTrackerBarrelCollection", {"Tracks from CLD"}),
                KeyValues("OuterTrackerEndcapCollection", {"Tracks from CLD"}),
                KeyValues("VertexBarrelCollection", {"Tracks from CLD"}),
                KeyValues("VertexEndcapCollection", {"Tracks from CLD"}),
                KeyValues("tracks", {"Tracks from CLD"}),
            },
            {   
                KeyValues("test", {"test"})      
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
    StatusCode initialize() {

        geoManager = new TGeoManager("Geometry", "CLD geometry");

        std::string geoPath = geoPath_m.value();
        geoManager->Import(geoPath.c_str());

        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(new genfit::ConstField(0., 0., 20.)); // 20 kGauss

        const auto detector = m_geoSvc->getDetector();
        const auto surfMan = detector->extension<dd4hep::rec::SurfaceManager>();

        surfaceMap_vertex = surfMan->map("Vertex");
        surfaceMap_InnerTrackers = surfMan->map("InnerTrackers");
        surfaceMap_OuterTrackers = surfMan->map("OuterTrackers");

        return StatusCode::SUCCESS;

    }

    
    std::tuple<IntColl> operator()( const VertexHitsColl& InnerTrackerBarrelCollection,
                                    const VertexHitsColl& InnerTrackerEndcapCollection,
                                    const VertexHitsColl& OuterTrackerBarrelCollection,
                                    const VertexHitsColl& OuterTrackerEndcapCollection,
                                    const VertexHitsColl& VertexBarrelCollection,
                                    const VertexHitsColl& VertexEndcapCollection,
                                    const CLD_trackColl& tracks) const override 
    {
        
       
        int num_hits = InnerTrackerBarrelCollection.size() +InnerTrackerEndcapCollection.size() +OuterTrackerBarrelCollection.size() +OuterTrackerEndcapCollection.size() +VertexBarrelCollection.size() +VertexEndcapCollection.size();

        int track_num = 0;
        for (auto track : tracks)
        {
            int hits_in_track = track.getTrackerHits().size();

            if (hits_in_track == num_hits)
            {   
                 genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

                // particle pdg code; muon hypothesis
                const int pdg = 211;

                // start values for the fit, e.g. from pattern recognition
                TVector3 pos(0, 0, 0);
                TVector3 mom(2.,0.,0.);

                // trackrep and create track
                genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
                genfit::Track fitTrack(rep, pos, mom);

                ++track_num;
                int vtx_idx(0);
                for (auto hit : InnerTrackerBarrelCollection)
                {
                    int detID(0);
                    auto cellID0 = hit.getCellID();
                    dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_InnerTrackers->find(cellID0);
                    const dd4hep::rec::ISurface* surf  = sI->second;
                    dd4hep::rec::Vector3D u = surf->u();
                    // dd4hep::rec::Vector3D v = surf->v();
                    // dd4hep::rec::Vector3D origin = surf->origin();
                    
                    // // convert 3D global position to 2D local position
                    // auto pos_hit = hit.getPosition();
                    // dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                    // dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                    // TVectorD rawHitCoords(2);
                    // rawHitCoords[0] = local_pos[0]/10; //cm
                    // rawHitCoords[1] = local_pos[1]/10; //cm

                    // TMatrixDSym rawHitCov(2);
                    // rawHitCov(0,0) = 5e-4;
                    // rawHitCov(0,1) = 0;
                    // rawHitCov(1,0) = 0;
                    // rawHitCov(1,1) = 5e-4;


                    // // create measurement 
                    // genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                    // // add plane
                    // TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                    // TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                    // TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                    // measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                    // fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                    // o.Print();
                    // u_.Print();
                    // v_.Print();

                    // rawHitCov.Print();
                    // rawHitCoords.Print();

                    // std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;

                }

                // for (auto hit : InnerTrackerEndcapCollection)
                // {

                //     int detID(1);
                //     auto cellID0 = hit.getCellID();
                //     dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_InnerTrackers->find(cellID0);
                //     const dd4hep::rec::ISurface* surf  = sI->second;
                //     dd4hep::rec::Vector3D u = surf->u();
                //     dd4hep::rec::Vector3D v = surf->v();
                //     dd4hep::rec::Vector3D origin = surf->origin();
                    
                //     // convert 3D global position to 2D local position
                //     auto pos_hit = hit.getPosition();
                //     dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                //     dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                //     TVectorD rawHitCoords(2);
                //     rawHitCoords[0] = local_pos[0]/10; //cm
                //     rawHitCoords[1] = local_pos[1]/10; //cm

                //     TMatrixDSym rawHitCov(2);
                //     rawHitCov(0,0) = 5e-4;
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = 5e-4;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                //     TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     o.Print();
                //     u_.Print();
                //     v_.Print();

                //     rawHitCov.Print();
                //     rawHitCoords.Print();

                //     std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;
                    
                // }

                // for (auto hit : OuterTrackerBarrelCollection)
                // {
                //     int detID(2);
                //     auto cellID0 = hit.getCellID();
                //     dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_OuterTrackers->find(cellID0);
                //     const dd4hep::rec::ISurface* surf  = sI->second;
                //     dd4hep::rec::Vector3D u = surf->u();
                //     dd4hep::rec::Vector3D v = surf->v();
                //     dd4hep::rec::Vector3D origin = surf->origin();
                    
                //     // convert 3D global position to 2D local position
                //     auto pos_hit = hit.getPosition();
                //     dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                //     dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                //     TVectorD rawHitCoords(2);
                //     rawHitCoords[0] = local_pos[0]/10; //cm
                //     rawHitCoords[1] = local_pos[1]/10; //cm

                //     TMatrixDSym rawHitCov(2);
                //     rawHitCov(0,0) = 5e-4;
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = 5e-4;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                //     TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     o.Print();
                //     u_.Print();
                //     v_.Print();

                //     rawHitCov.Print();
                //     rawHitCoords.Print();

                //     std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;
                // }

                // for (auto hit : OuterTrackerEndcapCollection)
                // {
                //     int detID(0);
                //     auto cellID0 = hit.getCellID();
                //     dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_OuterTrackers->find(cellID0);
                //     const dd4hep::rec::ISurface* surf  = sI->second;
                //     dd4hep::rec::Vector3D u = surf->u();
                //     dd4hep::rec::Vector3D v = surf->v();
                //     dd4hep::rec::Vector3D origin = surf->origin();
                    
                //     // convert 3D global position to 2D local position
                //     auto pos_hit = hit.getPosition();
                //     dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                //     dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                //     TVectorD rawHitCoords(2);
                //     rawHitCoords[0] = local_pos[0]/10; //cm
                //     rawHitCoords[1] = local_pos[1]/10; //cm

                //     TMatrixDSym rawHitCov(2);
                //     rawHitCov(0,0) = 5e-4;
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = 5e-4;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                //     TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     o.Print();
                //     u_.Print();
                //     v_.Print();

                //     rawHitCov.Print();
                //     rawHitCoords.Print();

                //     std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;
                // }

                // for (auto hit : VertexBarrelCollection)
                // {
                //     int detID(0);
                //     auto cellID0 = hit.getCellID();
                //     dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_vertex->find(cellID0);
                //     const dd4hep::rec::ISurface* surf  = sI->second;
                //     dd4hep::rec::Vector3D u = surf->u();
                //     dd4hep::rec::Vector3D v = surf->v();
                //     dd4hep::rec::Vector3D origin = surf->origin();
                    
                //     // convert 3D global position to 2D local position
                //     auto pos_hit = hit.getPosition();
                //     dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                //     dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                //     TVectorD rawHitCoords(2);
                //     rawHitCoords[0] = local_pos[0]/10; //cm
                //     rawHitCoords[1] = local_pos[1]/10; //cm

                //     TMatrixDSym rawHitCov(2);
                //     rawHitCov(0,0) = 5e-4;
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = 5e-4;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                //     TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     o.Print();
                //     u_.Print();
                //     v_.Print();

                //     rawHitCov.Print();
                //     rawHitCoords.Print();

                //     std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;
                // }

                // for (auto hit : VertexEndcapCollection)
                // {
                //     int detID(0);
                //     auto cellID0 = hit.getCellID();
                //     dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap_vertex->find(cellID0);
                //     const dd4hep::rec::ISurface* surf  = sI->second;
                //     dd4hep::rec::Vector3D u = surf->u();
                //     dd4hep::rec::Vector3D v = surf->v();
                //     dd4hep::rec::Vector3D origin = surf->origin();
                    
                //     // convert 3D global position to 2D local position
                //     auto pos_hit = hit.getPosition();
                //     dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                //     dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                //     TVectorD rawHitCoords(2);
                //     rawHitCoords[0] = local_pos[0]/10; //cm
                //     rawHitCoords[1] = local_pos[1]/10; //cm

                //     TMatrixDSym rawHitCov(2);
                //     rawHitCov(0,0) = 5e-4;
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = 5e-4;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0]/10,u[1]/10,u[2]/10);
                //     TVector3 v_(v[0]/10,v[1]/10,v[2]/10);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     o.Print();
                //     u_.Print();
                //     v_.Print();

                //     rawHitCov.Print();
                //     rawHitCoords.Print();

                //     std::cout << "Plane:" << cellID0 << " hit: " << vtx_idx-1 << " Det: " << detID << std::endl;
                // }

                // std::cout << "-------------" << std::endl;
                // std::cout << "-------------" << std::endl;
        
        
                // fitter->processTrack(&fitTrack);
                // fitTrack.getFittedState().Print();

            }
        }

        std::cout << "-------------------" << std::endl;
        IntColl test;
        return std::make_tuple(std::move(test));

    } 

    private:

        TGeoManager* geoManager;
        genfit::MaterialEffects* materialEffects;
        genfit::FieldManager* fieldManager;

        Gaudi::Property<std::string> geoPath_m{this, "geoPath_m", "/afs/cern.ch/user/a/adevita/public/workDir/test/fitting_validation/CLD_validation/TGeo_CLD.root", "geoPath_m"};
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

        SmartIF<IGeoSvc> m_geoSvc;

        const dd4hep::rec::SurfaceMap* surfaceMap_vertex;
        const dd4hep::rec::SurfaceMap* surfaceMap_InnerTrackers;
        const dd4hep::rec::SurfaceMap* surfaceMap_OuterTrackers;

        
        

};

DECLARE_COMPONENT(GGTF_fitter_IDEAv3)


