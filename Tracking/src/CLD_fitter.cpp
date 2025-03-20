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

#include "edm4hep/ReconstructedParticleCollection.h"
using ParticleColl = edm4hep::ReconstructedParticleCollection;

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

#include "podio/ROOTReader.h"
#include "podio/Frame.h"
#include <KalmanFitterInfo.h>


/** @struct GGTF_fitter_IDEAv3
 *
 *
 */
 

struct CLD_fitter final : 
        k4FWCore::MultiTransformer< std::tuple<DoubleColl,DoubleColl,DoubleColl,DoubleColl,DoubleColl,DoubleColl>(const ParticleColl&,const CLD_trackColl&)> 
            
                                                                                            
{
    CLD_fitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("MCparticles", {"MCparticles"}),
                KeyValues("tracks", {"Tracks from CLD"}),
            },
            {   
                KeyValues("costheta_genfit", {"test"}),
                KeyValues("costheta_CLD", {"test"}),
                KeyValues("pt_genfit", {"test"}),
                KeyValues("pt_CLD", {"test"}),
                KeyValues("chi2_genfit", {"test"}),
                KeyValues("chi2_CLD", {"test"})      
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
    StatusCode initialize() {

        geoManager = new TGeoManager("Geometry", "CLD geometry");
        geoManager->Import(geoPath_m.value().c_str());

        materialEffects = genfit::MaterialEffects::getInstance();
        materialEffects->init(new genfit::TGeoMaterialInterface());

        fieldManager = genfit::FieldManager::getInstance();
        fieldManager->init(new genfit::ConstField(0., 0., m_Bz.value())); // kGauss

        const auto detector = m_geoSvc->getDetector();
        const auto surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
        surfaceMap_vertex = surfMan->map("Vertex");
        surfaceMap_InnerTrackers = surfMan->map("InnerTrackers");
        surfaceMap_OuterTrackers = surfMan->map("OuterTrackers");

        return StatusCode::SUCCESS;

    }

    
    std::tuple<DoubleColl,DoubleColl,DoubleColl,DoubleColl,DoubleColl,DoubleColl> operator()( const ParticleColl& MCparticles,const CLD_trackColl& tracks) const override 
    {
        
        DoubleColl pt_genfit;
        DoubleColl pt_CLD;
        DoubleColl costheta_genfit;
        DoubleColl costheta_CLD;
        DoubleColl chi2_CLD;
        DoubleColl chi2_genfit;

        double pt = 0.;
        double costheta = 0.;
        edm4hep::Track track;
        bool isTrack = false;
        for (auto MCparticle : MCparticles)
        {
            auto idx = MCparticle.getPDG();

            
            if (idx == 211)
            {
                auto p = MCparticle.getMomentum();
                TVector3 mom(p.x, p.y, p.z);
                pt = mom.Pt();
                costheta = mom.CosTheta();
                auto pi_track = MCparticle.getTracks();

                if (pi_track.size())
                {
                    track = pi_track[0];
                    isTrack = true;
                }

            }

        }

        std::cout << "True pt: " << pt << std::endl; 
        std::cout << "True cos(theta): " << costheta << std::endl;
        std::cout << "" << std::endl;

        if(isTrack)
        {
            auto hits_in_track = track.getTrackerHits();

            if (1)
            {   

                ///////////////////////////////////////////////////
                //////////////// CLDreconstruction ////////////////
                ///////////////////////////////////////////////////

                double c_light = 2.99792458e8;
                // double mchp = 0.139570;
                double a = c_light * 1e3 * 1e-15;
                double Bz = 2.0;

                auto trackstate = track.getTrackStates()[0];
               
                double omega = trackstate.omega;
                double pt_CLD_val = a * Bz / abs(omega);
                double phi = trackstate.phi;
                double pz = trackstate.tanLambda * pt_CLD_val;
                double px = pt_CLD_val * std::cos(phi);
                double py = pt_CLD_val * std::sin(phi);
                double p_CLD = std::sqrt(px * px + py * py + pz * pz);

                // double energy = std::sqrt(p * p + mchp * mchp);
                // double theta = std::acos(pz / p);
                // std::cout << "Energy: " << energy << std::endl;
                // std::cout << "Theta: " << theta << std::endl;

                TVector3 mom_CLD = TVector3(px, py, pz);
                std::cout << "CLDreco Momentum: " << "[ " << mom_CLD.X() << " " << mom_CLD.Y() << " " << mom_CLD.Z() << " ] -> " << p_CLD << std::endl;
                std::cout << "CLDreco Pt: " << mom_CLD.Pt() << std::endl;
                std::cout << "CLDreco cosTheta: " << mom_CLD.CosTheta() << std::endl;
                std::cout << "CLDreco Phi: " << mom_CLD.Phi() << std::endl;
               
                float chi2CLD_val = track.getChi2();
                int ndfCLD_val = track.getNdf();
                chi2_CLD.push_back((ndfCLD_val > 0 ? chi2CLD_val / ndfCLD_val : -1));
                
                costheta_CLD.push_back(std::abs(pz / p_CLD));
                pt_CLD.push_back(pt_CLD_val);

                std::cout << "CLD Chi2: " << chi2CLD_val << ", CLD NDF: " << ndfCLD_val << ", CLD Chi2/NDF: " << (ndfCLD_val > 0 ? chi2CLD_val / ndfCLD_val : -1) << std::endl;
                std::cout << "Number of hits in CLD track: " << hits_in_track.size() << std::endl;
                std::cout << "" << std::endl;

                ////////////////////////////////////////
                //////////////// GENFIT ////////////////
                ////////////////////////////////////////

                genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
                const int pdg = 211;        // pion
                TVector3 pos(0, 0, 0);      // cm
                TVector3 mom(0.,0.,10.);    // GeV/c

                // trackrep and create track
                genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
                genfit::Track fitTrack(rep, pos, mom);
                
                int vtx_idx(0);
                for (auto hit : hits_in_track)
                {

                    auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
                    auto cellID0 = vtx_hit.getCellID();

                    int detID(0);
                    dd4hep::rec::SurfaceMap::const_iterator sI;
                    if (surfaceMap_vertex && (sI = surfaceMap_vertex->find(cellID0)) != surfaceMap_vertex->end()) {
                        detID= 0;
                        
                    } 
                    else if (surfaceMap_InnerTrackers && (sI = surfaceMap_InnerTrackers->find(cellID0)) != surfaceMap_InnerTrackers->end()) {
                        detID= 1;
                        
                    } 
                    else if (surfaceMap_OuterTrackers && (sI = surfaceMap_OuterTrackers->find(cellID0)) != surfaceMap_OuterTrackers->end()) {
                        detID= 2;
                        
                    } 
                    else {
                        std::cerr << "Error: Surface not found for cellID: " << cellID0 << std::endl;
                    }


                    const dd4hep::rec::ISurface* surf  = sI->second;
                    dd4hep::rec::Vector3D u = surf->u();
                    dd4hep::rec::Vector3D v = surf->v();
                    dd4hep::rec::Vector3D origin = surf->origin();
                    
                    // convert 3D global position to 2D local position
                    auto pos_hit = vtx_hit.getPosition();
                    dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);
                    dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
                    
                    TVectorD rawHitCoords(2);
                    rawHitCoords[0] = local_pos[0]/10; //cm
                    rawHitCoords[1] = local_pos[1]/10; //cm

                    TMatrixDSym rawHitCov(2);
                    rawHitCov(0,0) = std::pow(vtx_hit.getDu()/10,2); // cm^2
                    rawHitCov(0,1) = 0;
                    rawHitCov(1,0) = 0;
                    rawHitCov(1,1) = std::pow(vtx_hit.getDv()/10,2); // cm^2
            
                    // create measurement 
                    genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                    // add plane
                    TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                    TVector3 u_(u[0],u[1],u[2]);
                    TVector3 v_(v[0],v[1],v[2]);
                    measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);

                    // fill track with measurement
                    fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));


                }
                
                // genfit fit
                try {
                    fitter->processTrack(&fitTrack);
                
                    TVector3 position, momentum;
                    TMatrixDSym covariance(6);
                
                    fitTrack.getFittedState().getPosMomCov(position, momentum, covariance);
                    double genfit_pt_val = momentum.Pt();
                    double genfit_costheta_val = momentum.CosTheta();

                    std::cout << "Genfit momentum: " << "[ " << momentum[0] << " " << momentum[1] << " " << momentum[2] << " ] -> " << momentum.Mag() << std::endl;
                    std::cout << "Genfit Pt: " << genfit_pt_val << std::endl;
                    std::cout << "Genfit cosTheta: " << momentum.CosTheta() << std::endl;
                    std::cout << "Genfit Phi: " << momentum.Phi() << std::endl;

                    double chi2_val = fitTrack.getFitStatus()->getChi2();
                    int ndf = fitTrack.getFitStatus()->getNdf();

                    std::cout << "Chi2: " << chi2_val << ", NDF: " << ndf << ", Chi2/NDF: " << (ndf > 0 ? chi2_val / ndf : -1) << std::endl;


                    std::cout << "Total TrackPoints: " << fitTrack.getNumPoints() << std::endl;
                    int trackInfo(0);
                    for (unsigned int i = 0; i < fitTrack.getNumPoints(); ++i) {
                        auto* tp = fitTrack.getPoint(i);
                        const genfit::KalmanFitterInfo* fi = dynamic_cast<const genfit::KalmanFitterInfo*>(tp->getFitterInfo(fitTrack.getCardinalRep()));

                        if (fi) 
                        {
                            ++trackInfo;
                        }
                    }
                    std::cout << "Total TrackPoints with info: " << trackInfo << std::endl;

                    chi2_genfit.push_back((ndf > 0 ? chi2_val / ndf : -1));

                    if ((ndf > 0 ? chi2_val / ndf : -1)>0)
                    {
                        pt_genfit.push_back(genfit_pt_val);
                        costheta_genfit.push_back(std::abs(genfit_costheta_val));
                    }
                    else
                    {
                        pt_genfit.push_back(-1);
                        costheta_genfit.push_back(-1);
                    }
                    
                    
                } catch (const genfit::Exception& e) {
                
                    std::cerr << "Genfit failure" << std::endl;
                    chi2_genfit.push_back(-1);
                    pt_genfit.push_back(-1);
                    costheta_genfit.push_back(-1);

                } catch (const std::exception& e) {
                    
                    std::cerr << "Genfit failure" << std::endl;
                    chi2_genfit.push_back(-1);
                    pt_genfit.push_back(-1);
                    costheta_genfit.push_back(-1);

                } catch (...) {
                    
                    std::cerr << "Genfit failure" << std::endl;
                    chi2_genfit.push_back(-1);
                    pt_genfit.push_back(-1);
                    costheta_genfit.push_back(-1);

                }
            }
        }
        std::cout << "-------------------" << std::endl;

        return std::make_tuple(std::move(costheta_genfit),std::move(costheta_CLD),std::move(pt_genfit),std::move(pt_CLD),std::move(chi2_genfit),std::move(chi2_CLD));

    } 

    private:

        TGeoManager* geoManager;
        genfit::MaterialEffects* materialEffects;
        genfit::FieldManager* fieldManager;

        Gaudi::Property<std::string> geoPath_m{this, "geoPath_m", "/afs/cern.ch/user/a/adevita/public/workDir/test/fitting_validation/CLD_validation/TGeo_CLD.root", "geoPath_m"};
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        Gaudi::Property<double> m_Bz{this, "Bz", 2., "Bz value"};

        SmartIF<IGeoSvc> m_geoSvc;

        const dd4hep::rec::SurfaceMap* surfaceMap_vertex;
        const dd4hep::rec::SurfaceMap* surfaceMap_InnerTrackers;
        const dd4hep::rec::SurfaceMap* surfaceMap_OuterTrackers;  

};

DECLARE_COMPONENT(CLD_fitter)


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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                // }

                // vtx_idx = 0;
                // for (auto hit : VertexEndcapCollection)
                // {
                //     int detID(1);
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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                // }

                // vtx_idx = 0;
                // for (auto hit : InnerTrackerBarrelCollection)
                // {
                //     int detID(2);
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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                    


                // }

                // vtx_idx = 0;
                // for (auto hit : InnerTrackerEndcapCollection)
                // {

                //     int detID(3);
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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                    
                // }

                // vtx_idx = 0;
                // for (auto hit : OuterTrackerBarrelCollection)
                // {
                //     int detID(4);
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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;


                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);

                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                // }

                // vtx_idx = 0;
                // for (auto hit : OuterTrackerEndcapCollection)
                // {
                //     int detID(5);
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
                //     rawHitCov(0,0) = std::pow(hit.getDu()/10,2);
                //     rawHitCov(0,1) = 0;
                //     rawHitCov(1,0) = 0;
                //     rawHitCov(1,1) = std::pow(hit.getDv()/10,2);
                //     // rawHitCov.Print();

                //     // TMatrixDSym rawHitCov(2);
                //     // rawHitCov(0,0) = 5e-9;
                //     // rawHitCov(0,1) = 0;
                //     // rawHitCov(1,0) = 0;
                //     // rawHitCov(1,1) = 5e-9;



                //     // create measurement 
                //     genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, detID, ++vtx_idx, nullptr);

                //     // add plane
                //     TVector3 o(origin[0]/10,origin[1]/10,origin[2]/10);
                //     TVector3 u_(u[0],u[1],u[2]);
                //     TVector3 v_(v[0],v[1],v[2]);
                    
                //     measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);
                //     fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                //     // std::cout << "rawHitCoords" << std::endl;
                //     // rawHitCoords.Print();
                //     // std::cout << "rawHitCov" << std::endl;
                //     // rawHitCov.Print();

                //     // std::cout << "O" << std::endl;
                //     // o.Print();
                //     // std::cout << "U" << std::endl;
                //     // u_.Print();
                //     // std::cout << "V" << std::endl;
                //     // v_.Print();

                //     // std::cout << "DetID: " << detID << std::endl;
                //     // std::cout << "PlaneID: " << cellID0 << std::endl;
                //     // std::cout << "HitID: " << vtx_idx-1 << std::endl;
                // }


                // std::cout << "-------------" << std::endl;
                // std::cout << "-------------" << std::endl;