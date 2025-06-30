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

// //=== Standard Library ===
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>  
// #include <fstream>  
// #include <filesystem>  
// #include <iostream>
// #include <iterator> 
// #include <map>
// #include <memory>
// #include <numeric>
// #include <queue>
// #include <sstream>
// #include <string>
// #include <typeinfo>
// #include <vector>

// //=== ROOT / C++ ABI ===
// #include <TEveManager.h>
// #include <TGeoManager.h>
// #include <TVector3.h>
// #include <cxxabi.h>

// //=== GenFit ===
// #include <ConstField.h>
// #include <DAF.h>
// #include <EventDisplay.h>
// #include <Exception.h>
// #include <FieldManager.h>
// #include <KalmanFitterInfo.h>
// #include <KalmanFitterRefTrack.h>
// #include <MaterialEffects.h>
// #include <PlanarMeasurement.h>
// #include <RKTrackRep.h>
// #include <StateOnPlane.h>
// #include <Track.h>
// #include <TrackPoint.h>
// #include <TGeoMaterialInterface.h>


// //=== Gaudi / k4FWCore ===
// #include "Gaudi/Property.h"
// #include "k4FWCore/DataHandle.h"
// #include "k4FWCore/Transformer.h"

// //=== DD4hep / DDRec / DDSegmentation ===
// #include "DD4hep/DD4hepUnits.h"
// #include "DD4hep/DetType.h"
// #include "DD4hep/DetElement.h"
// #include "DD4hep/Detector.h"
// #include "DD4hep/DetectorSelector.h"
// #include "DD4hep/Readout.h"
// #include <DDRec/DetectorData.h>
// #include <DDRec/Vector3D.h>
// #include <DDSegmentation/BitFieldCoder.h>
// #include "DD4hep/Fields.h"
// #include "DDRec/SurfaceManager.h"
// #include "DetectorChecksum.hpp"
// #include <DD4hep/Printout.h>

// //=== podio / edm4hep ===
// #include "podio/Frame.h"
// #include "podio/ROOTReader.h"
// #include "edm4hep/TrackerHitPlaneCollection.h"
// #include "extension/SenseWireHitCollection.h"
// #include "extension/TrackCollection.h"
// #include "edm4hep/TrackCollection.h"

// //=== k4Interface ===
// #include "k4Interface/IGeoSvc.h"
// #include "k4Interface/IUniqueIDGenSvc.h"

// // Define collection types
// #include "podio/UserDataCollection.h"

// #include "utils.hpp"

// struct checkGeometry final : 
//         k4FWCore::MultiTransformer< std::tuple< extension::TrackCollection>(const edm4hep::TrackCollection&)>                                                                         
// {
//     checkGeometry(const std::string& name, ISvcLocator* svcLoc) : 
//         MultiTransformer ( name, svcLoc,
//             {   KeyValues("tracks_input", {"tracks_input"})  },
//             {   KeyValues("extensionTrackCollection", {"extensionTrackCollection"})   }) {}
         
//     StatusCode initialize() { 
        
//         std::string xml_file = "/eos/user/a/adevita/saveSpace/k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml";
//         dd4hep::Detector& description = Detector::getInstance();
//         description.fromCompact(xml_file); 

//         dd4hep::DetElement top = description.world();
//         dd4hep::detail::DetectorChecksum checksum(description);

//         checksum.analyzeDetector(top);
//         auto hash_entry = checksum.handleHeader();

//         std::cout << "Checksum: " << hash_entry.hash << std::endl;

        
        
        
        
//         return StatusCode::SUCCESS;}

    
//     std::tuple< extension::TrackCollection> operator()( const edm4hep::TrackCollection& tracks_input) const override                                                 
//     {
        
//         extension::TrackCollection extensionTrackCollection; 
//         return std::make_tuple( std::move(extensionTrackCollection));
        
//     } 
    
//     StatusCode finalize() {return StatusCode::SUCCESS;}

//     private:
       


// };

// DECLARE_COMPONENT(checkGeometry)
