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

//=== Standard Library ===
#include <algorithm>
#include <cmath>
#include <cstdlib>  
#include <fstream>  
#include <filesystem>  
#include <iostream>
#include <iterator> 
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

//=== ROOT / C++ ABI ===
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <cxxabi.h>

//=== GenFit ===
#include <ConstField.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <PlanarMeasurement.h>
#include <RKTrackRep.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <TGeoMaterialInterface.h>

//=== Gaudi / k4FWCore ===
#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

//=== DD4hep / DDRec / DDSegmentation ===
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>

//=== podio / edm4hep ===
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackState.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "Wire_measurement.hpp"
#include "GenfitTrack.hpp"
#include "Planar_measurement.hpp"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>
#include "utils.hpp"

// Define collection types
#include "podio/UserDataCollection.h"
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;

#include <DD4hep/Detector.h>
#include <DD4hep/Printout.h>
#include <DD4hep/DetectorTools.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Objects.h>

#include "DetectorCheckSum.hpp" 


/** @struct comparisonXmlRoot

*
*
*/

using namespace dd4hep;
using namespace dd4hep::detail;

struct comparisonXmlRoot final : 
        k4FWCore::MultiTransformer< std::tuple< extension::TrackCollection,IntColl,IntColl >(  const edm4hep::MCParticleCollection&,
                                                                                                                                        const extension::SenseWireHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
                                                                                                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> 
            
                                                                                            
{
    comparisonXmlRoot(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("mcParticles", {"mcParticles"}),
                KeyValues("DC_associations", {"DC_associations"}),
                KeyValues("VTXD_links", {"VTXD_links"}),
                KeyValues("VTXB_links", {"VTXIB_links"}),
                KeyValues("wrapperB_links", {"wrapperB_links"}),
                KeyValues("wrapperD_links", {"wrapperD_links"})
            },
            {   
                KeyValues("fittedTracks", {"fittedTracks"}),
                KeyValues("indexTrack", {"indexTrack"}),
                KeyValues("indexMCparticle", {"indexMCparticle"})
            
            }) {}
         
    StatusCode initialize() {

        
        std::string xml_file = "/eos/user/a/adevita/saveSpace/k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml";
        Detector& description = Detector::getInstance();
        description.fromCompact(xml_file); 

        DetElement top = description.world();
        dd4hep::detail::DetectorChecksum checksum(description);

        checksum.analyzeDetector(top);
        auto hash_entry = checksum.handleHeader();

        std::cout << "Checksum: " << hash_entry.hash << std::endl;

        return StatusCode::SUCCESS;

    }

    
    std::tuple<extension::TrackCollection,IntColl,IntColl> operator()(  const edm4hep::MCParticleCollection& mcParticles,
                                                        const extension::SenseWireHitSimTrackerHitLinkCollection& DC_asso,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXD_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& VTXB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperB_links,
                                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& wrapperD_links) const override 
    {
        
        IntColl indexTrack;
        IntColl indexMCparticle;
        extension::TrackCollection trackColl;

        return std::make_tuple(std::move(trackColl),std::move(indexTrack),std::move(indexMCparticle));
        
    } 

    public:


};

DECLARE_COMPONENT(comparisonXmlRoot)
