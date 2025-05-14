#ifndef IDEA_TRACK_H
#define IDEA_TRACK_H

// Standard library
#include <algorithm>
#include <cmath>
#include <cstdlib>  // For getenv
#include <fstream>  // For std::ifstream
#include <iostream>
#include <iterator> // For std::istreambuf_iterator
#include <map>
#include <memory> 
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>  // For std::runtime_error
#include <string>
#include <typeinfo>
#include <vector>
#include <filesystem>  // For std::filesystem::path

// ROOT
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <TVectorD.h>

// GenFit
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
#include <TGeoMaterialInterface.h>
#include <Track.h>
#include <TrackPoint.h>
#include "SI_measurement.hpp"
#include "DC_measurement.hpp"

// EDM4hep
#include "extension/MutableTrack.h"
#include "extension/TrackCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

// k4FWCore & k4Interface
#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// Torch
#include <torch/torch.h>

// podio
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "podio/UserDataCollection.h"

// Define collection types for podio
using DoubleColl = podio::UserDataCollection<double>;
using IntColl = podio::UserDataCollection<int>;
using FloatColl = podio::UserDataCollection<float>;

// Custom project headers
#include "utils.hpp"

namespace GENFIT {

class IDEAtrack {
public:
    IDEAtrack(const extension::Track& track, const dd4hep::rec::SurfaceManager* surfMan, const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder);
    ~IDEAtrack();

    void checkInitialization();
    void init(const extension::Track& track_init);
    void initPDG(int pdg) { _particle_hypothesis = pdg; }
    void initPOS(TVector3 pos) { _posInit = pos; }
    void initMOM(TVector3 mom) { _momInit = mom; }

    void createGenFitTrack();
    bool fit(double Beta_init, double Beta_final, double Beta_steps);

    genfit::Track* getTrack_genfit() { return genfitTrack_; }
    genfit::AbsTrackRep* getRep_genfit() { return genfitTrackRep_; }
    extension::MutableTrack& getTrack_edm4hep() { return edm4hepTrack_; }

private:

    int _particle_hypothesis = 211;
    TVector3 _posInit{0., 0., 0.};
    TVector3 _momInit{0., 0., 0.};

    genfit::AbsTrackRep* genfitTrackRep_;    
    genfit::AbsTrackRep* genfitTrackRepBack_;    

    genfit::Track* genfitTrack_;
    genfit::Track* backwardGenfitTrack_;
    
    extension::MutableTrack edm4hepTrack_;  

    const dd4hep::rec::SurfaceManager* _surfMan;
    const dd4hep::rec::SurfaceMap* surfaceMap_vertex;
    const dd4hep::rec::SurfaceMap* surfaceMap_wrapper_barrel;
    const dd4hep::rec::SurfaceMap* surfaceMap_wrapper_endcap;

    const dd4hep::rec::DCH_info* _dch_info;
    const dd4hep::DDSegmentation::BitFieldCoder* _dc_decoder;
};

}  // namespace GENFIT

#endif // IDEA_TRACK_H
