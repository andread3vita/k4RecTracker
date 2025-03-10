#ifndef IDEA_FITTER_H
#define IDEA_FITTER_H


#include <KalmanFitterRefTrack.h>
#include <Track.h>
#include <TrackPoint.h>
#include <RKTrackRep.h>

#include <PlanarMeasurement.h>
#include <WirePointMeasurement.h> 

#include <TVectorD.h>
#include <vector>

#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "extension/TrackCollection.h"

namespace IDEAtracking {

class IDEAtrack {
public:

    IDEAtrack(extension::Track track, TVector3 pos, TVector3 mom, int pdg);


    void insertPoint(genfit::AbsMeasurement* measurement);
    void processTrack();

    genfit::Track* getTrack_genfit() const;
    extension::Track* getTrack_edm4hep() const;

private:
   
    genfit::AbsKalmanFitter* genfitFitter_;      
    genfit::AbsTrackRep* genfitTrackRep_;  
          
    genfit::Track* genfitTrack_;  
    extension::Track* extensionTrack_;  


};

} 

#endif // IDEA_FITTER_H