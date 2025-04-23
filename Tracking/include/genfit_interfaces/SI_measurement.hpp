#ifndef SI_MEASUREMENT_H
#define SI_MEASUREMENT_H

#include <PlanarMeasurement.h>
#include <StateOnPlane.h>
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <iostream>
#include <TrackPoint.h>
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"
#include "DDRec/SurfaceManager.h"
#include "RectangularFinitePlane.h"

namespace GENFIT {

class SI_measurement {
public:

    // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
    SI_measurement(const edm4hep::TrackerHitPlane& hit, const dd4hep::rec::SurfaceMap* surfaceMap, const int det_idx, const int hit_idx);

    // Getter for genfit::WirePointMeasurement
    genfit::PlanarMeasurement* getGenFit() const;

private:
    genfit::PlanarMeasurement* genfitHit_; // Pointer to store the measurement
};

} 

#endif // SI_MEASUREMENT_H