ifndef IDEATRACKING_H
#define IDEATRACKING_H

#include <PlanarMeasurement.h>
#include <StateOnPlane.h>
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <iostream>
#include <TrackPoint.h>

namespace IDEAtracking {

class VTX_measurement {
public:

    // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
    VTX_measurement(const edm4hep::TrackerHit3D& hit, const int det_idx, const int plane_idx, const int hit_idx);

    // Getter for genfit::WirePointMeasurement
    genfit::PlanarMeasurement* getGenFit() const;

private:
    genfit::PlanarMeasurement* genfitHit_; // Pointer to store the measurement
};

} 

#endif // IDETRACKING_H