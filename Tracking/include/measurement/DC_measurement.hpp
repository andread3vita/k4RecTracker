#ifndef IDEATRACKING_H
#define IDEATRACKING_H

#include <WirePointMeasurement.h> 
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <iostream>
#include <TrackPoint.h>

namespace IDEAtracking {

class DC_measurement {
public:

    // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
    DC_measurement(const extension::SenseWireHit& hit, const int det_idx, const int hit_idx);

    // Getter for genfit::WirePointMeasurement
    genfit::WirePointMeasurement* getGenFit() const;

private:
    genfit::WirePointMeasurement* genfitHit_; // Pointer to store the measurement
};

} 

#endif // IDETRACKING_H
