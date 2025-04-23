#ifndef DC_MEASUREMENT_H
#define DC_MEASUREMENT_H

// Standard
#include <iostream>

// ROOT
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

// GenFit / reconstruction framework
#include <PlanarMeasurement.h>
#include <RectangularFinitePlane.h>
#include <StateOnPlane.h>
#include <TrackPoint.h>
#include <WirePointMeasurement.h>

// EDM / interface
#include "DDRec/SurfaceManager.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"
#include "DD4hep/Detector.h" 
#include "DDSegmentation/BitFieldCoder.h"
#include "DDRec/DCH_info.h"


namespace GENFIT {

class DC_measurement {
public:

    // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
    DC_measurement(const extension::SenseWireHit& hit, const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder, const int det_idx, const int hit_idx);

    // Getter for genfit::WirePointMeasurement
    genfit::WirePointMeasurement* getGenFit() const;

private:
    genfit::WirePointMeasurement* genfitHit_; // Pointer to store the measurement
};

} 

#endif // GENFIT


// #ifndef DC_MEASUREMENT_H
// #define DC_MEASUREMENT_H

// // Standard
// #include <iostream>

// // ROOT
// #include <TMatrixD.h>
// #include <TMatrixDSym.h>
// #include <TVectorD.h>

// // GenFit / reconstruction framework
// #include <PlanarMeasurement.h>
// #include <RectangularFinitePlane.h>
// #include <StateOnPlane.h>
// #include <TrackPoint.h>
// #include <WireMeasurementNew.h>

// // EDM / interface
// #include "DDRec/SurfaceManager.h"
// #include "extension/SenseWireHitCollection.h"
// #include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
// #include "k4Interface/IGeoSvc.h"
// #include "k4Interface/IUniqueIDGenSvc.h"


// namespace GENFIT {

// class DC_measurement {
// public:

//     // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
//     DC_measurement(const extension::SenseWireHit& hit, const int det_idx, const int hit_idx);

//     // Getter for genfit::WirePointMeasurement
//     genfit::WireMeasurementNew* getGenFit() const;

// private:
//     genfit::WireMeasurementNew* genfitHit_; // Pointer to store the measurement
// };

// } 

// #endif // GENFIT
