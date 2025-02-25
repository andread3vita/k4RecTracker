#include "VTX_measurement.hpp"
#include "TDecompSVD.h"

namespace IDEAtracking {

// Constructor implementation
VTX_measurement::VTX_measurement(const edm4hep::TrackerHit3D& hit, const int det_idx, const int plane_idx, const int hit_idx) {
    

    TVectorD rawHitCoords(2);
    auto position = hit.getPosition();
    // auto position_on_plane = convert_to_local(position);
    rawHitCoords[0] = position_on_plane[0];
    rawHitCoords[1] = position_on_plane[1];

    TMatrixDSym rawHitCov(2);
    auto cov_matrix = hit.getCovMatrix();
    // auto covariance_on_plane = convert_cov_to_local(cov_matrix);
    rawHitCov[0,0] = covariance_on_plane[0,0];
    rawHitCov[0,1] = covariance_on_plane[0,1];
    rawHitCov[1,0] = covariance_on_plane[1,0];
    rawHitCov[1,1] = covariance_on_plane[1,1];

    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(nullptr);
    genfitHit_ = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, trackPoint);

    // add plane
    // TVector3 o = hit.getOrigin();
    // TVector3 u = hit.getU();
    // TVector3 v = hit.getV();
    // genfitHit_->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u,v)), plane_idx);
}

// Get implementation
genfit::PlanarMeasurement* VTX_measurement::getGenFit() const {
    return genfitHit_;
}

} // namespace IDEtracking
