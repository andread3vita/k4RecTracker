#include "measurement/DC_measurement.hpp"
#include "TDecompSVD.h"

namespace IDEAtracking {

// Constructor implementation
DC_measurement::DC_measurement(const extension::SenseWireHit& hit, const int det_idx, const int hit_idx) {
    
    float wireStereoAngle = hit.getWireStereoAngle();
    float wireAzimuthalAngle = hit.getWireAzimuthalAngle();
    edm4hep::Vector3d position = hit.getPosition();
    double positionAlongWireError = hit.getPositionAlongWireError();
    float distanceToWire = hit.getDistanceToWire();
    float distanceToWireError = hit.getDistanceToWireError();


    // Wire direction
    float d_x = std::sin(wireStereoAngle) * std::cos(wireAzimuthalAngle);
    float d_y = std::sin(wireStereoAngle) * std::sin(wireAzimuthalAngle);
    float d_z = std::cos(wireStereoAngle);

    // wire extremities (   e.g. w1_x = x_0 + [(w1_z - z_0)/d_z]*d_x   )
    float w1_z = -2000.;
    float w1_x = position.x +(w1_z-position.z)/d_z*d_x;
    float w1_y = position.y +(w1_z-position.z)/d_z*d_y;
    float w2_z = 2000.;
    float w2_x = position.x +(w2_z-position.z)/d_z*d_x;
    float w2_y = position.y +(w2_z-position.z)/d_z*d_y;

    //Rdrift
    double Rdrift = distanceToWire;

    // wire position in local coordinates z_reco = d(left wire extremity , wire_pos)
    float zreco = std::sqrt(std::pow(w1_x-position.x,2)+std::pow(w1_y-position.y,2)+std::pow(w1_z-position.z,2));

    // genfit::WirePointMeasurement
    TVectorD rawHitCoords(8);
    rawHitCoords(0) = w1_x;             // wire1 X
    rawHitCoords(1) = w1_y;             // wire1 Y
    rawHitCoords(2) = w1_z;             // wire1 Z
    rawHitCoords(3) = w2_x;             // wire2 X
    rawHitCoords(4) = w2_y;             // wire2 Y
    rawHitCoords(5) = w2_z;             // wire2 Z
    rawHitCoords(6) = Rdrift;   // Rdrift
    rawHitCoords(7) = zreco;            // zreco

    // Covariance matrix
    // double w1_z_sigma = 0.;
    // double w1_x_sigma = std::sqrt(std::pow(positionAlongWireError,2)+std::pow(positionAlongWireError/d_z*d_x,2));
    // double w1_y_sigma = std::sqrt(std::pow(positionAlongWireError,2)+std::pow(positionAlongWireError/d_z*d_y,2));
    // double w2_z_sigma = 0.;
    // double w2_x_sigma = std::sqrt(std::pow(positionAlongWireError,2)+std::pow(positionAlongWireError/d_z*d_x,2));
    // double w2_y_sigma = std::sqrt(std::pow(positionAlongWireError,2)+std::pow(positionAlongWireError/d_z*d_y,2));

    double w1_z_sigma = 0.;
    double w1_x_sigma = 0.;
    double w1_y_sigma = 0.;
    double w2_z_sigma = 0.;
    double w2_x_sigma = 0.;
    double w2_y_sigma = 0.;

    double Rdrift_sigma = distanceToWireError;

    double dz_dx = (w1_x - position.x) / zreco;
    double dz_dy = (w1_y - position.y) / zreco;
    double dz_dz = (w1_z - position.z) / zreco;

    double zreco_sigma = std::sqrt(
        std::pow(dz_dx*positionAlongWireError, 2) +
        std::pow(dz_dy*positionAlongWireError, 2) +
        std::pow(dz_dz*positionAlongWireError, 2)
    );

    TMatrixDSym rawHitCov(8);

    rawHitCov(0, 0) = w1_x_sigma * w1_x_sigma;          // Variance for w1_x
    rawHitCov(1, 1) = w1_y_sigma * w1_y_sigma;          // Variance for w1_y
    rawHitCov(2, 2) = w1_z_sigma * w1_z_sigma;          // Variance for w1_z
    rawHitCov(3, 3) = w2_x_sigma * w2_x_sigma;          // Variance for w2_x
    rawHitCov(4, 4) = w2_y_sigma * w2_y_sigma;          // Variance for w2_y
    rawHitCov(5, 5) = w2_z_sigma * w2_z_sigma;          // Variance for w2_z
    rawHitCov(6, 6) = Rdrift_sigma * Rdrift_sigma;      // Variance for Rdrift
    rawHitCov(7, 7) = zreco_sigma * zreco_sigma;        // Variance for zreco

    rawHitCov(0, 1) = 0;                                // Covariance between w1_x and w1_y
    rawHitCov(0, 2) = 0;                                // Covariance between w1_x and w1_z
    rawHitCov(0, 3) = 0;                                // Covariance between w1_x and w2_x
    rawHitCov(0, 4) = 0;                                // Covariance between w1_x and w2_y
    rawHitCov(0, 5) = 0;                                // Covariance between w1_x and w2_z
    rawHitCov(0, 6) = 0;                                // Covariance between w1_x and Rdrift
    rawHitCov(0, 7) = 0;                                // Covariance between w1_x and zreco

    rawHitCov(1, 2) = 0;                                // Covariance between w1_y and w1_z
    rawHitCov(1, 3) = 0;                                // Covariance between w1_y and w2_x
    rawHitCov(1, 4) = 0;                                // Covariance between w1_y and w2_y
    rawHitCov(1, 5) = 0;                                // Covariance between w1_y and w2_z
    rawHitCov(1, 6) = 0;                                // Covariance between w1_y and Rdrift
    rawHitCov(1, 7) = 0;                                // Covariance between w1_y and zreco

    rawHitCov(2, 3) = 0;                                // Covariance between w1_z and w2_x
    rawHitCov(2, 4) = 0;                                // Covariance between w1_z and w2_y
    rawHitCov(2, 5) = 0;                                // Covariance between w1_z and w2_z
    rawHitCov(2, 6) = 0;                                // Covariance between w1_z and Rdrift
    rawHitCov(2, 7) = 0;                                // Covariance between w1_z and zreco

    rawHitCov(3, 4) = 0;                                // Covariance between w2_x and w2_y
    rawHitCov(3, 5) = 0;                                // Covariance between w2_x and w2_z
    rawHitCov(3, 6) = 0;                                // Covariance between w2_x and Rdrift
    rawHitCov(3, 7) = 0;                                // Covariance between w2_x and zreco

    rawHitCov(4, 5) = 0;                                // Covariance between w2_y and w2_z
    rawHitCov(4, 6) = 0;                                // Covariance between w2_y and Rdrift
    rawHitCov(4, 7) = 0;                                // Covariance between w2_y and zreco

    rawHitCov(5, 6) = 0;                                // Covariance between w2_z and Rdrift
    rawHitCov(5, 7) = 0;                                // Covariance between w2_z and zreco

    rawHitCov(6, 7) = 0;                                // Covariance between Rdrift and zreco

    rawHitCov(1, 0) = rawHitCov(0, 1);
    rawHitCov(2, 0) = rawHitCov(0, 2);
    rawHitCov(3, 0) = rawHitCov(0, 3);
    rawHitCov(4, 0) = rawHitCov(0, 4);
    rawHitCov(5, 0) = rawHitCov(0, 5);
    rawHitCov(6, 0) = rawHitCov(0, 6);
    rawHitCov(7, 0) = rawHitCov(0, 7);

    rawHitCov(2, 1) = rawHitCov(1, 2);
    rawHitCov(3, 2) = rawHitCov(2, 3);
    rawHitCov(4, 3) = rawHitCov(3, 4);
    rawHitCov(5, 4) = rawHitCov(4, 5);
    rawHitCov(6, 5) = rawHitCov(5, 6);
    rawHitCov(7, 6) = rawHitCov(6, 7);
    

    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(nullptr);
    genfitHit_ = new genfit::WirePointMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, trackPoint);

}

// Get implementation
genfit::WirePointMeasurement* DC_measurement::getGenFit() const {
    return genfitHit_;
}

} // namespace IDEtracking
