/*
 * Copyright (c) 2020-2024 Key4hep-Project.
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

#include "GenfitWireMeasurement.hpp"
#include "TDecompSVD.h"

namespace GenfitInterface {

    // Constructor implementation
    WireMeasurement::WireMeasurement(const edm4hep::SenseWireHit& hit, const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder,const int det_idx, const int hit_idx, const int debug_lvl) {
        
        int cellid = hit.getCellID();

        double wireStereoAngle = hit.getWireStereoAngle();
        double wireAzimuthalAngle = hit.getWireAzimuthalAngle();
        auto pos = hit.getPosition();

        TVector3 position(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);  // cm
        double positionAlongWireError = hit.getPositionAlongWireError() * dd4hep::mm;   // cm

        double distanceToWire = hit.getDistanceToWire() * dd4hep::mm;                    // cm
        double distanceToWireError = hit.getDistanceToWireError() * dd4hep::mm;          // cm

        // Wire direction
        TVector3 direction(0,0,1);
        direction.RotateX(wireStereoAngle);
        direction.RotateZ(wireAzimuthalAngle);

        // Compute wire extremities
        int ilayer          = dch_info->CalculateILayerFromCellIDFields(decoder->get(cellid, "layer"), decoder->get(cellid, "superlayer"));
        int nphi            = decoder->get(cellid, "nphi");
        auto& l             = dch_info->database.at(ilayer);
        int stereosign      = l.StereoSign();
        double rz0          = l.radius_sw_z0;
        double dphi         = dch_info->twist_angle;
        double kappa        = (1. / dch_info->Lhalf) * tan(dphi / 2);
        
        // point 1
        double x1 = rz0;                                          // cm
        double y1 = -stereosign * rz0 * kappa * dch_info->Lhalf;  // cm
        double z1 = -dch_info->Lhalf;                             // cm
        TVector3 p1(x1, y1, z1);
        
        // point 2
        double x2 = rz0;                                         // cm
        double y2 = stereosign * rz0 * kappa * dch_info->Lhalf;  // cm
        double z2 = dch_info->Lhalf;                             // cm
        TVector3 p2(x2, y2, z2);
        
        // calculate phi rotation of whole twisted tube, ie, rotation at z=0
        double phi_z0 = dch_info->Calculate_wire_phi_z0(ilayer, nphi);
        p1.RotateZ(phi_z0);
        p2.RotateZ(phi_z0);

        // Compute drift distance
        double Rdrift = distanceToWire;       // cm

        // Compute position along the wire in local coordinates: z_reco = distance(left wire extremity, wire_pos)
        double zreco = (position - p1).Mag();  //cm

        // Create the coordinate array (8D): p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, drift distance, position along the wire
        TVectorD rawHitCoords(8);
        rawHitCoords(0) = p1.x();               // wire1 X
        rawHitCoords(1) = p1.y();               // wire1 Y
        rawHitCoords(2) = p1.z();               // wire1 Z
        rawHitCoords(3) = p2.x();               // wire2 X
        rawHitCoords(4) = p2.y();               // wire2 Y
        rawHitCoords(5) = p2.z();               // wire2 Z
        rawHitCoords(6) = Rdrift;               // Rdrift
        rawHitCoords(7) = zreco;                // zreco

        // Create covariance matrix
        double sigma_wire = 1e-4;               // almost 0 because they are fixed positions
        double w1_z_sigma = sigma_wire;
        double w1_x_sigma = sigma_wire;
        double w1_y_sigma = sigma_wire;
        double w2_z_sigma = sigma_wire;
        double w2_x_sigma = sigma_wire;
        double w2_y_sigma = sigma_wire;

        double Rdrift_sigma = distanceToWireError;  // cm

        // double dz_dx = ( p1.x() - position.x() ) / zreco; // cm
        // double dz_dy = ( p1.y() - position.y() ) / zreco; // cm
        // double dz_dz = ( p1.z() - position.z() ) / zreco; // cms

        // double zreco_sigma = std::sqrt(
        //     std::pow(dz_dx * positionAlongWireError, 2) +
        //     std::pow(dz_dy * positionAlongWireError, 2) +
        //     std::pow(dz_dz * positionAlongWireError, 2)
        // ); // cm

        double zreco_sigma = positionAlongWireError; // cm

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

        // // // inflate errors
        // // double sigmaNoise = 1;
        // // for (int i = 6; i < 8; ++i) {
        // //     rawHitCov(i, i) += sigmaNoise * sigmaNoise;
        // // }
       
        // rawHitCov(6, 6) += Rdrift_sigma * 1.5 * Rdrift_sigma * 1.5;
        // rawHitCov(7, 7) += zreco_sigma * 1.5 * zreco_sigma * 1.5;

        // Create the genfit::WirePointMeasurement
        m_genfitHit = new genfit::WirePointMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, nullptr);

        if (debug_lvl > 0) {

            std::cout << "Wire measurement settings:" << std::endl;
            std::cout << "Position along the wire [cm]: " << position.X() << " " << position.Y() << " " << position.Z() << std::endl;
            std::cout << "Position along the wire error [cm]: " << positionAlongWireError << std::endl;
            std::cout << "Distance to the wire [cm]: " << distanceToWire << std::endl;
            std::cout << "Distance to the wire error [cm]: " << distanceToWireError << std::endl;

            std::cout << "rawHitCoords = [ ";
            for (int i = 0; i < rawHitCoords.GetNrows(); ++i) {
                std::cout << rawHitCoords[i] << " ,";
            }
            std::cout << "]" << std::endl;    
            std::cout << "Wire position error [cm]: " << sigma_wire << std::endl;
            std::cout << "Drift distance error [cm]: " << Rdrift_sigma << std::endl;
            std::cout << "Zreco error [cm]: " << zreco_sigma << std::endl;

            std::cout << "detID: " << det_idx << std::endl;
            std::cout << "hitID: " << hit_idx << std::endl;
            std::cout << "" << std::endl;
        }

    }

} // namespace GENFIT