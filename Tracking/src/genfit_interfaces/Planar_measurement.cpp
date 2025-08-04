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

#include "Planar_measurement.hpp"

namespace GenfitInterface {

    Planar_measurement::Planar_measurement(const edm4hep::TrackerHitPlane& hit, const int det_idx, const int hit_idx, const int debug_lvl=0) {

               
        // Create a 2D local measurement from a 3D global measurement
        auto pos_hit = hit.getPosition();                                                                           // mm
        dd4hep::rec::Vector3D global_pos(dd4hep::mm * pos_hit.x,dd4hep::mm * pos_hit.y,dd4hep::mm * pos_hit.z);     // cm
        dd4hep::rec::Vector2D local_pos(0.,0.);                                                                     // cm
        TVector3 Origin(global_pos[0],global_pos[1],global_pos[2]);

        TVectorD rawHitCoords(2);
        rawHitCoords[0] = local_pos[0]; //cm
        rawHitCoords[1] = local_pos[1]; //cm


        // Create the 2x2 Covariance Matrix
        auto U_theta = hit.getU()[0];
        auto U_phi = hit.getU()[1];
        TVector3 U(std::sin(U_theta) * std::cos(U_phi),std::sin(U_theta) * std::sin(U_phi),std::cos(U_theta)); U = U.Unit();

        auto V_theta = hit.getV()[0];
        auto V_phi = hit.getV()[1];
        TVector3 V(std::sin(V_theta) * std::cos(V_phi),std::sin(V_theta) * std::sin(V_phi),std::cos(V_theta)); V = V.Unit();

        TMatrixDSym rawHitCov(2);
        double sigma_u = hit.getDu();                       // mm
        double sigma_v = hit.getDv();                       // mm 
        rawHitCov(0,0) = std::pow(dd4hep::mm * sigma_u, 2); // cm^2
        rawHitCov(0,1) = 0;
        rawHitCov(1,0) = 0;
        rawHitCov(1,1) = std::pow(dd4hep::mm * sigma_v, 2); // cm^2
        
        
        // Create genfit::PlanarMeasurement
        genfitHit_ = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, nullptr);

        // Add plane
        auto cellID0 = hit.getCellID();
        genfit::SharedPlanePtr plane(new genfit::DetPlane(Origin,U,V));
        genfitHit_->setPlane(plane, cellID0);

        if (debug_lvl > 0) {
            std::cout << "Planar measurement created with the following parameters:" << std::endl;

            std::cout << "rawHitCoords: " << rawHitCoords[0] << "," << rawHitCoords[1] << std::endl;
            std::cout << "rawHitCov: " << std::pow(dd4hep::mm * sigma_u, 2) << "," << std::pow(dd4hep::mm * sigma_v, 2) << std::endl;
            std::cout << "detId: " << det_idx << std::endl;
            std::cout << "hitID: " << hit_idx << std::endl;
            std::cout << "O: " << Origin[0] <<" , " << Origin[1] <<" , " << Origin[2] << std::endl;
            std::cout << "U: " << U[0] <<" , " << U[1] <<" , " << U[2] << std::endl;
            std::cout << "V: " << V[0] <<" , " << V[1] <<" , " << V[2] << std::endl;
            std::cout << "PlaneID: " << cellID0 << std::endl;
            std::cout << "" << std::endl;
        }

    }

} // namespace GENFIT
