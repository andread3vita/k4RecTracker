#include "Planar_measurement.hpp"
#include "TDecompSVD.h"


namespace GENFIT {

    Planar_measurement::Planar_measurement(const edm4hep::TrackerHitPlane& hit, const int det_idx, const int hit_idx, const int debug_lvl=0) {

        auto cellID0 = hit.getCellID();
        
                        
        // convert 3D global position to 2D local position
        auto pos_hit = hit.getPosition();                                                                           // mm
        dd4hep::rec::Vector3D global_pos(dd4hep::mm * pos_hit.x,dd4hep::mm * pos_hit.y,dd4hep::mm * pos_hit.z);     // cm
        dd4hep::rec::Vector2D local_pos(0.,0.);                                                                     // cm
        dd4hep::rec::Vector3D origin = global_pos;

        auto U_theta = hit.getU()[0];
        auto U_phi = hit.getU()[1];
        dd4hep::rec::Vector3D u(std::sin(U_theta) * std::cos(U_phi),std::sin(U_theta) * std::sin(U_phi),std::cos(U_theta));

        auto V_theta = hit.getV()[0];
        auto V_phi = hit.getV()[1];
        dd4hep::rec::Vector3D v(std::sin(V_theta) * std::cos(V_phi),std::sin(V_theta) * std::sin(V_phi),std::cos(V_theta));
                    
        TVectorD rawHitCoords(2);
        rawHitCoords[0] = local_pos[0]; //cm
        rawHitCoords[1] = local_pos[1]; //cm

        TMatrixDSym rawHitCov(2);
        double sigma_u = hit.getDu();                       // mm
        double sigma_v = hit.getDv();                       // mm 
        rawHitCov(0,0) = std::pow(dd4hep::mm * sigma_u, 2); // cm^2
        rawHitCov(0,1) = 0;
        rawHitCov(1,0) = 0;
        rawHitCov(1,1) = std::pow(dd4hep::mm * sigma_v, 2); // cm^2
        
        
        // create measurement 
        genfitHit_ = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, nullptr);

        if (debug_lvl > 0) {
            std::cout << "Planar measurement created with the following parameters:" << std::endl;

            std::cout << "rawHitCoords: " << rawHitCoords[0] << "," << rawHitCoords[1] << std::endl;
            std::cout << "rawHitCov: " << std::pow(dd4hep::mm * sigma_u, 2) << "," << std::pow(dd4hep::mm * sigma_v, 2) << std::endl;
            std::cout << "detId: " << det_idx << std::endl;
            std::cout << "hitID: " << hit_idx << std::endl;
            std::cout << "O: " << origin[0] <<" , " << origin[1] <<" , " << origin[2] << std::endl;
            std::cout << "U: " << u[0] <<" , " << u[1] <<" , " << u[2] << std::endl;
            std::cout << "V: " << v[0] <<" , " << v[1] <<" , " << v[2] << std::endl;
            std::cout << "PlaneID: " << cellID0 << std::endl;
            std::cout << "" << std::endl;
        }
        

        // add plane
        TVector3 o(origin[0],origin[1],origin[2]); // cm
        TVector3 u_(u[0],u[1],u[2]); u_ = u_.Unit();
        TVector3 v_(v[0],v[1],v[2]); v_ = v_.Unit();
        genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u_,v_));
        genfitHit_->setPlane(plane, cellID0);

    }

    // Get implementation
    genfit::PlanarMeasurement* Planar_measurement::getGenFit() const {
        return genfitHit_;
    }

} // namespace GENFIT
