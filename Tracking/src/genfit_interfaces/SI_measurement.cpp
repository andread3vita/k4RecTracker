#include "SI_measurement.hpp"
#include "TDecompSVD.h"


namespace GENFIT {

    SI_measurement::SI_measurement(const edm4hep::TrackerHitPlane& hit, const dd4hep::rec::SurfaceMap* surfaceMap, const int det_idx, const int hit_idx) {

        auto cellID0 = hit.getCellID();
        dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(cellID0);

        const dd4hep::rec::ISurface* surf  = sI->second;
        dd4hep::rec::Vector3D u = surf->u();
        dd4hep::rec::Vector3D v = surf->v();
        dd4hep::rec::Vector3D origin = surf->origin();
                        
        // convert 3D global position to 2D local position
        auto pos_hit = hit.getPosition();                                               // mm
        dd4hep::rec::Vector3D global_pos(pos_hit.x,pos_hit.y,pos_hit.z);                // mm
        dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos); // cm
                    
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

        // add plane
        TVector3 o(origin[0],origin[1],origin[2]); // cm
        TVector3 u_(u[0],u[1],u[2]);
        TVector3 v_(v[0],v[1],v[2]);

        // std::cout << "rawHitCoords: " << rawHitCoords[0] << "," << rawHitCoords[1] << std::endl;
        // std::cout << "rawHitCov: " << std::pow(dd4hep::mm * sigma_u, 2) << "," << std::pow(dd4hep::mm * sigma_v, 2) << std::endl;
        // std::cout << "detId: " << det_idx << std::endl;
        // std::cout << "hitID: " << hit_idx << std::endl;
        // std::cout << "O: " << origin[0] <<" , " << origin[1] <<" , " << origin[2] << std::endl;
        // std::cout << "U: " << u[0] <<" , " << u[1] <<" , " << u[2] << std::endl;
        // std::cout << "V: " << v[0] <<" , " << v[1] <<" , " << v[2] << std::endl;
        // std::cout << "PlaneID: " << cellID0 << std::endl;
        // std::cout << "" << std::endl;

        // double lengthU = surf->length_along_u();
        // double lengthV = surf->length_along_v();
        // std::cout << "lengthU: " << lengthU << std::endl;
        // std::cout << "lengthV: " << lengthV << std::endl;   
        // std::cout << "" << std::endl;
        // genfit::RectangularFinitePlane* pixelPlane = new genfit::RectangularFinitePlane(-lengthU/2.,lengthU/2.,-lengthV/2.,lengthV/2.);
        genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u_,v_));

        // plane->setFinitePlane(pixelPlane);
        genfitHit_->setPlane(plane, cellID0);

    }

    // Get implementation
    genfit::PlanarMeasurement* SI_measurement::getGenFit() const {
        return genfitHit_;
    }

} // namespace GENFIT
