#include "VTX_measurement.hpp"
#include "TDecompSVD.h"

namespace IDEAtracking {

    VTX_measurement::VTX_measurement(const edm4hep::TrackerHitPlane& hit, const dd4hep::rec::SurfaceMap* surfaceMap, const int det_idx, const int hit_idx) {

        auto cellID0 = hit.getCellID();
        dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(cellID0);
        const dd4hep::rec::ISurface* surf  = sI->second;
        dd4hep::rec::Vector3D u = surf->u();
        dd4hep::rec::Vector3D v = surf->v();
        dd4hep::rec::Vector3D origin = surf->origin();
        dd4hep::rec::Vector3D normal = surf->normal();
        
        // convert 3D global position to 2D local position
        auto pos = hit.getPosition();
        dd4hep::rec::Vector3D global_pos(pos.x,pos.y,pos.z);
        dd4hep::rec::Vector2D local_pos = surf->globalToLocal(dd4hep::mm * global_pos);
        
        TVectorD rawHitCoords(2);
        rawHitCoords[0] = local_pos[0];
        rawHitCoords[1] = local_pos[1];

        TMatrixD R(3,3);
        for (int i = 0; i < 3; i++) {
            R(i, 0) = u[i];
            R(i, 1) = v[i];
            R(i, 2) = normal[i];
        }

        // Convert 3x3 global covariance matrix to 2x2 local covariance martrix
        auto cov_matrix = hit.getCovMatrix();
        auto values_cov = cov_matrix.values;

        TMatrixDSym sigma(3);
        sigma(0, 0) = values_cov[0];  // sigma_{xx}
        sigma(0, 1) = values_cov[1];  // sigma_{xy}
        sigma(0, 2) = values_cov[2];  // sigma_{xz}
        sigma(1, 1) = values_cov[3];  // sigma_{yy}
        sigma(1, 2) = values_cov[4];  // sigma_{yz}
        sigma(2, 2) = values_cov[5];  // sigma_{zz}

        TMatrixD Rt(TMatrixD::kTransposed, R);
        TMatrixD Sigma_prime = Rt * sigma * R;
        
        TMatrixD Sigma_local(2,2);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                Sigma_local(i, j) = Sigma_prime(i, j);
            }
        }

        TMatrixDSym rawHitCov(2);
        rawHitCov(0,0) = Sigma_local(0,0);
        rawHitCov(0,1) = Sigma_local(0,1);
        rawHitCov(1,0) = Sigma_local(1,0);
        rawHitCov(1,1) = Sigma_local(1,1);

        rawHitCov(0,0) = 1e-6;
        rawHitCov(0,1) = 0;
        rawHitCov(1,0) = 0;
        rawHitCov(1,1) = 1e-6;


        // create measurement 
        genfit::TrackPoint* trackPoint = new genfit::TrackPoint(nullptr);
        genfitHit_ = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, trackPoint);

        // add plane
        TVector3 o(origin[0],origin[1],origin[2]);
        TVector3 u_(u[0],u[1],u[2]);
        TVector3 v_(v[0],v[1],v[2]);
        
        genfitHit_->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(o,u_,v_)), cellID0);

    }

    // Get implementation
    genfit::PlanarMeasurement* VTX_measurement::getGenFit() const {
        return genfitHit_;
    }

} // namespace IDEtracking
