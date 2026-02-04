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

#include "GenfitTrack.hpp"
#include "KalmanFitterInfo.h"
#include "AbsFitterInfo.h"
#include "EventDisplay.h"
#include <TMatrixDSymEigen.h>
#include <TDecompChol.h>

void printJacobian(double J[6][5]) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nJacobian J (6 x 5):\n\n";

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(14) << J[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void printCcart(const std::array<double,36>& C_cart)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nC_cart (6 x 6):\n\n";

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            std::cout << std::setw(14) << C_cart[i*6 + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void printChelix(const std::array<double,25>& C_helix)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nC_helix (5 x 5):\n\n";

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(14) << C_helix[i*5 + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

void helixToCartesianAndCov(

    TVector3 gen_position, TVector3 gen_momentum, int charge,
    const std::array<double,25>& C_helix,    
    double Bz,
    double Prx_cm, double Pry_cm, 
    std::array<double,36>& C_cart)     
{   

    double x0_PCA = gen_position.X();
    double y0_PCA = gen_position.Y();
    double pt = gen_momentum.Perp(); 
    double px = gen_momentum.X();
    double py = gen_momentum.Y();
    double pz = gen_momentum.Z();    

    double phi0 = gen_momentum.Phi();  
    double tanLambda = pz / pt;
    double omega = charge * Bz / pt;
    double d0 = - (Prx_cm - x0_PCA) * sin(phi0) + (Pry_cm - y0_PCA) * cos(phi0);

    // std::cout << "\nHelix parameters:\n";
    // std::cout << "d0: " << d0 << " mm\n";
    // std::cout << "phi0: " << phi0 << " rad\n";
    // std::cout << "omega: " << omega << " 1/mm\n";
    // std::cout << "z0: " << gen_position.Z() / dd4hep::mm << " mm\n";
    // std::cout << "tanLambda: " << tanLambda << "\n";

               
    double s = std::sin(phi0);
    double c = std::cos(phi0);

    // Jacobian J[row][col] = d(cartesian variable) / d(helix parameter)
    // rows:    0=x, 1=y, 2=z, 3=px, 4=py, 5=pz
    // columns: 0=d0, 1=phi0, 2=omega, 3=z0, 4=tanLambda

    double J[6][5] = {0};

    // =====================
    // x
    // =====================
    J[0][0] =  1.0 / s;  
    // dx / dd0

    J[0][1] = c / (s * s) * d0 - (y0_PCA - Pry_cm) / (c * c * std::tan(phi0) * std::tan(phi0));
    // dx / dphi0

    J[0][2] =  0.0;  
    // dx / domega

    J[0][3] =  0.0;  
    // dx / dz0

    J[0][4] =  0.0;  
    // dx / dtheta


    // =====================
    // y
    // =====================
    J[1][0] = - 1.0 / c;  
    // dy / dd0

    J[1][1] =  -s / (c * c) * d0 - (Prx_cm - x0_PCA) / (c * c);
    // dy / dphi0

    J[1][2] =  0.0;  
    // dy / domega

    J[1][3] =  0.0;  
    // dy / dz0

    J[1][4] =  0.0;  
    // dy / dtheta


    // =====================
    // z
    // =====================
    J[2][0] =  0.0;  
    // dz / dd0

    J[2][1] =  0.0;  
    // dz / dphi0

    J[2][2] =  0.0;  
    // dz / domega

    J[2][3] =  1.0;  
    // dz / dz0

    J[2][4] =  0.0;  
    // dz / dtheta


    // =====================
    // px = p * cos(phi0) * sin(theta)
    // =====================
    J[3][0] =  0.0;  
    // dpx / dd0

    J[3][1] = - px * std::tan(phi0);
    // dpx / dphi0

    J[3][2] =  - px / std::abs(omega);
    // dpx / domega

    J[3][3] =  0.0;  
    // dpx / dz0

    J[3][4] = 0;
    // dpx / dtheta


    // =====================
    // py = p * sin(phi0) * sin(theta)
    // =====================
    J[4][0] =  0.0;  
    // dpy / dd0

    J[4][1] =  - py / std::tan(phi0);
    // dpy / dphi0

    J[4][2] =  - py / std::abs(omega);
    // dpy / domega

    J[4][3] =  0.0;  
    // dpy / dz0

    J[4][4] = 0.0;
    // dpy / dtheta

    // =====================
    // pz = p * cos(theta)
    // =====================
    J[5][0] =  0.0;  
    // dpz / dd0

    J[5][1] =  0.0;  
    // dpz / dphi0

    J[5][2] =  - py / std::abs(omega);
    // dpz / domega

    J[5][3] =  0.0;  
    // dpz / dz0

    J[5][4] = pz / tanLambda;
    // dpz / dtheta

    // std::cout << omega << std::endl;
    // printJacobian(J);

    C_cart.fill(0.0);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double sum = 0.0;
            for (int k = 0; k < 5; ++k) {
                for (int l = 0; l < 5; ++l) {
                   
                    sum += J[i][k] * C_helix[k*5 + l] * J[j][l];
                }
            }
            C_cart[i*6 + j] = sum;
        }
    }

    // printChelix(C_helix);
    // printCcart(C_cart);
}

bool isPositiveSemiDefinite(const TMatrixDSym& M, double tol = 1e-10)
{
    TMatrixDSymEigen eig(M);
    TVectorD eigenValues = eig.GetEigenValues();

    for (int i = 0; i < eigenValues.GetNrows(); ++i) {
        if (eigenValues[i] < -tol)
            return false;
    }
    return true;
}

std::pair<TVector3, double> computeD0(TVector3 position, TVector3 momentum, int charge, TVector3 refPoint, double Bz)
{
    // -------------------------------
    // INPUTS (assumed to be defined)
    // -------------------------------
    // TVector3 position;   // point on the trajectory
    // TVector3 momentum;   // momentum at that point
    // int charge;          // particle charge (+1 or -1)
    // TVector3 refPoint;   // reference point (e.g. beam line)
    // double Bz;           // magnetic field along +z (Bz > 0)

    // -------------------------------
    // 1. Transverse momentum
    // -------------------------------
    double px = momentum.X();
    double py = momentum.Y();
    double pt = momentum.Perp();

    if (pt == 0.0) {
        throw std::runtime_error("Transverse momentum is zero");
    }

    // -------------------------------
    // 2. Radius of curvature
    // -------------------------------
    double R = pt / (std::abs(charge) * Bz);

    // -------------------------------
    // 3. Curvature direction
    // -------------------------------
    // For Bz > 0:
    //  - positive charge → counter-clockwise motion
    //  - negative charge → clockwise motion
    double sign = (charge > 0 ? +1.0 : -1.0);

    // Unit vector perpendicular to momentum in the XY plane
    TVector3 normal(
        -sign * py / pt,
        sign * px / pt,
        0.0
    );

    // -------------------------------
    // 4. Circle center in the XY plane
    // -------------------------------
    TVector3 center = position + R * normal;

    // Circle equation:
    // (x - center.X())^2 + (y - center.Y())^2 = R^2

    // -------------------------------
    // 5. Closest point on the circle
    //    to the reference point
    // -------------------------------
    TVector3 v = refPoint - center;

    // Distance in the transverse plane
    double vxy = std::sqrt(v.X()*v.X() + v.Y()*v.Y());

    if (vxy == 0.0) {
        throw std::runtime_error("Reference point coincides with circle center");
    }

    // Radial unit vector from center to refPoint
    TVector3 u(v.X()/vxy, v.Y()/vxy, 0.0);

    // Point on the circle closest to refPoint
    TVector3 closestPoint = center + R * u;

    // -------------------------------
    // 6. Tangent to the circle
    //    at the closest point
    // -------------------------------
    // Radius vector at the closest point
    TVector3 r = closestPoint - center;

    // Tangent is perpendicular to the radius
    // Direction follows the particle motion
    TVector3 tangent(
        -sign * r.Y(),
        sign * r.X(),
        0.0
    );

    tangent = tangent.Unit();

    // -------------------------------
    // 7. Angle between the tangent
    //    and the x-axis
    // -------------------------------
    // Oriented angle in the range [-pi, pi]
    double angleToX = std::atan2(tangent.Y(), tangent.X());

    // -------------------------------
    // Outputs:
    // -------------------------------
    // TVector3 center        -> circle center
    // double   R             -> radius
    // TVector3 closestPoint  -> point of closest approach to refPoint
    // TVector3 tangent       -> tangent direction at that point
    // double   angleToX      -> angle between tangent and x-axis (rad)

    return std::make_pair(closestPoint, angleToX);

}

namespace GenfitInterface {

    GenfitTrack::GenfitTrack(   const edm4hep::Track& track, const int dir, 
                                const int particle_hypothesis, 
                                const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder,
                                std::optional<int> maxHitForLoopers,
                                std::optional<TVector3> initial_position, std::optional<TVector3> initial_momentum)

        :   m_particle_hypothesis(particle_hypothesis), 
            m_direction(dir),

            m_posInit(0., 0., 0.), 
            m_momInit(0., 0., 0.), 
            m_genfitTrackRep(nullptr), 
            m_genfitTrack(nullptr), 
            m_edm4hepTrack(),

            m_dch_info(dch_info),
            m_dc_decoder(decoder)

            
    {   

        checkInitialization();
        init(track, dir, maxHitForLoopers, initial_position, initial_momentum);
    }

    GenfitTrack::~GenfitTrack() {}

    /**
    * @brief Check if required Genfit components are properly initialized.
    * 
    * This method verifies whether the singleton instances of `genfit::FieldManager`
    * and `genfit::MaterialEffects` have been initialized. These components are essential
    * for tracking operations in the Genfit framework. 
    * 
    * If either component is not initialized, an error message is printed to `std::cerr` 
    * and the program exits with `EXIT_FAILURE`.
    * 
    * @note This method should be called before performing any tracking-related operations 
    * that depend on magnetic field or material effects.
    */
    void GenfitTrack::checkInitialization() {

        if (!genfit::FieldManager::getInstance()->isInitialized()) {
            std::cerr << "Error: FieldManager is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
            std::cerr << "Error: MaterialEffects is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

    /**
    * @brief Initialize the GenfitTrack from a given input track.
    * 
    * This method initializes the internal state of the `GenfitTrack` object using
    * an input track of type `edm4hep::Track`. The following operations are performed:
    * - Initializes the internal `m_edm4hepTrack` with a new `MutableTrack`.
    * - Sorts the tracker hits of the input track by their Euclidean distance from the origin.
    * - Fills the internal track (`m_edm4hepTrack`) with the sorted hits.
    * - Extracts the positions of the first two hits and the last hit to initialize:
    *   - `m_posInit`: the initial position (first hit).
    *   - `m_momInit`: the initial direction (unit vector from first to second hit).
    *   - `m_firstHit_referencePoint`: copy of the first hit position.
    *   - `m_lastHit_referencePoint`: position of the last hit.
    *
    * @param track_init The input track containing tracker hits used for initialization.
    * @param initial_position Optional initial position for the track. If not provided, the first hit position is used.
    * @param initial_pT Optional initial pT for the track. If not provided, the direction from the first to the second hit is used.
    * @param maxHitForLoopers Optional maximum number of hits to consider for looper tracks. If not provided, all hits are used.
    *
    */
    void GenfitTrack::init(const edm4hep::Track& track_init, int dir, std::optional<int> maxHitForLoopers, std::optional<TVector3> initial_position /*cm*/, std::optional<TVector3> initial_momentum /*GeV/c*/) {

        // Initialize the m_edm4hepTrack
        m_edm4hepTrack = edm4hep::MutableTrack();

        auto hits_in_track = track_init.getTrackerHits();
        double z_max = -1e9; double x_withZmax = 0.0; double y_withZmax = 0.0;
        double z_min = 1e9; double x_withZmin = 0.0; double y_withZmin = 0.0;

        auto index = 0;
        for (auto hit : hits_in_track) {

            const auto pos_Hit = hit.getPosition();
            if (pos_Hit.z > z_max) {

                z_max = pos_Hit.z;
                x_withZmax = pos_Hit.x;
                y_withZmax = pos_Hit.y;
            }
            if (pos_Hit.z < z_min) {
                z_min = pos_Hit.z;
                x_withZmin = pos_Hit.x;
                y_withZmin = pos_Hit.y;

            }
            index++;

        }

        double first_x = 0;
        double first_y = 0;
        double first_z = 0;

        //Take x y z correspoding to the Z which is closer to zero
        if (std::abs(z_min) < std::abs(z_max)) {

            if (dir > 0) {
                first_x = x_withZmin;
                first_y = y_withZmin;
                first_z = z_min;
            }
            else {
                first_x = x_withZmax;
                first_y = y_withZmax;
                first_z = z_max;
            }

        } else {

            if (dir > 0) {
                first_x = x_withZmax;
                first_y = y_withZmax;
                first_z = z_max;
            }
            else {

                first_x = x_withZmin;
                first_y = y_withZmin;
                first_z = z_min;
            }

        }

        double delta_y = std::abs(y_withZmax - y_withZmin);
        double delta_z = std::abs(z_max - z_min);
        double cos_theta = delta_z / std::sqrt(std::pow(delta_y,2) + std::pow(delta_z,2));

        if ( std::abs(cos_theta) < 0.01 ) {

            double z_min_R = std::sqrt( x_withZmin * x_withZmin + y_withZmin * y_withZmin + z_min * z_min);
            double z_max_R = std::sqrt( x_withZmax * x_withZmax + y_withZmax * y_withZmax + z_max * z_max);

            if (z_min_R < z_max_R)
            {   
                if (dir > 0) {
                    first_x = x_withZmin;
                    first_y = y_withZmin;
                    first_z = z_min;
                }
                else
                {
                    first_x = x_withZmax;
                    first_y = y_withZmax;
                    first_z = z_max;
                }
           
            }
            else
            {   
                if (dir > 0) {
                    first_x = x_withZmax;
                    first_y = y_withZmax;
                    first_z = z_max;
                }
                else
                {
                    first_x = x_withZmin;
                    first_y = y_withZmin;
                    first_z = z_min;
                }
               
            }
            
        }

        std::vector<std::pair<float, int>> hitDistIndices{};
        index = 0;
        for (auto hit : hits_in_track) {

            const auto pos_siHit = hit.getPosition();
            const auto distance = std::sqrt(std::pow(pos_siHit.x - first_x, 2) + std::pow(pos_siHit.y - first_y, 2) + std::pow(pos_siHit.z - first_z, 2));

            hitDistIndices.emplace_back(distance, index++);

        }
        
        // Fill the internal track with sorted hits
        std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);

        int idx_fill_track = 0;
        int maxHit = maxHitForLoopers.value_or(hits_in_track.size());
        for (const auto& [_, idx] : hitDistIndices) {

            if (idx_fill_track >= maxHit)
            {
                break;
            }

            m_edm4hepTrack.addToTrackerHits(hits_in_track[idx]);
            idx_fill_track++;
        }

        // initialize track
        int index_loopHit = 0;
        TVector3 secondHit_referencePoint(0,0,0);

        auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
        for (auto hit : hits_for_genfit) {
            if (index_loopHit == 0) {
     
                auto position = hit.getPosition();
                m_firstHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z); // cm
            }
            if (index_loopHit == 1) {

                auto position = hit.getPosition();
                secondHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z); // cm
            
            }
            if (static_cast<size_t>(index_loopHit) == hits_for_genfit.size() - 1) {

                auto position = hit.getPosition();

                if (dir > 0)
                m_lastHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z); // cm

            }
            index_loopHit++;
        }

        m_posInit = initial_position.value_or( m_firstHit_referencePoint );                                     // cm
        m_momInit = initial_momentum.value_or((secondHit_referencePoint - m_firstHit_referencePoint).Unit());   // GeV/c
       
    }

    /**
    * @brief Creates a Genfit track object from the initialized state and tracker hits.
    * 
    * This method constructs the `genfit::Track` and associated `genfit::TrackRep` using the
    * initial position and momentum previously set during initialization. It assigns a default
    * covariance matrix (10 cm position resolution and 1 GeV momentum resolution), constructs
    * the state vector, and creates a `genfit::RKTrackRep` using the specified particle hypothesis.
    * 
    * Tracker hits from the internal `m_edm4hepTrack` are processed and converted to Genfit-compatible
    * measurements. Based on the hit type, the method creates either planar or wire measurements.
    * These measurements are inserted into the `genfit::Track` as `genfit::TrackPoint` objects.
    * 
    * The supported hit types and their corresponding detector type IDs are:
    * - `edm4hep::TrackerHitPlane`: Planar detector hit
    * - `edm4hep::SenseWireHit`: Drift chamber wire hit
    * 
    * If an unsupported hit type is encountered, the method prints an error message and exits.
    * 
    * @param debug_lvl The debug verbosity level to pass to the measurement constructors.
    * 
    * @note The method deletes any existing Genfit track or track representation before creating new ones.
    * 
    * @warning The method will terminate the program with `std::exit(EXIT_FAILURE)` if it encounters
    * an unknown hit type.
    */
    void GenfitTrack::createGenFitTrack(bool preFitInizialization, int debug_lvl) {

        delete m_genfitTrackRep;
        delete m_genfitTrack;

        TMatrixDSym covState(6);
        if ( preFitInizialization ) {

            // Reference point
            double Prx = m_firstHit_referencePoint.X(); // cm
            double Pry = m_firstHit_referencePoint.Y(); // cm
            double Prz = m_firstHit_referencePoint.Z(); // cm

            std::array<double,25> C_helix{};
            C_helix.fill(0.0);

            // columns: 0=d0, 1=phi0, 2=omega, 3=z0, 4=tanLambda

            C_helix[0*5 + 0] = 0.05 * 0.05;                         // d0 : 500 um = 0.05 cm
            C_helix[1*5 + 1] = 0.1 * 0.1;                           // phi0 : 0.1 rad
            
            double Bz = 2.;     
            double pt = m_momInit.Perp();    
            
            int charge = 0;
            if (std::abs(m_particle_hypothesis) == 11 || std::abs(m_particle_hypothesis) == 13)
                charge = (m_particle_hypothesis > 0) ? -1 : +1;
            else
                charge = (m_particle_hypothesis > 0) ? +1 : -1;

            double omega = charge * Bz / pt;                        // q * Bz / pt -> CHECK SIGN OF THE CHARGE

            C_helix[2*5 + 2] = std::pow(0.5 * omega, 2);            // omega : 0.5*omega 
            C_helix[3*5 + 3] = std::pow(0.1 * Prz, 2);              // z0 : 0.1*z0 
            C_helix[4*5 + 4] = 0.1 * 0.1;                           // tanLambda : 0.1  

	        
            std::array<double,36> C_cart{};
            helixToCartesianAndCov(
                m_posInit, m_momInit, charge,
                C_helix, Bz,
                Prx, Pry,
                C_cart
            );

            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    covState(i,j) = C_cart[i*6 + j];
                }
            }

        }
        else
        {

            // Default covariance: 10 cm position resolution, 1 GeV momentum resolution
            covState(0, 0) = 1000.0 * 10.0;  // sigma_x^2
            covState(1, 1) = 100.0 * 100.0;  // sigma_y^2
            covState(2, 2) = 100.0 * 100.0;  // sigma_z^2
            covState(3, 3) = 10.0 * 10.0;    // sigma_px^2
            covState(4, 4) = 10.0 * 10.0;    // sigma_py^2
            covState(5, 5) = 10.0 * 10.0;    // sigma_pz^2

        }

        

        if (!isPositiveSemiDefinite(covState)) {

            std::cerr << "Track Covariance Matrix is not positive definite!" << std::endl;
            std::exit(EXIT_FAILURE);

        
        }

        // Create stateVec
        TVectorD stateVec(6);
        
        // pos
        stateVec[0] = m_posInit.X();
        stateVec[1] = m_posInit.Y();
        stateVec[2] = m_posInit.Z();

        stateVec[3] = m_momInit.X();
        stateVec[4] = m_momInit.Y();
        stateVec[5] = m_momInit.Z();


        m_genfitTrackRep = new genfit::RKTrackRep(m_particle_hypothesis);
        m_genfitTrack = new genfit::Track(m_genfitTrackRep, stateVec, covState);

        auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
 
        int hit_idx(0);
        int detID(-1);
        for (auto hit : hits_for_genfit)
        {
            // The idea is that we can use the hit.getType() to identify the subdetector. 
            // Here is the convention:
            // 0: Vertex
            // 1: Drift Chamber
            // 2: Wrapper Barrel
            // 3: Wrapper Endcap
            // 5: Inner Barrel
            // 6: Outer Barrel
            // 7: Inner Endcap
            // 8: Outer Endcap

            
            auto cellID0 = hit.getCellID();
            detID = hit.getType();

            if (hit.isA<edm4hep::TrackerHitPlane>())
            {
               
                auto planar_hit =  hit.as<edm4hep::TrackerHitPlane>();
                GenfitInterface::PlanarMeasurement measurement = GenfitInterface::PlanarMeasurement(planar_hit, detID, ++hit_idx, debug_lvl);
                m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
                
            }
            else if (hit.isA<edm4hep::SenseWireHit>()) 
            {
                
                auto wire_hit =  hit.as<edm4hep::SenseWireHit>();
                GenfitInterface::WireMeasurement measurement = GenfitInterface::WireMeasurement(wire_hit, m_dch_info, m_dc_decoder, detID, ++hit_idx, debug_lvl);
                m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));

                       
            } 
            else 
            {

                std::cerr << "Error: No hits with cellID: " << cellID0 << std::endl;
                std::exit(EXIT_FAILURE);
             
            }

        }

    }

    /**
    * @brief Fit the Genfit track using the Deterministic Annealing Filter (DAF).
    * 
    * This function performs a forward and backward track fit using the Genfit DAF algorithm,
    * and populates the EDM4hep track states at the IP, first hit, and last hit positions.
    * 
    * @param Beta_init  Initial annealing parameter (temperature).
    * @param Beta_final Final annealing parameter.
    * @param Beta_steps Number of steps in the annealing schedule.
    * @param Bz         Value of z-component of the magnetic field at the center of the detector
    * @param debug_lvl  debug level: output if > 0
    * 
    * @return true if the fit was successful, false otherwise.
    */
    bool GenfitTrack::fit(int charge, double Beta_init = 100., double Beta_final=0.1, double Beta_steps=10, double Bz = 2., int debug_lvl = 0) {

        edm4hep::Track Track_temp = m_edm4hepTrack;
        for (size_t i = 0; i < Track_temp.trackStates_size(); ++i) {

            Track_temp.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* m_genfitFitter = new genfit::DAF(true, 1e-3,1e-3);
            m_genfitFitter->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);
            m_genfitFitter->setProbCut(1e-5);
            m_genfitFitter->setConvergenceDeltaWeight(1e-2); 

            int debug_lvl_fit = 0;
            if (debug_lvl > 1)
            {
                debug_lvl_fit = 2;
            }

           

            m_genfitFitter->setDebugLvl(debug_lvl_fit);

            debug_lvl_fit = 2;

            // Process track
            genfit::Track genfitTrack = *m_genfitTrack;
            genfit::AbsTrackRep* trackRep = genfitTrack.getTrackRep(0);
            m_genfitFitter->processTrackWithRep(&genfitTrack, trackRep);

            // Update edm4hep track state
            genfit::MeasuredStateOnPlane fittedState;
            TVector3 gen_position, gen_momentum;
            TMatrixDSym covariancePosMom(6);
            
            double x0_PCA;
            double y0_PCA;
            double z0_PCA; 

            double pz;  
            double pt;

            double d0;
            double z0;

            double phi;
            double omega;  
            double tanLambda; 
            
            if (m_genfitFitter->isTrackFitted(&genfitTrack, trackRep))
            {

            
                // trackState First Hit
                fittedState = genfitTrack.getFittedState();
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecFirstHit = fittedState.getState();
                
                edm4hep::TrackState trackStateFirstHit;
                x0_PCA = gen_position.X();     // cm
                y0_PCA = gen_position.Y();     // cm
                z0_PCA = gen_position.Z();     // cm
                pz = gen_momentum.Z();         // gev
                pt = gen_momentum.Perp();      // gev
            
                auto infoComputeD0 = computeD0(TVector3(x0_PCA, y0_PCA, z0_PCA), gen_momentum, charge, m_IP_referencePoint, Bz);
                
                d0 = ( - (m_IP_referencePoint.X() - infoComputeD0.first.X() )*sin(infoComputeD0.second) + (m_IP_referencePoint.Y() - infoComputeD0.first.Y())*cos(infoComputeD0.second) ) / dd4hep::mm; // mm
                z0 = ( infoComputeD0.first.Z() - m_IP_referencePoint.Z() ) / dd4hep::mm;   
                phi = gen_momentum.Phi();      // rad                                                                                                             // mm
                tanLambda = pz / pt;
                omega =  charge * Bz / abs(pt); // a.u. because I am not multiplying by c

                trackStateFirstHit.D0 = d0;
                trackStateFirstHit.Z0 = z0;
                trackStateFirstHit.phi = phi;
                trackStateFirstHit.omega = omega;
                trackStateFirstHit.tanLambda = tanLambda;
                // trackStateFirstHit.time = 0.;

                trackStateFirstHit.referencePoint = edm4hep::Vector3f(x0_PCA / dd4hep::mm, y0_PCA / dd4hep::mm, z0_PCA / dd4hep::mm);
                trackStateFirstHit.location = edm4hep::TrackState::AtFirstHit;

                if (debug_lvl_fit > 1)
                {
                    std::cout << "GenfitTrackFitter    DEBUG : TrackState at First Hit: " << std::endl;
                    std::cout << "  D0: " << trackStateFirstHit.D0 << " mm" << std::endl;
                    std::cout << "  Z0: " << trackStateFirstHit.Z0 << " mm" << std::endl;
                    std::cout << "  phi: " << trackStateFirstHit.phi << " rad" << std::endl;
                    std::cout << "  omega: " << trackStateFirstHit.omega << " a.u." << std::endl;
                    std::cout << "  tanLambda: " << trackStateFirstHit.tanLambda << std::endl;
                    std::cout << "  location: " << trackStateFirstHit.location << std::endl;
                    std::cout << "  reference point: (" << trackStateFirstHit.referencePoint.x << ", " << trackStateFirstHit.referencePoint.y << ", " << trackStateFirstHit.referencePoint.z << ") mm" << std::endl;
                    std::cout << "  PCA point: (" << infoComputeD0.first.X() << ", " << infoComputeD0.first.Y() << ", " << infoComputeD0.first.Z() << ") mm" << std::endl;
                }

                // trackState lastHit
                fittedState = genfitTrack.getFittedState(genfitTrack.getNumPoints()-1);
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecLastHit = fittedState.getState();
                
                edm4hep::TrackState trackStateLastHit;
                x0_PCA = gen_position.X();     // cm
                y0_PCA = gen_position.Y();     // cm
                z0_PCA = gen_position.Z();     // cm
                pz = gen_momentum.Z();         // gev
                pt = gen_momentum.Perp();      // gev

                auto infoComputeD0_lastHit = computeD0(TVector3(x0_PCA, y0_PCA, z0_PCA), gen_momentum, charge, m_IP_referencePoint, Bz);
                d0 = ( - (m_IP_referencePoint.X() - infoComputeD0_lastHit.first.X() )*sin(infoComputeD0_lastHit.second) + (m_IP_referencePoint.Y() - infoComputeD0_lastHit.first.Y())*cos(infoComputeD0_lastHit.second) ) / dd4hep::mm; // mm
                z0 = ( infoComputeD0_lastHit.first.Z() - m_IP_referencePoint.Z() ) / dd4hep::mm;                                                                                                                                        // mm
                phi = gen_momentum.Phi();      // rad
                tanLambda = pz / pt;    
                omega =  charge * Bz / abs(pt); // a.u. because I am not multiplying by c

                trackStateLastHit.D0 = d0;
                trackStateLastHit.Z0 = z0;
                trackStateLastHit.phi = phi;
                trackStateLastHit.omega = omega;
                trackStateLastHit.tanLambda = tanLambda;
                // trackStateLastHit.time = 0.;
                
                trackStateLastHit.referencePoint = edm4hep::Vector3f(x0_PCA / dd4hep::mm, y0_PCA / dd4hep::mm, z0_PCA / dd4hep::mm);
                trackStateLastHit.location = edm4hep::TrackState::AtLastHit;
                
                if (debug_lvl_fit > 1)
                {
                    std::cout << "GenfitTrackFitter    DEBUG : TrackState at Last Hit: " << std::endl;
                    std::cout << "  D0: " << trackStateLastHit.D0 << " mm" << std::endl;
                    std::cout << "  Z0: " << trackStateLastHit.Z0 << " mm" << std::endl;
                    std::cout << "  phi: " << trackStateLastHit.phi << " rad" << std::endl;
                    std::cout << "  omega: " << trackStateLastHit.omega << " a.u." << std::endl;
                    std::cout << "  tanLambda: " << trackStateLastHit.tanLambda << std::endl;
                    std::cout << "  location: " << trackStateLastHit.location << std::endl;
                    std::cout << "  reference point: (" << trackStateLastHit.referencePoint.x << ", " << trackStateLastHit.referencePoint.y << ", " << trackStateLastHit.referencePoint.z << ") mm" << std::endl;
                    std::cout << "  PCA point: (" << infoComputeD0_lastHit.first.X() << ", " << infoComputeD0_lastHit.first.Y() << ", " << infoComputeD0_lastHit.first.Z() << ") mm" << std::endl;
                }

                if (m_direction < 0)
                {
                    trackRep->setPropDir(-1);
                    
                }
                
                //take first fitted point
                genfit::TrackPoint* tp = genfitTrack.getPointWithFitterInfo(0);
                if (tp == NULL) {std::cout << "Track has no TrackPoint with fitterInfo!(but fitstatus ok?)"<<std::endl;}
                auto* fi = static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(trackRep));

                //extrapolate rep to target plane
                try{

                    fittedState=fi->getFittedState(true);
                    trackRep->extrapolateToLine(fittedState,TVector3(0,0,0),TVector3(0,0,1));

                    fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                    auto stateVecIP = fittedState.getState();

                    gen_momentum.SetX(-gen_momentum.X());
                    gen_momentum.SetY(-gen_momentum.Y());
                    gen_momentum.SetZ(-gen_momentum.Z());
                    
                    edm4hep::TrackState trackStateIP;
                    x0_PCA = gen_position.X();     // cm
                    y0_PCA = gen_position.Y();     // cm
                    z0_PCA = gen_position.Z();     // cm
                    pz = gen_momentum.Z();         // gev
                    pt = gen_momentum.Perp();      // gev

                    auto infoComputeD0_IP = computeD0(TVector3(x0_PCA, y0_PCA, z0_PCA), gen_momentum, charge, m_IP_referencePoint, Bz);
                    d0 = ( - (m_IP_referencePoint.X() - infoComputeD0_IP.first.X() )*sin(infoComputeD0_IP.second) + (m_IP_referencePoint.Y() - infoComputeD0_IP.first.Y())*cos(infoComputeD0_IP.second) ) / dd4hep::mm; // mm
                    z0 = ( infoComputeD0_IP.first.Z() - m_IP_referencePoint.Z() ) / dd4hep::mm;                                                                                                                         // mm
                    phi = gen_momentum.Phi();      // rad
                    tanLambda = pz / pt;    
                    omega =  charge * Bz / abs(pt); // a.u. because I am not multiplying by c   

                    trackStateIP.D0 = d0;
                    trackStateIP.Z0 = z0;
                    trackStateIP.phi = phi;
                    trackStateIP.omega = omega;
                    trackStateIP.tanLambda = tanLambda;
                    // trackStateIP.time = 0.;
                    
                    trackStateIP.referencePoint = edm4hep::Vector3f(x0_PCA / dd4hep::mm, y0_PCA / dd4hep::mm, z0_PCA / dd4hep::mm);
                    trackStateIP.location = edm4hep::TrackState::AtIP;

                    if (debug_lvl_fit > 1)
                    {
                        std::cout << "GenfitTrackFitter    DEBUG : TrackState at IP: " << std::endl;
                        std::cout << "  D0: " << trackStateIP.D0 << " mm" << std::endl;
                        std::cout << "  Z0: " << trackStateIP.Z0 << " mm" << std::endl;
                        std::cout << "  phi: " << trackStateIP.phi << " rad" << std::endl;
                        std::cout << "  omega: " << trackStateIP.omega << " a.u." << std::endl;
                        std::cout << "  tanLambda: " << trackStateIP.tanLambda << std::endl;
                        std::cout << "  location: " << trackStateIP.location << std::endl;
                        std::cout << "  reference point: (" << trackStateIP.referencePoint.x << ", " << trackStateIP.referencePoint.y << ", " << trackStateIP.referencePoint.z << ") mm" << std::endl;
                        std::cout << "  PCA point: (" << infoComputeD0_IP.first.X() << ", " << infoComputeD0_IP.first.Y() << ", " << infoComputeD0_IP.first.Z() << ") mm" << std::endl;
                    }

                    m_edm4hepTrack.addToTrackStates(trackStateIP);
                    m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
                    m_edm4hepTrack.addToTrackStates(trackStateLastHit);

                }catch(...){

                    return false;
                }

                if (m_genfitFitter->isTrackFitted(&genfitTrack, trackRep))
                {
                    
                    m_edm4hepTrack.setChi2(genfitTrack.getFitStatus()->getChi2());
                    m_edm4hepTrack.setNdf(genfitTrack.getFitStatus()->getNdf());

                }
                else
                {

                    m_edm4hepTrack.setChi2(-1);
                    m_edm4hepTrack.setNdf(-1);
                    return false;

                }

            }
            else
            {

                m_edm4hepTrack.setChi2(-1);
                m_edm4hepTrack.setNdf(-1);
                return false;

            }

            return m_genfitFitter->isTrackFitted(&genfitTrack, trackRep);
            
        }
        catch(...)
        {
            return false;
        }
 
    }
}
