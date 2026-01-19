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
    double Prx_mm, double Pry_mm, 
    std::array<double,36>& C_cart)     
{   

    double x0_PCA = gen_position.X() / dd4hep::mm;
    double y0_PCA = gen_position.Y() / dd4hep::mm;
    double pt = gen_momentum.Perp(); 
    double px = gen_momentum.X();
    double py = gen_momentum.Y();
    double pz = gen_momentum.Z();    

    double phi0 = gen_momentum.Phi();  
    double tanLambda = pz / pt;
    double omega = charge * Bz / pt;
    double d0 = - (Prx_mm - x0_PCA) * sin(phi0) + (Pry_mm - y0_PCA) * cos(phi0);

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

    J[0][1] = c / (s * s) * d0 - (y0_PCA - Pry_mm) / (c * c * std::tan(phi0) * std::tan(phi0));
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

    J[1][1] =  -s / (c * c) * d0 - (Prx_mm - x0_PCA) / (c * c);
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

namespace GenfitInterface {

    GenfitTrack::GenfitTrack(const edm4hep::Track& track,
                             const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder, 
                             const int particle_hypothesis, std::optional<int> maxHitForLoopers, std::optional<int> maxHitForLoopers_back, std::optional<TVector3> initial_position, std::optional<TVector3> initial_momentum, std::optional<TVector3> initial_position_back, std::optional<TVector3> initial_momentum_back)

        :   m_particle_hypothesis(particle_hypothesis), 
            m_posInit(0., 0., 0.), 
            m_momInit(0., 0., 0.), 
            m_genfitTrackRep(nullptr), 
            m_genfitTrack(nullptr), 
            m_edm4hepTrack(), 
            m_dch_info(dch_info),
            m_dc_decoder(decoder),

            m_posInit_back(0., 0., 0.),
            m_momInit_back(0., 0., 0.),
            m_genfitTrackRep_back(nullptr),
            m_genfitTrack_back(nullptr),
            m_edm4hepTrack_back(),
            m_dch_info_back(dch_info),
            m_dc_decoder_back(decoder)

            
    {   

        checkInitialization();
        init(track, maxHitForLoopers, maxHitForLoopers_back, initial_position, initial_momentum, initial_position_back, initial_momentum_back);

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
    void GenfitTrack::init(const edm4hep::Track& track_init, std::optional<int> maxHitForLoopers, std::optional<int> maxHitForLoopers_back, std::optional<TVector3> initial_position, std::optional<TVector3> initial_momentum, std::optional<TVector3> initial_position_back, std::optional<TVector3> initial_momentum_back) {

        // Initialize the m_edm4hepTrack
        m_edm4hepTrack = edm4hep::MutableTrack();

        // // Sort the hits by distance from the origin
        // std::vector<std::pair<float, int>> hitDistIndices{};
        // int index = 0;

        // auto hits_in_track = track_init.getTrackerHits();
        // for (auto hit : hits_in_track) {

        //     const auto pos_siHit = hit.getPosition();
        //     const auto distance = std::sqrt(pos_siHit.x * pos_siHit.x + pos_siHit.y * pos_siHit.y + pos_siHit.z * pos_siHit.z);
        //     hitDistIndices.emplace_back(distance, index++);

        // }
        
        // // Fill the internal track with sorted hits
        // std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);
        // int idx_fill_track = 0;
        // int maxHit = maxHitForLoopers.value_or(hits_in_track.size());
        // for (const auto& [_, idx] : hitDistIndices) {

        //     if (idx_fill_track >= maxHit)
        //     {
        //         break;
        //     }

        //     m_edm4hepTrack.addToTrackerHits(hits_in_track[idx]);
        //     idx_fill_track++;
        // }

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

            first_x = x_withZmin;
            first_y = y_withZmin;
            first_z = z_min;

        } else {

            first_x = x_withZmax;
            first_y = y_withZmax;
            first_z = z_max;

        }

        double delta_y = std::abs(y_withZmax - y_withZmin);
        double delta_z = std::abs(z_max - z_min);
        double cos_theta = delta_z / std::sqrt(std::pow(delta_y,2) + std::pow(delta_z,2));

        if ( std::abs(cos_theta) < 0.01 ) {

            double z_min_R = std::sqrt( x_withZmin * x_withZmin + y_withZmin * y_withZmin + z_min * z_min);
            double z_max_R = std::sqrt( x_withZmax * x_withZmax + y_withZmax * y_withZmax + z_max * z_max);

            if (z_min_R < z_max_R)
            {
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

        std::vector<std::pair<float, int>> hitDistIndices{};
        index = 0;
        for (auto hit : hits_in_track) {

            const auto pos_siHit = hit.getPosition();
            const auto distance = std::sqrt(std::pow(pos_siHit.x - first_x, 2) + std::pow(pos_siHit.y - first_y, 2) + std::pow(pos_siHit.z - first_z, 2));

            std::cout << "Hit " << index << ": pos = " << pos_siHit.x << ", " << pos_siHit.y << ", " << pos_siHit.z << " | distance to first hit = " << distance << std::endl;
            hitDistIndices.emplace_back(distance, index++);

        }
        
        // Fill the internal track with sorted hits
        std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);
        // std::ranges::reverse(hitDistIndices);

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
        TVector3 m_secondHit_referencePoint(0,0,0);

        auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
        for (auto hit : hits_for_genfit) {
            if (index_loopHit == 0) {
     
                auto position = hit.getPosition();
                m_firstHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (index_loopHit == 1) {

                auto position = hit.getPosition();
                m_secondHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (static_cast<size_t>(index_loopHit) == hits_for_genfit.size() - 1) {

                auto position = hit.getPosition();
                m_lastHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);

            }
            index_loopHit++;
        }

        m_posInit = initial_position.value_or( m_firstHit_referencePoint );
        m_momInit = initial_momentum.value_or((m_secondHit_referencePoint - m_firstHit_referencePoint).Unit());

        // Do the same for the backward track
        m_edm4hepTrack_back = edm4hep::MutableTrack();

        if (first_x == x_withZmin)
        {
            first_x = x_withZmax;
            first_y = y_withZmax;
            first_z = z_max;
        }
        else if (first_x == x_withZmax)
        {
            first_x = x_withZmin;
            first_y = y_withZmin;
            first_z = z_min;
        }
        else
        {
            std::exit(EXIT_FAILURE);
        }

        if ( std::abs(cos_theta) < 0.01 ) {

            double z_min_R = std::sqrt( x_withZmin * x_withZmin + y_withZmin * y_withZmin + z_min * z_min);
            double z_max_R = std::sqrt( x_withZmax * x_withZmax + y_withZmax * y_withZmax + z_max * z_max);

            if (z_min_R > z_max_R)
            {
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

        std::vector<std::pair<float, int>> hitDistIndices_back{};
        index = 0;
        for (auto hit : hits_in_track) {

            const auto pos_siHit = hit.getPosition();
            const auto distance = std::sqrt(std::pow(pos_siHit.x - first_x, 2) + std::pow(pos_siHit.y - first_y, 2) + std::pow(pos_siHit.z - first_z, 2));
            hitDistIndices_back.emplace_back(distance, index++);

        }
        
        // Fill the internal track with sorted hits
        std::ranges::sort(hitDistIndices_back, {}, &std::pair<float, int>::first);
        // std::ranges::reverse(hitDistIndices);

        idx_fill_track = 0;
        maxHit = maxHitForLoopers_back.value_or(hits_in_track.size());
        for (const auto& [_, idx] : hitDistIndices_back) {

            if (idx_fill_track >= maxHit)
            {
                break;
            }

            m_edm4hepTrack_back.addToTrackerHits(hits_in_track[idx]);
            idx_fill_track++;
        }

        // initialize track
        index_loopHit = 0;
        TVector3 m_secondHit_referencePoint_back(0,0,0);

        hits_for_genfit = m_edm4hepTrack_back.getTrackerHits();
        for (auto hit : hits_for_genfit) {
            if (index_loopHit == 0) {
     
                auto position = hit.getPosition();
                m_firstHit_referencePoint_back = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);   // cm
            }
            if (index_loopHit == 1) {

                auto position = hit.getPosition();
                m_secondHit_referencePoint_back = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);  // cm
            }
            if (static_cast<size_t>(index_loopHit) == hits_for_genfit.size() - 1) {

                auto position = hit.getPosition();
                m_lastHit_referencePoint_back = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);    // cm

            }
            index_loopHit++;
        }

        m_posInit_back = initial_position_back.value_or( m_firstHit_referencePoint_back );
        m_momInit_back = initial_momentum_back.value_or((m_secondHit_referencePoint_back - m_firstHit_referencePoint_back).Unit());

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
    void GenfitTrack::createGenFitTrack(int dir, int debug_lvl, std::optional<double> pz_initial, std::optional<double> pT_initial) {

        if (dir > 0)
        {
            delete m_genfitTrackRep;
            delete m_genfitTrack;
        }
        else
        {
            delete m_genfitTrackRep_back;
            delete m_genfitTrack_back;
        }
        

        TMatrixDSym covState(6);
        if ( pT_initial.has_value() && pz_initial.has_value()) {

            // if (m_momInit.Mag() < 1.0) {
                
            //     TVector3 dir_init = m_momInit.Unit();
            //     pT_initial = dir_init.Perp();
            //     pz_initial = dir_init.Z();

            // }


	        // m_momInit.Print();

            // Reference point
            double Prx = m_firstHit_referencePoint.X() / dd4hep::mm;
            double Pry = m_firstHit_referencePoint.Y() / dd4hep::mm;
            double Prz = m_firstHit_referencePoint.Z() / dd4hep::mm;

            // Magnetic field constants
           

            std::array<double,25> C_helix{};
            C_helix.fill(0.0);

            // columns: 0=d0, 1=phi0, 2=omega, 3=z0, 4=tanLambda

            C_helix[0*5 + 0] = 0.05 * 0.05;                         // d0 : 500 um = 0.05 cm
            C_helix[1*5 + 1] = 0.1 * 0.1;                           // phi0 : 0.1 rad
            
            double Bz = 2.;     
            double pt = m_momInit.Perp();    
            int charge = getHypotesisCharge(m_particle_hypothesis);  
            double omega = charge * Bz / pt;                        // -> 1/mm 1/(mm**2 * s)

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
                    covState(i,j) = C_cart[i*6 + j]; // * 1000;
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
            // std::exit(EXIT_FAILURE);
	    
	        // covState(0, 0) = 1000.0 * 10.0;  // sigma_x^2
            // covState(1, 1) = 100.0 * 100.0;  // sigma_y^2
            // covState(2, 2) = 100.0 * 100.0;  // sigma_z^2
            // covState(3, 3) = 10.0 * 10.0;    // sigma_px^2
            // covState(4, 4) = 10.0 * 10.0;    // sigma_py^2
            // covState(5, 5) = 10.0 * 10.0;    // sigma_pz^2

        
        }

	    // covState.Print();

        // Create stateVec
        TVectorD stateVec(6);
        
        // pos
        stateVec[0] = m_posInit.X();
        stateVec[1] = m_posInit.Y();
        stateVec[2] = m_posInit.Z();

        // mom
        if (pT_initial.has_value() && pz_initial.has_value()) {

            // if (m_momInit.Mag() < 1.0) {

            //     TVector3 dir_init = m_momInit.Unit();
            //     pT_initial = dir_init.Perp();
            //     pz_initial = dir_init.Z();

            //     stateVec[3] = dir_init.X();
            //     stateVec[4] = dir_init.Y();
            //     stateVec[5] = dir_init.Z();

            //     m_momInit = dir_init;

            // }
            // else
            // {

                stateVec[3] = m_momInit.X();
                stateVec[4] = m_momInit.Y();
                stateVec[5] = m_momInit.Z();

            // }

        }
        else {

            stateVec[3] = m_momInit.X();
            stateVec[4] = m_momInit.Y();
            stateVec[5] = m_momInit.Z();

        }

        // stateVec.Print();

        if (dir > 0)
        {
            m_genfitTrackRep = new genfit::RKTrackRep(m_particle_hypothesis);
            m_genfitTrack = new genfit::Track(m_genfitTrackRep, stateVec, covState);
        }
        else
        {
            m_genfitTrackRep_back = new genfit::RKTrackRep(m_particle_hypothesis);
            m_genfitTrack_back = new genfit::Track(m_genfitTrackRep_back, stateVec, covState);
        }

        auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
        if (dir > 0)
        {
            hits_for_genfit = m_edm4hepTrack.getTrackerHits();
        }
        else
        {
            hits_for_genfit = m_edm4hepTrack_back.getTrackerHits();
        }
 
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

                if (dir > 0)
                {
                    m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
                }
                else
                {
                    m_genfitTrack_back->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack_back));
                }
                
            }
            else if (hit.isA<edm4hep::SenseWireHit>()) 
            {
                
                auto wire_hit =  hit.as<edm4hep::SenseWireHit>();


                // int cellid = wire_hit.getCellID();

                // float wireStereoAngle = wire_hit.getWireStereoAngle();
                // float wireAzimuthalAngle = wire_hit.getWireAzimuthalAngle();
                
                // // Wire hit
                // auto pos = wire_hit.getPosition();
                // TVector3 position(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm);   // cm

                // // Wire direction
                // TVector3 direction(0,0,1);
                // direction.RotateX(wireStereoAngle);
                // direction.RotateZ(wireAzimuthalAngle);

                // double Zplane_right = 200.0;       //  -20.0
                // double Zplane_left = -200.0;       //  20.0

                // if (direction.Z() == 0) {
                //     std::cout << "Direction is orthogonal to Z" << std::endl;
                // } else {

                //     double t = (Zplane_right - position.Z()) / direction.Z();
                //     TVector3 Pint = position + t * direction;
                //     std::cout << "Intersection at 200 cm: "
                //             << Pint.X() << " "
                //             << Pint.Y() << " "
                //             << Pint.Z() << std::endl;

                //     t = (Zplane_left - position.Z()) / direction.Z();
                //     Pint = position + t * direction;
                //     std::cout << "Intersection at -200 cm: "
                //             << Pint.X() << " "
                //             << Pint.Y() << " "
                //             << Pint.Z() << std::endl;
                // }

                // // Compute wire extremities
                // int ilayer          = m_dch_info->CalculateILayerFromCellIDFields(m_dc_decoder->get(cellid, "layer"), m_dc_decoder->get(cellid, "superlayer"));
                // int nphi            = m_dc_decoder->get(cellid, "nphi");
                // auto& l             = m_dch_info->database.at(ilayer);
                // int stereosign      = l.StereoSign();
                // double rz0          = l.radius_sw_z0;
                // double dphi         = m_dch_info->twist_angle;
                // double kappa        = (1. / m_dch_info->Lhalf) * tan(dphi / 2);
                
                // // point 1
                // double x1 = rz0;                                          // cm
                // double y1 = -stereosign * rz0 * kappa * m_dch_info->Lhalf;  // cm
                // double z1 = -m_dch_info->Lhalf;                             // cm
                // TVector3 p1(x1, y1, z1);
                
                // // point 2
                // double x2 = rz0;                                         // cm
                // double y2 = stereosign * rz0 * kappa * m_dch_info->Lhalf;  // cm
                // double z2 = m_dch_info->Lhalf;                             // cm
                // TVector3 p2(x2, y2, z2);
                
                // // calculate phi rotation of whole twisted tube, ie, rotation at z=0
                // double phi_z0 = m_dch_info->Calculate_wire_phi_z0(ilayer, nphi);
                // p1.RotateZ(phi_z0);
                // p2.RotateZ(phi_z0);

                // std::cout << "Wire extremities: " << std::endl;
                // std::cout << "Point 1: " << p1.X() << " " << p1.Y() << " " << p1.Z() << std::endl;
                // std::cout << "Point 2: " << p2.X() << " " << p2.Y() << " " << p2.Z() << std::endl;



                GenfitInterface::WireMeasurement measurement = GenfitInterface::WireMeasurement(wire_hit, m_dch_info, m_dc_decoder, detID, ++hit_idx, debug_lvl);
                if (dir > 0)
                {
                    m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
                }
                else
                {
                    m_genfitTrack_back->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack_back));
                }

                       
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
    bool GenfitTrack::fit(int dir = 1, double Beta_init = 100., double Beta_final=0.1, double Beta_steps=10, double Bz = 2., int debug_lvl = 0) {

        edm4hep::Track Track_temp = (dir > 0) ? m_edm4hepTrack : m_edm4hepTrack_back;
        for (size_t i = 0; i < Track_temp.trackStates_size(); ++i) {

            Track_temp.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* m_genfitFitter = new genfit::DAF(true, 1e-3,1e-3);
            m_genfitFitter->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);
            m_genfitFitter->setProbCut(1e-5);
            m_genfitFitter->setConvergenceDeltaWeight(1e-2); 
            m_genfitFitter->setDebugLvl(debug_lvl);
            
            // Process forward fit
            genfit::Track forwardTrack = *m_genfitTrack;
            genfit::AbsTrackRep* forwardRep = forwardTrack.getTrackRep(0);
            
            // auto list_of_points = forwardTrack.getPoints();
            // int idx2 = 0;
            // for (auto pt : list_of_points)
            // {   
            //     if (idx2 == 10)
            //     {   
            //         auto RAWmeas = pt->getRawMeasurement(0);

            //         std::cout << "Before fitting, TrackPoint n° " << idx2 << " measurements: " << std::endl;
            //         RAWmeas->Print();
            //     }
            //     ++idx2;
               
            // }

            m_genfitFitter->processTrackWithRep(&forwardTrack, forwardRep);

            // int numPoints = forwardTrack.getNumPoints();
            // for (int idx = 0; idx < numPoints; idx++)
            // {

            //     auto pt = forwardTrack.getPoint(idx);
            //     auto fitterInfo = pt->getFitterInfo(forwardRep);
            //     genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(fitterInfo);
                                    
            //     unsigned int nMeas = kfi->getNumMeasurements();
            //     if (1)
            //     {

            //         for(unsigned int j=0; j<nMeas; j++) {
                        
            //             if ( idx==10 )
            //             {
            //                 genfit::MeasurementOnPlane *mop = kfi->getMeasurementOnPlane(j);
            //                 std::cout << "Measurement n° " << idx << " weight: " << mop->getWeight() << std::endl;

            //                 // TVector3 posMeas = mop->getPos();
            //                 // // TMatrixDSym cov(6);
            //                 // // mop->getPosMomCov(pos, mom, cov);

            //                 // std::cout << "Measurement n° " << idx << " pos: " << posMeas.X() << " " << posMeas.Y() << " " << posMeas.Z() << std::endl;

            //                 mop->Print();

            //                 // if (mop->getWeight() < 10e-5)
            //                 // {
            //                 //     isWeightLow.push_back(-1);
            //                 // }
            //                 // else
            //                 // {
            //                 //     isWeightLow.push_back(1);
            //                 // }
            //             }
                        

            //         }
        


            //     }
            //     std::cout << "---------" << std::endl;

            // }

            // for (auto el : idx_to_remove)
            // {
            //     fitTrack.deletePoint(el);
            // }

         
            // if (!m_genfitFitter->isTrackFitted(&forwardTrack, forwardRep)) {
            //     return false;
            // }

            // int numPoints = forwardTrack.getNumPoints();
            // std::cout << numPoints << std::endl;
            // for (int idx = 0; idx < numPoints; idx++)
            // {
            //     std::cout << "here" << std::endl;
            //     auto pt = forwardTrack.getPoint(idx);
            //     auto rawMeas = pt->getRawMeasurements();
            //     const genfit::MeasuredStateOnPlane& measuredState = forwardTrack.getFittedState(idx);

            //     std::cout << typeid(rawMeas).name() << std::endl;
            //     for (size_t i = 0; i < rawMeas.size(); ++i) {
                    

            //         auto* meas = rawMeas[i];
            //         if (!meas) continue;

            //         std::cout << "rawMeas[" << i << "] type = "
            //                 << typeid(*meas).name()
            //                 << std::endl;
            //     }
                
            //     for (auto el : rawMeas)
            //     {

            //         std::vector<genfit::MeasurementOnPlane*> measPlanes = el->constructMeasurementsOnPlane(measuredState);

            //         for (auto el1: measPlanes)
            //         {
            //             std::cout << el1->getWeight() << std::endl;
            //         }

            //     }
            

                
            //     // std::cout << rawMeas.size() << std::endl;
               

                

            //     // for (size_t i = 0; i < rawMeas.size(); ++i) {
                    
            //     //     auto* meas = rawMeas[i];
            //     //     if (!meas) continue;
            //     //     auto test = meas->constructMeasurementsOnPlane(state);
            //     //     for (auto el : test)
            //     //     {
            //     //         auto weight = el->getWeight();
            //     //         std::cout << weight << std::endl;
            //     //     }
            //     //     std::cout << "--" << std::endl;
            //     // }
            // }


            // int numPoints = forwardTrack.getNumPoints();
            // for (int idx = 0; idx < numPoints; idx++)
            // {

            //     auto pt = forwardTrack.getPoint(idx);
            //     auto fitterInfo = pt->getFitterInfo(forwardRep);
            //     genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(fitterInfo);
            
            //     unsigned int nMeas = kfi->getNumMeasurements();

            //     for(unsigned int j=0; j<nMeas; j++) {

            //         genfit::MeasurementOnPlane *mop = kfi->getMeasurementOnPlane(j);
            //         double chi2 = 9999;
            //         double residual = 9999;

            //         std::cout << mop->getWeight() << std::endl;

            //     }
            //     std::cout << "---------" << std::endl;

            // }

            // // Process backward fit
            // genfit::Track backwardTrack = forwardTrack;
            // backwardTrack.reverseTrack();
            // m_genfitFitter->processTrack(&backwardTrack);
            // genfit::AbsTrackRep* backwardRep = backwardTrack.getTrackRep(0);
            // if (!m_genfitFitter->isTrackFitted(&backwardTrack, backwardRep)) {
            //     return false;
            // }
            
            // Update edm4hep track state
            genfit::MeasuredStateOnPlane fittedState;
            TVector3 gen_position, gen_momentum;
            TMatrixDSym covariancePosMom(6);
            
            double x0_PCA;
            double y0_PCA;
            double z0_PCA; 
            double pz;  

            double pt;

            double phi0;
            double tanLambda;
            double d0;
            double z0;
            double omega;  

            double c_light = 2.99792458e8;
            double a = c_light * 1e3 * 1e-15;
            
            int charge = getHypotesisCharge(m_particle_hypothesis);
            if (m_genfitFitter->isTrackFitted(&forwardTrack, forwardRep))
            {

                // trackState First Hit
                fittedState = forwardTrack.getFittedState();
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecFirstHit = fittedState.getState();
                
                edm4hep::TrackState trackStateFirstHit;
                x0_PCA = gen_position.X() / dd4hep::mm;
                y0_PCA = gen_position.Y() / dd4hep::mm;
                z0_PCA = gen_position.Z() / dd4hep::mm;
                pz = gen_momentum.Z();      
                pt = gen_momentum.Perp(); 

                phi0 = gen_momentum.Phi();  
                tanLambda = pz / pt;
                omega = (charge < 0) ? -a * Bz / abs(pt) : a * Bz / abs(pt);

                d0 = -(m_firstHit_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (m_firstHit_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                z0 = z0_PCA - m_firstHit_referencePoint.Z() / dd4hep::mm;
               
                
                trackStateFirstHit.D0 = d0;
                trackStateFirstHit.Z0 = z0;
                trackStateFirstHit.phi = phi0;
                trackStateFirstHit.omega = omega;
                trackStateFirstHit.tanLambda = tanLambda;
                // trackStateFirstHit.time = time;
                // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                // trackStateFirstHit.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, m_firstHit_referencePoint, timeError, covariancePosMom);
                trackStateFirstHit.referencePoint = edm4hep::Vector3f(m_firstHit_referencePoint.X() / dd4hep::mm, m_firstHit_referencePoint.Y() / dd4hep::mm, m_firstHit_referencePoint.Z() / dd4hep::mm);
                trackStateFirstHit.location = edm4hep::TrackState::AtFirstHit;

                // trackState lastHit
                fittedState = forwardTrack.getFittedState(forwardTrack.getNumPoints()-1);
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecLastHit = fittedState.getState();
                
                edm4hep::TrackState trackStateLastHit;
                x0_PCA = gen_position.X() / dd4hep::mm;
                y0_PCA = gen_position.Y() / dd4hep::mm;
                z0_PCA = gen_position.Z() / dd4hep::mm;
                pz = gen_momentum.Z();      
                pt = gen_momentum.Perp(); 

                phi0 = gen_momentum.Phi();  
                tanLambda = pz / pt;
                omega = (charge < 0) ? -a * Bz / abs(pt) : a * Bz / abs(pt);

                d0 = -(m_lastHit_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (m_lastHit_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                z0 = z0_PCA - m_lastHit_referencePoint.Z() / dd4hep::mm;

                trackStateLastHit.D0 = d0;
                trackStateLastHit.Z0 = z0;
                trackStateLastHit.phi = phi0;
                trackStateLastHit.omega = omega;
                trackStateLastHit.tanLambda = tanLambda;
                // trackStateLastHit.time = time;
                // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                // trackStateLastHit.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, m_lastHit_referencePoint, timeError, covariancePosMom);
                trackStateLastHit.referencePoint = edm4hep::Vector3f(m_lastHit_referencePoint.X() / dd4hep::mm, m_lastHit_referencePoint.Y() / dd4hep::mm, m_lastHit_referencePoint.Z() / dd4hep::mm);
                trackStateLastHit.location = edm4hep::TrackState::AtLastHit;

                // propagate to IP
                if(dir>0)
                {
                    forwardTrack.reverseTrack();
                }
                m_genfitFitter->processTrackWithRep(&forwardTrack, forwardRep);
                try{

                    TVector3 IP(0, 0, 0);
                    fittedState = forwardTrack.getFittedState(forwardTrack.getNumPoints()-1);
                    forwardRep->extrapolateToPoint(fittedState, IP);

                    fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                    auto stateVecIP = fittedState.getState();

                    gen_momentum.SetX(-gen_momentum.X());
                    gen_momentum.SetY(-gen_momentum.Y());
                    gen_momentum.SetZ(-gen_momentum.Z());
                    
                    edm4hep::TrackState trackStateIP;
                    x0_PCA = gen_position.X() / dd4hep::mm;
                    y0_PCA = gen_position.Y() / dd4hep::mm;
                    z0_PCA = gen_position.Z() / dd4hep::mm;
                    pz = gen_momentum.Z();      
                    pt = gen_momentum.Perp(); 

                    phi0 = gen_momentum.Phi();  
                    tanLambda = pz / pt;
                    omega = (charge < 0) ? -a * Bz / abs(pt) : a * Bz / abs(pt);

                    d0 = -(m_IP_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (m_IP_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                    z0 = z0_PCA - m_IP_referencePoint.Z() / dd4hep::mm;

                    trackStateIP.D0 = d0;
                    trackStateIP.Z0 = z0;
                    trackStateIP.phi = phi0;
                    trackStateIP.omega = omega;
                    trackStateIP.tanLambda = tanLambda;
                    // trackStateIP.time = time;
                    // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                    // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                    // trackStateIP.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, m_IP_referencePoint, timeError, covariancePosMom);
                    trackStateIP.referencePoint = edm4hep::Vector3f(m_IP_referencePoint.X() / dd4hep::mm, m_IP_referencePoint.Y() / dd4hep::mm, m_IP_referencePoint.Z() / dd4hep::mm);
                    trackStateIP.location = edm4hep::TrackState::AtIP;
                    
                    if (dir>0)
                    {
                        m_edm4hepTrack.addToTrackStates(trackStateIP);
                        m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
                        m_edm4hepTrack.addToTrackStates(trackStateLastHit);
                    }
                    else
                    {
                        m_edm4hepTrack_back.addToTrackStates(trackStateIP);
                        m_edm4hepTrack_back.addToTrackStates(trackStateFirstHit);
                        m_edm4hepTrack_back.addToTrackStates(trackStateLastHit);
                    }
                   

                }
                catch(...)
                {
                    return false;
                }

                if (!m_genfitFitter->isTrackFitted(&forwardTrack, forwardRep)) {
                    return false;
                }
                
                if (m_genfitFitter->isTrackFitted(&forwardTrack,forwardRep))
                {
                    
                    m_edm4hepTrack.setChi2(forwardTrack.getFitStatus()->getChi2());
                    m_edm4hepTrack.setNdf(forwardTrack.getFitStatus()->getNdf());

                }
                else
                {

                    m_edm4hepTrack.setChi2(-1);
                    m_edm4hepTrack.setNdf(-1);

                }

            }
            else
            {

                m_edm4hepTrack.setChi2(-1);
                m_edm4hepTrack.setNdf(-1);

            }

            return m_genfitFitter->isTrackFitted(&forwardTrack,forwardRep);
            
        }
        catch(...)
        {
            return false;
        }
 
    }
}
