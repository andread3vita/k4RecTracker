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


namespace GenfitInterface {

    GenfitTrack::GenfitTrack(const edm4hep::Track& track,
                             const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder, 
                             const int particle_hypothesis, std::optional<TVector3> initial_position, std::optional<TVector3> initial_momentum)
        :   m_particle_hypothesis(particle_hypothesis), 
            m_posInit(0., 0., 0.), 
            m_momInit(0., 0., 0.), 
            m_genfitTrackRep(nullptr), 
            m_genfitTrack(nullptr), 
            m_edm4hepTrack(), 
            m_dch_info(dch_info),
            m_dc_decoder(decoder)
    {   

        checkInitialization();
        init(track, initial_position, initial_momentum);
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
    * @param initial_momentum Optional initial momentum for the track. If not provided, the direction from the first to the second hit is used.
    *
    */
    void GenfitTrack::init(const edm4hep::Track& track_init, std::optional<TVector3> initial_position, std::optional<TVector3> initial_momentum) {

        // Initialize the m_edm4hepTrack
        m_edm4hepTrack = edm4hep::MutableTrack();

        // Sort the hits by distance from the origin
        std::vector<std::pair<float, int>> hitDistIndices{};
        int index = 0;

        auto hits_in_track = track_init.getTrackerHits();
        for (auto hit : hits_in_track) {

            const auto pos_siHit = hit.getPosition();
            const auto distance = std::sqrt(pos_siHit.x * pos_siHit.x + pos_siHit.y * pos_siHit.y + pos_siHit.z * pos_siHit.z);
            hitDistIndices.emplace_back(distance, index++);

        }

        std::ranges::sort(hitDistIndices, {}, &std::pair<float, int>::first);

        // fill the m_edm4hepTrack
        for (const auto& [_, idx] : hitDistIndices) {
            m_edm4hepTrack.addToTrackerHits(hits_in_track[idx]);
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
            if (static_cast<size_t>(index_loopHit) == hits_for_genfit.size()-1) {

                auto position = hit.getPosition();
                m_lastHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);

            }
            index_loopHit++;
        }

        m_posInit = initial_position.value_or(m_firstHit_referencePoint);
        m_momInit = initial_momentum.value_or((m_secondHit_referencePoint - m_firstHit_referencePoint).Unit());

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
    void GenfitTrack::createGenFitTrack(int debug_lvl) {

        delete m_genfitTrackRep;
        delete m_genfitTrack;

        // Create covState ( 1 GeV , 10 cm )      
        TMatrixDSym covState(6);
        covState.Zero();  

        for (int i = 0; i < 3; ++i) {
            covState(i,i) = 100.;  // (10 cm)^2
        }
        for (int i = 3; i < 6; ++i) {
            covState(i,i) = 1.;  // (1 GeV)^2
        }

        // Create stateVec
        TVectorD stateVec(6);
        
        // pos
        stateVec[0] = m_posInit.X();
        stateVec[1] = m_posInit.Y();
        stateVec[2] = m_posInit.Z();

        // mom
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

    // /**
    // * @brief Fit the Genfit track using the Deterministic Annealing Filter (DAF).
    // * 
    // * This function performs a forward and backward track fit using the Genfit DAF algorithm,
    // * and populates the EDM4hep track states at the IP, first hit, and last hit positions.
    // * 
    // * @param Beta_init  Initial annealing parameter (temperature).
    // * @param Beta_final Final annealing parameter.
    // * @param Beta_steps Number of steps in the annealing schedule.
    // * @param Bz         Value of z-component of the magnetic field at the center of the detector
    // * @param debug_lvl  debug level: output if > 0
    // * 
    // * @return true if the fit was successful, false otherwise.
    // */
    bool GenfitTrack::fit(double Beta_init = 10., double Beta_final=0.1, double Beta_steps=10, double Bz = 2., int debug_lvl = 0) {

        for (size_t i = 0; i < m_edm4hepTrack.trackStates_size(); ++i) {

            m_edm4hepTrack.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* m_genfitFitter = new genfit::DAF(true, 1e-3,1e-3);
            genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
            genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
            m_genfitFitter->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);
            m_genfitFitter->setDebugLvl(debug_lvl);
            
            // Process forward fit
            genfit::Track forwardTrack = *m_genfitTrack;
            m_genfitFitter->processTrack(&forwardTrack);
            genfit::AbsTrackRep* forwardRep = forwardTrack.getTrackRep(0);
            if (!m_genfitFitter->isTrackFitted(&forwardTrack, forwardRep)) {
                return false;
            }

            // Process backward fit
            genfit::Track backwardTrack = forwardTrack;
            backwardTrack.reverseTrack();
            m_genfitFitter->processTrack(&backwardTrack);
            genfit::AbsTrackRep* backwardRep = backwardTrack.getTrackRep(0);
            if (!m_genfitFitter->isTrackFitted(&backwardTrack, backwardRep)) {
                return false;
            }
            
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
            if (m_genfitFitter->isTrackFitted(&forwardTrack,forwardRep))
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
                try{

                    TVector3 IP(0, 0, 0);
                    fittedState = backwardTrack.getFittedState(backwardTrack.getNumPoints()-1);
                    backwardRep->extrapolateToPoint(fittedState, IP);

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
                    
                    m_edm4hepTrack.addToTrackStates(trackStateIP);
                    m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
                    m_edm4hepTrack.addToTrackStates(trackStateLastHit);

                }
                catch(...)
                {
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