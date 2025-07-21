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

    GenfitTrack::GenfitTrack(const extension::Track& track,
                             const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder, 
                             const int particle_hypothesis, const TVector3& initial_position, const TVector3& initial_momentum)
        :   _particle_hypothesis(particle_hypothesis), 
            _posInit(0., 0., 0.), 
            _momInit(0., 0., 0.), 
            genfitTrackRep_(nullptr), 
            genfitTrack_(nullptr), 
            edm4hepTrack_(), 
            _dch_info(dch_info),
            _dc_decoder(decoder)
    {   

        checkInitialization();
        init(track,initial_position,initial_momentum);
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
    * an input track of type `extension::Track`. The following operations are performed:
    * - Initializes the internal `edm4hepTrack_` with a new `MutableTrack`.
    * - Sorts the tracker hits of the input track by their Euclidean distance from the origin.
    * - Fills the internal track (`edm4hepTrack_`) with the sorted hits.
    * - Extracts the positions of the first two hits and the last hit to initialize:
    *   - `_posInit`: the initial position (first hit).
    *   - `_momInit`: the initial direction (unit vector from first to second hit).
    *   - `firstHit_referencePoint`: copy of the first hit position.
    *   - `lastHit_referencePoint`: position of the last hit.
    *
    * @param track_init The input track containing tracker hits used for initialization.
    * @param initial_position Optional initial position for the track. If not provided, the first hit position is used.
    * @param initial_momentum Optional initial momentum for the track. If not provided, the direction from the first to the second hit is used.
    *
    * @note If the initial position or momentum is provided, they are used to set `_posInit` and `_momInit`.
    * The default values for these two vectors are (-1,-1,-1) for initial_position and (0,0,0) for initial_momentum. 
    * If these default values are used, the method computes the initial position and momentum based on the first two hits.
    *
    */
    void GenfitTrack::init(const extension::Track& track_init, const TVector3& initial_position, const TVector3& initial_momentum) {

        // Initialize the edm4hepTrack_
        edm4hepTrack_ = extension::MutableTrack();

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

        // fill the edm4hepTrack_
        for (const auto& [_, idx] : hitDistIndices) {
            edm4hepTrack_.addToTrackerHits(hits_in_track[idx]);
        }

        // initialize track
        int index_loopHit = 0;
        TVector3 first_hit(0,0,0);
        TVector3 second_hit(0,0,0);

        auto hits_for_genfit = edm4hepTrack_.getTrackerHits();
        for (auto hit : hits_for_genfit) {
            if (index_loopHit == 0) {
     
                auto position = hit.getPosition();
                first_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (index_loopHit == 1) {

                auto position = hit.getPosition();
                second_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (static_cast<size_t>(index_loopHit) == hits_for_genfit.size()-1) {

                auto position = hit.getPosition();
                lastHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);

            }
            index_loopHit++;
        }

        if (initial_position != TVector3(-1,-1,-1)) {
            _posInit = initial_position;
        }
        else
        {
           _posInit = first_hit;
        }

        if (initial_momentum != TVector3(0,0,0)) {
            _momInit = initial_momentum;
        }
        else
        {
            _momInit = (second_hit - first_hit).Unit();
        }
       
        firstHit_referencePoint = first_hit;

    }

    
    /**
    * @brief Creates a Genfit track object from the initialized state and tracker hits.
    * 
    * This method constructs the `genfit::Track` and associated `genfit::TrackRep` using the
    * initial position and momentum previously set during initialization. It assigns a default
    * covariance matrix (10 cm position resolution and 1 GeV momentum resolution), constructs
    * the state vector, and creates a `genfit::RKTrackRep` using the specified particle hypothesis.
    * 
    * Tracker hits from the internal `edm4hepTrack_` are processed and converted to Genfit-compatible
    * measurements. Based on the hit type, the method creates either planar or wire measurements.
    * These measurements are inserted into the `genfit::Track` as `genfit::TrackPoint` objects.
    * 
    * The supported hit types and their corresponding detector type IDs are:
    * - `edm4hep::TrackerHitPlane`: Planar detector hit
    * - `extension::SenseWireHit`: Drift chamber wire hit
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

        delete genfitTrackRep_;
        delete genfitTrack_;

        // Create covState ( 1 GeV + 10 cm )      
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
        stateVec[0] = _posInit.X();
        stateVec[1] = _posInit.Y();
        stateVec[2] = _posInit.Z();

        // mom
        stateVec[3] = _momInit.X();
        stateVec[4] = _momInit.Y();
        stateVec[5] = _momInit.Z();


        genfitTrackRep_ = new genfit::RKTrackRep(_particle_hypothesis);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, stateVec, covState);
                
        auto hits_for_genfit = edm4hepTrack_.getTrackerHits();

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
                GenfitInterface::Planar_measurement measurement = GenfitInterface::Planar_measurement(planar_hit,detID,++hit_idx,debug_lvl);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                
            }
            else if (hit.isA<extension::SenseWireHit>()) 
            {

                auto wire_hit =  hit.as<extension::SenseWireHit>();
                GenfitInterface::Wire_measurement measurement = GenfitInterface::Wire_measurement(wire_hit,_dch_info,_dc_decoder,detID,++hit_idx,debug_lvl);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                       
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
    * 
    * @return true if the fit was successful, false otherwise.
    */
    bool GenfitTrack::fit(double Beta_init = 10., double Beta_final=0.1, double Beta_steps=10, double Bz = 2.) {

        for (size_t i = 0; i < edm4hepTrack_.trackStates_size(); ++i) {

            edm4hepTrack_.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* genfitFitter_ = new genfit::DAF(true, 1e-3,1e-3);
            genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
            genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
            genfitFitter_->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);
            // genfitFitter_->setDebugLvl(1);
            
            // Process forward fit
            genfit::Track forwardTrack = *genfitTrack_;
            genfitFitter_->processTrack(&forwardTrack);
            genfit::AbsTrackRep* forwardRep = forwardTrack.getTrackRep(0);
            if (!genfitFitter_->isTrackFitted(&forwardTrack, forwardRep)) {
                return false;
            }

            // Process backward fit
            genfit::Track backwardTrack = forwardTrack;
            backwardTrack.reverseTrack();
            genfitFitter_->processTrack(&backwardTrack);
            genfit::AbsTrackRep* backwardRep = backwardTrack.getTrackRep(0);
            if (!genfitFitter_->isTrackFitted(&backwardTrack, backwardRep)) {
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
            
            int charge = getHypotesisCharge(_particle_hypothesis);
            if (genfitFitter_->isTrackFitted(&forwardTrack,forwardRep))
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

                d0 = -(firstHit_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (firstHit_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                z0 = z0_PCA - firstHit_referencePoint.Z() / dd4hep::mm;
               
                
                trackStateFirstHit.D0 = d0;
                trackStateFirstHit.Z0 = z0;
                trackStateFirstHit.phi = phi0;
                trackStateFirstHit.omega = omega;
                trackStateFirstHit.tanLambda = tanLambda;
                // trackStateFirstHit.time = time;
                // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                // trackStateFirstHit.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, firstHit_referencePoint, timeError, covariancePosMom);
                trackStateFirstHit.referencePoint = edm4hep::Vector3f(firstHit_referencePoint.X() / dd4hep::mm, firstHit_referencePoint.Y() / dd4hep::mm, firstHit_referencePoint.Z() / dd4hep::mm);
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

                d0 = -(lastHit_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (lastHit_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                z0 = z0_PCA - lastHit_referencePoint.Z() / dd4hep::mm;

                trackStateLastHit.D0 = d0;
                trackStateLastHit.Z0 = z0;
                trackStateLastHit.phi = phi0;
                trackStateLastHit.omega = omega;
                trackStateLastHit.tanLambda = tanLambda;
                // trackStateLastHit.time = time;
                // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                // trackStateLastHit.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, lastHit_referencePoint, timeError, covariancePosMom);
                trackStateLastHit.referencePoint = edm4hep::Vector3f(lastHit_referencePoint.X() / dd4hep::mm, lastHit_referencePoint.Y() / dd4hep::mm, lastHit_referencePoint.Z() / dd4hep::mm);
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

                    d0 = -(IP_referencePoint.X() / dd4hep::mm - x0_PCA)*sin(phi0) + (IP_referencePoint.Y() / dd4hep::mm - y0_PCA)*cos(phi0);
                    z0 = z0_PCA - IP_referencePoint.Z() / dd4hep::mm;

                    trackStateIP.D0 = d0;
                    trackStateIP.Z0 = z0;
                    trackStateIP.phi = phi0;
                    trackStateIP.omega = omega;
                    trackStateIP.tanLambda = tanLambda;
                    // trackStateIP.time = time;
                    // TVectorD trackStateGenfit(x0_PCA,y0_PCA,Z0_PCA,px,py,pz):
                    // TVectorD params(omega, phi0,d0,z0,tanLambda,time);
                    // trackStateIP.covMatrix = computeTrackStateCovMatrix(trackStateGenfit, params, IP_referencePoint, timeError, covariancePosMom);
                    trackStateIP.referencePoint = edm4hep::Vector3f(IP_referencePoint.X() / dd4hep::mm, IP_referencePoint.Y() / dd4hep::mm, IP_referencePoint.Z() / dd4hep::mm);
                    trackStateIP.location = edm4hep::TrackState::AtIP;
                    
                    edm4hepTrack_.addToTrackStates(trackStateIP);
                    edm4hepTrack_.addToTrackStates(trackStateFirstHit);
                    edm4hepTrack_.addToTrackStates(trackStateLastHit);

                }
                catch(...)
                {
                    return false;
                }

                if (genfitFitter_->isTrackFitted(&forwardTrack,forwardRep))
                {
                    
                    edm4hepTrack_.setChi2(forwardTrack.getFitStatus()->getChi2());
                    edm4hepTrack_.setNdf(forwardTrack.getFitStatus()->getNdf());

                }
                else
                {

                    edm4hepTrack_.setChi2(-1);
                    edm4hepTrack_.setNdf(-1);

                }

                
            }
            else
            {

                edm4hepTrack_.setChi2(-1);
                edm4hepTrack_.setNdf(-1);

            }

            return genfitFitter_->isTrackFitted(&forwardTrack,forwardRep);
            
        }
        catch(...)
        {
            return false;
        }

 
    }
}
