#include "GenfitTrack.hpp"
#include<ranges>

#define CHECK_MAP(m) if (!(m)) std::cerr << #m " is null!\n";

namespace GENFIT {

    GenfitTrack::GenfitTrack(const extension::Track& track,const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder, const int particle_hypothesis)
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
        init(track);
    }

    GenfitTrack::~GenfitTrack() {}

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

    void GenfitTrack::init(const extension::Track& track_init) {

        // // Check if the track is empty
        // if (track_init.getTrackerHits().empty()) {
        //     std::cerr << "Error: Track has no hits!" << std::endl;
        //     std::exit(EXIT_FAILURE);
        // }

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
            index_loopHit++;
        }

        _posInit = first_hit;
        _momInit = (second_hit - first_hit).Unit();

        // _posInit.Print();
        // _momInit.Print();

    }

    void GenfitTrack::createGenFitTrack(int debug_lvl) {

        delete genfitTrackRep_;
        delete genfitTrack_;

        genfitTrackRep_ = new genfit::RKTrackRep(_particle_hypothesis);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, _posInit, _momInit);
                
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
                GENFIT::Planar_measurement measurement = GENFIT::Planar_measurement(planar_hit,detID,++hit_idx,debug_lvl);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                
            }
            else if (hit.isA<extension::SenseWireHit>()) 
            {

                auto wire_hit =  hit.as<extension::SenseWireHit>();
                GENFIT::Wire_measurement measurement = GENFIT::Wire_measurement(wire_hit,_dch_info,_dc_decoder,detID,++hit_idx,debug_lvl);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                       
            } 
            else 
            {

                std::cerr << "Error: No hits with cellID: " << cellID0 << std::endl;
                std::exit(EXIT_FAILURE);
             
            }

        }

    }

    bool GenfitTrack::fit(double Beta_init = 10., double Beta_final=0.1, double Beta_steps=10) {

        for (size_t i = 0; i < edm4hepTrack_.trackStates_size(); ++i) {

            edm4hepTrack_.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* genfitFitter_ = new genfit::DAF(true, 1e-3,1e-3);
            genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
            genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
            genfitFitter_->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);

            // genfit::KalmanFitterRefTrack* genfitFitter_ = new genfit::KalmanFitterRefTrack();
            // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
            // genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
           
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
            double x0;
            double y0;
            double z0_hit;
            double phi0;
            double pz;   
            double pt;
            double tanLambda;
            double c_light = 2.99792458e8;
            double a = c_light * 1e3 * 1e-15;
            double d0;
            double z0;
            double omega;  
            int charge = getHypotesisCharge(_particle_hypothesis);
            if (genfitFitter_->isTrackFitted(&forwardTrack,forwardRep))
            {

                // trackState First Hit
                fittedState = forwardTrack.getFittedState();
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecFirstHit = fittedState.getState();
                
                edm4hep::TrackState trackStateFirstHit;
                x0 = gen_position.X()*10.;
                y0 = gen_position.Y()*10.;
                z0_hit = gen_position.Z()*10.;  
                phi0 = gen_momentum.Phi();  
                pz = gen_momentum.Z();      
                pt = gen_momentum.Perp();  
                tanLambda = pz / pt;
                d0 = -x0 * sin(phi0) + y0 * cos(phi0);
                z0 = z0_hit - d0 * tanLambda;
                omega = a * 2.0 / abs(pt);
                if (charge < 0)
                {
                    omega = -omega;
                }
                
                trackStateFirstHit.D0 = d0;
                trackStateFirstHit.Z0 = z0;
                trackStateFirstHit.phi = phi0;
                trackStateFirstHit.omega = omega;
                trackStateFirstHit.tanLambda = tanLambda;
                // trackStateFirstHit.time = 
                // trackStateFirstHit.covMatrix = 
                trackStateFirstHit.referencePoint = edm4hep::Vector3f(gen_position[0]*10., gen_position[1]*10., gen_position[2]*10.);
                trackStateFirstHit.location = edm4hep::TrackState::AtFirstHit;

                // trackState lastHit
                fittedState = forwardTrack.getFittedState(forwardTrack.getNumPoints()-1);
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecLastHit = fittedState.getState();
                
                edm4hep::TrackState trackStateLastHit;
                x0 = gen_position.X()*10.;
                y0 = gen_position.Y()*10.;
                z0_hit = gen_position.Z()*10.; 
                phi0 = gen_momentum.Phi();  
                pz = gen_momentum.Z();      
                pt = gen_momentum.Perp();  
                tanLambda = pz / pt;
                d0 = -x0 * sin(phi0) + y0 * cos(phi0);
                z0 = z0_hit - d0 * tanLambda;
                omega = a * 2.0 / abs(pt);
                if (charge < 0)
                {
                    omega = -omega;
                }

                trackStateLastHit.D0 = d0;
                trackStateLastHit.Z0 = z0;
                trackStateLastHit.phi = phi0;
                trackStateLastHit.omega = omega;
                trackStateLastHit.tanLambda = tanLambda;
                // trackStateLastHit.time = 
                // trackStateLastHit.covMatrix = 
                trackStateLastHit.referencePoint = edm4hep::Vector3f(gen_position[0]*10., gen_position[1]*10., gen_position[2]*10.);
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
                    x0 = gen_position.X()*10.;
                    y0 = gen_position.Y()*10.;
                    z0_hit = gen_position.Z()*10.; 
                    phi0 = gen_momentum.Phi();  
                    pz = gen_momentum.Z();      
                    pt = gen_momentum.Perp();  
                    tanLambda = pz / pt;
                    d0 = -x0 * sin(phi0) + y0 * cos(phi0);
                    z0 = z0_hit - d0 * tanLambda;
                    omega = a * 2.0 / abs(pt);
                    if (charge < 0)
                    {
                        omega = -omega;
                    }

                    trackStateIP.D0 = d0;
                    trackStateIP.Z0 = z0;
                    trackStateIP.phi = phi0;
                    trackStateIP.omega = omega;
                    trackStateIP.tanLambda = tanLambda;
                    // trackStateLastHit.time = 
                    // trackStateLastHit.covMatrix = 
                    trackStateIP.referencePoint = edm4hep::Vector3f(gen_position[0]*10., gen_position[1]*10., gen_position[2]*10.);
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
