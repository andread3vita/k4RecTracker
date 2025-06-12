#include "IDEAtrack.hpp"
#include<ranges>

#define CHECK_MAP(m) if (!(m)) std::cerr << #m " is null!\n";

namespace GENFIT {

    // (d0, phi0, omega, z0, tanLambda, time)
    TVectorD computeTrackParams(const TVectorD& x) {

        double c_light = 2.99792458e8;
        double a = c_light * 1e3 * 1e-15;

        double x0 = x(0) * 10.;
        double y0 = x(1) * 10.;
        double z0_hit = x(2) * 10.;
        double px = x(3);
        double py = x(4);
        double pz = x(5);

        double pt = std::sqrt(px*px + py*py);
        double phi0 = std::atan2(py, px);
        double tanLambda = pz / pt;
        double d0 = -x0 * std::sin(phi0) + y0 * std::cos(phi0);
        double z0 = z0_hit - d0 * tanLambda;
        double omega = a * 2.0 / std::abs(pt);
        double time = 0.0;  // placeholder for time

        TVectorD y(6);
        y[0] = d0;
        y[1] = phi0;
        y[2] = omega;
        y[3] = z0;
        y[4] = tanLambda;
        y[5] = time;
        return y;
    }

    TMatrixD computeJacobian(const TVectorD& x, double epsilon = 1e-5) {

        TVectorD y0 = computeTrackParams(x);
        TMatrixD J(6, 6);

        for (int i = 0; i < 6; ++i) {
            TVectorD x_eps = x;
            x_eps[i] += epsilon;
            TVectorD y_eps = computeTrackParams(x_eps);
            for (int j = 0; j < 6; ++j) {
                J(j, i) = (y_eps[j] - y0[j]) / epsilon;
            }
        }

        return J;
    }

    IDEAtrack::IDEAtrack(const extension::Track& track, const dd4hep::rec::SurfaceManager* surfMan,const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder)
        :   _particle_hypothesis(211), 
            _posInit(0., 0., 0.), 
            _momInit(0., 0., 0.), 
            genfitTrackRep_(nullptr), 
            genfitTrack_(nullptr), 
            edm4hepTrack_(), 
            _surfMan(surfMan), 
            surfaceMap_vertex(nullptr), 
            surfaceMap_wrapper_barrel(nullptr),
            surfaceMap_wrapper_endcap(nullptr),
            _dch_info(dch_info),
            _dc_decoder(decoder)
    {   

        checkInitialization();

        if (!_surfMan) {
            throw std::runtime_error("Error: SurfaceManager is null!");
        }

        if (!_dch_info) {
            throw std::runtime_error("Error: DCH_INFO is null!");
        }

        surfaceMap_vertex = _surfMan->map("Vertex");
        surfaceMap_wrapper_barrel = _surfMan->map("SiWrB");
        surfaceMap_wrapper_endcap = _surfMan->map("SiWrD");

        CHECK_MAP(surfaceMap_vertex);
        CHECK_MAP(surfaceMap_wrapper_barrel);
        CHECK_MAP(surfaceMap_wrapper_endcap);

        init(track);
    }

    IDEAtrack::~IDEAtrack() {}

    void IDEAtrack::checkInitialization() {

        // if (!gGeoManager) {
        //     std::cerr << "Error: TGeoManager is not initialized!" << std::endl;
        //     std::exit(EXIT_FAILURE);
        // }

        // if (!genfit::FieldManager::getInstance()->isInitialized()) {
        //     std::cerr << "Error: FieldManager is not initialized!" << std::endl;
        //     std::exit(EXIT_FAILURE);
        // }

        // if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
        //     std::cerr << "Error: MaterialEffects is not initialized!" << std::endl;
        //     std::exit(EXIT_FAILURE);
        // }
    }

    void IDEAtrack::init(const extension::Track& track_init) {

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

        _posInit.Print();
        _momInit.Print();

    }

    void IDEAtrack::createGenFitTrack() {

        genfitTrackRep_ = new genfit::RKTrackRep(_particle_hypothesis);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, _posInit, _momInit);
        
        auto hits_for_genfit = edm4hepTrack_.getTrackerHits();

        int hit_idx(0);
        for (auto hit : hits_for_genfit)
        {

            int detID(-1);
            auto cellID0 = hit.getCellID();
            if (hit.isA<edm4hep::TrackerHitPlane>())
            {
               
               dd4hep::rec::SurfaceMap::const_iterator sI;    

                if (surfaceMap_vertex && (sI = surfaceMap_vertex->find(cellID0)) != surfaceMap_vertex->end()) {

                    detID = 0;
                    auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
                    GENFIT::SI_measurement measurement = GENFIT::SI_measurement(vtx_hit,surfaceMap_vertex,detID,++hit_idx);
                    genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                    
                            
                } 
                else if (surfaceMap_wrapper_barrel && (sI = surfaceMap_wrapper_barrel->find(cellID0)) != surfaceMap_wrapper_barrel->end()) {

                    detID = 2;
                    auto wrapper_hit =  hit.as<edm4hep::TrackerHitPlane>();
                    GENFIT::SI_measurement measurement = GENFIT::SI_measurement(wrapper_hit,surfaceMap_wrapper_barrel,detID,++hit_idx);
                    genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                    
                            
                } 
                else if (surfaceMap_wrapper_endcap && (sI = surfaceMap_wrapper_endcap->find(cellID0)) != surfaceMap_wrapper_endcap->end()) {
    
                    detID = 3;
                    auto wrapper_hit =  hit.as<edm4hep::TrackerHitPlane>();
                    GENFIT::SI_measurement measurement = GENFIT::SI_measurement(wrapper_hit,surfaceMap_wrapper_endcap,detID,++hit_idx);
                    genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_)); 
                            
                } 
            }
            else if (hit.isA<extension::SenseWireHit>()) {

                detID = 1;
                auto dc_hit =  hit.as<extension::SenseWireHit>();
                GENFIT::DC_measurement measurement = GENFIT::DC_measurement(dc_hit,_dch_info,_dc_decoder,detID,++hit_idx);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                
                        
            } 
            else {

                std::cerr << "Error: No hits with cellID: " << cellID0 << std::endl;
                std::exit(EXIT_FAILURE);
             
            }

        }

    }

    bool IDEAtrack::fit(double Beta_init = 10., double Beta_final=0.1, double Beta_steps=10) {

        try{

            genfit::DAF* genfitFitter_ = new genfit::DAF(true, 1e-3,1e-3);
            genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
            genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
            genfitFitter_->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);

            // genfitFitter_->setDebugLvl(1);

            genfit::Track forwardTrack = *genfitTrack_;
            genfitFitter_->processTrack(&forwardTrack);

            genfit::AbsTrackRep* forwardRep = forwardTrack.getTrackRep(0);
            if (!genfitFitter_->isTrackFitted(&forwardTrack, forwardRep)) {
                return false;
            }

            int nUsedHits = 0;
            for (unsigned int i = 0; i < forwardTrack.getNumPoints(); ++i) {
                genfit::TrackPoint* point = forwardTrack.getPoint(i);

                if (point->hasFitterInfo(forwardRep)) {
                    ++nUsedHits;
                }
            }

            // std::cout << "Number of used hits: " << nUsedHits << std::endl;
            
            genfit::Track backwardTrack = forwardTrack;
            backwardTrack.reverseTrack();
            genfitFitter_->processTrack(&backwardTrack);

            genfit::AbsTrackRep* backwardRep = backwardTrack.getTrackRep(0);
            if (!genfitFitter_->isTrackFitted(&backwardTrack, backwardRep)) {
                return false;
            }
            
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

                // TVectorD x = TVectorD(gen_position.X(),gen_position.Y(),gen_position.Z(),gen_momentum.X(),gen_momentum.Y(),gen_momentum.Z());
                // TVectorD y = computeTrackParams(x);

                // double omega = y(2);
                if (_particle_hypothesis < 0)
                {
                    omega = -omega; // _particle_hypothesis uses the pdg numbering scheme, so the positive particles have positive number encodings
                }
                // y(2) = omega;

                // TMatrixD J = computeJacobian(x);    
                // TMatrixD JT = TMatrixD(TMatrixD::kTransposed, J);   // J^T
                // TMatrixD temp(6, 6);
                // temp.Mult(J, covariancePosMom);                     // temp = J * Cx
                // TMatrixD Cy(6, 6);
                // Cy.Mult(temp, JT);                                  // Cy = temp * J^T

                // std::array<float, 21> covVals;
                // int index = 0;
                // for (int i = 0; i < 6; ++i) {
                //     for (int j = i; j < 6; ++j) {
                //         covVals[index++] = static_cast<float>(Cy(i, j));
                //     }
                // }
                
                // trackStateFirstHit.D0 = x(0);
                // trackStateFirstHit.phi = x(1);
                // trackStateFirstHit.omega = x(2);
                // trackStateFirstHit.Z0 = x(3);
                // trackStateFirstHit.tanLambda = x(4);
                // trackStateFirstHit.time = x(5);         // time is not used in the current implementation
                // trackStateFirstHit.covMatrix = edm4hep::CovMatrix6f(
                //     covVals[0], covVals[1], covVals[2], covVals[3], covVals[4], covVals[5],
                //     covVals[6], covVals[7], covVals[8], covVals[9], covVals[10],
                //     covVals[11], covVals[12], covVals[13], covVals[14], covVals[15],
                //     covVals[16], covVals[17], covVals[18], covVals[19], covVals[20]
                // );
                // trackStateFirstHit.referencePoint = edm4hep::Vector3f(gen_position[0]*10., gen_position[1]*10., gen_position[2]*10.);
                // trackStateFirstHit.location = edm4hep::TrackState::AtFirstHit;

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
                if (_particle_hypothesis < 0)
                {
                    omega = -omega; // _particle_hypothesis uses the pdg numbering scheme, so the positive particles have positive number encodings
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
                    if (_particle_hypothesis < 0)
                    {
                        omega = -omega; // _particle_hypothesis uses the pdg numbering scheme, so the positive particles have positive number encodings
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
                    
                    // std::cout << forwardTrack.getFitStatus()->getChi2() << std::endl;
                    // std::cout << forwardTrack.getFitStatus()->getNdf() << std::endl;
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
