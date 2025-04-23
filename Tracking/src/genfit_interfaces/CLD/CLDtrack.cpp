#include "CLDtrack.hpp"
#include<ranges>

namespace GENFIT {

    CLDtrack::CLDtrack(const edm4hep::Track& track, const dd4hep::rec::SurfaceManager* surfMan)
        :   _particle_hypothesis(211), 
            _posInit(0., 0., 0.), 
            _momInit(0., 0., 0.), 
            genfitTrackRep_(nullptr), 
            genfitTrack_(nullptr), 
            edm4hepTrack_(), 
            _surfMan(surfMan), 
            surfaceMap_vertex(nullptr), 
            surfaceMap_InnerTrackers(nullptr), 
            surfaceMap_OuterTrackers(nullptr)
    {   

        checkInitialization();

        if (!_surfMan) {
            throw std::runtime_error("Error: SurfaceManager is null!");
        }

        surfaceMap_vertex = surfMan->map("Vertex");
        surfaceMap_InnerTrackers = surfMan->map("InnerTrackers");
        surfaceMap_OuterTrackers = surfMan->map("OuterTrackers");

        init(track);
    }

    CLDtrack::~CLDtrack() {}

    void CLDtrack::checkInitialization() {

        if (!gGeoManager) {
            std::cerr << "Error: TGeoManager is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (!genfit::FieldManager::getInstance()->isInitialized()) {
            std::cerr << "Error: FieldManager is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
            std::cerr << "Error: MaterialEffects is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void CLDtrack::init(const edm4hep::Track& track_init) {

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
        TVector3 last_hit(0,0,0);
        TVector3 lastButOne_hit(0,0,0);

        auto hits_for_genfit = edm4hepTrack_.getTrackerHits();
        for (auto hit : hits_for_genfit) {
            if (index_loopHit == 0) {
                auto vtx_hit = hit.as<edm4hep::TrackerHitPlane>();
                auto position = vtx_hit.getPosition();
                first_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (index_loopHit == 1) {
                auto vtx_hit = hit.as<edm4hep::TrackerHitPlane>();
                auto position = vtx_hit.getPosition();
                second_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (index_loopHit == static_cast<int>(hits_for_genfit.size())-1) {
                auto vtx_hit = hit.as<edm4hep::TrackerHitPlane>();
                auto position = vtx_hit.getPosition();
                last_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            if (index_loopHit == static_cast<int>(hits_for_genfit.size())-2) {
                auto vtx_hit = hit.as<edm4hep::TrackerHitPlane>();
                auto position = vtx_hit.getPosition();
                lastButOne_hit = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z);
            }
            index_loopHit++;
        }

        _posInit = first_hit;
        _momInit = (second_hit - first_hit).Unit();
        _posInit_back = last_hit;
        _momInit_back = (lastButOne_hit - last_hit).Unit();


    }

    void CLDtrack::createGenFitTrack() {

        genfitTrackRep_ = new genfit::RKTrackRep(_particle_hypothesis);
        // genfitTrackRepBack_ = new genfit::RKTrackRep(_particle_hypothesis);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, _posInit, _momInit);
        

        auto hits_for_genfit = edm4hepTrack_.getTrackerHits();
        int vtx_idx(0);
        for (auto hit : hits_for_genfit)
        {

            auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
            auto cellID0 = vtx_hit.getCellID();     
                    
            int detID(-1);
            dd4hep::rec::SurfaceMap::const_iterator sI;
            if (surfaceMap_vertex && (sI = surfaceMap_vertex->find(cellID0)) != surfaceMap_vertex->end()) {

                detID= 0;
                GENFIT::SI_measurement measurement = GENFIT::SI_measurement(vtx_hit,surfaceMap_vertex,detID,++vtx_idx);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                
                        
            } 
            else if (surfaceMap_InnerTrackers && (sI = surfaceMap_InnerTrackers->find(cellID0)) != surfaceMap_InnerTrackers->end()) {

                detID= 1;
                GENFIT::SI_measurement measurement = GENFIT::SI_measurement(vtx_hit,surfaceMap_InnerTrackers,detID,++vtx_idx);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
            
            } 
            else if (surfaceMap_OuterTrackers && (sI = surfaceMap_OuterTrackers->find(cellID0)) != surfaceMap_OuterTrackers->end()) {

                detID= 2;
                GENFIT::SI_measurement measurement = GENFIT::SI_measurement(vtx_hit,surfaceMap_OuterTrackers,detID,++vtx_idx);
                genfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), genfitTrack_));
                        
            } 
            else {

                std::cerr << "Error: Surface not found for cellID: " << cellID0 << std::endl;
                std::exit(EXIT_FAILURE);
             
            }

        }

        // genfit::AbsKalmanFitter* genfitFitter_ = new genfit::KalmanFitterRefTrack();
        // genfitFitter_->processTrack(genfitTrack_);
        // genfit::MeasuredStateOnPlane fittedState = genfitTrack_->getFittedState(genfitTrack_->getNumPoints()-1);
        // TVector3 gen_position, gen_momentum;
        // TMatrixDSym covariancePosMom(6);
        // fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);

        // gen_momentum[0] = -gen_momentum[0];
        // gen_momentum[1] = -gen_momentum[1];
        // gen_momentum[2] = -gen_momentum[2];

        // _posInit_back = gen_position;
        // _momInit_back = gen_momentum;
        // backwardGenfitTrack_ = new genfit::Track(genfitTrackRepBack_, _posInit_back, _momInit_back);

        // int vtx_idx_back(0);
        // for (auto hit : std::views::reverse(hits_for_genfit))
        // {    

        //     auto vtx_hit =  hit.as<edm4hep::TrackerHitPlane>();
        //     auto cellID0 = vtx_hit.getCellID();     
                    
        //     int detID(-1);
        //     dd4hep::rec::SurfaceMap::const_iterator sI;
        //     if (surfaceMap_vertex && (sI = surfaceMap_vertex->find(cellID0)) != surfaceMap_vertex->end()) {

        //         detID= 2;
        //         CLDtracking::CLD_measurement measurement = CLDtracking::CLD_measurement(vtx_hit,surfaceMap_vertex,detID,++vtx_idx_back);
        //         backwardGenfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), backwardGenfitTrack_));
                
                        
        //     } 
        //     else if (surfaceMap_InnerTrackers && (sI = surfaceMap_InnerTrackers->find(cellID0)) != surfaceMap_InnerTrackers->end()) {

        //         detID= 1;
        //         CLDtracking::CLD_measurement measurement = CLDtracking::CLD_measurement(vtx_hit,surfaceMap_InnerTrackers,detID,++vtx_idx_back);
        //         backwardGenfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), backwardGenfitTrack_));
            
        //     } 
        //     else if (surfaceMap_OuterTrackers && (sI = surfaceMap_OuterTrackers->find(cellID0)) != surfaceMap_OuterTrackers->end()) {

        //         detID= 0;
        //         CLDtracking::CLD_measurement measurement = CLDtracking::CLD_measurement(vtx_hit,surfaceMap_OuterTrackers,detID,++vtx_idx_back);
        //         backwardGenfitTrack_->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), backwardGenfitTrack_));
                        
        //     } 
        //     else {

        //         std::cerr << "Error: Surface not found for cellID: " << cellID0 << std::endl;
        //         std::exit(EXIT_FAILURE);
             
        //     }

        // }

       
        // _posInit.Print();
        // _momInit.Print();
        // _posInit_back.Print();
        // _momInit_back.Print();
        // std::cout << "" << std::endl;




    }

    bool CLDtrack::fit() {

        try{

            genfit::AbsKalmanFitter* genfitFitter_ = new genfit::KalmanFitterRefTrack();
            // genfit::AbsKalmanFitter* genfitFitter_back = new genfit::KalmanFitterRefTrack();

            genfitTrack_->checkConsistency();
            // backwardGenfitTrack_->checkConsistency();

            genfit::Track forwardTrack = *genfitTrack_;
            
            genfitFitter_->processTrack(&forwardTrack);
            // std::cout << "Chi2 forward: " << forwardTrack.getFitStatus()->getChi2() << std::endl;

            genfit::Track backwardTrack = forwardTrack;
            backwardTrack.reverseTrack();
            genfitFitter_->processTrack(&backwardTrack);
            // std::cout << "Chi2 backward: " << backwardTrack.getFitStatus()->getChi2() << std::endl;

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

            genfit::AbsTrackRep* forwardRep = forwardTrack.getTrackRep(0);
            genfit::AbsTrackRep* backwardRep = backwardTrack.getTrackRep(0);

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
                if (_particle_hypothesis < 0)
                {
                    omega = -omega; // _particle_hypothesis uses the pdg numbering scheme, so the positive particles have positive number encodings
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
