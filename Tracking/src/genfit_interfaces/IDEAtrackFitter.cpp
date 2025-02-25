
#include "IDEAtrackFitter.hpp"

namespace IDEAtracking {


    IDEAtrackFitter::IDEAtrackFitter(TVector3 pos, TVector3 mom, int pdg)
    {

        genfitTrackRep_ = new genfit::RKTrackRep(pdg);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, pos, mom);


    }


    void IDEAtrackFitter::insertPoint(genfit::AbsMeasurement* measurement) {
        if (!genfitTrack_) {
            throw std::runtime_error("genfitTrack_ is null in insertPoint()");
        }
        
        if (!measurement) {
            throw std::invalid_argument("measurement is null in insertPoint()");
        }
    
        genfitTrack_->insertPoint(new genfit::TrackPoint(measurement, genfitTrack_));
    }
    

    extension::Track* IDEAtrackFitter::getTrack_edm4hep() const {
        return extensionTrack_;
    }

    genfit::Track* IDEAtrackFitter::getTrack_genfit() const {
        return genfitTrack_;
    }

    void IDEAtrackFitter::processTrack()
    {
        genfitTrack_->checkConsistency();
        genfitFitter_ = new genfit::KalmanFitterRefTrack();

        genfitFitter_->processTrack(genfitTrack_);

        genfitTrack_->checkConsistency();

    }
    

}