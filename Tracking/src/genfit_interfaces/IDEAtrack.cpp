
#include "IDEAtrack.hpp"

namespace IDEAtracking {


    IDEAtrack::IDEAtrack(extension::Track track, TVector3 pos, TVector3 mom, int pdg)
    {

        extensionTrack_ = &track;
        genfitTrackRep_ = new genfit::RKTrackRep(pdg);
        genfitTrack_ = new genfit::Track(genfitTrackRep_, pos, mom);


    }


    void IDEAtrack::insertPoint(genfit::AbsMeasurement* measurement) {
        if (!genfitTrack_) {
            throw std::runtime_error("genfitTrack_ is null in insertPoint()");
        }
        
        if (!measurement) {
            throw std::invalid_argument("measurement is null in insertPoint()");
        }
    
        genfitTrack_->insertPoint(new genfit::TrackPoint(measurement, genfitTrack_));
    }
    

    genfit::Track* IDEAtrack::getTrack_genfit() const {
        return genfitTrack_;
    }

    extension::Track* IDEAtrack::getTrack_edm4hep() const {
        return extensionTrack_;
    }

    void IDEAtrack::processTrack()
    {
        genfitTrack_->checkConsistency();
        genfitFitter_ = new genfit::KalmanFitterRefTrack();

        genfitFitter_->processTrack(genfitTrack_);

        genfitTrack_->checkConsistency();

    }
    

}