//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit
/// A genfit hit can be created
///
/// In this file, including:
///   a genfit hit class
///
///   Units are following GenfitUnit
#ifndef RECGENFITALG_GENFITHIT_H
#define RECGENFITALG_GENFITHIT_H
#include "TVector3.h"

namespace extension{
    class TrackerHit;
}

namespace edm4hep{
    class SimTrackerHit;
}

class GenfitHit{
    public:
        GenfitHit(const extension::TrackerHit* trackerHit,
                const edm4hep::SimTrackerHit* simTrackerHit,
                double driftDistance,
                double driftDistanceErr);

        ~GenfitHit(){;}
        double getDriftDistance()const{return m_driftDistance;}
        double getDriftDistanceErr()const{return m_driftDistanceErr;}
        const edm4hep::SimTrackerHit* getSimTrackerHit()const{return m_simTrackerHit;}
        const extension::TrackerHit* getTrackerHit()const{return m_trackerHit;}
        TVector3 getEnd0()const;
        TVector3 getEnd1()const;
        TVector3 getTruthPos()const;
        TVector3 getTruthMom()const;
        double getMaxDistance()const{return 0.6*1.4;}//FIXME
        // int getLeftRightAmbig()const;

        void setDriftDistance(double d){m_driftDistance=d;}
        void setDriftDistanceErr(double de){m_driftDistanceErr=de;}
        void print()const;

    private:
        
        const extension::TrackerHit* m_trackerHit;
        const edm4hep::SimTrackerHit* m_simTrackerHit;
        double m_driftDistance;
        double m_driftDistanceErr;
};

#endif