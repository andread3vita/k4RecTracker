//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit
/// A genfit track can be created, fitted and extrapolated
/// Track with only track representation(s)
///
/// In this file, including:
///   a genfit track class
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Zhang Yao (zhangyao@ihep.ac.cn)
///
//////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITTRACK_H
#define RECGENFITALG_GENFITTRACK_H

#include "GenfitFitter.hpp"
#include "GenfitHit.hpp"

//ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TLorentzVector.h"

//Gaudi
#include "GaudiKernel/SmartIF.h"

//STL
#include <vector>

#include "WireMeasurementDC.hpp"
#include "PlanarMeasurement.hpp"

class TLorentzVector;
class IGeomSvc;
class WireMeasurementDC;

namespace genfit{
    class Track;
    class FitStatus;
    class AbsTrackRep;
    class RKTrackRep;
    class KalmanFittedStateOnPlane;
    class PlanarMeasurement;
}

namespace extension{
    class TrackerHit;
    class Track;
    class MutableTrack;
    class TrackCollection;
    class TrackerHitCollection;

}

namespace edm4hep{
    class MCParticle;
    class SimTrackerHitCollection;
    class ReconstructedParticle;
    class MutableReconstructedParticle;
    class MCRecoTrackerAssociationCollection;
    class SimTrackerHit;
    class Vector3d;
    class Vector3f;
    class MCParticleCollection;
}

class GenfitTrack {
    friend int GenfitFitter::processTrack(
            GenfitTrack* track, bool resort);
    
    friend int GenfitFitter::processTrack(
            genfit::Track* track, bool resort);

    friend int GenfitFitter::processTrackWithRep(
            GenfitTrack* track, int repID, bool resort);

    public:
    GenfitTrack(SmartIF<IGeomSvc> geom, const char* name="GenfitTrack");
    virtual ~GenfitTrack();

    ///Create genfit track from MCParticle
    bool createGenfitTrackFromEDM4HepTrack(int pidType,const edm4hep::Track& track,
            double eventStartTime,bool isUseCovTrack);

    double checkGenfitTrack(int& dcFitNum);

    /// ---------Add measurements---------
    ///Add one WireMeasurements
    bool addWireMeasurement(const edm4hep::TrackerHit& trackerHit,
            edm4hep::SimTrackerHit simTrackerHitAsso,int dcID,float sigma);
  
    ///Add one silicon hits
    bool addSiliconMeasurement(edm4hep::TrackerHit* hit,TVector3 pos,TVector3 mom,
            float sigmaU,float sigmaV,int cellID,int hitID);

    //return dd4hep unit
    double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca,edm4hep::MCParticle mcParticle,
            unsigned long long cellID, int repID=1,
            bool stopAtBoundary=false,
            bool calcJacobianNoise=true);
    //return dd4hep unit
    double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca,edm4hep::MCParticle mcParticle,
            HelixClass helix,unsigned long long cellID, int repID=1,
            bool stopAtBoundary=false,
            bool calcJacobianNoise=true);
    //return dd4hep unit
    double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca,TVector3 startPos,TVector3 startMom,
            int chargeId, unsigned long long cellID, int repID=1, bool stopAtBoundary=false,
            bool calcJacobianNoise=true)const;
    double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca,TVector3 startPos,TVector3 startMom,
            unsigned long long cellID,int pdg, bool stopAtBoundary=false,
            bool calcJacobianNoise=true)const;
    /// Extrapolate the track to the point
    /// Output: pos and mom of POCA point to point
    /// Input: genfitTrack,point,repID,stopAtBoundary and calcAverageState
    /// repID same with pidType
    double extrapolateToPoint(TVector3& pos, TVector3& mom, TMatrixDSym& cov,
            const TVector3& point, int repID=0,
            bool stopAtBoundary = false, bool calcJacobianNoise = true) const;

    /// Extrapolate the track to the cyliner at fixed raidus
    /// Output: pos and mom at fixed radius
    /// Input: genfitTrack, radius of cylinder at center of the origin,
    ///        repID, stopAtBoundary and calcAverageState
    double extrapolateToCylinder(TVector3& pos, TVector3& mom,
            double radius, const TVector3 linePoint,
            const TVector3 lineDirection, int hitID =0, int repID=0,
            bool stopAtBoundary=false, bool calcJacobianNoise=true);

    double extrapolateToPlane(TVector3 startPos, TVector3 startMom,
            genfit::SharedPlanePtr plane,int pdg,
            bool stopAtBoundary = false, bool calcJacobianNoise = true) const;

    bool getPosMomCovMOP(int hitID, TLorentzVector& pos, TVector3& mom,
            TMatrixDSym& cov, int repID=0) const;

    /// get the seed position and momentum from track
    const TLorentzVector getSeedStatePos() const;
    const TVector3 getSeedStateMom() const;
    void getSeedStateMom(TLorentzVector& pos, TVector3& mom) const;
    unsigned int getNumPoints() const;
    unsigned int getNumPointsWithMeasurement() const;
    unsigned int getNumPointsDet(int cellID) const;
    void getSeedStateCov(TVectorD& seedState,TMatrixDSym& covSeed);

    /// get the fitted track status
    const genfit::FitStatus* getFitStatus(int repID=0) const;
    int getFittedState(TLorentzVector& pos, TVector3& mom, TMatrixDSym& cov,
            int trackPointId=0, int repID=0, bool biased=true) const;
    int getNumPointsWithFittedInfo(int repID=0) const;
    bool getFirstPointWithFittedInfo(int repID=0) const;
    bool fitSuccess(int repID) const;

    /// get the wire infomation
    int getDetIDWithFitterInfo(int hitID, int idRaw=0) const;

    int getPDG(int id=0) const;
    int getPDGCharge(int id=0) const;

    /// print genfit track
    void printSeed() const;//print seed
    void printFitted(int repID=0) const;//print seed
    void print(TLorentzVector pos, TVector3 mom, const char* comment="") const;
    void print() const;


    bool sortedHitOnTrack(int pdg);
    void checkTrackPoint();

    /// get name of this object
    const char* getName() const {return m_name;}
    genfit::Track* getTrack() const{return m_track;}

    /// Add a track representation
    genfit::RKTrackRep* addTrackRep(int pdgType,int charge);

    /// Get a hit according to index
    GenfitHit* GetHit(long unsigned int i) const;

    private:

    /// ---------Add a Genfit track-------
    bool createGenfitTrack(int pdgType,int charge,
            TVectorD trackParam, TMatrixDSym covMInit_6);

    int getDetTypeID(unsigned long long cellID) const;
    const char* m_name;

    ///Note! private functions are using genfit unit, cm and MeV

    genfit::AbsTrackRep* getRep(int id=0) const;
    bool getMOP(int hitID, genfit::MeasuredStateOnPlane& mop,
            genfit::AbsTrackRep* trackRep=nullptr) const;
    const dd4hep::rec::ISurface* getISurface(edm4hep::TrackerHit* hit);
    void getSeedCov(TMatrixDSym& cov);
    void getAssoSimTrackerHit(
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            edm4hep::TrackerHit* trackerHit,
            edm4hep::SimTrackerHit& simTrackerHit) const;
    void getAssoSimTrackerHit2(
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            const edm4hep::TrackerHit& trackerHit,
            edm4hep::SimTrackerHit& simTrackerHit) const;
    void getEndPointsOfWire(unsigned long long cellID,TVector3& end0,TVector3& end1)const;
    void getTrackFromEDMTrack(const edm4hep::Track& edm4HepTrack,
            double& charge, TVectorD& trackParam, TMatrixDSym& cov) const;
    void getTrackFromMCPartile(const edm4hep::MCParticle mcParticle,
            TVectorD& trackParam, TMatrixDSym& cov) const;
    void getPosMomFromMCPartile(const edm4hep::MCParticle mcParticle,
            TVector3& pos,TVector3& mom) const;
    void clearGenfitHitVec();
    void getISurfaceOUV(const dd4hep::rec::ISurface* iSurface,TVector3& o,
            TVector3& u,TVector3& v,double& lengthU,double& lengthV);
    void getMeasurementAndCov(edm4hep::TrackerHit* hit,TVector3& pos,TMatrixDSym& cov);
    int getSigmas(int cellID,std::vector<float> sigmaUVec,
        std::vector<float> sigmaVVec,float& sigmaU,float& sigmaV)const;
    bool isCDCHit(edm4hep::TrackerHit* hit);
    GenfitHit* makeAGenfitHit(edm4hep::TrackerHit* trackerHit,
            edm4hep::SimTrackerHit* simTrackerHitAsso,
            double sigma,bool truthAmbig,double skipCorner,double skipNear);
    void getSortedTrackerHits(std::vector<edm4hep::TrackerHit*>& hits,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            std::vector<edm4hep::TrackerHit*>& sortedDCTrackerHits,
            int sortMethod);
    void getSortedTrackerHitsTrF(std::vector<edm4hep::TrackerHit*> trackerHits,
            std::vector<edm4hep::TrackerHit*>& sortedDCTrackerHits,
            int sortMethod=1);

    genfit::Track* m_track;/// track
    int m_debug;/// debug level
    int m_debugLocal;/// debug level local

    SmartIF<IGeomSvc> m_geomSvc;
    const GenfitField* m_genfitField;//pointer to genfit field
    const dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
    const dd4hep::DDSegmentation::BitFieldCoder* m_decoderDC;

    static const int s_PDG[2][5];

    std::vector<GenfitHit*> m_genfitHitVec;

};
#endif