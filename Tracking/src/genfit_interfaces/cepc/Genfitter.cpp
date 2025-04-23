#include "GenfitFitter.hpp"
#include "GenfitTrack.hpp" //TODO

//Gaudi
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ISvcLocator.h"

//External
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

//genfit
#include <Track.h>
#include <Exception.h>
#include <FieldManager.h>
#include <TGeoMaterialInterface.h>
#include <TGeoManager.h>
#include <MaterialEffects.h>
#include <MeasuredStateOnPlane.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <DAF.h>
#include <AbsKalmanFitter.h>
#include <KalmanFitterInfo.h>

//ROOT
#include <TVector3.h>
#include <TGeoManager.h>

//STL
#include <iostream>
#include <string>
#include <string.h>


GenfitFitter::~GenfitFitter(){
    delete m_absKalman;
}

GenfitFitter::GenfitFitter(const char* type):
    m_absKalman(nullptr)
    ,m_fitterType(type)
    ,m_minIterations(4)
    ,m_maxIterations(10)
    ,m_deltaPval(1e-3)
    ,m_relChi2Change(0.2)
    ,m_blowUpFactor(500)
    ,m_resetOffDiagonals(true)
    ,m_blowUpMaxVal(1.e6)
    ,m_multipleMeasurementHandling(genfit::unweightedClosestToPredictionWire)
    ,m_maxFailedHits(-1)
    ,m_deltaWeight(1e-3)
    ,m_annealingBetaStart(100)
    ,m_annealingBetaStop(0.01)
    ,m_annealingNSteps(0.01)
    ,m_noEffects(false)
    ,m_energyLossBetheBloch(true)
    ,m_noiseBetheBloch(true)
    ,m_noiseCoulomb(true)
    ,m_energyLossBrems(false)
    ,m_noiseBrems(false)
    ,m_ignoreBoundariesBetweenEqualMaterials(true)
{
    
    init();
}

/// initialize genfit fitter, old fitter will be deleted
int GenfitFitter::init(bool deleteOldFitter)
{
    if(deleteOldFitter && m_absKalman) delete m_absKalman;

    if (m_fitterType=="DAFRef") {

        m_absKalman = new genfit::DAF(true,getDeltaPval(),getConvergenceDeltaWeight());
    }
    else if (m_fitterType=="DAF") {
        m_absKalman = new genfit::DAF(false,getDeltaPval(),getConvergenceDeltaWeight());
    }
    else if (m_fitterType=="KalmanFitter") {
        m_absKalman = new genfit::KalmanFitter(getMaxIterations());
    }
    else if (m_fitterType=="KalmanFitterRefTrack") {
        m_absKalman = new genfit::KalmanFitterRefTrack(getMaxIterations());
    }
    else {
        m_absKalman = nullptr;
        return -1;
    }

    return 0;
}

/// Fit a track from a candidate track
int GenfitFitter::processTrackWithRep(GenfitTrack* track,int repID,bool resort)
{

    if(track->getNumPoints()<=0){
        return false;
    }
    try{
        m_absKalman->processTrackWithRep(track->getTrack(), track->getRep(repID),resort);
    }catch(genfit::Exception& e){
      
        return false;
    }
    return true;
}

/// Fit a track from a candidate track
int GenfitFitter::processTrack(GenfitTrack* track, bool resort)
{
  
    if(track->getNumPoints()<=0){
        return false;
    }

    try{
        m_absKalman->processTrack(track->getTrack(),resort);
    }catch(genfit::Exception& e){
        std::cout<<"Genfit exception caught "<<std::endl;
        return false;
    }
    
    return true;
}

/// Fit a genfit::track from a candidate track
int GenfitFitter::processTrack(genfit::Track* track, bool resort)
{
    if(track->getNumPoints()<=0){
        return false;
    }

    try{
        m_absKalman->processTrack(track,resort);
    }catch(genfit::Exception& e){
        std::cout<<"Genfit exception caught "<<std::endl;
        return false;
    }
    
    return true;
} // End of Process genfit::Track


void GenfitFitter::setFitterType(const char* val)
{
    std::string oldFitterType=m_fitterType;
    m_fitterType = val;
    init(oldFitterType==val);
}

GenfitFitter& GenfitFitter::operator=(
        const GenfitFitter& r)
{
    m_fitterType          = r.m_fitterType;
    m_minIterations       = r.m_minIterations;
    m_maxIterations       = r.m_maxIterations;
    m_deltaPval           = r.m_deltaPval;
    m_relChi2Change       = r.m_relChi2Change;
    m_blowUpFactor        = r.m_blowUpFactor;
    m_resetOffDiagonals   = r.m_resetOffDiagonals;
    m_blowUpMaxVal        = r.m_blowUpMaxVal;
    m_multipleMeasurementHandling = r.m_multipleMeasurementHandling;
    m_maxFailedHits       = r.m_maxFailedHits;
    m_annealingNSteps     = r.m_annealingNSteps;
    m_deltaWeight         = r.m_deltaWeight;
    m_annealingBetaStart  = r.m_annealingBetaStart;
    m_annealingBetaStop   = r.m_annealingBetaStop;
    m_noEffects           = r.m_noEffects;
    m_energyLossBetheBloch= r.m_energyLossBetheBloch;
    m_noiseBetheBloch     = r.m_noiseBetheBloch;
    m_noiseCoulomb        = r.m_noiseCoulomb;
    m_energyLossBrems     = r.m_energyLossBrems;
    m_noiseBrems          = r.m_noiseBrems;
    m_ignoreBoundariesBetweenEqualMaterials
        = r.m_ignoreBoundariesBetweenEqualMaterials;
    m_mscModelName        = r.m_mscModelName;
    return *this;
}


///Setters of AbsKalmanFitter
void GenfitFitter::setMinIterations(unsigned int val)
{
    m_absKalman->setMinIterations(val);
    m_minIterations = val;
}

void GenfitFitter::setMaxIterations(unsigned int val)
{
    m_absKalman->setMaxIterations(val);
    m_maxIterations = val;
}

void GenfitFitter::setMaxIterationsBetas(double bStart,double bFinal,
        unsigned int val)
{
    m_absKalman->setMaxIterations(val);
    m_maxIterations = val;
    genfit::DAF* daf = dynamic_cast<genfit::DAF*> (m_absKalman);
    daf->setAnnealingScheme(bStart,bFinal,val);
}

void GenfitFitter::setDeltaPval(double val)
{
    m_absKalman->setDeltaPval(val);
    m_deltaPval = val;
}

void GenfitFitter::setRelChi2Change(double val)
{
    m_absKalman->setRelChi2Change(val);
    m_relChi2Change = val;
}

void GenfitFitter::setBlowUpFactor(double val)
{
    m_absKalman->setBlowUpFactor(val);
    m_blowUpFactor = val;
    if (m_fitterType=="DAFRef" || m_fitterType=="DAF") {
        getDAF()->getKalman()->setBlowUpFactor(m_blowUpFactor);
    }
}

void GenfitFitter::setResetOffDiagonals(bool val)
{
    m_absKalman->setResetOffDiagonals(val);
    m_resetOffDiagonals = val;
    if (m_fitterType=="DAFRef" || m_fitterType=="DAF") {
        getDAF()->getKalman()->setResetOffDiagonals(m_resetOffDiagonals);
    }
}

void GenfitFitter::setBlowUpMaxVal(double val)
{
    m_absKalman->setBlowUpMaxVal(val);
    m_blowUpMaxVal = val;
    if (m_fitterType=="DAFRef" || m_fitterType=="DAF") {
        getDAF()->getKalman()->setBlowUpMaxVal(m_blowUpMaxVal);
    }
}

void GenfitFitter::setMultipleMeasurementHandling(
        genfit::eMultipleMeasurementHandling val)
{
    m_absKalman->setMultipleMeasurementHandling(val);
    m_multipleMeasurementHandling = val;
}

void GenfitFitter::setMaxFailedHits(int val)
{
    m_absKalman->setMaxFailedHits(val);
    m_maxFailedHits = val;
}

///DAF setters
void GenfitFitter::setConvergenceDeltaWeight(double val)
{
    genfit::DAF* daf = getDAF();
    if(nullptr != daf) daf->setConvergenceDeltaWeight(val);
    m_deltaWeight = val;
}
void GenfitFitter::setAnnealingScheme(
        double bStart, double bFinal, unsigned int nSteps)
{
    genfit::DAF* daf = getDAF();
    if(nullptr != daf) daf->setAnnealingScheme(bStart, bFinal, nSteps);
    m_annealingBetaStart = bStart;
    m_annealingBetaStop = bFinal;
    m_annealingNSteps = nSteps;
}

///Material effects
void GenfitFitter::setNoEffects(bool val)
{
    genfit::MaterialEffects::getInstance()->setNoEffects();
    m_noEffects = val;
}

void GenfitFitter::setEnergyLossBetheBloch(bool val)
{
    genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(val);
    m_energyLossBetheBloch = val;
}

void GenfitFitter::setNoiseBetheBloch(bool val)
{
    genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(val);
    m_noiseBetheBloch = val;
}

void GenfitFitter::setNoiseCoulomb(bool val)
{
    genfit::MaterialEffects::getInstance()->setNoiseCoulomb(val);
    m_noiseCoulomb = val;
}

void GenfitFitter::setEnergyLossBrems(bool val)
{
    genfit::MaterialEffects::getInstance()->setEnergyLossBrems(val);
    m_energyLossBrems = val;
}

void GenfitFitter::setNoiseBrems(bool val)
{
    genfit::MaterialEffects::getInstance()->setNoiseBrems(val);
    m_noiseBrems = val;
}

void GenfitFitter::setIgnoreBoundariesBetweenEqualMaterials(bool val)
{
    genfit::MaterialEffects::getInstance()->
        ignoreBoundariesBetweenEqualMaterials(val);
    m_ignoreBoundariesBetweenEqualMaterials = val;
}

void GenfitFitter::setMscModelName(std::string val)
{
    genfit::MaterialEffects::getInstance()->setMscModel(val);
    m_mscModelName = val;
}

genfit::DAF* GenfitFitter::getDAF()
{
    genfit::DAF* daf = nullptr;
    try{
        daf = dynamic_cast<genfit::DAF*> (m_absKalman);
    }catch(...){

        return nullptr;
    }
    return daf;
}

genfit::KalmanFitterRefTrack* GenfitFitter::getKalRef()
{
    genfit::KalmanFitterRefTrack* ref=nullptr;
    try{
        ref = dynamic_cast<genfit::KalmanFitterRefTrack*> (m_absKalman);
    }catch(...){
        return nullptr;
    }
    return ref;
}

void GenfitFitter::SetRunEvent(int event){
    m_absKalman->SetRunEvent(event);
}