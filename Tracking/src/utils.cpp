#include "utils.hpp"

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag) {

    dd4hep::rec::LayeredCalorimeterData * theExtension = 0;
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    if (detector.detectors().empty()) {
        throw std::runtime_error("Detector is not initialized: no DetElements found.");
    }

    const std::vector< dd4hep::DetElement>& theDetectors = dd4hep::DetectorSelector(detector).detectors(  includeFlag, excludeFlag );
    if( theDetectors.size()  != 1 ){

      std::stringstream es ;
      es << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
        << " --- found detectors : " ;
      for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
        es << theDetectors.at(i).name() << ", " ;
      }
      throw std::runtime_error( es.str() ) ;
    }

    theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

    return theExtension;
}

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit,double Bz) {

    edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

    double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
   
    // get extrapolated momentum from the helix with ref point at last hit
    double momAtCalorimeter[] = {0.,0.,0.};
    helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

    // produce new helix at calorimeter position
    auto helixAtCalorimeter = HelixClass_double();
    helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), Bz);

    // fill the TrackState parameters
    trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
    trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
    trackState_AtCalorimeter.phi = std::atan2(momAtCalorimeter[1], momAtCalorimeter[0]);
    trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
    trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
    trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
    trackState_AtCalorimeter.referencePoint = edm4hep::Vector3f(posAtCalorimeter[0],
                                                                posAtCalorimeter[1],
                                                                posAtCalorimeter[2]);
    return trackState_AtCalorimeter;
}

int getHypotesisCharge(int pdg) {
    if (pdg == 11)
    {
      return -1;
    }
    else if (pdg == -11)
    {
      return 1;
    }
    else if (pdg == 13)
    {
      return -1;
    }
    else if (pdg == -13)
    {
      return 1;
    }
    else if (pdg == 211)
    {
      return 1;
    }
    else if (pdg == -211)
    {
      return -1;
    }
    else if (pdg == 321)
    {
      return 1;
    }
    else if (pdg == -321)
    {
      return -1;
    }
    else if (pdg == 2212)
    {
      return 1;
    }
    else if (pdg == -2212)
    {
      return -1;
    }
  
  return 0; // Default case, should not happen
}


edm4hep::TrackState extrapolateToCalorimeter(   double posAtLastHit[3],
                                                double momAtLastHit[3],
                                                double m_eCalBarrelInnerR,
                                                double m_eCalBarrelMaxZ,
                                                double m_eCalEndCapInnerR,
                                                double m_eCalEndCapOuterR,
                                                double m_eCalEndCapInnerZ,
                                                double m_eCalEndCapOuterZ,
                                                double Bz,
                                                double charge) {
  
  double c_light = 2.99792458e8;
  double a = c_light * 1e3 * 1e-15;

  auto helixAtLastHit = HelixClass_double();
  helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, charge, Bz);

  // Propagation to Endcap
  if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
    float minGenericTime(std::numeric_limits<float>::max());
                        
    // create helix to project
    // rather than using parameters at production, better to use those from
    // last hit
    pandora::CartesianVector pos_lasthit(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
    pandora::CartesianVector mom_lasthit(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);

    const pandora::Helix helix(pos_lasthit, mom_lasthit, charge ,Bz);
    const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
    const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
                        
    // First project to endcap
    pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
    bool hasEndCapProjection(false);
    if (m_eCalEndCapInnerR>0) {
      float genericTime(std::numeric_limits<float>::max());
      const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint, endCapProjection, genericTime));
      float x = endCapProjection.GetX();
      float y = endCapProjection.GetY();
      float r = std::sqrt(x*x+y*y);
      if (
          (pandora::STATUS_CODE_SUCCESS == statusCode) &&
          (genericTime < minGenericTime) &&
          (r >= m_eCalEndCapInnerR) &&
          (r <= m_eCalEndCapOuterR)
      ) {
          minGenericTime = genericTime;
          bestECalProjection = endCapProjection;
          hasEndCapProjection = true;
        }
    }
                                    
                                    
    // Then project to barrel surface(s), and keep projection
    // if extrapolation is within the z acceptance of the detector
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
    bool hasBarrelProjection = false;
    if (m_eCalBarrelInnerR>0) {

      float genericTime(std::numeric_limits<float>::max());
      const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
                            
      if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (std::fabs(barrelProjection.GetZ())<= m_eCalBarrelMaxZ)) {
        hasBarrelProjection = true;
        if (genericTime < minGenericTime) {
          minGenericTime = genericTime;
          secondBestECalProjection = bestECalProjection;
          bestECalProjection = barrelProjection;
        }
        else {
         secondBestECalProjection = barrelProjection;
        }
      }
    }
           
    // store extrapolation to calo
    // by default, store extrapolation with lower arrival time
    // get extrapolated position
    edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit,Bz);
    // double omega_lastHit = trackState_AtCalorimeter.omega;
    // double pt_lasthit = a * Bz / abs(omega_lastHit);
    // double phi_lasthit = trackState_AtCalorimeter.phi;
    // double pz_lasthit = trackState_AtCalorimeter.tanLambda * pt_lasthit;
    // double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
    // double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
    // auto ref_lastHit = trackState_AtCalorimeter.referencePoint;
    // attach the TrackState to the track

    return trackState_AtCalorimeter;      
  }
  else
  {
    std::cerr << "Error: Calorimeter parameters are not properly defined!" << std::endl;
    std::exit(EXIT_FAILURE);
  }                                  

}