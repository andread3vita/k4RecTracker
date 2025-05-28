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

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit,double m_Bz) {

    edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

    double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
   
    // get extrapolated momentum from the helix with ref point at last hit
    double momAtCalorimeter[] = {0.,0.,0.};
    helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

    // produce new helix at calorimeter position
    auto helixAtCalorimeter = HelixClass_double();
    helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), m_Bz);

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