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



void FillTrackWithCalorimeterExtrapolation(
    extension::MutableTrack& edm4hep_track,
    double m_Bz,
    int charge,
    double a,
    double m_eCalBarrelInnerR,
    double m_eCalBarrelMaxZ,
    double m_eCalEndCapInnerR,
    double m_eCalEndCapOuterR,
    double m_eCalEndCapInnerZ
) {

  auto trackStateLastHit = edm4hep_track.getTrackStates()[2];
  double omega_lastHit = trackStateLastHit.omega;
  double pt_lasthit = a * m_Bz / abs(omega_lastHit);
  double phi_lasthit = trackStateLastHit.phi;
  double pz_lasthit = trackStateLastHit.tanLambda * pt_lasthit;
  double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
  double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
  auto ref_lastHit = trackStateLastHit.referencePoint;

  // produce new helix at last hit position
  double posAtLastHit[] = {ref_lastHit[0], ref_lastHit[1], ref_lastHit[2]};
  double momAtLastHit[] = {px_lasthit, py_lasthit, pz_lasthit};
  auto helixAtLastHit = HelixClass_double();
  helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, charge, m_Bz);

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

    const pandora::Helix helix(pos_lasthit, mom_lasthit, charge ,m_Bz);
    const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
    const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
                        
    // First project to endcap
    pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
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
        }
    }
                                    
                                    
    // Then project to barrel surface(s), and keep projection
    // if extrapolation is within the z acceptance of the detector
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
    if (m_eCalBarrelInnerR>0) {

      float genericTime(std::numeric_limits<float>::max());
      const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
      
      if (
          (pandora::STATUS_CODE_SUCCESS == statusCode) &&
          (std::fabs(barrelProjection.GetZ())<= m_eCalBarrelMaxZ)
      ) {
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
    edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit,m_Bz);
    omega_lastHit = trackState_AtCalorimeter.omega;
    pt_lasthit = a * m_Bz / abs(omega_lastHit);
    phi_lasthit = trackState_AtCalorimeter.phi;
    pz_lasthit = trackState_AtCalorimeter.tanLambda * pt_lasthit;
    px_lasthit = pt_lasthit * std::cos(phi_lasthit);
    py_lasthit = pt_lasthit * std::sin(phi_lasthit);
    ref_lastHit = trackState_AtCalorimeter.referencePoint;
    // attach the TrackState to the track
    edm4hep_track.addToTrackStates(trackState_AtCalorimeter);
  
  }
}
torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor unassigned, float tbeta) {
    int n_points = unassigned.size(0);
    int size_b = betas.size(0);
    auto select_condpoints = betas.gt(tbeta);
    auto mask_unassigned = torch::zeros({size_b}, torch::dtype(torch::kBool));
    for (int i = 0; i < n_points-1; ++i) {
        auto ii = unassigned.index({i});
        mask_unassigned.index_put_({ii.squeeze()}, true);
    }
    select_condpoints = mask_unassigned * select_condpoints;
    auto indices_condpoints = select_condpoints.nonzero();
    auto betas_condpoints = -betas.index({indices_condpoints});
    auto sorted_indices = torch::argsort(betas_condpoints, /*dim=*/0, /*descending=*/false);
    indices_condpoints = indices_condpoints.index({sorted_indices});
    return indices_condpoints;
}

// Main clustering function
torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta, float td) {

    torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {num_rows,4}, torch::kFloat32);
    auto rows_output_model = torch::arange(0, output_model_tensor.size(0), torch::kLong);
    auto coord_ind_output_model = torch::arange(0, 3, torch::kLong);
    auto betas = output_model_tensor.index({rows_output_model,3});
    auto X = output_model_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(0, 3)});
    int n_points = betas.size(0);
    auto select_condpoints = betas.gt(tbeta);
    auto indices_condpoints = find_condpoints(betas, torch::arange(n_points), tbeta);
    
    auto clustering = torch::zeros({n_points}, torch::kLong);
    int index_assignation = 1;
    auto unassigned = torch::arange(n_points);
    while (indices_condpoints.size(0) > 0 && unassigned.size(0) > 0) {
        auto index_condpoint = indices_condpoints.index({0});
        auto d = (X.index({unassigned}) - X.index({index_condpoint})).norm(2, -1).squeeze();
        auto mask_distance = d.lt(td);
        auto  mask_distance_ind = mask_distance.nonzero();
        auto assigned_to_this_condpoint = unassigned.index({mask_distance_ind});
        clustering.index_put_({assigned_to_this_condpoint}, index_assignation);
        auto mask_distance_out = d.ge(td);
        auto  mask_distance_ind_out = mask_distance_out.squeeze(0).nonzero();
        unassigned = unassigned.index({mask_distance_ind_out});
        indices_condpoints = find_condpoints(betas, unassigned, tbeta);
        index_assignation += 1;
    }
    
    return clustering;
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