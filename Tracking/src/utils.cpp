#include "utils.hpp"

// Function to compute cos(theta) and transverse momentum (pt) from momentum components
std::vector<double> computeCosThetaAndPt(double px, double py, double pz) {    
    // Calculate the magnitude of the momentum vector P
    double P = std::sqrt(px * px + py * py + pz * pz);    
    
    // Calculate the absolute value of cos(theta), where theta is the angle between the 
    // momentum vector and the z-axis
    double cosTheta = std::abs(pz / P);    
    
    // Calculate the transverse momentum pt (momentum in the xy-plane)
    double pt = std::sqrt(px * px + py * py);        
    
    // Return cosTheta and pt as a vector of doubles
    return {cosTheta, pt};
}

// Function to compute cos(theta), transverse momentum (pt), and azimuthal angle (phi) from momentum components
std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz) {
    // Calculate the magnitude of the momentum vector P
    double P = std::sqrt(px * px + py * py + pz * pz);
    
    // Calculate the absolute value of cos(theta)
    double cosTheta = std::abs(pz / P);
    
    // Calculate the transverse momentum pt
    double pt = std::sqrt(px * px + py * py);
    
    // Calculate the azimuthal angle phi (angle in the xy-plane)
    double phi = std::atan2(py, px);

    // Return cosTheta, pt, and phi as a vector of doubles
    return {cosTheta, pt, phi};
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
    // auto clustering = -1 * torch::ones({n_points}, torch::kLong);
    // int index_assignation = 0;
    // auto unassigned = torch::arange(n_points);
    // while (indices_condpoints.size(0) > 0 && unassigned.size(0) > 0) {
    //     auto index_condpoint = indices_condpoints.index({0});
    //     auto d = (X.index({unassigned}) - X.index({index_condpoint})).norm(2, -1).squeeze();
    //     auto mask_distance = d.lt(td);
    //     auto  mask_distance_ind = mask_distance.nonzero();
    //     auto assigned_to_this_condpoint = unassigned.index({mask_distance_ind});
    //     clustering.index_put_({assigned_to_this_condpoint}, index_assignation);
    //     auto mask_distance_out = d.ge(td);
    //     auto  mask_distance_ind_out = mask_distance_out.squeeze(0).nonzero();
    //     unassigned = unassigned.index({mask_distance_ind_out});
    //     indices_condpoints = find_condpoints(betas, unassigned, tbeta);
    //     index_assignation += 1;
    // }
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

std::vector<float> local_to_global(std::vector<float> local_pos,std::vector<float> x_prime,std::vector<float> y_prime,std::vector<float> z_prime,std::vector<float> wire_pos)
{

   
    std::vector<float> global_pos(3, 0.0f);
    global_pos[0] = x_prime[0]*local_pos[0] +x_prime[1]*local_pos[1] + x_prime[2]*local_pos[2] + wire_pos[0];
    global_pos[1] = y_prime[0]*local_pos[0] +y_prime[1]*local_pos[1] + y_prime[2]*local_pos[2] + wire_pos[1];
    global_pos[2] = z_prime[0]*local_pos[0] +z_prime[1]*local_pos[1] + z_prime[2]*local_pos[2] + wire_pos[2];

    return global_pos;

}

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit) {

    edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

    double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
    
    // get extrapolated momentum from the helix with ref point at last hit
    double momAtCalorimeter[] = {0.,0.,0.};
    helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

    // produce new helix at calorimeter position
    auto helixAtCalorimeter = HelixClass_double();
    helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), 2.);

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

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag, dd4hep::Detector* mainDetector_pntr) {
  
    dd4hep::Detector& mainDetector = *mainDetector_pntr;

    dd4hep::rec::LayeredCalorimeterData * theExtension = 0;

    // dd4hep::Detector & mainDetector = *(mGeo->getDetector());
    const std::vector< dd4hep::DetElement>& theDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );
    theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

    return theExtension;
}
