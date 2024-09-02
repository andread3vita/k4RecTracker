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
    
    int n_points = unassigned.size(0); // Number of unassigned points
    int size_b = betas.size(0); // Total number of points
    auto select_condpoints = betas.gt(tbeta); // Create a mask of points where beta > tbeta
    auto mask_unassigned = torch::zeros({size_b}, torch::dtype(torch::kBool)); // Mask for unassigned points

    // Mark the unassigned points in the mask
    for (int i = 0; i < n_points-1; ++i) {
        auto ii = unassigned.index({i});
        mask_unassigned.index_put_({ii.squeeze()}, true);
    }

    // Update select_condpoints to only include unassigned points
    select_condpoints = mask_unassigned * select_condpoints;

    // Get the indices of the conditional points
    auto indices_condpoints = select_condpoints.nonzero();
    
    // Get the negative beta values at the conditional points (for sorting)
    auto betas_condpoints = -betas.index({indices_condpoints});
    
    // Sort the conditional points based on their beta values in ascending order
    auto sorted_indices = torch::argsort(betas_condpoints, /*dim=*/0, /*descending=*/false);
    indices_condpoints = indices_condpoints.index({sorted_indices});

    // Return the sorted indices of the conditional points
    return indices_condpoints;
}

torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta, float td) {
    // Convert the input vector to a Torch tensor with shape (num_rows, 4)
    torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {num_rows, 4}, torch::kFloat32);

    // Create index tensors for rows and coordinates
    auto rows_output_model = torch::arange(0, output_model_tensor.size(0), torch::kLong);
    auto coord_ind_output_model = torch::arange(0, 3, torch::kLong);

    // Extract the beta values from the last column of the tensor
    auto betas = output_model_tensor.index({rows_output_model, 3});

    // Extract the 3D coordinates (x, y, z) from the first three columns
    auto X = output_model_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(0, 3)});

    int n_points = betas.size(0); // Number of points
    auto select_condpoints = betas.gt(tbeta); // Initial selection of conditional points
    auto indices_condpoints = find_condpoints(betas, torch::arange(n_points), tbeta); // Find initial conditional points

    // Initialize clustering tensor with -1 (indicating unassigned points)
    auto clustering = -1 * torch::ones({n_points}, torch::kLong);
    int index_assignation = 0; // Index for cluster assignment
    auto unassigned = torch::arange(n_points); // List of unassigned points

    // Main clustering loop
    while (indices_condpoints.size(0) > 0 && unassigned.size(0) > 0) {
        auto index_condpoint = indices_condpoints.index({0}); // Get the first conditional point
        auto d = (X.index({unassigned}) - X.index({index_condpoint})).norm(2, -1).squeeze(); // Compute distances to the conditional point

        // Create a mask for points within the distance threshold td
        auto mask_distance = d.lt(td);
        auto mask_distance_ind = mask_distance.nonzero();

        // Assign points within the threshold to the current cluster
        auto assigned_to_this_condpoint = unassigned.index({mask_distance_ind});
        clustering.index_put_({assigned_to_this_condpoint}, index_assignation);

        // Update the list of unassigned points
        auto mask_distance_out = d.ge(td);
        auto mask_distance_ind_out = mask_distance_out.squeeze(0).nonzero();
        unassigned = unassigned.index({mask_distance_ind_out});

        // Find the next set of conditional points
        indices_condpoints = find_condpoints(betas, unassigned, tbeta);
        index_assignation += 1; // Increment the cluster index
    }

    // Return the clustering result
    return clustering;
}

std::string print_shape1(const std::vector<std::int64_t>& v) {
    std::stringstream ss("");
    for (std::size_t i = 0; i < v.size() - 1; i++) ss << v[i] << "x";
    ss << v[v.size() - 1];
    return ss.str();
}
