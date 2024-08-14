#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cmath>
#include <torch/torch.h>

// Function to compute cosine of the theta angle and transverse momentum (pt)
std::vector<double> computeCosThetaAndPt(double px, double py, double pz);

std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz);


/**
 * @find_condpoints Finds and returns the indices of conditional points based on a threshold value.
 *
 * This function identifies points from the given `betas` tensor that satisfy the following criteria:
 * 1. Their beta value exceeds the specified threshold `tbeta`.
 * 2. They are among the points listed in the `unassigned` tensor.
 * 3. The resulting indices are sorted in ascending order based on their beta values.
 *
 * @param betas A tensor containing beta values for all points.
 * @param unassigned A tensor containing indices of unassigned points.
 * @param tbeta A float representing the threshold value that beta values must exceed.
 *
 * @return A tensor containing the sorted indices of the conditional points that are both
 *         unassigned and exceed the threshold value `tbeta`.
 *
 * The function performs the following steps:
 * 1. Creates a mask to identify points where beta values are greater than `tbeta`.
 * 2. Updates the mask to reflect which of these points are unassigned.
 * 3. Filters the conditional points to include only those that are unassigned.
 * 4. Extracts the indices of the filtered points.
 * 5. Sorts these indices based on their beta values in ascending order.
 * 6. Returns the sorted indices of the conditional points.
 */
torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor unassigned, float tbeta);

/**
 * @get_clustering Performs clustering based on beta values and spatial distances.
 *
 * This function performs clustering of points based on their beta values and spatial distances in 3D space.
 * The points are clustered based on a specified beta threshold and a distance threshold. The function 
 * assigns cluster indices to the points and returns a tensor representing these cluster assignments.
 *
 * @param output_vector A vector of floats containing the data for the points. Each point has 4 values: 
 *                      x, y, z (coordinates), and beta (value).
 * @param num_rows The number of rows (points) in the `output_vector`.
 * @param tbeta A float representing the threshold value for beta. Points with beta values greater than this 
 *               threshold are considered for clustering.
 * @param td A float representing the distance threshold. Points within this distance from a conditional point 
 *           are assigned to the same cluster.
 *
 * @return A tensor of size `(num_rows)` containing the cluster assignments for each point. Points that are not 
 *         assigned to any cluster will have a value of `-1`.
 *
 * The function performs the following steps:
 * 1. Converts the input vector into a 2D Torch tensor with shape `(num_rows, 4)`, where the last column contains 
 *    beta values and the first three columns contain 3D coordinates.
 * 2. Extracts the beta values and the 3D coordinates from the tensor.
 * 3. Initializes the clustering tensor to `-1` (indicating unassigned points) and sets up a list of unassigned points.
 * 4. Iteratively assigns points to clusters based on their distances to the current conditional point and the beta 
 *    value threshold.
 * 5. Updates the list of unassigned points and finds new conditional points for the next cluster.
 * 6. Increments the cluster index for each new cluster found.
 * 7. Returns the tensor with cluster assignments for all points.
 */
torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows, float tbeta, float td);

#endif // UTILS_HPP
