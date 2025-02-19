#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cmath>
#include <torch/torch.h>

// Function to compute cosine of the theta angle and transverse momentum (pt)
std::vector<double> computeCosThetaAndPt(double px, double py, double pz);

std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz);

torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor unassigned, float tbeta);

torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta, float td);

std::vector<float> local_to_global(std::vector<float> local_pos,std::vector<float> x_prime,std::vector<float> y_prime,std::vector<float> z_prime,std::vector<float> wire_pos);

#endif // UTILS_HPP
