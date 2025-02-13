#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cmath>
#include <torch/torch.h>

// Function to compute cosine of the theta angle and transverse momentum (pt)
std::vector<double> computeCosThetaAndPt(double px, double py, double pz);

std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz);

#endif // UTILS_HPP
