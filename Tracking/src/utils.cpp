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

