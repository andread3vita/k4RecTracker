#include "utils.hpp"


std::vector<double> computeCosThetaAndPt(double px, double py, double pz) {    
    
    double P = std::sqrt(px * px + py * py + pz * pz);    
    double cosTheta = std::abs(pz / P);    
    double pt = std::sqrt(px * px + py * py);        
    
    return {cosTheta, pt};
}

std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz) {
    double P = std::sqrt(px * px + py * py + pz * pz);
    double cosTheta = std::abs(pz / P);
    double pt = std::sqrt(px * px + py * py);
    double phi = std::atan2(py, px);

    return {cosTheta, pt, phi};
}
