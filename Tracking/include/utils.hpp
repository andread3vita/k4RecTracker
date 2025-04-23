#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <cmath>
#include <torch/torch.h>
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>
#include "utils.hpp"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>

#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"

// Function to compute cosine of the theta angle and transverse momentum (pt)
std::vector<double> computeCosThetaAndPt(double px, double py, double pz);

std::vector<double> computeCosThetaPtAndPhi(double px, double py, double pz);

torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor unassigned, float tbeta);

torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta, float td);

std::vector<float> local_to_global(std::vector<float> local_pos,std::vector<float> x_prime,std::vector<float> y_prime,std::vector<float> z_prime,std::vector<float> wire_pos);

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag, dd4hep::Detector* mainDetector_pntr);

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit);


#endif // UTILS_HPP
