#ifndef UTILS_HPP
#define UTILS_HPP

//=== Standard Library ===
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <limits>

//=== DD4hep / DDRec ===
#include "DD4hep/DetType.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include <DDRec/DetectorData.h>

//=== edm4hep ===
#include "edm4hep/TrackState.h"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>
#include <torch/torch.h>

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag);

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit,double m_Bz);

int getHypotesisCharge(int pdg);

edm4hep::TrackState extrapolateToCalorimeter(   double posAtLastHit[3],
                                                double momAtLastHit[3],
                                                double m_eCalBarrelInnerR,
                                                double m_eCalBarrelMaxZ,
                                                double m_eCalEndCapInnerR,
                                                double m_eCalEndCapOuterR,
                                                double m_eCalEndCapInnerZ,
                                                double m_eCalEndCapOuterZ,
                                                double Bz,
                                                double charge);

torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor unassigned, float tbeta);

torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta, float td);


#endif // UTILS_HPP