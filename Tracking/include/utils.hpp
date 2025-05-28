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

dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag);

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit,double m_Bz);

int getHypotesisCharge(int pdg);

#endif // UTILS_HPP