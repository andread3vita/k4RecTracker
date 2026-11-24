// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2DATA_H
#define EXTENSION_DriftChamberDigiV2DATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiV2Data
 *  Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).
 *  @author: A. Tolosa-Delgado, B. Francois, CERN
 */
class DriftChamberDigiV2Data {
public:
  std::uint64_t cellID{}; ///< ID of the sensor that created this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  float time{}; ///< time of the hit [ns]
  float eDep{}; ///< energy deposited on the hit [GeV]
  float eDepError{}; ///< error measured on eDep [GeV]
  ::edm4hep::Vector3d position{}; ///< point on the sensitive wire (SW) which is closest to the simhit [mm]
  ::edm4hep::Vector3d directionSW{}; ///< direction of SW
  float distanceToWire{}; ///< distance hit-wire [mm]
  std::uint32_t nCluster{}; ///< number of clusters associated to this hit

  unsigned int nElectrons_begin{};
  unsigned int nElectrons_end{};
};


namespace v1 {
using DriftChamberDigiV2Data = extension::DriftChamberDigiV2Data;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
