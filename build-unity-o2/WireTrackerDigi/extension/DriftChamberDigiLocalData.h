// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocalDATA_H
#define EXTENSION_DriftChamberDigiLocalDATA_H

#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiLocalData
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiLocalData {
public:
  std::uint64_t cellID{}; ///< ID of the wire that created this hit
  float distanceToWire{}; ///< smeared distance of closest approach between the wire and the hit [mm]
  float zPositionAlongWire{}; ///< smeared z position in the local wire coordinate system [mm]
  float time{}; ///< time of the hit [ns].
  float eDep{}; ///< energy deposited on the hit [GeV].
  float eDepError{}; ///< error measured on eDep [GeV].
  std::uint32_t clusterCount{}; ///< number of clusters associated to this hit

};


namespace v1 {
using DriftChamberDigiLocalData = extension::DriftChamberDigiLocalData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
