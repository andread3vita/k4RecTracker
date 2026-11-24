// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiDATA_H
#define EXTENSION_DriftChamberDigiDATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class DriftChamberDigiData
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiData {
public:
  std::uint64_t cellID{}; ///< ID of the wire that created this hit
  ::edm4hep::Vector3d leftPosition{}; ///< position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  ::edm4hep::Vector3d rightPosition{}; ///< position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  float time{}; ///< time of the hit [ns].
  float eDep{}; ///< energy deposited on the hit [GeV].
  float eDepError{}; ///< error measured on eDep [GeV].
  std::uint32_t clusterCount{}; ///< number of clusters associated to this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  ::edm4hep::Vector3d position{}; ///< hit position [mm]

};


namespace v1 {
using DriftChamberDigiData = extension::DriftChamberDigiData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
