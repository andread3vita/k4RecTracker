// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitDATA_H
#define EXTENSION_SenseWireHitDATA_H

#include "edm4hep/Vector3d.h"
#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class SenseWireHitData
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class SenseWireHitData {
public:
  std::uint64_t cellID{}; ///< ID of the sensor that created this hit
  std::int32_t type{}; ///< type of the raw data hit
  std::int32_t quality{}; ///< quality bit flag of the hit
  float time{}; ///< time of the hit [ns]
  float eDep{}; ///< energy deposited by the hit [GeV]
  float eDepError{}; ///< error on eDep [GeV]
  float wireStereoAngle{}; ///< angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float wireAzimuthalAngle{}; ///< azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  ::edm4hep::Vector3d position{}; ///< point on the sense wire which is closest to the hit (center of the circle) [mm]
  double positionAlongWireError{}; ///< error on the hit position along the wire direction [mm]
  float distanceToWire{}; ///< distance between the hit and the wire (radius of the circle) [mm]
  float distanceToWireError{}; ///< error on distanceToWire [mm]

  unsigned int nElectrons_begin{};
  unsigned int nElectrons_end{};
};


namespace v1 {
using SenseWireHitData = extension::SenseWireHitData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
