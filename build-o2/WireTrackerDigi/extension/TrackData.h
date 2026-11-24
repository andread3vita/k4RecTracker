// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackDATA_H
#define EXTENSION_TrackDATA_H

#include <cstdint>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


/** @class TrackData
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class TrackData {
public:
  std::int32_t type{}; ///< flagword that defines the type of track.Bits 16-31 are used internally
  float chi2{}; ///< Chi^2 of the track fit
  std::int32_t ndf{}; ///< number of degrees of freedom of the track fit
  float dEdx{}; ///< dEdx of the track
  float dEdxError{}; ///< error of dEdx
  float radiusOfInnermostHit{}; ///< radius of the innermost hit that has been used in the track fit

  unsigned int subdetectorHitNumbers_begin{};
  unsigned int subdetectorHitNumbers_end{};
  unsigned int trackStates_begin{};
  unsigned int trackStates_end{};
  unsigned int dxQuantities_begin{};
  unsigned int dxQuantities_end{};
  unsigned int trackerHits_begin{};
  unsigned int trackerHits_end{};
  unsigned int tracks_begin{};
  unsigned int tracks_end{};
};


namespace v1 {
using TrackData = extension::TrackData;
} // namespace v1


/// @}
/// \endcond
} // namespace extension


#endif
