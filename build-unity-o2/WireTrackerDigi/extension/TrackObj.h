// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackOBJ_H
#define EXTENSION_TrackOBJ_H

// data model specific includes
#include "extension/TrackData.h"
#include "extension/TrackerHit.h"
#include "edm4hep/Quantity.h"
#include "edm4hep/TrackState.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class Track;

class TrackObj {
public:
  /// constructor
  TrackObj();
  /// copy constructor (does a deep-copy of relation containers)
  TrackObj(const TrackObj&);
  /// constructor from ObjectID and TrackData
  /// does not initialize the internal relation containers
  TrackObj(const podio::ObjectID& id, const TrackData& data);
  /// No assignment operator
  TrackObj& operator=(const TrackObj&) = delete;
  virtual ~TrackObj();

public:
  podio::ObjectID id{};
  TrackData data;
  std::vector<extension::TrackerHit>* m_trackerHits{nullptr};
  std::vector<extension::Track>* m_tracks{nullptr};
  std::vector<std::int32_t>* m_subdetectorHitNumbers{nullptr};
  std::vector<edm4hep::TrackState>* m_trackStates{nullptr};
  std::vector<edm4hep::Quantity>* m_dxQuantities{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
