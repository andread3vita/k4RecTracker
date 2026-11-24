// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkOBJ_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkOBJ_H

// data model specific includes
#include "extension/SenseWireHitSimTrackerHitLinkData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class SenseWireHit;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class SenseWireHitSimTrackerHitLink;

class SenseWireHitSimTrackerHitLinkObj {
public:
  /// constructor
  SenseWireHitSimTrackerHitLinkObj();
  /// copy constructor (does a deep-copy of relation containers)
  SenseWireHitSimTrackerHitLinkObj(const SenseWireHitSimTrackerHitLinkObj&);
  /// constructor from ObjectID and SenseWireHitSimTrackerHitLinkData
  /// does not initialize the internal relation containers
  SenseWireHitSimTrackerHitLinkObj(const podio::ObjectID& id, const SenseWireHitSimTrackerHitLinkData& data);
  /// No assignment operator
  SenseWireHitSimTrackerHitLinkObj& operator=(const SenseWireHitSimTrackerHitLinkObj&) = delete;
  virtual ~SenseWireHitSimTrackerHitLinkObj();

public:
  podio::ObjectID id{};
  SenseWireHitSimTrackerHitLinkData data;
  std::unique_ptr<extension::SenseWireHit> m_from{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_to{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
