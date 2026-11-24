// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2AssociationOBJ_H
#define EXTENSION_MCRecoDriftChamberDigiV2AssociationOBJ_H

// data model specific includes
#include "extension/MCRecoDriftChamberDigiV2AssociationData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class DriftChamberDigiV2;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class MCRecoDriftChamberDigiV2Association;

class MCRecoDriftChamberDigiV2AssociationObj {
public:
  /// constructor
  MCRecoDriftChamberDigiV2AssociationObj();
  /// copy constructor (does a deep-copy of relation containers)
  MCRecoDriftChamberDigiV2AssociationObj(const MCRecoDriftChamberDigiV2AssociationObj&);
  /// constructor from ObjectID and MCRecoDriftChamberDigiV2AssociationData
  /// does not initialize the internal relation containers
  MCRecoDriftChamberDigiV2AssociationObj(const podio::ObjectID& id, const MCRecoDriftChamberDigiV2AssociationData& data);
  /// No assignment operator
  MCRecoDriftChamberDigiV2AssociationObj& operator=(const MCRecoDriftChamberDigiV2AssociationObj&) = delete;
  virtual ~MCRecoDriftChamberDigiV2AssociationObj();

public:
  podio::ObjectID id{};
  MCRecoDriftChamberDigiV2AssociationData data;
  std::unique_ptr<extension::DriftChamberDigiV2> m_digi{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_sim{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
