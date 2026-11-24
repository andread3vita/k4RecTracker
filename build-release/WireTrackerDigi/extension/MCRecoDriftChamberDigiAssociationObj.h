// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationOBJ_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationOBJ_H

// data model specific includes
#include "extension/MCRecoDriftChamberDigiAssociationData.h"

#include "podio/ObjectID.h"
#include <memory>
// forward declarations
namespace extension {
class DriftChamberDigi;
}
namespace edm4hep {
class SimTrackerHit;
}


namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class MCRecoDriftChamberDigiAssociation;

class MCRecoDriftChamberDigiAssociationObj {
public:
  /// constructor
  MCRecoDriftChamberDigiAssociationObj();
  /// copy constructor (does a deep-copy of relation containers)
  MCRecoDriftChamberDigiAssociationObj(const MCRecoDriftChamberDigiAssociationObj&);
  /// constructor from ObjectID and MCRecoDriftChamberDigiAssociationData
  /// does not initialize the internal relation containers
  MCRecoDriftChamberDigiAssociationObj(const podio::ObjectID& id, const MCRecoDriftChamberDigiAssociationData& data);
  /// No assignment operator
  MCRecoDriftChamberDigiAssociationObj& operator=(const MCRecoDriftChamberDigiAssociationObj&) = delete;
  virtual ~MCRecoDriftChamberDigiAssociationObj();

public:
  podio::ObjectID id{};
  MCRecoDriftChamberDigiAssociationData data;
  std::unique_ptr<extension::DriftChamberDigi> m_digi{nullptr};
  std::unique_ptr<edm4hep::SimTrackerHit> m_sim{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
