// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiOBJ_H
#define EXTENSION_DriftChamberDigiOBJ_H

// data model specific includes
#include "extension/DriftChamberDigiData.h"

#include "podio/ObjectID.h"



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigi;

class DriftChamberDigiObj {
public:
  /// constructor
  DriftChamberDigiObj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiObj(const DriftChamberDigiObj&);
  /// constructor from ObjectID and DriftChamberDigiData
  /// does not initialize the internal relation containers
  DriftChamberDigiObj(const podio::ObjectID& id, const DriftChamberDigiData& data);
  /// No assignment operator
  DriftChamberDigiObj& operator=(const DriftChamberDigiObj&) = delete;
  virtual ~DriftChamberDigiObj() = default;

public:
  podio::ObjectID id{};
  DriftChamberDigiData data;
};

/// @}
/// \endcond
} // namespace extension


#endif
