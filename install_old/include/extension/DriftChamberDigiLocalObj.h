// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocalOBJ_H
#define EXTENSION_DriftChamberDigiLocalOBJ_H

// data model specific includes
#include "extension/DriftChamberDigiLocalData.h"

#include "podio/ObjectID.h"



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigiLocal;

class DriftChamberDigiLocalObj {
public:
  /// constructor
  DriftChamberDigiLocalObj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiLocalObj(const DriftChamberDigiLocalObj&);
  /// constructor from ObjectID and DriftChamberDigiLocalData
  /// does not initialize the internal relation containers
  DriftChamberDigiLocalObj(const podio::ObjectID& id, const DriftChamberDigiLocalData& data);
  /// No assignment operator
  DriftChamberDigiLocalObj& operator=(const DriftChamberDigiLocalObj&) = delete;
  virtual ~DriftChamberDigiLocalObj() = default;

public:
  podio::ObjectID id{};
  DriftChamberDigiLocalData data;
};

/// @}
/// \endcond
} // namespace extension


#endif
