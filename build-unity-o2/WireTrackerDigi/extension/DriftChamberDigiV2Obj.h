// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2OBJ_H
#define EXTENSION_DriftChamberDigiV2OBJ_H

// data model specific includes
#include "extension/DriftChamberDigiV2Data.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class DriftChamberDigiV2;

class DriftChamberDigiV2Obj {
public:
  /// constructor
  DriftChamberDigiV2Obj();
  /// copy constructor (does a deep-copy of relation containers)
  DriftChamberDigiV2Obj(const DriftChamberDigiV2Obj&);
  /// constructor from ObjectID and DriftChamberDigiV2Data
  /// does not initialize the internal relation containers
  DriftChamberDigiV2Obj(const podio::ObjectID& id, const DriftChamberDigiV2Data& data);
  /// No assignment operator
  DriftChamberDigiV2Obj& operator=(const DriftChamberDigiV2Obj&) = delete;
  virtual ~DriftChamberDigiV2Obj();

public:
  podio::ObjectID id{};
  DriftChamberDigiV2Data data;
  std::vector<std::uint16_t>* m_nElectrons{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
