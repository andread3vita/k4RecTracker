// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitOBJ_H
#define EXTENSION_SenseWireHitOBJ_H

// data model specific includes
#include "extension/SenseWireHitData.h"
#include <vector>

#include "podio/ObjectID.h"
#include <vector>



namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{

class SenseWireHit;

class SenseWireHitObj {
public:
  /// constructor
  SenseWireHitObj();
  /// copy constructor (does a deep-copy of relation containers)
  SenseWireHitObj(const SenseWireHitObj&);
  /// constructor from ObjectID and SenseWireHitData
  /// does not initialize the internal relation containers
  SenseWireHitObj(const podio::ObjectID& id, const SenseWireHitData& data);
  /// No assignment operator
  SenseWireHitObj& operator=(const SenseWireHitObj&) = delete;
  virtual ~SenseWireHitObj();

public:
  podio::ObjectID id{};
  SenseWireHitData data;
  std::vector<std::uint16_t>* m_nElectrons{nullptr};
};

/// @}
/// \endcond
} // namespace extension


#endif
