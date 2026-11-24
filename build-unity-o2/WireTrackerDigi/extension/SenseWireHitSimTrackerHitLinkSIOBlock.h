// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkSIOBlock_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkSIOBlock_H

#include "extension/DatamodelDefinition.h"

#include "podio/SIOBlock.h"

#include <sio/api.h>
#include <sio/io_device.h>
#include <sio/version.h>

#include <typeindex>
#include <string>

namespace podio {
  struct CollectionReadBuffers;
}

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


class SenseWireHitSimTrackerHitLinkSIOBlock: public podio::SIOBlock {
public:
  SenseWireHitSimTrackerHitLinkSIOBlock() :
  SIOBlock("SenseWireHitSimTrackerHitLink", sio::version::encode_version(extension::meta::schemaVersion, 0)) {
    podio::SIOBlockFactory::instance().registerBlockForCollection("extension::SenseWireHitSimTrackerHitLink", this);
  }

  SenseWireHitSimTrackerHitLinkSIOBlock(const std::string& name) :
  SIOBlock(name, sio::version::encode_version(extension::meta::schemaVersion, 0)) {}

  // Read the collection data from the device
  void read(sio::read_device& device, sio::version_type version) override;

  // Write the collection data to the device
  void write(sio::write_device& device) override;

  SIOBlock* create(const std::string& name) const override { return new SenseWireHitSimTrackerHitLinkSIOBlock(name); }
};

static SenseWireHitSimTrackerHitLinkSIOBlock _dummySenseWireHitSimTrackerHitLinkSIOBlock;

/// @}
/// \endcond
} // namespace extension


#endif
