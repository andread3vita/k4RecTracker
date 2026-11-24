// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackSIOBlock_H
#define EXTENSION_TrackSIOBlock_H

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


class TrackSIOBlock: public podio::SIOBlock {
public:
  TrackSIOBlock() :
  SIOBlock("Track", sio::version::encode_version(extension::meta::schemaVersion, 0)) {
    podio::SIOBlockFactory::instance().registerBlockForCollection("extension::Track", this);
  }

  TrackSIOBlock(const std::string& name) :
  SIOBlock(name, sio::version::encode_version(extension::meta::schemaVersion, 0)) {}

  // Read the collection data from the device
  void read(sio::read_device& device, sio::version_type version) override;

  // Write the collection data to the device
  void write(sio::write_device& device) override;

  SIOBlock* create(const std::string& name) const override { return new TrackSIOBlock(name); }
};

static TrackSIOBlock _dummyTrackSIOBlock;

/// @}
/// \endcond
} // namespace extension


#endif
