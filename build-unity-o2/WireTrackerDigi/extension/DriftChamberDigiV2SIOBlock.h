// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2SIOBlock_H
#define EXTENSION_DriftChamberDigiV2SIOBlock_H

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


class DriftChamberDigiV2SIOBlock: public podio::SIOBlock {
public:
  DriftChamberDigiV2SIOBlock() :
  SIOBlock("DriftChamberDigiV2", sio::version::encode_version(extension::meta::schemaVersion, 0)) {
    podio::SIOBlockFactory::instance().registerBlockForCollection("extension::DriftChamberDigiV2", this);
  }

  DriftChamberDigiV2SIOBlock(const std::string& name) :
  SIOBlock(name, sio::version::encode_version(extension::meta::schemaVersion, 0)) {}

  // Read the collection data from the device
  void read(sio::read_device& device, sio::version_type version) override;

  // Write the collection data to the device
  void write(sio::write_device& device) override;

  SIOBlock* create(const std::string& name) const override { return new DriftChamberDigiV2SIOBlock(name); }
};

static DriftChamberDigiV2SIOBlock _dummyDriftChamberDigiV2SIOBlock;

/// @}
/// \endcond
} // namespace extension


#endif
