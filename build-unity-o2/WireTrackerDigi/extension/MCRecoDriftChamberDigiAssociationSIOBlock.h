// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationSIOBlock_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationSIOBlock_H

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


class MCRecoDriftChamberDigiAssociationSIOBlock: public podio::SIOBlock {
public:
  MCRecoDriftChamberDigiAssociationSIOBlock() :
  SIOBlock("MCRecoDriftChamberDigiAssociation", sio::version::encode_version(extension::meta::schemaVersion, 0)) {
    podio::SIOBlockFactory::instance().registerBlockForCollection("extension::MCRecoDriftChamberDigiAssociation", this);
  }

  MCRecoDriftChamberDigiAssociationSIOBlock(const std::string& name) :
  SIOBlock(name, sio::version::encode_version(extension::meta::schemaVersion, 0)) {}

  // Read the collection data from the device
  void read(sio::read_device& device, sio::version_type version) override;

  // Write the collection data to the device
  void write(sio::write_device& device) override;

  SIOBlock* create(const std::string& name) const override { return new MCRecoDriftChamberDigiAssociationSIOBlock(name); }
};

static MCRecoDriftChamberDigiAssociationSIOBlock _dummyMCRecoDriftChamberDigiAssociationSIOBlock;

/// @}
/// \endcond
} // namespace extension


#endif
