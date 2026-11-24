// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/MCRecoDriftChamberDigiV2AssociationSIOBlock.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationCollection.h"

#include "podio/CollectionBuffers.h"
#include "podio/CollectionBufferFactory.h"

#include <sio/block.h>
#include <sio/io_device.h>
#include <sio/version.h>

namespace extension {

/// \cond PODIO_INTERNAL


void MCRecoDriftChamberDigiV2AssociationSIOBlock::read(sio::read_device& device, sio::version_type version) {
  const auto& bufferFactory = podio::CollectionBufferFactory::instance();
  // TODO:
  // - Error handling of empty optional
  auto maybeBuffers = bufferFactory.createBuffers("extension::MCRecoDriftChamberDigiV2AssociationCollection", sio::version::major_version(version), m_subsetColl);
  m_readBuffers = std::move(maybeBuffers).value_or(podio::CollectionReadBuffers{});

  if (not m_subsetColl) {
    unsigned size(0);
    device.data( size );
    auto* dataVec = m_readBuffers.dataAsVector<extension::MCRecoDriftChamberDigiV2AssociationData>();
    dataVec->resize(size);
    podio::handlePODDataSIO(device, dataVec->data(), size);
  }

  //---- read ref collections -----
  const auto* refCols = m_readBuffers.references;
  for( auto& refC : *refCols ){
    unsigned size{0};
    device.data( size ) ;
    refC->resize(size) ;
    podio::handlePODDataSIO( device ,  refC->data(), size ) ;
  }

}

void MCRecoDriftChamberDigiV2AssociationSIOBlock::write(sio::write_device& device) {
  if (not m_subsetColl) {
    auto* dataVec = m_writeBuffers.dataAsVector<extension::MCRecoDriftChamberDigiV2AssociationData>();
    unsigned size = dataVec->size() ;
    device.data( size ) ;
    podio::handlePODDataSIO( device ,  dataVec->data(), size ) ;
  }

  //---- write ref collections -----
  const auto* refCols = m_writeBuffers.references;
  for( auto& refC : *refCols ){
    unsigned size = refC->size() ;
    device.data( size ) ;
    podio::handlePODDataSIO( device ,  refC->data(), size ) ;
  }

}

namespace {
  // Create one instance of the type in order to ensure that the SioBlock
  // library actually needs linking to the core library. Otherwise it is
  // possible that the registry is not populated when the SioBlock library is
  // loaded, e.g. when using the python bindings.
  const auto elemMCRecoDriftChamberDigiV2Association = extension::MCRecoDriftChamberDigiV2Association{};
}

/// \endcond
} // namespace extension

