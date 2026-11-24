// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/TrackSIOBlock.h"
#include "extension/TrackCollection.h"

#include "podio/CollectionBuffers.h"
#include "podio/CollectionBufferFactory.h"

#include <sio/block.h>
#include <sio/io_device.h>
#include <sio/version.h>

namespace extension {

/// \cond PODIO_INTERNAL


void TrackSIOBlock::read(sio::read_device& device, sio::version_type version) {
  const auto& bufferFactory = podio::CollectionBufferFactory::instance();
  // TODO:
  // - Error handling of empty optional
  auto maybeBuffers = bufferFactory.createBuffers("extension::TrackCollection", sio::version::major_version(version), m_subsetColl);
  m_readBuffers = std::move(maybeBuffers).value_or(podio::CollectionReadBuffers{});

  if (not m_subsetColl) {
    unsigned size(0);
    device.data( size );
    auto* dataVec = m_readBuffers.dataAsVector<extension::TrackData>();
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

  if (not m_subsetColl) {
    //---- read vector members
    const auto* vecMemInfo = m_readBuffers.vectorMembers;
    unsigned size{0};

    auto* vec0 = reinterpret_cast<std::vector<std::int32_t>*>(vecMemInfo->at(0).second);
    size = 0u;
    device.data(size);
    vec0->resize(size);
    podio::handlePODDataSIO(device, vec0->data(), size);

    auto* vec1 = reinterpret_cast<std::vector<edm4hep::TrackState>*>(vecMemInfo->at(1).second);
    size = 0u;
    device.data(size);
    vec1->resize(size);
    podio::handlePODDataSIO(device, vec1->data(), size);

    auto* vec2 = reinterpret_cast<std::vector<edm4hep::Quantity>*>(vecMemInfo->at(2).second);
    size = 0u;
    device.data(size);
    vec2->resize(size);
    podio::handlePODDataSIO(device, vec2->data(), size);

  }
}

void TrackSIOBlock::write(sio::write_device& device) {
  if (not m_subsetColl) {
    auto* dataVec = m_writeBuffers.dataAsVector<extension::TrackData>();
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

  if (not m_subsetColl) {
    //---- write vector members
    const auto* vecMemInfo = m_writeBuffers.vectorMembers;
    unsigned size{0};

    auto* vec0 = *reinterpret_cast<std::vector<std::int32_t>**>(vecMemInfo->at(0).second);
    size = vec0->size();
    device.data(size);
    podio::handlePODDataSIO(device, &(*vec0)[0], size);


    auto* vec1 = *reinterpret_cast<std::vector<edm4hep::TrackState>**>(vecMemInfo->at(1).second);
    size = vec1->size();
    device.data(size);
    podio::handlePODDataSIO(device, &(*vec1)[0], size);


    auto* vec2 = *reinterpret_cast<std::vector<edm4hep::Quantity>**>(vecMemInfo->at(2).second);
    size = vec2->size();
    device.data(size);
    podio::handlePODDataSIO(device, &(*vec2)[0], size);


  }
}

namespace {
  // Create one instance of the type in order to ensure that the SioBlock
  // library actually needs linking to the core library. Otherwise it is
  // possible that the registry is not populated when the SioBlock library is
  // loaded, e.g. when using the python bindings.
  const auto elemTrack = extension::Track{};
}

/// \endcond
} // namespace extension

