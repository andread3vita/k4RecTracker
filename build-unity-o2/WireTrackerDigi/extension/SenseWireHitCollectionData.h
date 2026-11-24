// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHit_CollectionData_H
#define EXTENSION_SenseWireHit_CollectionData_H

// datamodel specific includes
#include "extension/SenseWireHitData.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHit.h"

// schema evolution specific includes

// podio specific includes
#include "podio/CollectionBuffers.h"
#include "podio/ICollectionProvider.h"

#include <deque>
#include <memory>

namespace extension {

/// \cond PODIO_INTERNAL
/// \addtogroup extension_podio_internal
/// @{


using SenseWireHitObjPointerContainer = std::deque<SenseWireHitObj*>;
using SenseWireHitDataContainer = std::vector<SenseWireHitData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class SenseWireHitCollectionData {
public:
  /**
   * The Objs of this collection
   */
  SenseWireHitObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  SenseWireHitCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  SenseWireHitCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  SenseWireHitCollectionData(const SenseWireHitCollectionData&) = delete;
  SenseWireHitCollectionData& operator=(const SenseWireHitCollectionData&) = delete;
  SenseWireHitCollectionData(SenseWireHitCollectionData&& other) = default;
  SenseWireHitCollectionData& operator=(SenseWireHitCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~SenseWireHitCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();

  void createRelations(SenseWireHitObj* obj);

  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members
  podio::UVecPtr<std::uint16_t> m_vec_nElectrons{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<std::uint16_t>> m_vecs_nElectrons{}; /// pointers to individual member vectors

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<SenseWireHitDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
