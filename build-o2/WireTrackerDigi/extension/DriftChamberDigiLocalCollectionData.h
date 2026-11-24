// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocal_CollectionData_H
#define EXTENSION_DriftChamberDigiLocal_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiLocalData.h"
#include "extension/DriftChamberDigiLocalObj.h"
#include "extension/DriftChamberDigiLocal.h"

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


using DriftChamberDigiLocalObjPointerContainer = std::deque<DriftChamberDigiLocalObj*>;
using DriftChamberDigiLocalDataContainer = std::vector<DriftChamberDigiLocalData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiLocalCollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiLocalObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiLocalCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiLocalCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiLocalCollectionData(const DriftChamberDigiLocalCollectionData&) = delete;
  DriftChamberDigiLocalCollectionData& operator=(const DriftChamberDigiLocalCollectionData&) = delete;
  DriftChamberDigiLocalCollectionData(DriftChamberDigiLocalCollectionData&& other) = default;
  DriftChamberDigiLocalCollectionData& operator=(DriftChamberDigiLocalCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiLocalCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<DriftChamberDigiLocalDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
