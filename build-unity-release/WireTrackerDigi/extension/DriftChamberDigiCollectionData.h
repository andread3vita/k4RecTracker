// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigi_CollectionData_H
#define EXTENSION_DriftChamberDigi_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiData.h"
#include "extension/DriftChamberDigiObj.h"
#include "extension/DriftChamberDigi.h"

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


using DriftChamberDigiObjPointerContainer = std::deque<DriftChamberDigiObj*>;
using DriftChamberDigiDataContainer = std::vector<DriftChamberDigiData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiCollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiCollectionData(const DriftChamberDigiCollectionData&) = delete;
  DriftChamberDigiCollectionData& operator=(const DriftChamberDigiCollectionData&) = delete;
  DriftChamberDigiCollectionData(DriftChamberDigiCollectionData&& other) = default;
  DriftChamberDigiCollectionData& operator=(DriftChamberDigiCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiCollectionData() = default;

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
  std::unique_ptr<DriftChamberDigiDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
