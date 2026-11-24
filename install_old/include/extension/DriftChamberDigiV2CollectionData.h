// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2_CollectionData_H
#define EXTENSION_DriftChamberDigiV2_CollectionData_H

// datamodel specific includes
#include "extension/DriftChamberDigiV2Data.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2.h"

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


using DriftChamberDigiV2ObjPointerContainer = std::deque<DriftChamberDigiV2Obj*>;
using DriftChamberDigiV2DataContainer = std::vector<DriftChamberDigiV2Data>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class DriftChamberDigiV2CollectionData {
public:
  /**
   * The Objs of this collection
   */
  DriftChamberDigiV2ObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  DriftChamberDigiV2CollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  DriftChamberDigiV2CollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  DriftChamberDigiV2CollectionData(const DriftChamberDigiV2CollectionData&) = delete;
  DriftChamberDigiV2CollectionData& operator=(const DriftChamberDigiV2CollectionData&) = delete;
  DriftChamberDigiV2CollectionData(DriftChamberDigiV2CollectionData&& other) = default;
  DriftChamberDigiV2CollectionData& operator=(DriftChamberDigiV2CollectionData&& other) = default;

  /**
   * Destructor
   */
  ~DriftChamberDigiV2CollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();

  void createRelations(DriftChamberDigiV2Obj* obj);

  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations

  // members to handle vector members
  podio::UVecPtr<std::uint16_t> m_vec_nElectrons{nullptr}; /// combined vector of all objects in collection
  std::vector<podio::UVecPtr<std::uint16_t>> m_vecs_nElectrons{}; /// pointers to individual member vectors

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<DriftChamberDigiV2DataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
