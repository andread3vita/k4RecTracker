// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2Association_CollectionData_H
#define EXTENSION_MCRecoDriftChamberDigiV2Association_CollectionData_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiV2AssociationData.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/MCRecoDriftChamberDigiV2Association.h"
#include "extension/DriftChamberDigiV2.h"
#include "edm4hep/SimTrackerHit.h"

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


using MCRecoDriftChamberDigiV2AssociationObjPointerContainer = std::deque<MCRecoDriftChamberDigiV2AssociationObj*>;
using MCRecoDriftChamberDigiV2AssociationDataContainer = std::vector<MCRecoDriftChamberDigiV2AssociationData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class MCRecoDriftChamberDigiV2AssociationCollectionData {
public:
  /**
   * The Objs of this collection
   */
  MCRecoDriftChamberDigiV2AssociationObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  MCRecoDriftChamberDigiV2AssociationCollectionData(const MCRecoDriftChamberDigiV2AssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiV2AssociationCollectionData& operator=(const MCRecoDriftChamberDigiV2AssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiV2AssociationCollectionData(MCRecoDriftChamberDigiV2AssociationCollectionData&& other) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionData& operator=(MCRecoDriftChamberDigiV2AssociationCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~MCRecoDriftChamberDigiV2AssociationCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::DriftChamberDigiV2> m_rel_digi{std::make_unique<std::vector<extension::DriftChamberDigiV2>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_sim{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<MCRecoDriftChamberDigiV2AssociationDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
