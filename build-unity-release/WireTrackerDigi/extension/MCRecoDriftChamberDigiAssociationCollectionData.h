// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociation_CollectionData_H
#define EXTENSION_MCRecoDriftChamberDigiAssociation_CollectionData_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiAssociationData.h"
#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/MCRecoDriftChamberDigiAssociation.h"
#include "extension/DriftChamberDigi.h"
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


using MCRecoDriftChamberDigiAssociationObjPointerContainer = std::deque<MCRecoDriftChamberDigiAssociationObj*>;
using MCRecoDriftChamberDigiAssociationDataContainer = std::vector<MCRecoDriftChamberDigiAssociationData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class MCRecoDriftChamberDigiAssociationCollectionData {
public:
  /**
   * The Objs of this collection
   */
  MCRecoDriftChamberDigiAssociationObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  MCRecoDriftChamberDigiAssociationCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  MCRecoDriftChamberDigiAssociationCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  MCRecoDriftChamberDigiAssociationCollectionData(const MCRecoDriftChamberDigiAssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiAssociationCollectionData& operator=(const MCRecoDriftChamberDigiAssociationCollectionData&) = delete;
  MCRecoDriftChamberDigiAssociationCollectionData(MCRecoDriftChamberDigiAssociationCollectionData&& other) = default;
  MCRecoDriftChamberDigiAssociationCollectionData& operator=(MCRecoDriftChamberDigiAssociationCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~MCRecoDriftChamberDigiAssociationCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::DriftChamberDigi> m_rel_digi{std::make_unique<std::vector<extension::DriftChamberDigi>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_sim{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<MCRecoDriftChamberDigiAssociationDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
