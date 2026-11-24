// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLink_CollectionData_H
#define EXTENSION_SenseWireHitSimTrackerHitLink_CollectionData_H

// datamodel specific includes
#include "extension/SenseWireHitSimTrackerHitLinkData.h"
#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHitSimTrackerHitLink.h"
#include "extension/SenseWireHit.h"
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


using SenseWireHitSimTrackerHitLinkObjPointerContainer = std::deque<SenseWireHitSimTrackerHitLinkObj*>;
using SenseWireHitSimTrackerHitLinkDataContainer = std::vector<SenseWireHitSimTrackerHitLinkData>;


/**
 * Class encapsulating everything related to storage of data that is needed by a
 * collection.
 */
class SenseWireHitSimTrackerHitLinkCollectionData {
public:
  /**
   * The Objs of this collection
   */
  SenseWireHitSimTrackerHitLinkObjPointerContainer entries{};

  /**
   * Default constructor setting up the necessary buffers
   */
  SenseWireHitSimTrackerHitLinkCollectionData();

  /**
   * Constructor from existing I/O buffers
   */
  SenseWireHitSimTrackerHitLinkCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl);

  /**
   * Non copy-able, move-only class
   */
  SenseWireHitSimTrackerHitLinkCollectionData(const SenseWireHitSimTrackerHitLinkCollectionData&) = delete;
  SenseWireHitSimTrackerHitLinkCollectionData& operator=(const SenseWireHitSimTrackerHitLinkCollectionData&) = delete;
  SenseWireHitSimTrackerHitLinkCollectionData(SenseWireHitSimTrackerHitLinkCollectionData&& other) = default;
  SenseWireHitSimTrackerHitLinkCollectionData& operator=(SenseWireHitSimTrackerHitLinkCollectionData&& other) = default;

  /**
   * Destructor
   */
  ~SenseWireHitSimTrackerHitLinkCollectionData() = default;

  void clear(bool isSubsetColl);

  podio::CollectionWriteBuffers getCollectionBuffers(bool isSubsetColl);

  void prepareForWrite(bool isSubsetColl);

  void prepareAfterRead(uint32_t collectionID);

  void makeSubsetCollection();


  bool setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl);

private:
  // members to handle 1-to-N-relations
  podio::UVecPtr<extension::SenseWireHit> m_rel_from{std::make_unique<std::vector<extension::SenseWireHit>>()}; ///< Relation buffer for read / write
  podio::UVecPtr<edm4hep::SimTrackerHit> m_rel_to{std::make_unique<std::vector<edm4hep::SimTrackerHit>>()}; ///< Relation buffer for read / write

  // members to handle vector members

  // I/O related buffers
  podio::CollRefCollection m_refCollections{};
  podio::VectorMembersInfo m_vecmem_info{};
  std::unique_ptr<SenseWireHitSimTrackerHitLinkDataContainer> m_data{nullptr};
};


/// @}
/// \endcond
} // namespace extension


#endif
