// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/MCRecoDriftChamberDigiAssociationCollectionData.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"

#include "extension/DriftChamberDigiCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include <podio/detail/RelationIOHelpers.h>

namespace extension {

/// \cond PODIO_INTERNAL


MCRecoDriftChamberDigiAssociationCollectionData::MCRecoDriftChamberDigiAssociationCollectionData() :  m_data(new MCRecoDriftChamberDigiAssociationDataContainer()) {
  m_refCollections.emplace_back(std::make_unique<std::vector<podio::ObjectID>>());
  m_refCollections.emplace_back(std::make_unique<std::vector<podio::ObjectID>>());
}

MCRecoDriftChamberDigiAssociationCollectionData::MCRecoDriftChamberDigiAssociationCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl) :
  m_rel_digi(new std::vector<extension::DriftChamberDigi>()),
  m_rel_sim(new std::vector<edm4hep::SimTrackerHit>()),
  m_refCollections(std::move(*buffers.references)),
  m_vecmem_info(std::move(*buffers.vectorMembers)) {
  // For subset collections we are done, for proper collections we still have to
  // populate the data and vector members
  if (!isSubsetColl) {
    m_data.reset(buffers.dataAsVector<extension::MCRecoDriftChamberDigiAssociationData>());
    buffers.data = nullptr;

  // The following is ugly code for the case when reading garbage data from ROOT
  // after calling dataAsVector, which can trigger infinite loops that consume
  // all the available memory and this works at least for GCC 15 and Clang 20 in
  // {Debug,RelWithDebInfo,Release} modes.
  // https://github.com/AIDASoft/podio/pull/817#issuecomment-3266748609 and
  // https://github.com/AIDASoft/podio/pull/842
  volatile std::uint64_t s = m_data->size();
  if (s > 1e15) {
    throw std::runtime_error("Bad data after reading: a collection is too big (extension::MCRecoDriftChamberDigiAssociation)");
  } else {
    if (s == 0) {
      for ([[maybe_unused]] const auto& _ : *m_data.get()) {
        throw std::runtime_error("Bad data after reading: zero-sized collection with data (extension::MCRecoDriftChamberDigiAssociation)");
      }
    }
  }
  // end of ugly


  }

  // Cleanup these to avoid leaking them
  delete buffers.references;
  buffers.references = nullptr;
  delete buffers.vectorMembers;
  buffers.vectorMembers = nullptr;
}

void MCRecoDriftChamberDigiAssociationCollectionData::clear(bool isSubsetColl) {
  if (isSubsetColl) {
    // We don't own the objects so no cleanup to do here
    entries.clear();
    // Clear the ObjectID I/O buffer
    for (const auto& pointer : m_refCollections) { pointer->clear(); }
    return;
  }

  // Normal collections manage a bit more and have to clean up a bit more
  if (m_data) {
    m_data->clear();
  }
  for (const auto& pointer : m_refCollections) { pointer->clear(); }
  if (m_rel_digi) {
    for (auto& item : (*m_rel_digi)) { item.unlink(); }
    m_rel_digi->clear();
  }

  if (m_rel_sim) {
    for (auto& item : (*m_rel_sim)) { item.unlink(); }
    m_rel_sim->clear();
  }

  for (auto& obj : entries) { delete obj; }
  entries.clear();
}

podio::CollectionWriteBuffers MCRecoDriftChamberDigiAssociationCollectionData::getCollectionBuffers(bool isSubsetColl) {
return {
    isSubsetColl ? nullptr : static_cast<void*>(m_data.get()),
    &m_refCollections, // only need to store the ObjectIDs of the referenced objects
    &m_vecmem_info
  };
}

void MCRecoDriftChamberDigiAssociationCollectionData::prepareForWrite(bool isSubsetColl) {
  for (const auto& pointer : m_refCollections) { pointer->clear(); }

  // If this is a subset collection use the relation storing mechanism to
  // store the ObjectIDs of all referenced objects and nothing else
  if (isSubsetColl) {
    for (const auto* obj : entries) {
      m_refCollections[0]->emplace_back(obj->id);
    }
    return;
  }

  // Normal collections have to store the data and all the relations
  m_data->reserve(entries.size());
  for (const auto& obj : entries) { m_data->push_back(obj->data); }


  for (auto& obj : entries) {
    if (obj->m_digi) {
      m_refCollections[0]->emplace_back(obj->m_digi->getObjectID());
    } else {
      m_refCollections[0]->push_back({podio::ObjectID::invalid, 0});
    }
  }

  for (auto& obj : entries) {
    if (obj->m_sim) {
      m_refCollections[1]->emplace_back(obj->m_sim->getObjectID());
    } else {
      m_refCollections[1]->push_back({podio::ObjectID::invalid, 0});
    }
  }

}

void MCRecoDriftChamberDigiAssociationCollectionData::prepareAfterRead(uint32_t collectionID) {
  int index = 0;
  for (const auto& data : *m_data) {
    auto obj = new MCRecoDriftChamberDigiAssociationObj({index, collectionID}, data);

    entries.emplace_back(obj);
    ++index;
  }

  // at this point we could clear the I/O data buffer, but we keep them intact
  // because then we can save a call to prepareForWrite
}



bool MCRecoDriftChamberDigiAssociationCollectionData::setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl) {
  if (isSubsetColl) {
    for (const auto& id : *m_refCollections[0]) {
      podio::CollectionBase* coll = nullptr;
      extension::MCRecoDriftChamberDigiAssociationObj* obj = nullptr;
      if (collectionProvider->get(id.collectionID, coll)) {
        auto* tmp_coll = static_cast<extension::MCRecoDriftChamberDigiAssociationCollection*>(coll);
        obj = tmp_coll->m_storage.entries[id.index];
      }
      entries.push_back(obj);
    }
    return true; // TODO: check success, how?
  }

  // Normal collections have to resolve all relations
  for (unsigned int i = 0, size = entries.size(); i != size; ++i) {
    const auto id = (*m_refCollections[0])[i];
    if (id.index != podio::ObjectID::invalid) {
      podio::CollectionBase* coll = nullptr;
      if (!collectionProvider->get(id.collectionID, coll)) {
        entries[i]->m_digi = nullptr;
        continue;
      }
      podio::detail::addSingleRelation(entries[i]->m_digi, coll, id);
    } else {
      entries[i]->m_digi = nullptr;
    }
  }

  for (unsigned int i = 0, size = entries.size(); i != size; ++i) {
    const auto id = (*m_refCollections[1])[i];
    if (id.index != podio::ObjectID::invalid) {
      podio::CollectionBase* coll = nullptr;
      if (!collectionProvider->get(id.collectionID, coll)) {
        entries[i]->m_sim = nullptr;
        continue;
      }
      podio::detail::addSingleRelation(entries[i]->m_sim, coll, id);
    } else {
      entries[i]->m_sim = nullptr;
    }
  }


  return true; // TODO: check success, how?
}

void MCRecoDriftChamberDigiAssociationCollectionData::makeSubsetCollection() {
  // Subset collections do not need all the data buffers that normal
  // collections need, so we can free them here
  m_vecmem_info.clear();

  m_data.reset(nullptr);

  m_rel_digi.reset(nullptr);
  m_rel_sim.reset(nullptr);

  // Subset collections need one vector of ObjectIDs for I/O purposes.
  m_refCollections.resize(1);
  m_refCollections[0] = std::make_unique<std::vector<podio::ObjectID>>();
}


/// \endcond
} // namespace extension

