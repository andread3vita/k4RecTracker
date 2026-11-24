// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/DriftChamberDigiV2CollectionData.h"
#include "extension/DriftChamberDigiV2Collection.h"

#include <numeric>

#include <podio/detail/RelationIOHelpers.h>

namespace extension {

/// \cond PODIO_INTERNAL


DriftChamberDigiV2CollectionData::DriftChamberDigiV2CollectionData() :  m_vec_nElectrons(new std::vector<std::uint16_t>()),
  m_data(new DriftChamberDigiV2DataContainer()) {
  m_vecmem_info.emplace_back("std::uint16_t", &m_vec_nElectrons);
}

DriftChamberDigiV2CollectionData::DriftChamberDigiV2CollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl) :
  m_refCollections(std::move(*buffers.references)),
  m_vecmem_info(std::move(*buffers.vectorMembers)) {
  // For subset collections we are done, for proper collections we still have to
  // populate the data and vector members
  if (!isSubsetColl) {
    m_data.reset(buffers.dataAsVector<extension::DriftChamberDigiV2Data>());
    buffers.data = nullptr;

  // The following is ugly code for the case when reading garbage data from ROOT
  // after calling dataAsVector, which can trigger infinite loops that consume
  // all the available memory and this works at least for GCC 15 and Clang 20 in
  // {Debug,RelWithDebInfo,Release} modes.
  // https://github.com/AIDASoft/podio/pull/817#issuecomment-3266748609 and
  // https://github.com/AIDASoft/podio/pull/842
  volatile std::uint64_t s = m_data->size();
  if (s > 1e15) {
    throw std::runtime_error("Bad data after reading: a collection is too big (extension::DriftChamberDigiV2)");
  } else {
    if (s == 0) {
      for ([[maybe_unused]] const auto& _ : *m_data.get()) {
        throw std::runtime_error("Bad data after reading: zero-sized collection with data (extension::DriftChamberDigiV2)");
      }
    }
  }
  // end of ugly


  m_vec_nElectrons.reset(podio::CollectionReadBuffers::asVector<std::uint16_t>(m_vecmem_info[0].second));
  }

  // Cleanup these to avoid leaking them
  delete buffers.references;
  buffers.references = nullptr;
  delete buffers.vectorMembers;
  buffers.vectorMembers = nullptr;
}

void DriftChamberDigiV2CollectionData::clear(bool isSubsetColl) {
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
  if (m_vec_nElectrons) {
    m_vec_nElectrons->clear();
  }
  m_vecs_nElectrons.clear();

  for (auto& obj : entries) { delete obj; }
  entries.clear();
}

podio::CollectionWriteBuffers DriftChamberDigiV2CollectionData::getCollectionBuffers(bool isSubsetColl) {
  if (!isSubsetColl) {
    // Make sure these point to the right place, even if a collection has been
    // moved since it has been created
    m_vecmem_info[0].second = &m_vec_nElectrons;
  }
return {
    isSubsetColl ? nullptr : static_cast<void*>(m_data.get()),
    &m_refCollections, // only need to store the ObjectIDs of the referenced objects
    &m_vecmem_info
  };
}

void DriftChamberDigiV2CollectionData::prepareForWrite(bool isSubsetColl) {
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

const auto nElectrons_size = std::accumulate(entries.begin(), entries.end(), 0,
  [](size_t sum, const DriftChamberDigiV2Obj* obj) { return sum + obj->m_nElectrons->size(); });
  m_vec_nElectrons->reserve(nElectrons_size);
  int nElectrons_index = 0;
  for (size_t i = 0, size = m_data->size(); i != size; ++i) {
    (*m_data)[i].nElectrons_begin = nElectrons_index;
    (*m_data)[i].nElectrons_end += nElectrons_index;
    nElectrons_index = (*m_data)[i].nElectrons_end;
    for (const auto& it : (*m_vecs_nElectrons[i])) { m_vec_nElectrons->push_back(it); }

  }
}

void DriftChamberDigiV2CollectionData::prepareAfterRead(uint32_t collectionID) {
  int index = 0;
  for (const auto& data : *m_data) {
    auto obj = new DriftChamberDigiV2Obj({index, collectionID}, data);

    obj->m_nElectrons = m_vec_nElectrons.get();
    entries.emplace_back(obj);
    ++index;
  }

  // at this point we could clear the I/O data buffer, but we keep them intact
  // because then we can save a call to prepareForWrite
}


void DriftChamberDigiV2CollectionData::createRelations(DriftChamberDigiV2Obj* obj) {
  // We take ownership of these here
  m_vecs_nElectrons.emplace_back(obj->m_nElectrons);
}

bool DriftChamberDigiV2CollectionData::setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl) {
  if (isSubsetColl) {
    for (const auto& id : *m_refCollections[0]) {
      podio::CollectionBase* coll = nullptr;
      extension::DriftChamberDigiV2Obj* obj = nullptr;
      if (collectionProvider->get(id.collectionID, coll)) {
        auto* tmp_coll = static_cast<extension::DriftChamberDigiV2Collection*>(coll);
        obj = tmp_coll->m_storage.entries[id.index];
      }
      entries.push_back(obj);
    }
    return true; // TODO: check success, how?
  }

  // Normal collections have to resolve all relations

  return true; // TODO: check success, how?
}

void DriftChamberDigiV2CollectionData::makeSubsetCollection() {
  // Subset collections do not need all the data buffers that normal
  // collections need, so we can free them here
  m_vecmem_info.clear();

  m_data.reset(nullptr);

  m_vec_nElectrons.reset(nullptr);

  // Subset collections need one vector of ObjectIDs for I/O purposes.
  m_refCollections.resize(1);
  m_refCollections[0] = std::make_unique<std::vector<podio::ObjectID>>();
}


/// \endcond
} // namespace extension

