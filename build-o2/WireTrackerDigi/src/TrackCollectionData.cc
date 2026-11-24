// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/TrackCollectionData.h"
#include "extension/TrackCollection.h"

#include "extension/DriftChamberDigiCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackerHit.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include <numeric>

#include <podio/detail/RelationIOHelpers.h>

namespace extension {

/// \cond PODIO_INTERNAL


TrackCollectionData::TrackCollectionData() :  m_vec_subdetectorHitNumbers(new std::vector<std::int32_t>()),
  m_vec_trackStates(new std::vector<edm4hep::TrackState>()),
  m_vec_dxQuantities(new std::vector<edm4hep::Quantity>()),
  m_data(new TrackDataContainer()) {
  m_refCollections.emplace_back(std::make_unique<std::vector<podio::ObjectID>>());
  m_refCollections.emplace_back(std::make_unique<std::vector<podio::ObjectID>>());
  m_vecmem_info.emplace_back("std::int32_t", &m_vec_subdetectorHitNumbers);
  m_vecmem_info.emplace_back("edm4hep::TrackState", &m_vec_trackStates);
  m_vecmem_info.emplace_back("edm4hep::Quantity", &m_vec_dxQuantities);
}

TrackCollectionData::TrackCollectionData(podio::CollectionReadBuffers&& buffers, bool isSubsetColl) :
  m_rel_trackerHits(new std::vector<extension::TrackerHit>()),
  m_rel_tracks(new std::vector<extension::Track>()),
  m_refCollections(std::move(*buffers.references)),
  m_vecmem_info(std::move(*buffers.vectorMembers)) {
  // For subset collections we are done, for proper collections we still have to
  // populate the data and vector members
  if (!isSubsetColl) {
    m_data.reset(buffers.dataAsVector<extension::TrackData>());
    buffers.data = nullptr;

  // The following is ugly code for the case when reading garbage data from ROOT
  // after calling dataAsVector, which can trigger infinite loops that consume
  // all the available memory and this works at least for GCC 15 and Clang 20 in
  // {Debug,RelWithDebInfo,Release} modes.
  // https://github.com/AIDASoft/podio/pull/817#issuecomment-3266748609 and
  // https://github.com/AIDASoft/podio/pull/842
  volatile std::uint64_t s = m_data->size();
  if (s > 1e15) {
    throw std::runtime_error("Bad data after reading: a collection is too big (extension::Track)");
  } else {
    if (s == 0) {
      for ([[maybe_unused]] const auto& _ : *m_data.get()) {
        throw std::runtime_error("Bad data after reading: zero-sized collection with data (extension::Track)");
      }
    }
  }
  // end of ugly


  m_vec_subdetectorHitNumbers.reset(podio::CollectionReadBuffers::asVector<std::int32_t>(m_vecmem_info[0].second));
  m_vec_trackStates.reset(podio::CollectionReadBuffers::asVector<edm4hep::TrackState>(m_vecmem_info[1].second));
  m_vec_dxQuantities.reset(podio::CollectionReadBuffers::asVector<edm4hep::Quantity>(m_vecmem_info[2].second));
  }

  // Cleanup these to avoid leaking them
  delete buffers.references;
  buffers.references = nullptr;
  delete buffers.vectorMembers;
  buffers.vectorMembers = nullptr;
}

void TrackCollectionData::clear(bool isSubsetColl) {
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
  // clear relations to trackerHits. Make sure to unlink() the reference data as they may be gone already.
  for (const auto& pointer : m_rel_trackerHits_tmp) {
    for (auto& item : *pointer) {
      item.unlink();
    }
  }
  m_rel_trackerHits_tmp.clear();
  if (m_rel_trackerHits) {
    for (auto& item : (*m_rel_trackerHits)) { item.unlink(); }
    m_rel_trackerHits->clear();
  }

  // clear relations to tracks. Make sure to unlink() the reference data as they may be gone already.
  for (const auto& pointer : m_rel_tracks_tmp) {
    for (auto& item : *pointer) {
      item.unlink();
    }
  }
  m_rel_tracks_tmp.clear();
  if (m_rel_tracks) {
    for (auto& item : (*m_rel_tracks)) { item.unlink(); }
    m_rel_tracks->clear();
  }

  if (m_vec_subdetectorHitNumbers) {
    m_vec_subdetectorHitNumbers->clear();
  }
  m_vecs_subdetectorHitNumbers.clear();

  if (m_vec_trackStates) {
    m_vec_trackStates->clear();
  }
  m_vecs_trackStates.clear();

  if (m_vec_dxQuantities) {
    m_vec_dxQuantities->clear();
  }
  m_vecs_dxQuantities.clear();

  for (auto& obj : entries) { delete obj; }
  entries.clear();
}

podio::CollectionWriteBuffers TrackCollectionData::getCollectionBuffers(bool isSubsetColl) {
  if (!isSubsetColl) {
    // Make sure these point to the right place, even if a collection has been
    // moved since it has been created
    m_vecmem_info[0].second = &m_vec_subdetectorHitNumbers;
    m_vecmem_info[1].second = &m_vec_trackStates;
    m_vecmem_info[2].second = &m_vec_dxQuantities;
  }
return {
    isSubsetColl ? nullptr : static_cast<void*>(m_data.get()),
    &m_refCollections, // only need to store the ObjectIDs of the referenced objects
    &m_vecmem_info
  };
}

void TrackCollectionData::prepareForWrite(bool isSubsetColl) {
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

  int trackerHits_index = 0;
  int tracks_index = 0;
const auto subdetectorHitNumbers_size = std::accumulate(entries.begin(), entries.end(), 0,
  [](size_t sum, const TrackObj* obj) { return sum + obj->m_subdetectorHitNumbers->size(); });
  m_vec_subdetectorHitNumbers->reserve(subdetectorHitNumbers_size);
  int subdetectorHitNumbers_index = 0;
const auto trackStates_size = std::accumulate(entries.begin(), entries.end(), 0,
  [](size_t sum, const TrackObj* obj) { return sum + obj->m_trackStates->size(); });
  m_vec_trackStates->reserve(trackStates_size);
  int trackStates_index = 0;
const auto dxQuantities_size = std::accumulate(entries.begin(), entries.end(), 0,
  [](size_t sum, const TrackObj* obj) { return sum + obj->m_dxQuantities->size(); });
  m_vec_dxQuantities->reserve(dxQuantities_size);
  int dxQuantities_index = 0;
  for (size_t i = 0, size = m_data->size(); i != size; ++i) {
    (*m_data)[i].trackerHits_begin = trackerHits_index;
    (*m_data)[i].trackerHits_end += trackerHits_index;
    trackerHits_index = (*m_data)[i].trackerHits_end;
    for (const auto& it : (*m_rel_trackerHits_tmp[i])) {
      if (it.getObjectID().index == podio::ObjectID::untracked) {
        throw std::runtime_error("Trying to persistify relation \"trackerHits\" that points to untracked object");
      }
      m_refCollections[0]->emplace_back(it.getObjectID());
    }

    (*m_data)[i].tracks_begin = tracks_index;
    (*m_data)[i].tracks_end += tracks_index;
    tracks_index = (*m_data)[i].tracks_end;
    for (const auto& it : (*m_rel_tracks_tmp[i])) {
      if (it.getObjectID().index == podio::ObjectID::untracked) {
        throw std::runtime_error("Trying to persistify relation \"tracks\" that points to untracked object");
      }
      m_refCollections[1]->emplace_back(it.getObjectID());
    }

    (*m_data)[i].subdetectorHitNumbers_begin = subdetectorHitNumbers_index;
    (*m_data)[i].subdetectorHitNumbers_end += subdetectorHitNumbers_index;
    subdetectorHitNumbers_index = (*m_data)[i].subdetectorHitNumbers_end;
    for (const auto& it : (*m_vecs_subdetectorHitNumbers[i])) { m_vec_subdetectorHitNumbers->push_back(it); }

    (*m_data)[i].trackStates_begin = trackStates_index;
    (*m_data)[i].trackStates_end += trackStates_index;
    trackStates_index = (*m_data)[i].trackStates_end;
    for (const auto& it : (*m_vecs_trackStates[i])) { m_vec_trackStates->push_back(it); }

    (*m_data)[i].dxQuantities_begin = dxQuantities_index;
    (*m_data)[i].dxQuantities_end += dxQuantities_index;
    dxQuantities_index = (*m_data)[i].dxQuantities_end;
    for (const auto& it : (*m_vecs_dxQuantities[i])) { m_vec_dxQuantities->push_back(it); }

  }
}

void TrackCollectionData::prepareAfterRead(uint32_t collectionID) {
  int index = 0;
  for (const auto& data : *m_data) {
    auto obj = new TrackObj({index, collectionID}, data);

    obj->m_trackerHits = m_rel_trackerHits.get();
    obj->m_tracks = m_rel_tracks.get();
    obj->m_subdetectorHitNumbers = m_vec_subdetectorHitNumbers.get();
    obj->m_trackStates = m_vec_trackStates.get();
    obj->m_dxQuantities = m_vec_dxQuantities.get();
    entries.emplace_back(obj);
    ++index;
  }

  // at this point we could clear the I/O data buffer, but we keep them intact
  // because then we can save a call to prepareForWrite
}


void TrackCollectionData::createRelations(TrackObj* obj) {
  // We take ownership of these here
  m_rel_trackerHits_tmp.emplace_back(obj->m_trackerHits);
  // We take ownership of these here
  m_rel_tracks_tmp.emplace_back(obj->m_tracks);
  // We take ownership of these here
  m_vecs_subdetectorHitNumbers.emplace_back(obj->m_subdetectorHitNumbers);
  // We take ownership of these here
  m_vecs_trackStates.emplace_back(obj->m_trackStates);
  // We take ownership of these here
  m_vecs_dxQuantities.emplace_back(obj->m_dxQuantities);
}

bool TrackCollectionData::setReferences(const podio::ICollectionProvider* collectionProvider, bool isSubsetColl) {
  if (isSubsetColl) {
    for (const auto& id : *m_refCollections[0]) {
      podio::CollectionBase* coll = nullptr;
      extension::TrackObj* obj = nullptr;
      if (collectionProvider->get(id.collectionID, coll)) {
        auto* tmp_coll = static_cast<extension::TrackCollection*>(coll);
        obj = tmp_coll->m_storage.entries[id.index];
      }
      entries.push_back(obj);
    }
    return true; // TODO: check success, how?
  }

  // Normal collections have to resolve all relations
  for (unsigned int i = 0, size = m_refCollections[0]->size(); i != size; ++i) {
    const auto id = (*m_refCollections[0])[i];
    if (id.index != podio::ObjectID::invalid) {
      podio::CollectionBase* coll = nullptr;
      if (!collectionProvider->get(id.collectionID, coll)) {
        m_rel_trackerHits->emplace_back(extension::TrackerHit::makeEmpty());
        continue;
      }
      podio::detail::addMultiRelation(*m_rel_trackerHits, coll, id);
    } else {
      m_rel_trackerHits->emplace_back(extension::TrackerHit::makeEmpty());
    }
  }

  for (unsigned int i = 0, size = m_refCollections[1]->size(); i != size; ++i) {
    const auto id = (*m_refCollections[1])[i];
    if (id.index != podio::ObjectID::invalid) {
      podio::CollectionBase* coll = nullptr;
      if (!collectionProvider->get(id.collectionID, coll)) {
        m_rel_tracks->emplace_back(extension::Track::makeEmpty());
        continue;
      }
      podio::detail::addMultiRelation(*m_rel_tracks, coll, id);
    } else {
      m_rel_tracks->emplace_back(extension::Track::makeEmpty());
    }
  }


  return true; // TODO: check success, how?
}

void TrackCollectionData::makeSubsetCollection() {
  // Subset collections do not need all the data buffers that normal
  // collections need, so we can free them here
  m_vecmem_info.clear();

  m_data.reset(nullptr);

  m_rel_trackerHits.reset(nullptr);
  m_rel_tracks.reset(nullptr);
  m_vec_subdetectorHitNumbers.reset(nullptr);
  m_vec_trackStates.reset(nullptr);
  m_vec_dxQuantities.reset(nullptr);

  // Subset collections need one vector of ObjectIDs for I/O purposes.
  m_refCollections.resize(1);
  m_refCollections[0] = std::make_unique<std::vector<podio::ObjectID>>();
}


/// \endcond
} // namespace extension

