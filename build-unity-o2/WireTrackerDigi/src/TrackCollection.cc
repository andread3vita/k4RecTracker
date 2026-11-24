// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/TrackCollection.h"
#include "extension/DatamodelDefinition.h"


#include "extension/DriftChamberDigiCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/TrackerHit.h"
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include <numeric>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


TrackCollection::TrackCollection(TrackCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

TrackCollection::~TrackCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

Track TrackCollection::operator[](std::size_t index) const {
  return Track(m_storage.entries[index]);
}

Track TrackCollection::at(std::size_t index) const {
  return Track(m_storage.entries.at(index));
}

MutableTrack TrackCollection::operator[](std::size_t index) {
  return MutableTrack(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableTrack TrackCollection::at(std::size_t index) {
  return MutableTrack(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t TrackCollection::size() const {
  return m_storage.entries.size();
}

std::size_t TrackCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool TrackCollection::empty() const {
  return m_storage.entries.empty();
}

void TrackCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableTrack TrackCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new TrackObj());
  m_storage.createRelations(obj);

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableTrack(podio::utils::MaybeSharedPtr(obj));
}

void TrackCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void TrackCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void TrackCollection::prepareAfterRead() {
  // No need to go through this again if we have already done it
  if (m_isPrepared) {
    return;
  }

  if (!m_isSubsetColl) {
    // Subset collections do not store any data that would require post-processing
    m_storage.prepareAfterRead(m_collectionID);
  }
  // Preparing a collection doesn't affect the underlying I/O buffers, so this
  // collection is still prepared
  m_isPrepared = true;
}

bool TrackCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void TrackCollection::push_back(const MutableTrack& object) {
  // We have to do different things here depending on whether this is a
  // subset collection or not. A normal collection cannot collect objects
  // that are already part of another collection, while a subset collection
  // can only collect such objects
  if (!m_isSubsetColl) {
    auto obj = object.m_obj;
    if (obj->id.index == podio::ObjectID::untracked) {
      obj->id = {static_cast<int>(m_storage.entries.size()), m_collectionID};
      m_storage.entries.push_back(obj.release());
      m_storage.createRelations(obj.get());
    } else {
      throw std::invalid_argument("Object already in a collection. Cannot add it to a second collection");
    }
  } else {

    push_back(Track(object));
  }
}

void TrackCollection::push_back(const Track& object) {
  if (!m_isSubsetColl) {
    throw std::invalid_argument("Immutable objects can only be added to subset collections");
  }
  auto obj = object.m_obj;
  if (obj->id.index < 0) {
    // This path is only possible if we arrive here from an untracked Mutable object
    throw std::invalid_argument("Object needs to be tracked by another collection in order for it to be storable in a subset collection");
  }
  m_storage.entries.push_back(obj.release());
}

podio::CollectionWriteBuffers TrackCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<std::int32_t> TrackCollection::type(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.type);
  }
  return tmp;
}

std::vector<float> TrackCollection::chi2(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.chi2);
  }
  return tmp;
}

std::vector<std::int32_t> TrackCollection::ndf(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.ndf);
  }
  return tmp;
}

std::vector<float> TrackCollection::dEdx(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.dEdx);
  }
  return tmp;
}

std::vector<float> TrackCollection::dEdxError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.dEdxError);
  }
  return tmp;
}

std::vector<float> TrackCollection::radiusOfInnermostHit(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.radiusOfInnermostHit);
  }
  return tmp;
}


size_t TrackCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT TrackCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createTrackBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::TrackCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new TrackDataContainer;
  // The number of ObjectID vectors is either 1 or the sum of OneToMany and
  // OneToOne relations
  const auto nRefs = isSubset ? 1 : 2 + 0;
  readBuffers.references = new podio::CollRefCollection(nRefs);
  for (auto& ref : *readBuffers.references) {
    // Make sure to place usable buffer pointers here
    ref = std::make_unique<std::vector<podio::ObjectID>>();
  }

  readBuffers.vectorMembers = new podio::VectorMembersInfo();
  if (!isSubset) {
    readBuffers.vectorMembers->reserve(3);
    readBuffers.vectorMembers->emplace_back("std::int32_t", new std::vector<std::int32_t>);
    readBuffers.vectorMembers->emplace_back("edm4hep::TrackState", new std::vector<edm4hep::TrackState>);
    readBuffers.vectorMembers->emplace_back("edm4hep::Quantity", new std::vector<edm4hep::Quantity>);
  }

  readBuffers.createCollection = [](podio::CollectionReadBuffers&& buffers, bool isSubsetColl) {
    TrackCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<TrackCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::TrackDataContainer*>(buffers.data);
      delete static_cast<std::vector<std::int32_t>*>((*buffers.vectorMembers)[0].second);
      delete static_cast<std::vector<edm4hep::TrackState>*>((*buffers.vectorMembers)[1].second);
      delete static_cast<std::vector<edm4hep::Quantity>*>((*buffers.vectorMembers)[2].second);
      buffers.data = nullptr;

    }
    delete buffers.references;
    buffers.references = nullptr;
    delete buffers.vectorMembers;
    buffers.vectorMembers = nullptr;
  };

  return readBuffers;
}



// The usual trick with an IIFE and a static variable inside a function and then
// making sure to call that function during shared library loading
bool registerTrackCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::TrackCollection", extension::meta::schemaVersion, createTrackBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::TrackCollection",
        schemaVersion,
        extension::meta::schemaVersion,
        podio::SchemaEvolution::noOpSchemaEvolution,
        podio::SchemaEvolution::Priority::AutoGenerated
      );
    }


    return true;
  }();
  return reg;
}

const auto registeredTrackCollection = registerTrackCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const TrackCollection& collection) {
  j = nlohmann::json::array();
  const auto subset = collection.isSubsetCollection();
  for (auto&& elem : collection) {
    if (subset) {
      j.emplace_back(elem.id());
    } else {
      j.emplace_back(elem);
    }
  }
}
#endif


TrackCollectionIterator::reference TrackCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index]) };
}

TrackCollectionIterator::pointer TrackCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index]);
  return &m_object;
}

TrackCollectionIterator& TrackCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

TrackCollectionIterator TrackCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

TrackCollectionIterator& TrackCollectionIterator::operator--() {
  --m_index;
  return *this;
}

TrackCollectionIterator TrackCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

TrackCollectionIterator& TrackCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

TrackCollectionIterator TrackCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

TrackCollectionIterator operator+(TrackCollectionIterator::difference_type n, const TrackCollectionIterator& it) {
    return it + n;
}

TrackCollectionIterator& TrackCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

TrackCollectionIterator TrackCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

TrackCollectionIterator::reference TrackCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index + n]) };
}

TrackCollectionIterator::difference_type TrackCollectionIterator::operator-(const TrackCollectionIterator& other) const {
    return m_index - other.m_index;
}



TrackMutableCollectionIterator::reference TrackMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index]) };
}

TrackMutableCollectionIterator::pointer TrackMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index]);
  return &m_object;
}

TrackMutableCollectionIterator& TrackMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

TrackMutableCollectionIterator TrackMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

TrackMutableCollectionIterator& TrackMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

TrackMutableCollectionIterator TrackMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

TrackMutableCollectionIterator& TrackMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

TrackMutableCollectionIterator TrackMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

TrackMutableCollectionIterator operator+(TrackMutableCollectionIterator::difference_type n, const TrackMutableCollectionIterator& it) {
    return it + n;
}

TrackMutableCollectionIterator& TrackMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

TrackMutableCollectionIterator TrackMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

TrackMutableCollectionIterator::reference TrackMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<TrackObj>((*m_collection)[m_index + n]) };
}

TrackMutableCollectionIterator::difference_type TrackMutableCollectionIterator::operator-(const TrackMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const TrackCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:        type:        chi2:         ndf:        dEdx:   dEdxError:radiusOfInnermostHit:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getType() << " "
      << std::setw(12) << el.getChi2() << " "
      << std::setw(12) << el.getNdf() << " "
      << std::setw(12) << el.getDEdx() << " "
      << std::setw(12) << el.getDEdxError() << " "
      << std::setw(12) << el.getRadiusOfInnermostHit() << " "
      << std::endl;

    o << "      trackerHits : ";
    for (unsigned j = 0, N = el.trackerHits_size(); j < N; ++j) {
      o << el.getTrackerHits(j).id() << " ";
    }
    o << std::endl;
    o << "      tracks : ";
    for (unsigned j = 0, N = el.tracks_size(); j < N; ++j) {
      o << el.getTracks(j).id() << " ";
    }
    o << std::endl;


    o << "      subdetectorHitNumbers : ";
    for (unsigned j = 0, N = el.subdetectorHitNumbers_size(); j < N; ++j) {
      o << el.getSubdetectorHitNumbers(j) << " ";
    }
    o << std::endl;
    o << "      trackStates : ";
    for (unsigned j = 0, N = el.trackStates_size(); j < N; ++j) {
      o << el.getTrackStates(j) << " ";
    }
    o << std::endl;
    o << "      dxQuantities : ";
    for (unsigned j = 0, N = el.dxQuantities_size(); j < N; ++j) {
      o << el.getDxQuantities(j) << " ";
    }
    o << std::endl;

  }

  o.flags(old_flags);
  return o;
}

void TrackCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

