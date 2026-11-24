// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/DatamodelDefinition.h"



#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


DriftChamberDigiLocalCollection::DriftChamberDigiLocalCollection(DriftChamberDigiLocalCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

DriftChamberDigiLocalCollection::~DriftChamberDigiLocalCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

DriftChamberDigiLocal DriftChamberDigiLocalCollection::operator[](std::size_t index) const {
  return DriftChamberDigiLocal(m_storage.entries[index]);
}

DriftChamberDigiLocal DriftChamberDigiLocalCollection::at(std::size_t index) const {
  return DriftChamberDigiLocal(m_storage.entries.at(index));
}

MutableDriftChamberDigiLocal DriftChamberDigiLocalCollection::operator[](std::size_t index) {
  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableDriftChamberDigiLocal DriftChamberDigiLocalCollection::at(std::size_t index) {
  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t DriftChamberDigiLocalCollection::size() const {
  return m_storage.entries.size();
}

std::size_t DriftChamberDigiLocalCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool DriftChamberDigiLocalCollection::empty() const {
  return m_storage.entries.empty();
}

void DriftChamberDigiLocalCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableDriftChamberDigiLocal DriftChamberDigiLocalCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new DriftChamberDigiLocalObj());

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(obj));
}

void DriftChamberDigiLocalCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void DriftChamberDigiLocalCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void DriftChamberDigiLocalCollection::prepareAfterRead() {
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

bool DriftChamberDigiLocalCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void DriftChamberDigiLocalCollection::push_back(const MutableDriftChamberDigiLocal& object) {
  // We have to do different things here depending on whether this is a
  // subset collection or not. A normal collection cannot collect objects
  // that are already part of another collection, while a subset collection
  // can only collect such objects
  if (!m_isSubsetColl) {
    auto obj = object.m_obj;
    if (obj->id.index == podio::ObjectID::untracked) {
      obj->id = {static_cast<int>(m_storage.entries.size()), m_collectionID};
      m_storage.entries.push_back(obj.release());
    } else {
      throw std::invalid_argument("Object already in a collection. Cannot add it to a second collection");
    }
  } else {

    push_back(DriftChamberDigiLocal(object));
  }
}

void DriftChamberDigiLocalCollection::push_back(const DriftChamberDigiLocal& object) {
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

podio::CollectionWriteBuffers DriftChamberDigiLocalCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<std::uint64_t> DriftChamberDigiLocalCollection::cellID(const size_t nElem) const {
  std::vector<std::uint64_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.cellID);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiLocalCollection::distanceToWire(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.distanceToWire);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiLocalCollection::zPositionAlongWire(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.zPositionAlongWire);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiLocalCollection::time(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.time);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiLocalCollection::eDep(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDep);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiLocalCollection::eDepError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDepError);
  }
  return tmp;
}

std::vector<std::uint32_t> DriftChamberDigiLocalCollection::clusterCount(const size_t nElem) const {
  std::vector<std::uint32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.clusterCount);
  }
  return tmp;
}


size_t DriftChamberDigiLocalCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT DriftChamberDigiLocalCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createDriftChamberDigiLocalBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::DriftChamberDigiLocalCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new DriftChamberDigiLocalDataContainer;
  // The number of ObjectID vectors is either 1 or the sum of OneToMany and
  // OneToOne relations
  const auto nRefs = isSubset ? 1 : 0 + 0;
  readBuffers.references = new podio::CollRefCollection(nRefs);
  for (auto& ref : *readBuffers.references) {
    // Make sure to place usable buffer pointers here
    ref = std::make_unique<std::vector<podio::ObjectID>>();
  }

  readBuffers.vectorMembers = new podio::VectorMembersInfo();
  if (!isSubset) {
    readBuffers.vectorMembers->reserve(0);
  }

  readBuffers.createCollection = [](podio::CollectionReadBuffers&& buffers, bool isSubsetColl) {
    DriftChamberDigiLocalCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<DriftChamberDigiLocalCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::DriftChamberDigiLocalDataContainer*>(buffers.data);
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
bool registerDriftChamberDigiLocalCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::DriftChamberDigiLocalCollection", extension::meta::schemaVersion, createDriftChamberDigiLocalBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::DriftChamberDigiLocalCollection",
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

const auto registeredDriftChamberDigiLocalCollection = registerDriftChamberDigiLocalCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiLocalCollection& collection) {
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


DriftChamberDigiLocalCollectionIterator::reference DriftChamberDigiLocalCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index]) };
}

DriftChamberDigiLocalCollectionIterator::pointer DriftChamberDigiLocalCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiLocalCollectionIterator& DriftChamberDigiLocalCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiLocalCollectionIterator DriftChamberDigiLocalCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiLocalCollectionIterator& DriftChamberDigiLocalCollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiLocalCollectionIterator DriftChamberDigiLocalCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiLocalCollectionIterator& DriftChamberDigiLocalCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiLocalCollectionIterator DriftChamberDigiLocalCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiLocalCollectionIterator operator+(DriftChamberDigiLocalCollectionIterator::difference_type n, const DriftChamberDigiLocalCollectionIterator& it) {
    return it + n;
}

DriftChamberDigiLocalCollectionIterator& DriftChamberDigiLocalCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiLocalCollectionIterator DriftChamberDigiLocalCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiLocalCollectionIterator::reference DriftChamberDigiLocalCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiLocalCollectionIterator::difference_type DriftChamberDigiLocalCollectionIterator::operator-(const DriftChamberDigiLocalCollectionIterator& other) const {
    return m_index - other.m_index;
}



DriftChamberDigiLocalMutableCollectionIterator::reference DriftChamberDigiLocalMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index]) };
}

DriftChamberDigiLocalMutableCollectionIterator::pointer DriftChamberDigiLocalMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiLocalMutableCollectionIterator& DriftChamberDigiLocalMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiLocalMutableCollectionIterator DriftChamberDigiLocalMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiLocalMutableCollectionIterator& DriftChamberDigiLocalMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiLocalMutableCollectionIterator DriftChamberDigiLocalMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiLocalMutableCollectionIterator& DriftChamberDigiLocalMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiLocalMutableCollectionIterator DriftChamberDigiLocalMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiLocalMutableCollectionIterator operator+(DriftChamberDigiLocalMutableCollectionIterator::difference_type n, const DriftChamberDigiLocalMutableCollectionIterator& it) {
    return it + n;
}

DriftChamberDigiLocalMutableCollectionIterator& DriftChamberDigiLocalMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiLocalMutableCollectionIterator DriftChamberDigiLocalMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiLocalMutableCollectionIterator::reference DriftChamberDigiLocalMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiLocalMutableCollectionIterator::difference_type DriftChamberDigiLocalMutableCollectionIterator::operator-(const DriftChamberDigiLocalMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const DriftChamberDigiLocalCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:      cellID:distanceToWire:zPositionAlongWire:        time:        eDep:   eDepError:clusterCount:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getCellID() << " "
      << std::setw(12) << el.getDistanceToWire() << " "
      << std::setw(12) << el.getZPositionAlongWire() << " "
      << std::setw(12) << el.getTime() << " "
      << std::setw(12) << el.getEDep() << " "
      << std::setw(12) << el.getEDepError() << " "
      << std::setw(12) << el.getClusterCount() << " "
      << std::endl;




  }

  o.flags(old_flags);
  return o;
}

void DriftChamberDigiLocalCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

