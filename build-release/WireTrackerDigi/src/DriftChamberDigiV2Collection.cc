// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/DriftChamberDigiV2Collection.h"
#include "extension/DatamodelDefinition.h"


#include <numeric>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


DriftChamberDigiV2Collection::DriftChamberDigiV2Collection(DriftChamberDigiV2CollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

DriftChamberDigiV2Collection::~DriftChamberDigiV2Collection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

DriftChamberDigiV2 DriftChamberDigiV2Collection::operator[](std::size_t index) const {
  return DriftChamberDigiV2(m_storage.entries[index]);
}

DriftChamberDigiV2 DriftChamberDigiV2Collection::at(std::size_t index) const {
  return DriftChamberDigiV2(m_storage.entries.at(index));
}

MutableDriftChamberDigiV2 DriftChamberDigiV2Collection::operator[](std::size_t index) {
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableDriftChamberDigiV2 DriftChamberDigiV2Collection::at(std::size_t index) {
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t DriftChamberDigiV2Collection::size() const {
  return m_storage.entries.size();
}

std::size_t DriftChamberDigiV2Collection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool DriftChamberDigiV2Collection::empty() const {
  return m_storage.entries.empty();
}

void DriftChamberDigiV2Collection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableDriftChamberDigiV2 DriftChamberDigiV2Collection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new DriftChamberDigiV2Obj());
  m_storage.createRelations(obj);

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(obj));
}

void DriftChamberDigiV2Collection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void DriftChamberDigiV2Collection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void DriftChamberDigiV2Collection::prepareAfterRead() {
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

bool DriftChamberDigiV2Collection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void DriftChamberDigiV2Collection::push_back(const MutableDriftChamberDigiV2& object) {
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

    push_back(DriftChamberDigiV2(object));
  }
}

void DriftChamberDigiV2Collection::push_back(const DriftChamberDigiV2& object) {
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

podio::CollectionWriteBuffers DriftChamberDigiV2Collection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<std::uint64_t> DriftChamberDigiV2Collection::cellID(const size_t nElem) const {
  std::vector<std::uint64_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.cellID);
  }
  return tmp;
}

std::vector<std::int32_t> DriftChamberDigiV2Collection::type(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.type);
  }
  return tmp;
}

std::vector<std::int32_t> DriftChamberDigiV2Collection::quality(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.quality);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiV2Collection::time(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.time);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiV2Collection::eDep(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDep);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiV2Collection::eDepError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDepError);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> DriftChamberDigiV2Collection::position(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.position);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> DriftChamberDigiV2Collection::directionSW(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.directionSW);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiV2Collection::distanceToWire(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.distanceToWire);
  }
  return tmp;
}

std::vector<std::uint32_t> DriftChamberDigiV2Collection::nCluster(const size_t nElem) const {
  std::vector<std::uint32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.nCluster);
  }
  return tmp;
}


size_t DriftChamberDigiV2Collection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT DriftChamberDigiV2Collection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createDriftChamberDigiV2Buffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::DriftChamberDigiV2Collection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new DriftChamberDigiV2DataContainer;
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
    readBuffers.vectorMembers->reserve(1);
    readBuffers.vectorMembers->emplace_back("std::uint16_t", new std::vector<std::uint16_t>);
  }

  readBuffers.createCollection = [](podio::CollectionReadBuffers&& buffers, bool isSubsetColl) {
    DriftChamberDigiV2CollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<DriftChamberDigiV2Collection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::DriftChamberDigiV2DataContainer*>(buffers.data);
      delete static_cast<std::vector<std::uint16_t>*>((*buffers.vectorMembers)[0].second);
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
bool registerDriftChamberDigiV2Collection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::DriftChamberDigiV2Collection", extension::meta::schemaVersion, createDriftChamberDigiV2Buffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::DriftChamberDigiV2Collection",
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

const auto registeredDriftChamberDigiV2Collection = registerDriftChamberDigiV2Collection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiV2Collection& collection) {
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


DriftChamberDigiV2CollectionIterator::reference DriftChamberDigiV2CollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index]) };
}

DriftChamberDigiV2CollectionIterator::pointer DriftChamberDigiV2CollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiV2CollectionIterator& DriftChamberDigiV2CollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiV2CollectionIterator DriftChamberDigiV2CollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiV2CollectionIterator& DriftChamberDigiV2CollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiV2CollectionIterator DriftChamberDigiV2CollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiV2CollectionIterator& DriftChamberDigiV2CollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiV2CollectionIterator DriftChamberDigiV2CollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiV2CollectionIterator operator+(DriftChamberDigiV2CollectionIterator::difference_type n, const DriftChamberDigiV2CollectionIterator& it) {
    return it + n;
}

DriftChamberDigiV2CollectionIterator& DriftChamberDigiV2CollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiV2CollectionIterator DriftChamberDigiV2CollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiV2CollectionIterator::reference DriftChamberDigiV2CollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiV2CollectionIterator::difference_type DriftChamberDigiV2CollectionIterator::operator-(const DriftChamberDigiV2CollectionIterator& other) const {
    return m_index - other.m_index;
}



DriftChamberDigiV2MutableCollectionIterator::reference DriftChamberDigiV2MutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index]) };
}

DriftChamberDigiV2MutableCollectionIterator::pointer DriftChamberDigiV2MutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiV2MutableCollectionIterator& DriftChamberDigiV2MutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiV2MutableCollectionIterator DriftChamberDigiV2MutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiV2MutableCollectionIterator& DriftChamberDigiV2MutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiV2MutableCollectionIterator DriftChamberDigiV2MutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiV2MutableCollectionIterator& DriftChamberDigiV2MutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiV2MutableCollectionIterator DriftChamberDigiV2MutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiV2MutableCollectionIterator operator+(DriftChamberDigiV2MutableCollectionIterator::difference_type n, const DriftChamberDigiV2MutableCollectionIterator& it) {
    return it + n;
}

DriftChamberDigiV2MutableCollectionIterator& DriftChamberDigiV2MutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiV2MutableCollectionIterator DriftChamberDigiV2MutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiV2MutableCollectionIterator::reference DriftChamberDigiV2MutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiV2MutableCollectionIterator::difference_type DriftChamberDigiV2MutableCollectionIterator::operator-(const DriftChamberDigiV2MutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const DriftChamberDigiV2Collection& v) {
  const auto old_flags = o.flags();
  o << "          id:      cellID:        type:     quality:        time:        eDep:   eDepError:    position: directionSW:distanceToWire:    nCluster:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getCellID() << " "
      << std::setw(12) << el.getType() << " "
      << std::setw(12) << el.getQuality() << " "
      << std::setw(12) << el.getTime() << " "
      << std::setw(12) << el.getEDep() << " "
      << std::setw(12) << el.getEDepError() << " "
      << std::setw(12) << el.getPosition() << " "
      << std::setw(12) << el.getDirectionSW() << " "
      << std::setw(12) << el.getDistanceToWire() << " "
      << std::setw(12) << el.getNCluster() << " "
      << std::endl;



    o << "      nElectrons : ";
    for (unsigned j = 0, N = el.nElectrons_size(); j < N; ++j) {
      o << el.getNElectrons(j) << " ";
    }
    o << std::endl;

  }

  o.flags(old_flags);
  return o;
}

void DriftChamberDigiV2Collection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

