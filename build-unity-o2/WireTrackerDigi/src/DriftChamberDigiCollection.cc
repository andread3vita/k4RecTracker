// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/DriftChamberDigiCollection.h"
#include "extension/DatamodelDefinition.h"



#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


DriftChamberDigiCollection::DriftChamberDigiCollection(DriftChamberDigiCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

DriftChamberDigiCollection::~DriftChamberDigiCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

DriftChamberDigi DriftChamberDigiCollection::operator[](std::size_t index) const {
  return DriftChamberDigi(m_storage.entries[index]);
}

DriftChamberDigi DriftChamberDigiCollection::at(std::size_t index) const {
  return DriftChamberDigi(m_storage.entries.at(index));
}

MutableDriftChamberDigi DriftChamberDigiCollection::operator[](std::size_t index) {
  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableDriftChamberDigi DriftChamberDigiCollection::at(std::size_t index) {
  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t DriftChamberDigiCollection::size() const {
  return m_storage.entries.size();
}

std::size_t DriftChamberDigiCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool DriftChamberDigiCollection::empty() const {
  return m_storage.entries.empty();
}

void DriftChamberDigiCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableDriftChamberDigi DriftChamberDigiCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new DriftChamberDigiObj());

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(obj));
}

void DriftChamberDigiCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void DriftChamberDigiCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void DriftChamberDigiCollection::prepareAfterRead() {
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

bool DriftChamberDigiCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void DriftChamberDigiCollection::push_back(const MutableDriftChamberDigi& object) {
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

    push_back(DriftChamberDigi(object));
  }
}

void DriftChamberDigiCollection::push_back(const DriftChamberDigi& object) {
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

podio::CollectionWriteBuffers DriftChamberDigiCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<std::uint64_t> DriftChamberDigiCollection::cellID(const size_t nElem) const {
  std::vector<std::uint64_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.cellID);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> DriftChamberDigiCollection::leftPosition(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.leftPosition);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> DriftChamberDigiCollection::rightPosition(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.rightPosition);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiCollection::time(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.time);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiCollection::eDep(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDep);
  }
  return tmp;
}

std::vector<float> DriftChamberDigiCollection::eDepError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDepError);
  }
  return tmp;
}

std::vector<std::uint32_t> DriftChamberDigiCollection::clusterCount(const size_t nElem) const {
  std::vector<std::uint32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.clusterCount);
  }
  return tmp;
}

std::vector<std::int32_t> DriftChamberDigiCollection::type(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.type);
  }
  return tmp;
}

std::vector<std::int32_t> DriftChamberDigiCollection::quality(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.quality);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> DriftChamberDigiCollection::position(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.position);
  }
  return tmp;
}


size_t DriftChamberDigiCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT DriftChamberDigiCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createDriftChamberDigiBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::DriftChamberDigiCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new DriftChamberDigiDataContainer;
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
    DriftChamberDigiCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<DriftChamberDigiCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::DriftChamberDigiDataContainer*>(buffers.data);
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
bool registerDriftChamberDigiCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::DriftChamberDigiCollection", extension::meta::schemaVersion, createDriftChamberDigiBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::DriftChamberDigiCollection",
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

const auto registeredDriftChamberDigiCollection = registerDriftChamberDigiCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiCollection& collection) {
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


DriftChamberDigiCollectionIterator::reference DriftChamberDigiCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index]) };
}

DriftChamberDigiCollectionIterator::pointer DriftChamberDigiCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiCollectionIterator& DriftChamberDigiCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiCollectionIterator DriftChamberDigiCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiCollectionIterator& DriftChamberDigiCollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiCollectionIterator DriftChamberDigiCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiCollectionIterator& DriftChamberDigiCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiCollectionIterator DriftChamberDigiCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiCollectionIterator operator+(DriftChamberDigiCollectionIterator::difference_type n, const DriftChamberDigiCollectionIterator& it) {
    return it + n;
}

DriftChamberDigiCollectionIterator& DriftChamberDigiCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiCollectionIterator DriftChamberDigiCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiCollectionIterator::reference DriftChamberDigiCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiCollectionIterator::difference_type DriftChamberDigiCollectionIterator::operator-(const DriftChamberDigiCollectionIterator& other) const {
    return m_index - other.m_index;
}



DriftChamberDigiMutableCollectionIterator::reference DriftChamberDigiMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index]) };
}

DriftChamberDigiMutableCollectionIterator::pointer DriftChamberDigiMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index]);
  return &m_object;
}

DriftChamberDigiMutableCollectionIterator& DriftChamberDigiMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

DriftChamberDigiMutableCollectionIterator DriftChamberDigiMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

DriftChamberDigiMutableCollectionIterator& DriftChamberDigiMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

DriftChamberDigiMutableCollectionIterator DriftChamberDigiMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

DriftChamberDigiMutableCollectionIterator& DriftChamberDigiMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

DriftChamberDigiMutableCollectionIterator DriftChamberDigiMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

DriftChamberDigiMutableCollectionIterator operator+(DriftChamberDigiMutableCollectionIterator::difference_type n, const DriftChamberDigiMutableCollectionIterator& it) {
    return it + n;
}

DriftChamberDigiMutableCollectionIterator& DriftChamberDigiMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

DriftChamberDigiMutableCollectionIterator DriftChamberDigiMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

DriftChamberDigiMutableCollectionIterator::reference DriftChamberDigiMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<DriftChamberDigiObj>((*m_collection)[m_index + n]) };
}

DriftChamberDigiMutableCollectionIterator::difference_type DriftChamberDigiMutableCollectionIterator::operator-(const DriftChamberDigiMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const DriftChamberDigiCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:      cellID:leftPosition:rightPosition:        time:        eDep:   eDepError:clusterCount:        type:     quality:    position:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getCellID() << " "
      << std::setw(12) << el.getLeftPosition() << " "
      << std::setw(12) << el.getRightPosition() << " "
      << std::setw(12) << el.getTime() << " "
      << std::setw(12) << el.getEDep() << " "
      << std::setw(12) << el.getEDepError() << " "
      << std::setw(12) << el.getClusterCount() << " "
      << std::setw(12) << el.getType() << " "
      << std::setw(12) << el.getQuality() << " "
      << std::setw(12) << el.getPosition() << " "
      << std::endl;




  }

  o.flags(old_flags);
  return o;
}

void DriftChamberDigiCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

