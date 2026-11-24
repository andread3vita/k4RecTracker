// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/SenseWireHitCollection.h"
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


SenseWireHitCollection::SenseWireHitCollection(SenseWireHitCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

SenseWireHitCollection::~SenseWireHitCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

SenseWireHit SenseWireHitCollection::operator[](std::size_t index) const {
  return SenseWireHit(m_storage.entries[index]);
}

SenseWireHit SenseWireHitCollection::at(std::size_t index) const {
  return SenseWireHit(m_storage.entries.at(index));
}

MutableSenseWireHit SenseWireHitCollection::operator[](std::size_t index) {
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableSenseWireHit SenseWireHitCollection::at(std::size_t index) {
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t SenseWireHitCollection::size() const {
  return m_storage.entries.size();
}

std::size_t SenseWireHitCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool SenseWireHitCollection::empty() const {
  return m_storage.entries.empty();
}

void SenseWireHitCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableSenseWireHit SenseWireHitCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new SenseWireHitObj());
  m_storage.createRelations(obj);

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(obj));
}

void SenseWireHitCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void SenseWireHitCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void SenseWireHitCollection::prepareAfterRead() {
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

bool SenseWireHitCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void SenseWireHitCollection::push_back(const MutableSenseWireHit& object) {
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

    push_back(SenseWireHit(object));
  }
}

void SenseWireHitCollection::push_back(const SenseWireHit& object) {
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

podio::CollectionWriteBuffers SenseWireHitCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<std::uint64_t> SenseWireHitCollection::cellID(const size_t nElem) const {
  std::vector<std::uint64_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.cellID);
  }
  return tmp;
}

std::vector<std::int32_t> SenseWireHitCollection::type(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.type);
  }
  return tmp;
}

std::vector<std::int32_t> SenseWireHitCollection::quality(const size_t nElem) const {
  std::vector<std::int32_t> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.quality);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::time(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.time);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::eDep(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDep);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::eDepError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.eDepError);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::wireStereoAngle(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.wireStereoAngle);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::wireAzimuthalAngle(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.wireAzimuthalAngle);
  }
  return tmp;
}

std::vector<edm4hep::Vector3d> SenseWireHitCollection::position(const size_t nElem) const {
  std::vector<edm4hep::Vector3d> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.position);
  }
  return tmp;
}

std::vector<double> SenseWireHitCollection::positionAlongWireError(const size_t nElem) const {
  std::vector<double> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.positionAlongWireError);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::distanceToWire(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.distanceToWire);
  }
  return tmp;
}

std::vector<float> SenseWireHitCollection::distanceToWireError(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.distanceToWireError);
  }
  return tmp;
}


size_t SenseWireHitCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT SenseWireHitCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createSenseWireHitBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::SenseWireHitCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new SenseWireHitDataContainer;
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
    SenseWireHitCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<SenseWireHitCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::SenseWireHitDataContainer*>(buffers.data);
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
bool registerSenseWireHitCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::SenseWireHitCollection", extension::meta::schemaVersion, createSenseWireHitBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::SenseWireHitCollection",
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

const auto registeredSenseWireHitCollection = registerSenseWireHitCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitCollection& collection) {
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


SenseWireHitCollectionIterator::reference SenseWireHitCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index]) };
}

SenseWireHitCollectionIterator::pointer SenseWireHitCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index]);
  return &m_object;
}

SenseWireHitCollectionIterator& SenseWireHitCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

SenseWireHitCollectionIterator SenseWireHitCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

SenseWireHitCollectionIterator& SenseWireHitCollectionIterator::operator--() {
  --m_index;
  return *this;
}

SenseWireHitCollectionIterator SenseWireHitCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

SenseWireHitCollectionIterator& SenseWireHitCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

SenseWireHitCollectionIterator SenseWireHitCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

SenseWireHitCollectionIterator operator+(SenseWireHitCollectionIterator::difference_type n, const SenseWireHitCollectionIterator& it) {
    return it + n;
}

SenseWireHitCollectionIterator& SenseWireHitCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

SenseWireHitCollectionIterator SenseWireHitCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

SenseWireHitCollectionIterator::reference SenseWireHitCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index + n]) };
}

SenseWireHitCollectionIterator::difference_type SenseWireHitCollectionIterator::operator-(const SenseWireHitCollectionIterator& other) const {
    return m_index - other.m_index;
}



SenseWireHitMutableCollectionIterator::reference SenseWireHitMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index]) };
}

SenseWireHitMutableCollectionIterator::pointer SenseWireHitMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index]);
  return &m_object;
}

SenseWireHitMutableCollectionIterator& SenseWireHitMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

SenseWireHitMutableCollectionIterator SenseWireHitMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

SenseWireHitMutableCollectionIterator& SenseWireHitMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

SenseWireHitMutableCollectionIterator SenseWireHitMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

SenseWireHitMutableCollectionIterator& SenseWireHitMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

SenseWireHitMutableCollectionIterator SenseWireHitMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

SenseWireHitMutableCollectionIterator operator+(SenseWireHitMutableCollectionIterator::difference_type n, const SenseWireHitMutableCollectionIterator& it) {
    return it + n;
}

SenseWireHitMutableCollectionIterator& SenseWireHitMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

SenseWireHitMutableCollectionIterator SenseWireHitMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

SenseWireHitMutableCollectionIterator::reference SenseWireHitMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<SenseWireHitObj>((*m_collection)[m_index + n]) };
}

SenseWireHitMutableCollectionIterator::difference_type SenseWireHitMutableCollectionIterator::operator-(const SenseWireHitMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const SenseWireHitCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:      cellID:        type:     quality:        time:        eDep:   eDepError:wireStereoAngle:wireAzimuthalAngle:    position:positionAlongWireError:distanceToWire:distanceToWireError:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getCellID() << " "
      << std::setw(12) << el.getType() << " "
      << std::setw(12) << el.getQuality() << " "
      << std::setw(12) << el.getTime() << " "
      << std::setw(12) << el.getEDep() << " "
      << std::setw(12) << el.getEDepError() << " "
      << std::setw(12) << el.getWireStereoAngle() << " "
      << std::setw(12) << el.getWireAzimuthalAngle() << " "
      << std::setw(12) << el.getPosition() << " "
      << std::setw(12) << el.getPositionAlongWireError() << " "
      << std::setw(12) << el.getDistanceToWire() << " "
      << std::setw(12) << el.getDistanceToWireError() << " "
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

void SenseWireHitCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

