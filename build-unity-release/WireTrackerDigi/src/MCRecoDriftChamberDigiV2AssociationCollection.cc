// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/MCRecoDriftChamberDigiV2AssociationCollection.h"
#include "extension/DatamodelDefinition.h"


#include "extension/DriftChamberDigiV2Collection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


MCRecoDriftChamberDigiV2AssociationCollection::MCRecoDriftChamberDigiV2AssociationCollection(MCRecoDriftChamberDigiV2AssociationCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

MCRecoDriftChamberDigiV2AssociationCollection::~MCRecoDriftChamberDigiV2AssociationCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

MCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::operator[](std::size_t index) const {
  return MCRecoDriftChamberDigiV2Association(m_storage.entries[index]);
}

MCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::at(std::size_t index) const {
  return MCRecoDriftChamberDigiV2Association(m_storage.entries.at(index));
}

MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::operator[](std::size_t index) {
  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::at(std::size_t index) {
  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t MCRecoDriftChamberDigiV2AssociationCollection::size() const {
  return m_storage.entries.size();
}

std::size_t MCRecoDriftChamberDigiV2AssociationCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool MCRecoDriftChamberDigiV2AssociationCollection::empty() const {
  return m_storage.entries.empty();
}

void MCRecoDriftChamberDigiV2AssociationCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new MCRecoDriftChamberDigiV2AssociationObj());

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(obj));
}

void MCRecoDriftChamberDigiV2AssociationCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void MCRecoDriftChamberDigiV2AssociationCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void MCRecoDriftChamberDigiV2AssociationCollection::prepareAfterRead() {
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

bool MCRecoDriftChamberDigiV2AssociationCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void MCRecoDriftChamberDigiV2AssociationCollection::push_back(const MutableMCRecoDriftChamberDigiV2Association& object) {
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

    push_back(MCRecoDriftChamberDigiV2Association(object));
  }
}

void MCRecoDriftChamberDigiV2AssociationCollection::push_back(const MCRecoDriftChamberDigiV2Association& object) {
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

podio::CollectionWriteBuffers MCRecoDriftChamberDigiV2AssociationCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<float> MCRecoDriftChamberDigiV2AssociationCollection::weight(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.weight);
  }
  return tmp;
}


size_t MCRecoDriftChamberDigiV2AssociationCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT MCRecoDriftChamberDigiV2AssociationCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createMCRecoDriftChamberDigiV2AssociationBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::MCRecoDriftChamberDigiV2AssociationCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new MCRecoDriftChamberDigiV2AssociationDataContainer;
  // The number of ObjectID vectors is either 1 or the sum of OneToMany and
  // OneToOne relations
  const auto nRefs = isSubset ? 1 : 0 + 2;
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
    MCRecoDriftChamberDigiV2AssociationCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<MCRecoDriftChamberDigiV2AssociationCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::MCRecoDriftChamberDigiV2AssociationDataContainer*>(buffers.data);
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
bool registerMCRecoDriftChamberDigiV2AssociationCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::MCRecoDriftChamberDigiV2AssociationCollection", extension::meta::schemaVersion, createMCRecoDriftChamberDigiV2AssociationBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::MCRecoDriftChamberDigiV2AssociationCollection",
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

const auto registeredMCRecoDriftChamberDigiV2AssociationCollection = registerMCRecoDriftChamberDigiV2AssociationCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2AssociationCollection& collection) {
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


MCRecoDriftChamberDigiV2AssociationCollectionIterator::reference MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index]) };
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator::pointer MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index]);
  return &m_object;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator& MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator& MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator--() {
  --m_index;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator& MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator operator+(MCRecoDriftChamberDigiV2AssociationCollectionIterator::difference_type n, const MCRecoDriftChamberDigiV2AssociationCollectionIterator& it) {
    return it + n;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator& MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator::reference MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index + n]) };
}

MCRecoDriftChamberDigiV2AssociationCollectionIterator::difference_type MCRecoDriftChamberDigiV2AssociationCollectionIterator::operator-(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& other) const {
    return m_index - other.m_index;
}



MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::reference MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index]) };
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::pointer MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index]);
  return &m_object;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator+(MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::difference_type n, const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& it) {
    return it + n;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::reference MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>((*m_collection)[m_index + n]) };
}

MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::difference_type MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator::operator-(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2AssociationCollection& v) {
  const auto old_flags = o.flags();
  o << "          id:      weight:" << '\n';

  for (const auto&& el : v) {
    o << std::scientific << std::showpos << std::setw(12) << el.id() << " "
      << std::setw(12) << el.getWeight() << " "
      << std::endl;


    o << "      digi : ";
    o << el.getDigi().id() << std::endl;
    o << "      sim : ";
    o << el.getSim().id() << std::endl;


  }

  o.flags(old_flags);
  return o;
}

void MCRecoDriftChamberDigiV2AssociationCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

