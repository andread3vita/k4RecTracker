// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "podio/CollectionBufferFactory.h"
#include "podio/SchemaEvolution.h"

#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/DatamodelDefinition.h"


#include "extension/DriftChamberDigiCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

// standard includes
#include <iomanip>
#include <mutex>
#include <stdexcept>

namespace extension {


MCRecoDriftChamberDigiAssociationCollection::MCRecoDriftChamberDigiAssociationCollection(MCRecoDriftChamberDigiAssociationCollectionData&& data, bool isSubsetColl) :
  m_isSubsetColl(isSubsetColl), m_storage(std::move(data)) {}

MCRecoDriftChamberDigiAssociationCollection::~MCRecoDriftChamberDigiAssociationCollection() {
  // Need to tell the storage how to clean-up
  m_storage.clear(m_isSubsetColl);
}

MCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::operator[](std::size_t index) const {
  return MCRecoDriftChamberDigiAssociation(m_storage.entries[index]);
}

MCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::at(std::size_t index) const {
  return MCRecoDriftChamberDigiAssociation(m_storage.entries.at(index));
}

MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::operator[](std::size_t index) {
  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(m_storage.entries[index]));
}

MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::at(std::size_t index) {
  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(m_storage.entries.at(index)));
}

std::size_t MCRecoDriftChamberDigiAssociationCollection::size() const {
  return m_storage.entries.size();
}

std::size_t MCRecoDriftChamberDigiAssociationCollection::max_size() const {
  const auto maxStorage = m_storage.entries.max_size();
  if (!m_isSubsetColl) {
    // non-subset collections shouldn't have more elements than the maximum index of ObjectID
    const auto maxIndex = std::numeric_limits<decltype(podio::ObjectID::index)>::max();
    return std::min<size_t>(maxIndex, maxStorage);
  }
  return maxStorage;
}

bool MCRecoDriftChamberDigiAssociationCollection::empty() const {
  return m_storage.entries.empty();
}

void MCRecoDriftChamberDigiAssociationCollection::setSubsetCollection(bool setSubset) {
  if (m_isSubsetColl != setSubset && !m_storage.entries.empty()) {
    throw std::logic_error("Cannot change the character of a collection that already contains elements");
  }

  if (setSubset) {
    m_storage.makeSubsetCollection();
  }
  m_isSubsetColl = setSubset;
}

MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::create() {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }

  auto obj = m_storage.entries.emplace_back(new MCRecoDriftChamberDigiAssociationObj());

  obj->id = {int(m_storage.entries.size() - 1), m_collectionID};
  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(obj));
}

void MCRecoDriftChamberDigiAssociationCollection::clear() {
  m_storage.clear(m_isSubsetColl);
  m_isPrepared = false;
}

void MCRecoDriftChamberDigiAssociationCollection::prepareForWrite() const {
  std::lock_guard lock{*m_storageMtx};
  if (m_isPrepared) {
    return;
  }
  m_storage.prepareForWrite(m_isSubsetColl);
  m_isPrepared = true;
}

void MCRecoDriftChamberDigiAssociationCollection::prepareAfterRead() {
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

bool MCRecoDriftChamberDigiAssociationCollection::setReferences(const podio::ICollectionProvider* collectionProvider) {
  return m_storage.setReferences(collectionProvider, m_isSubsetColl);
}

void MCRecoDriftChamberDigiAssociationCollection::push_back(const MutableMCRecoDriftChamberDigiAssociation& object) {
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

    push_back(MCRecoDriftChamberDigiAssociation(object));
  }
}

void MCRecoDriftChamberDigiAssociationCollection::push_back(const MCRecoDriftChamberDigiAssociation& object) {
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

podio::CollectionWriteBuffers MCRecoDriftChamberDigiAssociationCollection::getBuffers() {
  return m_storage.getCollectionBuffers(m_isSubsetColl);
}

std::vector<float> MCRecoDriftChamberDigiAssociationCollection::weight(const size_t nElem) const {
  std::vector<float> tmp;
  const auto valid_size = nElem != 0 ? std::min(nElem, m_storage.entries.size()) : m_storage.entries.size();
  tmp.reserve(valid_size);
  for (size_t i = 0; i < valid_size; ++i) {
    tmp.emplace_back(m_storage.entries[i]->data.weight);
  }
  return tmp;
}


size_t MCRecoDriftChamberDigiAssociationCollection::getDatamodelRegistryIndex() const {
  return extension::meta::DatamodelRegistryIndex::value();
}

podio::SchemaVersionT MCRecoDriftChamberDigiAssociationCollection::getSchemaVersion() const {
  return extension::meta::schemaVersion;
}

// anonymous namespace for registration with the CollectionBufferFactory. This
// ensures that we don't have to make up arbitrary namespace names here, since
// none of this is publicly visible
namespace {
 
podio::CollectionReadBuffers createMCRecoDriftChamberDigiAssociationBuffers(bool isSubset) {
  auto readBuffers = podio::CollectionReadBuffers{};
  readBuffers.type = "extension::MCRecoDriftChamberDigiAssociationCollection";
  readBuffers.schemaVersion = extension::meta::schemaVersion;
  readBuffers.data = isSubset ? nullptr : new MCRecoDriftChamberDigiAssociationDataContainer;
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
    MCRecoDriftChamberDigiAssociationCollectionData data(std::move(buffers), isSubsetColl);
    return std::make_unique<MCRecoDriftChamberDigiAssociationCollection>(std::move(data), isSubsetColl);
  };

  readBuffers.deleteBuffers = [](podio::CollectionReadBuffers& buffers) {
    if (buffers.data) {
      // If we have data then we are not a subset collection and we have to
      // clean up all type erased buffers by casting them back to something that
      // we can delete
      delete static_cast<extension::MCRecoDriftChamberDigiAssociationDataContainer*>(buffers.data);
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
bool registerMCRecoDriftChamberDigiAssociationCollection() {
  const static auto reg = []() {
    auto& factory = podio::CollectionBufferFactory::mutInstance();
    factory.registerCreationFunc("extension::MCRecoDriftChamberDigiAssociationCollection", extension::meta::schemaVersion, createMCRecoDriftChamberDigiAssociationBuffers);

    // Make the SchemaEvolution aware of the current version by
    // registering a no-op function for this and all preceding versions
    // will be overridden whenever an explicit action is required
    for (unsigned int schemaVersion=1; schemaVersion< extension::meta::schemaVersion+1; ++schemaVersion) {
      podio::SchemaEvolution::mutInstance().registerEvolutionFunc(
        "extension::MCRecoDriftChamberDigiAssociationCollection",
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

const auto registeredMCRecoDriftChamberDigiAssociationCollection = registerMCRecoDriftChamberDigiAssociationCollection();
} // namespace


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociationCollection& collection) {
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


MCRecoDriftChamberDigiAssociationCollectionIterator::reference MCRecoDriftChamberDigiAssociationCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index]) };
}

MCRecoDriftChamberDigiAssociationCollectionIterator::pointer MCRecoDriftChamberDigiAssociationCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index]);
  return &m_object;
}

MCRecoDriftChamberDigiAssociationCollectionIterator& MCRecoDriftChamberDigiAssociationCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

MCRecoDriftChamberDigiAssociationCollectionIterator MCRecoDriftChamberDigiAssociationCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

MCRecoDriftChamberDigiAssociationCollectionIterator& MCRecoDriftChamberDigiAssociationCollectionIterator::operator--() {
  --m_index;
  return *this;
}

MCRecoDriftChamberDigiAssociationCollectionIterator MCRecoDriftChamberDigiAssociationCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

MCRecoDriftChamberDigiAssociationCollectionIterator& MCRecoDriftChamberDigiAssociationCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

MCRecoDriftChamberDigiAssociationCollectionIterator MCRecoDriftChamberDigiAssociationCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

MCRecoDriftChamberDigiAssociationCollectionIterator operator+(MCRecoDriftChamberDigiAssociationCollectionIterator::difference_type n, const MCRecoDriftChamberDigiAssociationCollectionIterator& it) {
    return it + n;
}

MCRecoDriftChamberDigiAssociationCollectionIterator& MCRecoDriftChamberDigiAssociationCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

MCRecoDriftChamberDigiAssociationCollectionIterator MCRecoDriftChamberDigiAssociationCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

MCRecoDriftChamberDigiAssociationCollectionIterator::reference MCRecoDriftChamberDigiAssociationCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index + n]) };
}

MCRecoDriftChamberDigiAssociationCollectionIterator::difference_type MCRecoDriftChamberDigiAssociationCollectionIterator::operator-(const MCRecoDriftChamberDigiAssociationCollectionIterator& other) const {
    return m_index - other.m_index;
}



MCRecoDriftChamberDigiAssociationMutableCollectionIterator::reference MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator*() const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index]) };
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator::pointer MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator->() {
  m_object.m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index]);
  return &m_object;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator& MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator++() {
  ++m_index;
  return *this;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator++(int) {
  auto copy = *this;
  ++m_index;
  return copy;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator& MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator--() {
  --m_index;
  return *this;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator--(int) {
  auto copy = *this;
  --m_index;
  return copy;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator& MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator+=(difference_type n) {
  m_index += n;
  return *this;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator+(difference_type n) const {
  auto copy = *this;
  copy += n;
  return copy;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator+(MCRecoDriftChamberDigiAssociationMutableCollectionIterator::difference_type n, const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& it) {
    return it + n;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator& MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator-=(difference_type n) {
  m_index -= n;
  return *this;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator-(difference_type n) const {
  auto copy = *this;
  copy -= n;
  return copy;
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator::reference MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator[](difference_type n) const {
  return reference{ podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>((*m_collection)[m_index + n]) };
}

MCRecoDriftChamberDigiAssociationMutableCollectionIterator::difference_type MCRecoDriftChamberDigiAssociationMutableCollectionIterator::operator-(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& other) const {
    return m_index - other.m_index;
}



std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociationCollection& v) {
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

void MCRecoDriftChamberDigiAssociationCollection::print(std::ostream& os, bool flush) const {
  os << *this;
  if (flush) {
    os.flush();
  }
}


} // namespace extension

