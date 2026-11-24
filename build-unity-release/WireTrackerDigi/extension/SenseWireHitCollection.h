// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitCollection_H
#define EXTENSION_SenseWireHitCollection_H

// datamodel specific includes
#include "extension/SenseWireHit.h"
#include "extension/MutableSenseWireHit.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHitCollectionData.h"

// podio specific includes
#include "podio/ICollectionProvider.h"
#include "podio/CollectionBase.h"
#include "podio/detail/Pythonizations.h"
#include "podio/utilities/TypeHelpers.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

#include <string_view>
#include <vector>
#include <algorithm>
#include <ostream>
#include <mutex>
#include <memory>
#include <cstddef>

namespace podio {
  struct RelationNames;
}

namespace extension {

/// \addtogroup extension_user_facing
/// @{

class SenseWireHitCollectionIterator {
public:
  using value_type = SenseWireHit;
  using difference_type = ptrdiff_t;
  using reference = SenseWireHit;
  using pointer = SenseWireHit*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitCollectionIterator(size_t index, const SenseWireHitObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}), m_collection(collection) {}
  SenseWireHitCollectionIterator() = default;

  SenseWireHitCollectionIterator(const SenseWireHitCollectionIterator&) = default;
  SenseWireHitCollectionIterator(SenseWireHitCollectionIterator&&) = default;
  SenseWireHitCollectionIterator& operator=(const SenseWireHitCollectionIterator&) = default;
  SenseWireHitCollectionIterator& operator=(SenseWireHitCollectionIterator&&) = default;
  ~SenseWireHitCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitCollectionIterator& operator++();
  SenseWireHitCollectionIterator operator++(int);
  SenseWireHitCollectionIterator& operator--();
  SenseWireHitCollectionIterator operator--(int);
  SenseWireHitCollectionIterator& operator+=(difference_type n);
  SenseWireHitCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitCollectionIterator operator+(difference_type n, const SenseWireHitCollectionIterator& it);
  SenseWireHitCollectionIterator& operator-=(difference_type n);
  SenseWireHitCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitCollectionIterator& other) const;

private:
  size_t m_index{0};
  SenseWireHit m_object { podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr} };
  const SenseWireHitObjPointerContainer* m_collection{nullptr};
};


class SenseWireHitMutableCollectionIterator {
public:
  using value_type = SenseWireHit;
  using difference_type = ptrdiff_t;
  using reference = MutableSenseWireHit;
  using pointer = MutableSenseWireHit*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitMutableCollectionIterator(size_t index, const SenseWireHitObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}), m_collection(collection) {}
  SenseWireHitMutableCollectionIterator() = default;

  SenseWireHitMutableCollectionIterator(const SenseWireHitMutableCollectionIterator&) = default;
  SenseWireHitMutableCollectionIterator(SenseWireHitMutableCollectionIterator&&) = default;
  SenseWireHitMutableCollectionIterator& operator=(const SenseWireHitMutableCollectionIterator&) = default;
  SenseWireHitMutableCollectionIterator& operator=(SenseWireHitMutableCollectionIterator&&) = default;
  ~SenseWireHitMutableCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitMutableCollectionIterator& operator++();
  SenseWireHitMutableCollectionIterator operator++(int);
  SenseWireHitMutableCollectionIterator& operator--();
  SenseWireHitMutableCollectionIterator operator--(int);
  SenseWireHitMutableCollectionIterator& operator+=(difference_type n);
  SenseWireHitMutableCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitMutableCollectionIterator operator+(difference_type n, const SenseWireHitMutableCollectionIterator& it);
  SenseWireHitMutableCollectionIterator& operator-=(difference_type n);
  SenseWireHitMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableSenseWireHit m_object { podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr} };
  const SenseWireHitObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class SenseWireHitCollection : public podio::CollectionBase {
public:
  using value_type = SenseWireHit;
  using mutable_type = MutableSenseWireHit;
  using const_iterator = SenseWireHitCollectionIterator;
  using iterator = SenseWireHitMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  SenseWireHitCollection() = default;
  SenseWireHitCollection(SenseWireHitCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  SenseWireHitCollection(const SenseWireHitCollection& ) = delete;
  SenseWireHitCollection& operator=(const SenseWireHitCollection& ) = delete;
  SenseWireHitCollection(SenseWireHitCollection&&) = default;
  SenseWireHitCollection& operator=(SenseWireHitCollection&&) = default;

  ~SenseWireHitCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  SenseWireHitCollection(std::from_range_t, R&& range);
#endif

  /// Construct a SenseWireHit collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static SenseWireHitCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::SenseWireHitCollection";
  constexpr static std::string_view valueTypeName = "extension::SenseWireHit";
  constexpr static std::string_view dataTypeName = "extension::SenseWireHitData";

  /// Invalidates all iterators.
  void clear() final;

  /// Cppyy protocol to setup the pythonizations for this class. Not to be called directly.
  static void __cppjit_pythonize__(PyObject* klass, const std::string& name) {
          podio::detail::pythonizations::pythonize_subscript(klass, name);
  }
  /// Pythonizations for ROOT <= 6.40 before the switch to cppjit
  static void __cppyy_pythonize__(PyObject* klass, const std::string& name) {
        __cppjit_pythonize__(klass, name);
  }

  /// Print this collection to the passed stream
  void print(std::ostream& os=std::cout, bool flush=true) const final;

  /// Append a new object to the collection, and return this object.
  /// Does not invalidate iterators.
  MutableSenseWireHit create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableSenseWireHit create(Args&&... args);

  /// number of elements in the collection
  std::size_t size() const final;

  /// maximal number of elements in the collection
  std::size_t max_size() const final;

  /// Is the collection empty
  bool empty() const final;

  /// fully qualified type name
  const std::string_view getTypeName() const final { return typeName; }
  /// fully qualified type name of elements - with namespace
  const std::string_view getValueTypeName() const final { return valueTypeName; }
  /// fully qualified type name of stored POD elements - with namespace
  const std::string_view getDataTypeName() const final { return dataTypeName; }
  /// schema version
  podio::SchemaVersionT getSchemaVersion() const final;

  bool isSubsetCollection() const final {
    return m_isSubsetColl;
  }

  void setSubsetCollection(bool setSubset=true) final;

  /// Returns the const object of given index
  SenseWireHit operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableSenseWireHit operator[](std::size_t index);
  /// Returns the const object of given index
  SenseWireHit at(std::size_t index) const;
  /// Returns the object of given index
  MutableSenseWireHit at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableSenseWireHit& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const SenseWireHit& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (SenseWireHitObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
      );
    }
  }

  uint32_t getID() const final {
    return m_collectionID;
  }

  /// check if the collection has a valid ID
  bool hasID() const final {
    return getID() != static_cast<uint32_t>(podio::ObjectID::untracked) &&
        getID() != static_cast<uint32_t>(podio::ObjectID::invalid);
  }

  size_t getDatamodelRegistryIndex() const final;

  // support for the iterator protocol
  iterator begin() {
    return iterator(0, &m_storage.entries);
  }
  const_iterator begin() const {
    return const_iterator(0, &m_storage.entries);
  }
  const_iterator cbegin() const {
    return begin();
  }
  iterator end() {
    return iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator end() const {
    return const_iterator(m_storage.entries.size(), &m_storage.entries);
  }
  const_iterator cend() const {
    return end();
  }
  // reverse iterators
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const {
    return rbegin();
  }
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const {
    return rend();
  }


  std::vector<std::uint64_t> cellID(const size_t nElem = 0) const;
  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<std::int32_t> quality(const size_t nElem = 0) const;
  std::vector<float> time(const size_t nElem = 0) const;
  std::vector<float> eDep(const size_t nElem = 0) const;
  std::vector<float> eDepError(const size_t nElem = 0) const;
  std::vector<float> wireStereoAngle(const size_t nElem = 0) const;
  std::vector<float> wireAzimuthalAngle(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> position(const size_t nElem = 0) const;
  std::vector<double> positionAlongWireError(const size_t nElem = 0) const;
  std::vector<float> distanceToWire(const size_t nElem = 0) const;
  std::vector<float> distanceToWireError(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class SenseWireHitCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable SenseWireHitCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitCollection& v);

template<typename... Args>
MutableSenseWireHit SenseWireHitCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new SenseWireHitObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, SenseWireHitData} constructor
  obj->m_nElectrons = new std::vector<std::uint16_t>();
  m_storage.createRelations(obj);
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<SenseWireHitCollection::value_type> R>
SenseWireHitCollection::SenseWireHitCollection(std::from_range_t, R&& range)
  : SenseWireHitCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<SenseWireHitCollection::value_type> R>
SenseWireHitCollection SenseWireHitCollection::from(R&& range) {
  SenseWireHitCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitCollection& collection);
#endif

/// @}
} // namespace extension


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitCollection::typeName;
constexpr std::string_view extension::SenseWireHitCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitCollection::typeName;
constexpr std::string_view extension::SenseWireHitCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
