// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiV2Collection_H
#define EXTENSION_DriftChamberDigiV2Collection_H

// datamodel specific includes
#include "extension/DriftChamberDigiV2.h"
#include "extension/MutableDriftChamberDigiV2.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2CollectionData.h"

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

class DriftChamberDigiV2CollectionIterator {
public:
  using value_type = DriftChamberDigiV2;
  using difference_type = ptrdiff_t;
  using reference = DriftChamberDigiV2;
  using pointer = DriftChamberDigiV2*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiV2CollectionIterator(size_t index, const DriftChamberDigiV2ObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiV2CollectionIterator() = default;

  DriftChamberDigiV2CollectionIterator(const DriftChamberDigiV2CollectionIterator&) = default;
  DriftChamberDigiV2CollectionIterator(DriftChamberDigiV2CollectionIterator&&) = default;
  DriftChamberDigiV2CollectionIterator& operator=(const DriftChamberDigiV2CollectionIterator&) = default;
  DriftChamberDigiV2CollectionIterator& operator=(DriftChamberDigiV2CollectionIterator&&) = default;
  ~DriftChamberDigiV2CollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiV2CollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiV2CollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiV2CollectionIterator& operator++();
  DriftChamberDigiV2CollectionIterator operator++(int);
  DriftChamberDigiV2CollectionIterator& operator--();
  DriftChamberDigiV2CollectionIterator operator--(int);
  DriftChamberDigiV2CollectionIterator& operator+=(difference_type n);
  DriftChamberDigiV2CollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiV2CollectionIterator operator+(difference_type n, const DriftChamberDigiV2CollectionIterator& it);
  DriftChamberDigiV2CollectionIterator& operator-=(difference_type n);
  DriftChamberDigiV2CollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiV2CollectionIterator& other) const;

private:
  size_t m_index{0};
  DriftChamberDigiV2 m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr} };
  const DriftChamberDigiV2ObjPointerContainer* m_collection{nullptr};
};


class DriftChamberDigiV2MutableCollectionIterator {
public:
  using value_type = DriftChamberDigiV2;
  using difference_type = ptrdiff_t;
  using reference = MutableDriftChamberDigiV2;
  using pointer = MutableDriftChamberDigiV2*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  DriftChamberDigiV2MutableCollectionIterator(size_t index, const DriftChamberDigiV2ObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}), m_collection(collection) {}
  DriftChamberDigiV2MutableCollectionIterator() = default;

  DriftChamberDigiV2MutableCollectionIterator(const DriftChamberDigiV2MutableCollectionIterator&) = default;
  DriftChamberDigiV2MutableCollectionIterator(DriftChamberDigiV2MutableCollectionIterator&&) = default;
  DriftChamberDigiV2MutableCollectionIterator& operator=(const DriftChamberDigiV2MutableCollectionIterator&) = default;
  DriftChamberDigiV2MutableCollectionIterator& operator=(DriftChamberDigiV2MutableCollectionIterator&&) = default;
  ~DriftChamberDigiV2MutableCollectionIterator() = default;

  constexpr auto operator<=>(const DriftChamberDigiV2MutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const DriftChamberDigiV2MutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  DriftChamberDigiV2MutableCollectionIterator& operator++();
  DriftChamberDigiV2MutableCollectionIterator operator++(int);
  DriftChamberDigiV2MutableCollectionIterator& operator--();
  DriftChamberDigiV2MutableCollectionIterator operator--(int);
  DriftChamberDigiV2MutableCollectionIterator& operator+=(difference_type n);
  DriftChamberDigiV2MutableCollectionIterator operator+(difference_type n) const;
  friend DriftChamberDigiV2MutableCollectionIterator operator+(difference_type n, const DriftChamberDigiV2MutableCollectionIterator& it);
  DriftChamberDigiV2MutableCollectionIterator& operator-=(difference_type n);
  DriftChamberDigiV2MutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const DriftChamberDigiV2MutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableDriftChamberDigiV2 m_object { podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr} };
  const DriftChamberDigiV2ObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class DriftChamberDigiV2Collection : public podio::CollectionBase {
public:
  using value_type = DriftChamberDigiV2;
  using mutable_type = MutableDriftChamberDigiV2;
  using const_iterator = DriftChamberDigiV2CollectionIterator;
  using iterator = DriftChamberDigiV2MutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  DriftChamberDigiV2Collection() = default;
  DriftChamberDigiV2Collection(DriftChamberDigiV2CollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  DriftChamberDigiV2Collection(const DriftChamberDigiV2Collection& ) = delete;
  DriftChamberDigiV2Collection& operator=(const DriftChamberDigiV2Collection& ) = delete;
  DriftChamberDigiV2Collection(DriftChamberDigiV2Collection&&) = default;
  DriftChamberDigiV2Collection& operator=(DriftChamberDigiV2Collection&&) = default;

  ~DriftChamberDigiV2Collection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  DriftChamberDigiV2Collection(std::from_range_t, R&& range);
#endif

  /// Construct a DriftChamberDigiV2 collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static DriftChamberDigiV2Collection from(R&& range);

  constexpr static std::string_view typeName = "extension::DriftChamberDigiV2Collection";
  constexpr static std::string_view valueTypeName = "extension::DriftChamberDigiV2";
  constexpr static std::string_view dataTypeName = "extension::DriftChamberDigiV2Data";

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
  MutableDriftChamberDigiV2 create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableDriftChamberDigiV2 create(Args&&... args);

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
  DriftChamberDigiV2 operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableDriftChamberDigiV2 operator[](std::size_t index);
  /// Returns the const object of given index
  DriftChamberDigiV2 at(std::size_t index) const;
  /// Returns the object of given index
  MutableDriftChamberDigiV2 at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableDriftChamberDigiV2& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const DriftChamberDigiV2& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (DriftChamberDigiV2Obj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
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
  std::vector<edm4hep::Vector3d> position(const size_t nElem = 0) const;
  std::vector<edm4hep::Vector3d> directionSW(const size_t nElem = 0) const;
  std::vector<float> distanceToWire(const size_t nElem = 0) const;
  std::vector<std::uint32_t> nCluster(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class DriftChamberDigiV2CollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable DriftChamberDigiV2CollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiV2Collection& v);

template<typename... Args>
MutableDriftChamberDigiV2 DriftChamberDigiV2Collection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new DriftChamberDigiV2Obj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, DriftChamberDigiV2Data} constructor
  obj->m_nElectrons = new std::vector<std::uint16_t>();
  m_storage.createRelations(obj);
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<DriftChamberDigiV2Collection::value_type> R>
DriftChamberDigiV2Collection::DriftChamberDigiV2Collection(std::from_range_t, R&& range)
  : DriftChamberDigiV2Collection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<DriftChamberDigiV2Collection::value_type> R>
DriftChamberDigiV2Collection DriftChamberDigiV2Collection::from(R&& range) {
  DriftChamberDigiV2Collection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiV2Collection& collection);
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
constexpr std::string_view extension::DriftChamberDigiV2Collection::typeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiV2Collection::typeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::valueTypeName;
constexpr std::string_view extension::DriftChamberDigiV2Collection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
