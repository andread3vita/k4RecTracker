// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2AssociationCollection_H
#define EXTENSION_MCRecoDriftChamberDigiV2AssociationCollection_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiV2Association.h"
#include "extension/MutableMCRecoDriftChamberDigiV2Association.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationCollectionData.h"

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

class MCRecoDriftChamberDigiV2AssociationCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using difference_type = ptrdiff_t;
  using reference = MCRecoDriftChamberDigiV2Association;
  using pointer = MCRecoDriftChamberDigiV2Association*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiV2AssociationCollectionIterator(size_t index, const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiV2AssociationCollectionIterator() = default;

  MCRecoDriftChamberDigiV2AssociationCollectionIterator(const MCRecoDriftChamberDigiV2AssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator(MCRecoDriftChamberDigiV2AssociationCollectionIterator&&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator=(const MCRecoDriftChamberDigiV2AssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator=(MCRecoDriftChamberDigiV2AssociationCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiV2AssociationCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator++();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator++(int);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator--();
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator--(int);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiV2AssociationCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiV2AssociationCollectionIterator& it);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiV2AssociationCollectionIterator& other) const;

private:
  size_t m_index{0};
  MCRecoDriftChamberDigiV2Association m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* m_collection{nullptr};
};


class MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using difference_type = ptrdiff_t;
  using reference = MutableMCRecoDriftChamberDigiV2Association;
  using pointer = MutableMCRecoDriftChamberDigiV2Association*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(size_t index, const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator() = default;

  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator(MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator=(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator=(MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator++();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator++(int);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator--();
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator--(int);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& it);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableMCRecoDriftChamberDigiV2Association m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiV2AssociationObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class MCRecoDriftChamberDigiV2AssociationCollection : public podio::CollectionBase {
public:
  using value_type = MCRecoDriftChamberDigiV2Association;
  using mutable_type = MutableMCRecoDriftChamberDigiV2Association;
  using const_iterator = MCRecoDriftChamberDigiV2AssociationCollectionIterator;
  using iterator = MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  MCRecoDriftChamberDigiV2AssociationCollection() = default;
  MCRecoDriftChamberDigiV2AssociationCollection(MCRecoDriftChamberDigiV2AssociationCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  MCRecoDriftChamberDigiV2AssociationCollection(const MCRecoDriftChamberDigiV2AssociationCollection& ) = delete;
  MCRecoDriftChamberDigiV2AssociationCollection& operator=(const MCRecoDriftChamberDigiV2AssociationCollection& ) = delete;
  MCRecoDriftChamberDigiV2AssociationCollection(MCRecoDriftChamberDigiV2AssociationCollection&&) = default;
  MCRecoDriftChamberDigiV2AssociationCollection& operator=(MCRecoDriftChamberDigiV2AssociationCollection&&) = default;

  ~MCRecoDriftChamberDigiV2AssociationCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  MCRecoDriftChamberDigiV2AssociationCollection(std::from_range_t, R&& range);
#endif

  /// Construct a MCRecoDriftChamberDigiV2Association collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static MCRecoDriftChamberDigiV2AssociationCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::MCRecoDriftChamberDigiV2AssociationCollection";
  constexpr static std::string_view valueTypeName = "extension::MCRecoDriftChamberDigiV2Association";
  constexpr static std::string_view dataTypeName = "extension::MCRecoDriftChamberDigiV2AssociationData";

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
  MutableMCRecoDriftChamberDigiV2Association create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableMCRecoDriftChamberDigiV2Association create(Args&&... args);

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
  MCRecoDriftChamberDigiV2Association operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableMCRecoDriftChamberDigiV2Association operator[](std::size_t index);
  /// Returns the const object of given index
  MCRecoDriftChamberDigiV2Association at(std::size_t index) const;
  /// Returns the object of given index
  MutableMCRecoDriftChamberDigiV2Association at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableMCRecoDriftChamberDigiV2Association& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const MCRecoDriftChamberDigiV2Association& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (MCRecoDriftChamberDigiV2AssociationObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
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


  std::vector<float> weight(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class MCRecoDriftChamberDigiV2AssociationCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable MCRecoDriftChamberDigiV2AssociationCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2AssociationCollection& v);

template<typename... Args>
MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2AssociationCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new MCRecoDriftChamberDigiV2AssociationObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiV2AssociationCollection::value_type> R>
MCRecoDriftChamberDigiV2AssociationCollection::MCRecoDriftChamberDigiV2AssociationCollection(std::from_range_t, R&& range)
  : MCRecoDriftChamberDigiV2AssociationCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiV2AssociationCollection::value_type> R>
MCRecoDriftChamberDigiV2AssociationCollection MCRecoDriftChamberDigiV2AssociationCollection::from(R&& range) {
  MCRecoDriftChamberDigiV2AssociationCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2AssociationCollection& collection);
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
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiV2AssociationCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
