// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociationCollection_H
#define EXTENSION_MCRecoDriftChamberDigiAssociationCollection_H

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiAssociation.h"
#include "extension/MutableMCRecoDriftChamberDigiAssociation.h"
#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollectionData.h"

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

class MCRecoDriftChamberDigiAssociationCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using difference_type = ptrdiff_t;
  using reference = MCRecoDriftChamberDigiAssociation;
  using pointer = MCRecoDriftChamberDigiAssociation*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiAssociationCollectionIterator(size_t index, const MCRecoDriftChamberDigiAssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiAssociationCollectionIterator() = default;

  MCRecoDriftChamberDigiAssociationCollectionIterator(const MCRecoDriftChamberDigiAssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator(MCRecoDriftChamberDigiAssociationCollectionIterator&&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator=(const MCRecoDriftChamberDigiAssociationCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator=(MCRecoDriftChamberDigiAssociationCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiAssociationCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiAssociationCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociationCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator++();
  MCRecoDriftChamberDigiAssociationCollectionIterator operator++(int);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator--();
  MCRecoDriftChamberDigiAssociationCollectionIterator operator--(int);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiAssociationCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiAssociationCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiAssociationCollectionIterator& it);
  MCRecoDriftChamberDigiAssociationCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiAssociationCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiAssociationCollectionIterator& other) const;

private:
  size_t m_index{0};
  MCRecoDriftChamberDigiAssociation m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiAssociationObjPointerContainer* m_collection{nullptr};
};


class MCRecoDriftChamberDigiAssociationMutableCollectionIterator {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using difference_type = ptrdiff_t;
  using reference = MutableMCRecoDriftChamberDigiAssociation;
  using pointer = MutableMCRecoDriftChamberDigiAssociation*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(size_t index, const MCRecoDriftChamberDigiAssociationObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}), m_collection(collection) {}
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator() = default;

  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator(MCRecoDriftChamberDigiAssociationMutableCollectionIterator&&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator=(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator&) = default;
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator=(MCRecoDriftChamberDigiAssociationMutableCollectionIterator&&) = default;
  ~MCRecoDriftChamberDigiAssociationMutableCollectionIterator() = default;

  constexpr auto operator<=>(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator++();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator++(int);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator--();
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator--(int);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator+=(difference_type n);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator+(difference_type n) const;
  friend MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator+(difference_type n, const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& it);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator& operator-=(difference_type n);
  MCRecoDriftChamberDigiAssociationMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const MCRecoDriftChamberDigiAssociationMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableMCRecoDriftChamberDigiAssociation m_object { podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr} };
  const MCRecoDriftChamberDigiAssociationObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class MCRecoDriftChamberDigiAssociationCollection : public podio::CollectionBase {
public:
  using value_type = MCRecoDriftChamberDigiAssociation;
  using mutable_type = MutableMCRecoDriftChamberDigiAssociation;
  using const_iterator = MCRecoDriftChamberDigiAssociationCollectionIterator;
  using iterator = MCRecoDriftChamberDigiAssociationMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  MCRecoDriftChamberDigiAssociationCollection() = default;
  MCRecoDriftChamberDigiAssociationCollection(MCRecoDriftChamberDigiAssociationCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  MCRecoDriftChamberDigiAssociationCollection(const MCRecoDriftChamberDigiAssociationCollection& ) = delete;
  MCRecoDriftChamberDigiAssociationCollection& operator=(const MCRecoDriftChamberDigiAssociationCollection& ) = delete;
  MCRecoDriftChamberDigiAssociationCollection(MCRecoDriftChamberDigiAssociationCollection&&) = default;
  MCRecoDriftChamberDigiAssociationCollection& operator=(MCRecoDriftChamberDigiAssociationCollection&&) = default;

  ~MCRecoDriftChamberDigiAssociationCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  MCRecoDriftChamberDigiAssociationCollection(std::from_range_t, R&& range);
#endif

  /// Construct a MCRecoDriftChamberDigiAssociation collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static MCRecoDriftChamberDigiAssociationCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::MCRecoDriftChamberDigiAssociationCollection";
  constexpr static std::string_view valueTypeName = "extension::MCRecoDriftChamberDigiAssociation";
  constexpr static std::string_view dataTypeName = "extension::MCRecoDriftChamberDigiAssociationData";

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
  MutableMCRecoDriftChamberDigiAssociation create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableMCRecoDriftChamberDigiAssociation create(Args&&... args);

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
  MCRecoDriftChamberDigiAssociation operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableMCRecoDriftChamberDigiAssociation operator[](std::size_t index);
  /// Returns the const object of given index
  MCRecoDriftChamberDigiAssociation at(std::size_t index) const;
  /// Returns the object of given index
  MutableMCRecoDriftChamberDigiAssociation at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableMCRecoDriftChamberDigiAssociation& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const MCRecoDriftChamberDigiAssociation& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (MCRecoDriftChamberDigiAssociationObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
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
  friend class MCRecoDriftChamberDigiAssociationCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable MCRecoDriftChamberDigiAssociationCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociationCollection& v);

template<typename... Args>
MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociationCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new MCRecoDriftChamberDigiAssociationObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiAssociationCollection::value_type> R>
MCRecoDriftChamberDigiAssociationCollection::MCRecoDriftChamberDigiAssociationCollection(std::from_range_t, R&& range)
  : MCRecoDriftChamberDigiAssociationCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<MCRecoDriftChamberDigiAssociationCollection::value_type> R>
MCRecoDriftChamberDigiAssociationCollection MCRecoDriftChamberDigiAssociationCollection::from(R&& range) {
  MCRecoDriftChamberDigiAssociationCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociationCollection& collection);
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
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::typeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::valueTypeName;
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociationCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
