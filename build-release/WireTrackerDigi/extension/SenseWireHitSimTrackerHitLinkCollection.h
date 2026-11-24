// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLinkCollection_H
#define EXTENSION_SenseWireHitSimTrackerHitLinkCollection_H

// datamodel specific includes
#include "extension/SenseWireHitSimTrackerHitLink.h"
#include "extension/MutableSenseWireHitSimTrackerHitLink.h"
#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollectionData.h"

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

class SenseWireHitSimTrackerHitLinkCollectionIterator {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using difference_type = ptrdiff_t;
  using reference = SenseWireHitSimTrackerHitLink;
  using pointer = SenseWireHitSimTrackerHitLink*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitSimTrackerHitLinkCollectionIterator(size_t index, const SenseWireHitSimTrackerHitLinkObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}), m_collection(collection) {}
  SenseWireHitSimTrackerHitLinkCollectionIterator() = default;

  SenseWireHitSimTrackerHitLinkCollectionIterator(const SenseWireHitSimTrackerHitLinkCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator(SenseWireHitSimTrackerHitLinkCollectionIterator&&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator=(const SenseWireHitSimTrackerHitLinkCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator=(SenseWireHitSimTrackerHitLinkCollectionIterator&&) = default;
  ~SenseWireHitSimTrackerHitLinkCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitSimTrackerHitLinkCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLinkCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator++();
  SenseWireHitSimTrackerHitLinkCollectionIterator operator++(int);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator--();
  SenseWireHitSimTrackerHitLinkCollectionIterator operator--(int);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator+=(difference_type n);
  SenseWireHitSimTrackerHitLinkCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitSimTrackerHitLinkCollectionIterator operator+(difference_type n, const SenseWireHitSimTrackerHitLinkCollectionIterator& it);
  SenseWireHitSimTrackerHitLinkCollectionIterator& operator-=(difference_type n);
  SenseWireHitSimTrackerHitLinkCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitSimTrackerHitLinkCollectionIterator& other) const;

private:
  size_t m_index{0};
  SenseWireHitSimTrackerHitLink m_object { podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr} };
  const SenseWireHitSimTrackerHitLinkObjPointerContainer* m_collection{nullptr};
};


class SenseWireHitSimTrackerHitLinkMutableCollectionIterator {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using difference_type = ptrdiff_t;
  using reference = MutableSenseWireHitSimTrackerHitLink;
  using pointer = MutableSenseWireHitSimTrackerHitLink*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(size_t index, const SenseWireHitSimTrackerHitLinkObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}), m_collection(collection) {}
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator() = default;

  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator(SenseWireHitSimTrackerHitLinkMutableCollectionIterator&&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator=(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator&) = default;
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator=(SenseWireHitSimTrackerHitLinkMutableCollectionIterator&&) = default;
  ~SenseWireHitSimTrackerHitLinkMutableCollectionIterator() = default;

  constexpr auto operator<=>(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator++();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator++(int);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator--();
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator--(int);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator+=(difference_type n);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator+(difference_type n) const;
  friend SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator+(difference_type n, const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& it);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator& operator-=(difference_type n);
  SenseWireHitSimTrackerHitLinkMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const SenseWireHitSimTrackerHitLinkMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableSenseWireHitSimTrackerHitLink m_object { podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr} };
  const SenseWireHitSimTrackerHitLinkObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class SenseWireHitSimTrackerHitLinkCollection : public podio::CollectionBase {
public:
  using value_type = SenseWireHitSimTrackerHitLink;
  using mutable_type = MutableSenseWireHitSimTrackerHitLink;
  using const_iterator = SenseWireHitSimTrackerHitLinkCollectionIterator;
  using iterator = SenseWireHitSimTrackerHitLinkMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  SenseWireHitSimTrackerHitLinkCollection() = default;
  SenseWireHitSimTrackerHitLinkCollection(SenseWireHitSimTrackerHitLinkCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  SenseWireHitSimTrackerHitLinkCollection(const SenseWireHitSimTrackerHitLinkCollection& ) = delete;
  SenseWireHitSimTrackerHitLinkCollection& operator=(const SenseWireHitSimTrackerHitLinkCollection& ) = delete;
  SenseWireHitSimTrackerHitLinkCollection(SenseWireHitSimTrackerHitLinkCollection&&) = default;
  SenseWireHitSimTrackerHitLinkCollection& operator=(SenseWireHitSimTrackerHitLinkCollection&&) = default;

  ~SenseWireHitSimTrackerHitLinkCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  SenseWireHitSimTrackerHitLinkCollection(std::from_range_t, R&& range);
#endif

  /// Construct a SenseWireHitSimTrackerHitLink collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static SenseWireHitSimTrackerHitLinkCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::SenseWireHitSimTrackerHitLinkCollection";
  constexpr static std::string_view valueTypeName = "extension::SenseWireHitSimTrackerHitLink";
  constexpr static std::string_view dataTypeName = "extension::SenseWireHitSimTrackerHitLinkData";

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
  MutableSenseWireHitSimTrackerHitLink create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableSenseWireHitSimTrackerHitLink create(Args&&... args);

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
  SenseWireHitSimTrackerHitLink operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableSenseWireHitSimTrackerHitLink operator[](std::size_t index);
  /// Returns the const object of given index
  SenseWireHitSimTrackerHitLink at(std::size_t index) const;
  /// Returns the object of given index
  MutableSenseWireHitSimTrackerHitLink at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableSenseWireHitSimTrackerHitLink& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const SenseWireHitSimTrackerHitLink& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (SenseWireHitSimTrackerHitLinkObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
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
  friend class SenseWireHitSimTrackerHitLinkCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable SenseWireHitSimTrackerHitLinkCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitSimTrackerHitLinkCollection& v);

template<typename... Args>
MutableSenseWireHitSimTrackerHitLink SenseWireHitSimTrackerHitLinkCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new SenseWireHitSimTrackerHitLinkObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  return MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<SenseWireHitSimTrackerHitLinkCollection::value_type> R>
SenseWireHitSimTrackerHitLinkCollection::SenseWireHitSimTrackerHitLinkCollection(std::from_range_t, R&& range)
  : SenseWireHitSimTrackerHitLinkCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<SenseWireHitSimTrackerHitLinkCollection::value_type> R>
SenseWireHitSimTrackerHitLinkCollection SenseWireHitSimTrackerHitLinkCollection::from(R&& range) {
  SenseWireHitSimTrackerHitLinkCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitSimTrackerHitLinkCollection& collection);
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
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::typeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::typeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::valueTypeName;
constexpr std::string_view extension::SenseWireHitSimTrackerHitLinkCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
