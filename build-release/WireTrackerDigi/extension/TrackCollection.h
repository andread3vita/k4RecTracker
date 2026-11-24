// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_TrackCollection_H
#define EXTENSION_TrackCollection_H

// datamodel specific includes
#include "extension/Track.h"
#include "extension/MutableTrack.h"
#include "extension/TrackObj.h"
#include "extension/TrackCollectionData.h"

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

class TrackCollectionIterator {
public:
  using value_type = Track;
  using difference_type = ptrdiff_t;
  using reference = Track;
  using pointer = Track*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  TrackCollectionIterator(size_t index, const TrackObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<TrackObj>{nullptr}), m_collection(collection) {}
  TrackCollectionIterator() = default;

  TrackCollectionIterator(const TrackCollectionIterator&) = default;
  TrackCollectionIterator(TrackCollectionIterator&&) = default;
  TrackCollectionIterator& operator=(const TrackCollectionIterator&) = default;
  TrackCollectionIterator& operator=(TrackCollectionIterator&&) = default;
  ~TrackCollectionIterator() = default;

  constexpr auto operator<=>(const TrackCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const TrackCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  TrackCollectionIterator& operator++();
  TrackCollectionIterator operator++(int);
  TrackCollectionIterator& operator--();
  TrackCollectionIterator operator--(int);
  TrackCollectionIterator& operator+=(difference_type n);
  TrackCollectionIterator operator+(difference_type n) const;
  friend TrackCollectionIterator operator+(difference_type n, const TrackCollectionIterator& it);
  TrackCollectionIterator& operator-=(difference_type n);
  TrackCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const TrackCollectionIterator& other) const;

private:
  size_t m_index{0};
  Track m_object { podio::utils::MaybeSharedPtr<TrackObj>{nullptr} };
  const TrackObjPointerContainer* m_collection{nullptr};
};


class TrackMutableCollectionIterator {
public:
  using value_type = Track;
  using difference_type = ptrdiff_t;
  using reference = MutableTrack;
  using pointer = MutableTrack*;
  using iterator_category = std::input_iterator_tag;
  // `std::forward_iterator` is supported except that the pointers obtained with `operator->()`
  // remain valid as long as the iterator is valid, not as long as the range is valid.
  using iterator_concept = std::random_access_iterator_tag;

  TrackMutableCollectionIterator(size_t index, const TrackObjPointerContainer* collection) : m_index(index), m_object(podio::utils::MaybeSharedPtr<TrackObj>{nullptr}), m_collection(collection) {}
  TrackMutableCollectionIterator() = default;

  TrackMutableCollectionIterator(const TrackMutableCollectionIterator&) = default;
  TrackMutableCollectionIterator(TrackMutableCollectionIterator&&) = default;
  TrackMutableCollectionIterator& operator=(const TrackMutableCollectionIterator&) = default;
  TrackMutableCollectionIterator& operator=(TrackMutableCollectionIterator&&) = default;
  ~TrackMutableCollectionIterator() = default;

  constexpr auto operator<=>(const TrackMutableCollectionIterator& other) const {
    return m_index <=> other.m_index;
  }

  constexpr bool operator==(const TrackMutableCollectionIterator& x) const {
    return m_index ==  x.m_index;
  }

  reference operator*() const;
  pointer operator->();
  TrackMutableCollectionIterator& operator++();
  TrackMutableCollectionIterator operator++(int);
  TrackMutableCollectionIterator& operator--();
  TrackMutableCollectionIterator operator--(int);
  TrackMutableCollectionIterator& operator+=(difference_type n);
  TrackMutableCollectionIterator operator+(difference_type n) const;
  friend TrackMutableCollectionIterator operator+(difference_type n, const TrackMutableCollectionIterator& it);
  TrackMutableCollectionIterator& operator-=(difference_type n);
  TrackMutableCollectionIterator operator-(difference_type n) const;
  reference operator[](difference_type n) const;
  difference_type operator-(const TrackMutableCollectionIterator& other) const;

private:
  size_t m_index{0};
  MutableTrack m_object { podio::utils::MaybeSharedPtr<TrackObj>{nullptr} };
  const TrackObjPointerContainer* m_collection{nullptr};
};


/**
A Collection is identified by an ID.
*/
class TrackCollection : public podio::CollectionBase {
public:
  using value_type = Track;
  using mutable_type = MutableTrack;
  using const_iterator = TrackCollectionIterator;
  using iterator = TrackMutableCollectionIterator;
  using difference_type = ptrdiff_t;
  using size_type = size_t;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = std::reverse_iterator<iterator>;

  TrackCollection() = default;
  TrackCollection(TrackCollectionData&& data, bool isSubsetColl);
  // This is a move-only type
  TrackCollection(const TrackCollection& ) = delete;
  TrackCollection& operator=(const TrackCollection& ) = delete;
  TrackCollection(TrackCollection&&) = default;
  TrackCollection& operator=(TrackCollection&&) = default;

  ~TrackCollection() override;
#if defined(__cpp_lib_containers_ranges)
  template<podio::detail::RangeConvertibleTo<value_type> R>
  TrackCollection(std::from_range_t, R&& range);
#endif

  /// Construct a Track collection from a range of mutable or
  /// immutable handles
  template<podio::detail::RangeConvertibleTo<value_type> R>
  static TrackCollection from(R&& range);

  constexpr static std::string_view typeName = "extension::TrackCollection";
  constexpr static std::string_view valueTypeName = "extension::Track";
  constexpr static std::string_view dataTypeName = "extension::TrackData";

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
  MutableTrack create();

  /// Append a new object to the collection, and return this object.
  /// Initialized with the parameters given.
  /// Does not invalidate iterators.
  template<typename... Args>
  MutableTrack create(Args&&... args);

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
  Track operator[](std::size_t index) const;
  /// Returns the object of a given index
  MutableTrack operator[](std::size_t index);
  /// Returns the const object of given index
  Track at(std::size_t index) const;
  /// Returns the object of given index
  MutableTrack at(std::size_t index);


  /// Append object to the collection.
  /// Does not invalidate iterators.
  void push_back(const MutableTrack& object);
  /// Append an object to the (subset) collection.
  /// Does not invalidate iterators.
  void push_back(const Track& object);

  void prepareForWrite() const final;
  void prepareAfterRead() final;
  bool setReferences(const podio::ICollectionProvider* collectionProvider) final;

  /// Get the collection buffers for this collection
  podio::CollectionWriteBuffers getBuffers() final;

  void setID(uint32_t ID) final {
    m_collectionID = ID;
    if (!m_isSubsetColl) {
      std::for_each(m_storage.entries.begin(), m_storage.entries.end(),
                  [ID] (TrackObj* obj) { obj->id = {obj->id.index, static_cast<uint32_t>(ID)}; }
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


  std::vector<std::int32_t> type(const size_t nElem = 0) const;
  std::vector<float> chi2(const size_t nElem = 0) const;
  std::vector<std::int32_t> ndf(const size_t nElem = 0) const;
  std::vector<float> dEdx(const size_t nElem = 0) const;
  std::vector<float> dEdxError(const size_t nElem = 0) const;
  std::vector<float> radiusOfInnermostHit(const size_t nElem = 0) const;

private:
  // For setReferences, we need to give our own CollectionData access to our
  // private entries. Otherwise we would need to expose a public member function
  // that gives access to the Obj* which is definitely not what we want
  friend class TrackCollectionData;

  mutable bool m_isPrepared{false};
  bool m_isSubsetColl{false};
  uint32_t m_collectionID{static_cast<uint32_t>(podio::ObjectID::untracked)};
  mutable std::unique_ptr<std::mutex> m_storageMtx{std::make_unique<std::mutex>()};
  mutable TrackCollectionData m_storage{};
};

std::ostream& operator<<(std::ostream& o, const TrackCollection& v);

template<typename... Args>
MutableTrack TrackCollection::create(Args&&... args) {
  if (m_isSubsetColl) {
    throw std::logic_error("Cannot create new elements on a subset collection");
  }
  auto obj = new TrackObj({static_cast<int>(m_storage.entries.size()), m_collectionID}, {std::forward<Args>(args)...});
  m_storage.entries.push_back(obj);

  // Need to initialize the relation vectors manually for the {ObjectID, TrackData} constructor
  obj->m_trackerHits = new std::vector<extension::TrackerHit>();
  obj->m_tracks = new std::vector<extension::Track>();
  obj->m_subdetectorHitNumbers = new std::vector<std::int32_t>();
  obj->m_trackStates = new std::vector<edm4hep::TrackState>();
  obj->m_dxQuantities = new std::vector<edm4hep::Quantity>();
  m_storage.createRelations(obj);
  return MutableTrack(podio::utils::MaybeSharedPtr(obj));
}


#if defined(__cpp_lib_containers_ranges)
template<podio::detail::RangeConvertibleTo<TrackCollection::value_type> R>
TrackCollection::TrackCollection(std::from_range_t, R&& range)
  : TrackCollection() {
  setSubsetCollection(podio::detail::RangeOf<R, value_type>);

  for (auto&& elem : range) {
    push_back(std::forward<decltype(elem)>(elem));
  }
}
#endif

template<podio::detail::RangeConvertibleTo<TrackCollection::value_type> R>
TrackCollection TrackCollection::from(R&& range) {
  TrackCollection coll;
  if constexpr (podio::detail::RangeOf<R, value_type>) {
    coll.setSubsetCollection();
  }

  for (auto&& elem : range) {
    coll.push_back(std::forward<decltype(elem)>(elem));
  }

  return coll;
}


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const TrackCollection& collection);
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
constexpr std::string_view extension::TrackCollection::typeName;
constexpr std::string_view extension::TrackCollection::valueTypeName;
constexpr std::string_view extension::TrackCollection::dataTypeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::TrackCollection::typeName;
constexpr std::string_view extension::TrackCollection::valueTypeName;
constexpr std::string_view extension::TrackCollection::dataTypeName;
  #pragma GCC diagnostic pop
#endif


#endif
