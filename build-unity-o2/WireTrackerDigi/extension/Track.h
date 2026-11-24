// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_Track_H
#define EXTENSION_Track_H

#include "extension/TrackObj.h"

#include "extension/TrackerHit.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class TrackCollection;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::Track& obj);
};

namespace extension {

class MutableTrack;
class TrackCollection;
class TrackCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class Track
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class Track {

  friend class MutableTrack;
  friend class TrackCollection;
  friend class extension::TrackCollectionData;
  friend class TrackCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const Track & obj);

public:
  using mutable_type = MutableTrack;
  using collection_type = TrackCollection;

  /// default constructor
  Track() = default;

  /// Constructor initializing all members
  Track(const std::int32_t type, const float chi2, const std::int32_t ndf, const float dEdx, const float dEdxError, const float radiusOfInnermostHit);

  /// copy constructor
  Track(const Track& other) noexcept = default;

  /// move constructor
  Track(Track&& other) noexcept = default;

  /// copy-assignment operator
  Track& operator=(Track other) &; // Rebind this to other's internal object
  Track& operator=(Track other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableTrack clone(bool cloneRelations=true) const;

  /// destructor
  ~Track() = default;

  /// converting constructor from mutable object
  Track(const MutableTrack& other);

  static Track makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::Track";

  /// Access the flagword that defines the type of track.Bits 16-31 are used internally
  std::int32_t getType() const;

  /// Access the Chi^2 of the track fit
  float getChi2() const;

  /// Access the number of degrees of freedom of the track fit
  std::int32_t getNdf() const;

  /// Access the dEdx of the track
  float getDEdx() const;

  /// Access the error of dEdx
  float getDEdxError() const;

  /// Access the radius of the innermost hit that has been used in the track fit
  float getRadiusOfInnermostHit() const;



  std::size_t trackerHits_size() const;
  extension::TrackerHit getTrackerHits(std::size_t) const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_begin() const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_end() const;
  podio::RelationRange<extension::TrackerHit> getTrackerHits() const;
  std::size_t tracks_size() const;
  extension::Track getTracks(std::size_t) const;
  std::vector<extension::Track>::const_iterator tracks_begin() const;
  std::vector<extension::Track>::const_iterator tracks_end() const;
  podio::RelationRange<extension::Track> getTracks() const;
  std::size_t subdetectorHitNumbers_size() const;
  std::int32_t getSubdetectorHitNumbers(std::size_t) const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_begin() const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_end() const;
  podio::RelationRange<std::int32_t> getSubdetectorHitNumbers() const;
  std::size_t trackStates_size() const;
  edm4hep::TrackState getTrackStates(std::size_t) const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_begin() const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_end() const;
  podio::RelationRange<edm4hep::TrackState> getTrackStates() const;
  std::size_t dxQuantities_size() const;
  edm4hep::Quantity getDxQuantities(std::size_t) const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_begin() const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_end() const;
  podio::RelationRange<edm4hep::Quantity> getDxQuantities() const;


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from TrackObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<TrackObj>{nullptr}; }

  constexpr bool operator==(const Track& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableTrack& other) const;

  constexpr bool operator!=(const Track& other) const { return !(*this == other); }
  bool operator!=(const MutableTrack& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const Track& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<Track>;

  friend void swap(Track& a, Track& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing TrackObj
  explicit Track(podio::utils::MaybeSharedPtr<TrackObj> obj);
  Track(TrackObj* obj);

  podio::utils::MaybeSharedPtr<TrackObj> m_obj{new TrackObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const Track& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const Track& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::Track> {
  std::size_t operator()(const extension::Track& obj) const {
    return std::hash<extension::TrackObj*>{}(obj.m_obj.get());
  }
};


// This is needed to avoid triggering opening every library in LD_LIBRARY_PATH
// until it's fixed in ROOT. See https://github.com/root-project/root/issues/18489
// and https://github.com/AIDASoft/podio/issues/770
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunknown-warning-option"
  #pragma clang diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
  #pragma clang diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::Track::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::Track::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
