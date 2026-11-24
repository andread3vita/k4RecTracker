// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableTrack_H
#define EXTENSION_MutableTrack_H

#include "extension/TrackObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/Track.h"

#include "extension/TrackerHit.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class TrackCollection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableTrack
 *  Reconstructed track
 *  @author: EDM4hep authors
 */
class MutableTrack {

  friend class TrackCollection;
  friend class TrackMutableCollectionIterator;
  friend class Track;

public:
  using object_type = Track;
  using collection_type = TrackCollection;

  /// default constructor
  MutableTrack() = default;

  /// Constructor initializing all members
  MutableTrack(const std::int32_t type, const float chi2, const std::int32_t ndf, const float dEdx, const float dEdxError, const float radiusOfInnermostHit);

  /// copy constructor
  MutableTrack(const MutableTrack& other) noexcept = default;

  /// move constructor
  MutableTrack(MutableTrack&& other) noexcept = default;

  /// copy-assignment operator
  MutableTrack& operator=(MutableTrack other) &; // Rebind this to other's internal object
  MutableTrack& operator=(MutableTrack other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableTrack clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableTrack() = default;


public:

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



  /// Set the flagword that defines the type of track.Bits 16-31 are used internally
  void setType(const std::int32_t type);
  /// Get mutable reference to flagword that defines the type of track.Bits 16-31 are used internally
  std::int32_t& getType();

  /// Set the Chi^2 of the track fit
  void setChi2(const float chi2);
  /// Get mutable reference to Chi^2 of the track fit
  float& getChi2();

  /// Set the number of degrees of freedom of the track fit
  void setNdf(const std::int32_t ndf);
  /// Get mutable reference to number of degrees of freedom of the track fit
  std::int32_t& getNdf();

  /// Set the dEdx of the track
  void setDEdx(const float dEdx);
  /// Get mutable reference to dEdx of the track
  float& getDEdx();

  /// Set the error of dEdx
  void setDEdxError(const float dEdxError);
  /// Get mutable reference to error of dEdx
  float& getDEdxError();

  /// Set the radius of the innermost hit that has been used in the track fit
  void setRadiusOfInnermostHit(const float radiusOfInnermostHit);
  /// Get mutable reference to radius of the innermost hit that has been used in the track fit
  float& getRadiusOfInnermostHit();



  void addToTrackerHits(const extension::TrackerHit&);
  std::size_t trackerHits_size() const;
  extension::TrackerHit getTrackerHits(std::size_t) const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_begin() const;
  std::vector<extension::TrackerHit>::const_iterator trackerHits_end() const;
  podio::RelationRange<extension::TrackerHit> getTrackerHits() const;
  void addToTracks(const extension::Track&);
  std::size_t tracks_size() const;
  extension::Track getTracks(std::size_t) const;
  std::vector<extension::Track>::const_iterator tracks_begin() const;
  std::vector<extension::Track>::const_iterator tracks_end() const;
  podio::RelationRange<extension::Track> getTracks() const;
  void addToSubdetectorHitNumbers(const std::int32_t&);
  std::size_t subdetectorHitNumbers_size() const;
  std::int32_t getSubdetectorHitNumbers(std::size_t) const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_begin() const;
  std::vector<std::int32_t>::const_iterator subdetectorHitNumbers_end() const;
  podio::RelationRange<std::int32_t> getSubdetectorHitNumbers() const;
  void addToTrackStates(const edm4hep::TrackState&);
  std::size_t trackStates_size() const;
  edm4hep::TrackState getTrackStates(std::size_t) const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_begin() const;
  std::vector<edm4hep::TrackState>::const_iterator trackStates_end() const;
  podio::RelationRange<edm4hep::TrackState> getTrackStates() const;
  void addToDxQuantities(const edm4hep::Quantity&);
  std::size_t dxQuantities_size() const;
  edm4hep::Quantity getDxQuantities(std::size_t) const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_begin() const;
  std::vector<edm4hep::Quantity>::const_iterator dxQuantities_end() const;
  podio::RelationRange<edm4hep::Quantity> getDxQuantities() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from TrackObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<TrackObj>{nullptr}; }

  constexpr bool operator==(const MutableTrack& other) const { return m_obj == other.m_obj; }
  bool operator==(const Track& other) const;

  constexpr bool operator!=(const MutableTrack& other) const { return !(*this == other); }
  bool operator!=(const Track& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableTrack& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableTrack>;

  friend void swap(MutableTrack& a, MutableTrack& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing TrackObj
  explicit MutableTrack(podio::utils::MaybeSharedPtr<TrackObj> obj);

  podio::utils::MaybeSharedPtr<TrackObj> m_obj{new TrackObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableTrack& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableTrack> {
  std::size_t operator()(const extension::MutableTrack& obj) const {
    return std::hash<extension::TrackObj*>{}(obj.m_obj.get());
  }
};


#endif
