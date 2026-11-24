// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableSenseWireHitSimTrackerHitLink_H
#define EXTENSION_MutableSenseWireHitSimTrackerHitLink_H

#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/SenseWireHitSimTrackerHitLink.h"


#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class SenseWireHitSimTrackerHitLinkCollection;
class SenseWireHit;
class MutableSenseWireHit;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableSenseWireHitSimTrackerHitLink
 *  Link between a SenseWireHit and a SimTrackerHit
 *  @author: EDM4hep authors
 */
class MutableSenseWireHitSimTrackerHitLink {

  friend class SenseWireHitSimTrackerHitLinkCollection;
  friend class SenseWireHitSimTrackerHitLinkMutableCollectionIterator;
  friend class SenseWireHitSimTrackerHitLink;

public:
  using object_type = SenseWireHitSimTrackerHitLink;
  using collection_type = SenseWireHitSimTrackerHitLinkCollection;

  /// default constructor
  MutableSenseWireHitSimTrackerHitLink() = default;

  /// Constructor initializing all members
  MutableSenseWireHitSimTrackerHitLink(const float weight);

  /// copy constructor
  MutableSenseWireHitSimTrackerHitLink(const MutableSenseWireHitSimTrackerHitLink& other) noexcept = default;

  /// move constructor
  MutableSenseWireHitSimTrackerHitLink(MutableSenseWireHitSimTrackerHitLink&& other) noexcept = default;

  /// copy-assignment operator
  MutableSenseWireHitSimTrackerHitLink& operator=(MutableSenseWireHitSimTrackerHitLink other) &; // Rebind this to other's internal object
  MutableSenseWireHitSimTrackerHitLink& operator=(MutableSenseWireHitSimTrackerHitLink other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHitSimTrackerHitLink clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableSenseWireHitSimTrackerHitLink() = default;


public:

  /// Access the weight of this link
  float getWeight() const;


  /// Access the reference to the SenseWireHit
  const extension::SenseWireHit getFrom() const;
  /// Access the reference to the SimTrackerHit
  const edm4hep::SimTrackerHit getTo() const;

  /// Set the weight of this link
  void setWeight(const float weight);
  /// Get mutable reference to weight of this link
  float& getWeight();


  /// Set the reference to the SenseWireHit
  void setFrom(const extension::SenseWireHit& value);
  /// Set the reference to the SimTrackerHit
  void setTo(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitSimTrackerHitLinkObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}; }

  constexpr bool operator==(const MutableSenseWireHitSimTrackerHitLink& other) const { return m_obj == other.m_obj; }
  bool operator==(const SenseWireHitSimTrackerHitLink& other) const;

  constexpr bool operator!=(const MutableSenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }
  bool operator!=(const SenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableSenseWireHitSimTrackerHitLink& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableSenseWireHitSimTrackerHitLink>;

  friend void swap(MutableSenseWireHitSimTrackerHitLink& a, MutableSenseWireHitSimTrackerHitLink& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitSimTrackerHitLinkObj
  explicit MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj);

  podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> m_obj{new SenseWireHitSimTrackerHitLinkObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHitSimTrackerHitLink& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableSenseWireHitSimTrackerHitLink> {
  std::size_t operator()(const extension::MutableSenseWireHitSimTrackerHitLink& obj) const {
    return std::hash<extension::SenseWireHitSimTrackerHitLinkObj*>{}(obj.m_obj.get());
  }
};


#endif
