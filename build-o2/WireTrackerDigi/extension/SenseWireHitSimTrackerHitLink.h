// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHitSimTrackerHitLink_H
#define EXTENSION_SenseWireHitSimTrackerHitLink_H

#include "extension/SenseWireHitSimTrackerHitLinkObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
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


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::SenseWireHitSimTrackerHitLink& obj);
};

namespace extension {

class MutableSenseWireHitSimTrackerHitLink;
class SenseWireHitSimTrackerHitLinkCollection;
class SenseWireHitSimTrackerHitLinkCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class SenseWireHitSimTrackerHitLink
 *  Link between a SenseWireHit and a SimTrackerHit
 *  @author: EDM4hep authors
 */
class SenseWireHitSimTrackerHitLink {

  friend class MutableSenseWireHitSimTrackerHitLink;
  friend class SenseWireHitSimTrackerHitLinkCollection;
  friend class extension::SenseWireHitSimTrackerHitLinkCollectionData;
  friend class SenseWireHitSimTrackerHitLinkCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const SenseWireHitSimTrackerHitLink & obj);

public:
  using mutable_type = MutableSenseWireHitSimTrackerHitLink;
  using collection_type = SenseWireHitSimTrackerHitLinkCollection;

  /// default constructor
  SenseWireHitSimTrackerHitLink() = default;

  /// Constructor initializing all members
  SenseWireHitSimTrackerHitLink(const float weight);

  /// copy constructor
  SenseWireHitSimTrackerHitLink(const SenseWireHitSimTrackerHitLink& other) noexcept = default;

  /// move constructor
  SenseWireHitSimTrackerHitLink(SenseWireHitSimTrackerHitLink&& other) noexcept = default;

  /// copy-assignment operator
  SenseWireHitSimTrackerHitLink& operator=(SenseWireHitSimTrackerHitLink other) &; // Rebind this to other's internal object
  SenseWireHitSimTrackerHitLink& operator=(SenseWireHitSimTrackerHitLink other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHitSimTrackerHitLink clone(bool cloneRelations=true) const;

  /// destructor
  ~SenseWireHitSimTrackerHitLink() = default;

  /// converting constructor from mutable object
  SenseWireHitSimTrackerHitLink(const MutableSenseWireHitSimTrackerHitLink& other);

  static SenseWireHitSimTrackerHitLink makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::SenseWireHitSimTrackerHitLink";

  /// Access the weight of this link
  float getWeight() const;


  /// Access the reference to the SenseWireHit
  const extension::SenseWireHit getFrom() const;
  /// Access the reference to the SimTrackerHit
  const edm4hep::SimTrackerHit getTo() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitSimTrackerHitLinkObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>{nullptr}; }

  constexpr bool operator==(const SenseWireHitSimTrackerHitLink& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableSenseWireHitSimTrackerHitLink& other) const;

  constexpr bool operator!=(const SenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }
  bool operator!=(const MutableSenseWireHitSimTrackerHitLink& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const SenseWireHitSimTrackerHitLink& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<SenseWireHitSimTrackerHitLink>;

  friend void swap(SenseWireHitSimTrackerHitLink& a, SenseWireHitSimTrackerHitLink& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitSimTrackerHitLinkObj
  explicit SenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj);
  SenseWireHitSimTrackerHitLink(SenseWireHitSimTrackerHitLinkObj* obj);

  podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> m_obj{new SenseWireHitSimTrackerHitLinkObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHitSimTrackerHitLink& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitSimTrackerHitLink& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::SenseWireHitSimTrackerHitLink> {
  std::size_t operator()(const extension::SenseWireHitSimTrackerHitLink& obj) const {
    return std::hash<extension::SenseWireHitSimTrackerHitLinkObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::SenseWireHitSimTrackerHitLink::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHitSimTrackerHitLink::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
