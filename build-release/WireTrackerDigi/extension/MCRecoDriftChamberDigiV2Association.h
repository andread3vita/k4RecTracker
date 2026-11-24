// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiV2Association_H
#define EXTENSION_MCRecoDriftChamberDigiV2Association_H

#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiV2AssociationCollection;
class DriftChamberDigiV2;
class MutableDriftChamberDigiV2;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::MCRecoDriftChamberDigiV2Association& obj);
};

namespace extension {

class MutableMCRecoDriftChamberDigiV2Association;
class MCRecoDriftChamberDigiV2AssociationCollection;
class MCRecoDriftChamberDigiV2AssociationCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class MCRecoDriftChamberDigiV2Association
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiV2Association {

  friend class MutableMCRecoDriftChamberDigiV2Association;
  friend class MCRecoDriftChamberDigiV2AssociationCollection;
  friend class extension::MCRecoDriftChamberDigiV2AssociationCollectionData;
  friend class MCRecoDriftChamberDigiV2AssociationCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const MCRecoDriftChamberDigiV2Association & obj);

public:
  using mutable_type = MutableMCRecoDriftChamberDigiV2Association;
  using collection_type = MCRecoDriftChamberDigiV2AssociationCollection;

  /// default constructor
  MCRecoDriftChamberDigiV2Association() = default;

  /// Constructor initializing all members
  MCRecoDriftChamberDigiV2Association(const float weight);

  /// copy constructor
  MCRecoDriftChamberDigiV2Association(const MCRecoDriftChamberDigiV2Association& other) noexcept = default;

  /// move constructor
  MCRecoDriftChamberDigiV2Association(MCRecoDriftChamberDigiV2Association&& other) noexcept = default;

  /// copy-assignment operator
  MCRecoDriftChamberDigiV2Association& operator=(MCRecoDriftChamberDigiV2Association other) &; // Rebind this to other's internal object
  MCRecoDriftChamberDigiV2Association& operator=(MCRecoDriftChamberDigiV2Association other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiV2Association clone(bool cloneRelations=true) const;

  /// destructor
  ~MCRecoDriftChamberDigiV2Association() = default;

  /// converting constructor from mutable object
  MCRecoDriftChamberDigiV2Association(const MutableMCRecoDriftChamberDigiV2Association& other);

  static MCRecoDriftChamberDigiV2Association makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::MCRecoDriftChamberDigiV2Association";

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigiV2 getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiV2AssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}; }

  constexpr bool operator==(const MCRecoDriftChamberDigiV2Association& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableMCRecoDriftChamberDigiV2Association& other) const;

  constexpr bool operator!=(const MCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }
  bool operator!=(const MutableMCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MCRecoDriftChamberDigiV2Association& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MCRecoDriftChamberDigiV2Association>;

  friend void swap(MCRecoDriftChamberDigiV2Association& a, MCRecoDriftChamberDigiV2Association& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiV2AssociationObj
  explicit MCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj);
  MCRecoDriftChamberDigiV2Association(MCRecoDriftChamberDigiV2AssociationObj* obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> m_obj{new MCRecoDriftChamberDigiV2AssociationObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2Association& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2Association& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MCRecoDriftChamberDigiV2Association> {
  std::size_t operator()(const extension::MCRecoDriftChamberDigiV2Association& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiV2AssociationObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::MCRecoDriftChamberDigiV2Association::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiV2Association::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
