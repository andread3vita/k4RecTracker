// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MCRecoDriftChamberDigiAssociation_H
#define EXTENSION_MCRecoDriftChamberDigiAssociation_H

#include "extension/MCRecoDriftChamberDigiAssociationObj.h"


#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class MCRecoDriftChamberDigiAssociationCollection;
class DriftChamberDigi;
class MutableDriftChamberDigi;
}
namespace edm4hep {
class SimTrackerHit;
class MutableSimTrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::MCRecoDriftChamberDigiAssociation& obj);
};

namespace extension {

class MutableMCRecoDriftChamberDigiAssociation;
class MCRecoDriftChamberDigiAssociationCollection;
class MCRecoDriftChamberDigiAssociationCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class MCRecoDriftChamberDigiAssociation
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MCRecoDriftChamberDigiAssociation {

  friend class MutableMCRecoDriftChamberDigiAssociation;
  friend class MCRecoDriftChamberDigiAssociationCollection;
  friend class extension::MCRecoDriftChamberDigiAssociationCollectionData;
  friend class MCRecoDriftChamberDigiAssociationCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const MCRecoDriftChamberDigiAssociation & obj);

public:
  using mutable_type = MutableMCRecoDriftChamberDigiAssociation;
  using collection_type = MCRecoDriftChamberDigiAssociationCollection;

  /// default constructor
  MCRecoDriftChamberDigiAssociation() = default;

  /// Constructor initializing all members
  MCRecoDriftChamberDigiAssociation(const float weight);

  /// copy constructor
  MCRecoDriftChamberDigiAssociation(const MCRecoDriftChamberDigiAssociation& other) noexcept = default;

  /// move constructor
  MCRecoDriftChamberDigiAssociation(MCRecoDriftChamberDigiAssociation&& other) noexcept = default;

  /// copy-assignment operator
  MCRecoDriftChamberDigiAssociation& operator=(MCRecoDriftChamberDigiAssociation other) &; // Rebind this to other's internal object
  MCRecoDriftChamberDigiAssociation& operator=(MCRecoDriftChamberDigiAssociation other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiAssociation clone(bool cloneRelations=true) const;

  /// destructor
  ~MCRecoDriftChamberDigiAssociation() = default;

  /// converting constructor from mutable object
  MCRecoDriftChamberDigiAssociation(const MutableMCRecoDriftChamberDigiAssociation& other);

  static MCRecoDriftChamberDigiAssociation makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::MCRecoDriftChamberDigiAssociation";

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigi getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiAssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}; }

  constexpr bool operator==(const MCRecoDriftChamberDigiAssociation& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableMCRecoDriftChamberDigiAssociation& other) const;

  constexpr bool operator!=(const MCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }
  bool operator!=(const MutableMCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MCRecoDriftChamberDigiAssociation& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MCRecoDriftChamberDigiAssociation>;

  friend void swap(MCRecoDriftChamberDigiAssociation& a, MCRecoDriftChamberDigiAssociation& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiAssociationObj
  explicit MCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj);
  MCRecoDriftChamberDigiAssociation(MCRecoDriftChamberDigiAssociationObj* obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> m_obj{new MCRecoDriftChamberDigiAssociationObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociation& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociation& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MCRecoDriftChamberDigiAssociation> {
  std::size_t operator()(const extension::MCRecoDriftChamberDigiAssociation& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiAssociationObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociation::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::MCRecoDriftChamberDigiAssociation::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
