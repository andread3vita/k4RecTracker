// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableMCRecoDriftChamberDigiAssociation_H
#define EXTENSION_MutableMCRecoDriftChamberDigiAssociation_H

#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/MCRecoDriftChamberDigiAssociation.h"


#include "podio/utilities/MaybeSharedPtr.h"

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


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableMCRecoDriftChamberDigiAssociation
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MutableMCRecoDriftChamberDigiAssociation {

  friend class MCRecoDriftChamberDigiAssociationCollection;
  friend class MCRecoDriftChamberDigiAssociationMutableCollectionIterator;
  friend class MCRecoDriftChamberDigiAssociation;

public:
  using object_type = MCRecoDriftChamberDigiAssociation;
  using collection_type = MCRecoDriftChamberDigiAssociationCollection;

  /// default constructor
  MutableMCRecoDriftChamberDigiAssociation() = default;

  /// Constructor initializing all members
  MutableMCRecoDriftChamberDigiAssociation(const float weight);

  /// copy constructor
  MutableMCRecoDriftChamberDigiAssociation(const MutableMCRecoDriftChamberDigiAssociation& other) noexcept = default;

  /// move constructor
  MutableMCRecoDriftChamberDigiAssociation(MutableMCRecoDriftChamberDigiAssociation&& other) noexcept = default;

  /// copy-assignment operator
  MutableMCRecoDriftChamberDigiAssociation& operator=(MutableMCRecoDriftChamberDigiAssociation other) &; // Rebind this to other's internal object
  MutableMCRecoDriftChamberDigiAssociation& operator=(MutableMCRecoDriftChamberDigiAssociation other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiAssociation clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableMCRecoDriftChamberDigiAssociation() = default;


public:

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigi getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;

  /// Set the weight of this association
  void setWeight(const float weight);
  /// Get mutable reference to weight of this association
  float& getWeight();


  /// Set the reference to the digitized hit
  void setDigi(const extension::DriftChamberDigi& value);
  /// Set the reference to the simulated hit
  void setSim(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiAssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>{nullptr}; }

  constexpr bool operator==(const MutableMCRecoDriftChamberDigiAssociation& other) const { return m_obj == other.m_obj; }
  bool operator==(const MCRecoDriftChamberDigiAssociation& other) const;

  constexpr bool operator!=(const MutableMCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }
  bool operator!=(const MCRecoDriftChamberDigiAssociation& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableMCRecoDriftChamberDigiAssociation& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableMCRecoDriftChamberDigiAssociation>;

  friend void swap(MutableMCRecoDriftChamberDigiAssociation& a, MutableMCRecoDriftChamberDigiAssociation& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiAssociationObj
  explicit MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> m_obj{new MCRecoDriftChamberDigiAssociationObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiAssociation& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableMCRecoDriftChamberDigiAssociation> {
  std::size_t operator()(const extension::MutableMCRecoDriftChamberDigiAssociation& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiAssociationObj*>{}(obj.m_obj.get());
  }
};


#endif
