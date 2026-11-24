// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableMCRecoDriftChamberDigiV2Association_H
#define EXTENSION_MutableMCRecoDriftChamberDigiV2Association_H

#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/MCRecoDriftChamberDigiV2Association.h"


#include "podio/utilities/MaybeSharedPtr.h"

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


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableMCRecoDriftChamberDigiV2Association
 *  Association between a DriftChamberDigi and the corresponding simulated hit
 *  @author: B. Francois, CERN
 */
class MutableMCRecoDriftChamberDigiV2Association {

  friend class MCRecoDriftChamberDigiV2AssociationCollection;
  friend class MCRecoDriftChamberDigiV2AssociationMutableCollectionIterator;
  friend class MCRecoDriftChamberDigiV2Association;

public:
  using object_type = MCRecoDriftChamberDigiV2Association;
  using collection_type = MCRecoDriftChamberDigiV2AssociationCollection;

  /// default constructor
  MutableMCRecoDriftChamberDigiV2Association() = default;

  /// Constructor initializing all members
  MutableMCRecoDriftChamberDigiV2Association(const float weight);

  /// copy constructor
  MutableMCRecoDriftChamberDigiV2Association(const MutableMCRecoDriftChamberDigiV2Association& other) noexcept = default;

  /// move constructor
  MutableMCRecoDriftChamberDigiV2Association(MutableMCRecoDriftChamberDigiV2Association&& other) noexcept = default;

  /// copy-assignment operator
  MutableMCRecoDriftChamberDigiV2Association& operator=(MutableMCRecoDriftChamberDigiV2Association other) &; // Rebind this to other's internal object
  MutableMCRecoDriftChamberDigiV2Association& operator=(MutableMCRecoDriftChamberDigiV2Association other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableMCRecoDriftChamberDigiV2Association clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableMCRecoDriftChamberDigiV2Association() = default;


public:

  /// Access the weight of this association
  float getWeight() const;


  /// Access the reference to the digitized hit
  const extension::DriftChamberDigiV2 getDigi() const;
  /// Access the reference to the simulated hit
  const edm4hep::SimTrackerHit getSim() const;

  /// Set the weight of this association
  void setWeight(const float weight);
  /// Get mutable reference to weight of this association
  float& getWeight();


  /// Set the reference to the digitized hit
  void setDigi(const extension::DriftChamberDigiV2& value);
  /// Set the reference to the simulated hit
  void setSim(const edm4hep::SimTrackerHit& value);




  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from MCRecoDriftChamberDigiV2AssociationObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>{nullptr}; }

  constexpr bool operator==(const MutableMCRecoDriftChamberDigiV2Association& other) const { return m_obj == other.m_obj; }
  bool operator==(const MCRecoDriftChamberDigiV2Association& other) const;

  constexpr bool operator!=(const MutableMCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }
  bool operator!=(const MCRecoDriftChamberDigiV2Association& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableMCRecoDriftChamberDigiV2Association& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableMCRecoDriftChamberDigiV2Association>;

  friend void swap(MutableMCRecoDriftChamberDigiV2Association& a, MutableMCRecoDriftChamberDigiV2Association& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing MCRecoDriftChamberDigiV2AssociationObj
  explicit MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj);

  podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> m_obj{new MCRecoDriftChamberDigiV2AssociationObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiV2Association& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableMCRecoDriftChamberDigiV2Association> {
  std::size_t operator()(const extension::MutableMCRecoDriftChamberDigiV2Association& obj) const {
    return std::hash<extension::MCRecoDriftChamberDigiV2AssociationObj*>{}(obj.m_obj.get());
  }
};


#endif
