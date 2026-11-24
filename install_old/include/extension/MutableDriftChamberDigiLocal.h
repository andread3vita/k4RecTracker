// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigiLocal_H
#define EXTENSION_MutableDriftChamberDigiLocal_H

#include "extension/DriftChamberDigiLocalObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigiLocal.h"

#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiLocalCollection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigiLocal
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class MutableDriftChamberDigiLocal {

  friend class DriftChamberDigiLocalCollection;
  friend class DriftChamberDigiLocalMutableCollectionIterator;
  friend class DriftChamberDigiLocal;

public:
  using object_type = DriftChamberDigiLocal;
  using collection_type = DriftChamberDigiLocalCollection;

  /// default constructor
  MutableDriftChamberDigiLocal() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount);

  /// copy constructor
  MutableDriftChamberDigiLocal(const MutableDriftChamberDigiLocal& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigiLocal(MutableDriftChamberDigiLocal&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigiLocal& operator=(MutableDriftChamberDigiLocal other) &; // Rebind this to other's internal object
  MutableDriftChamberDigiLocal& operator=(MutableDriftChamberDigiLocal other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiLocal clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigiLocal() = default;


public:

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the smeared distance of closest approach between the wire and the hit [mm]
  float getDistanceToWire() const;

  /// Access the smeared z position in the local wire coordinate system [mm]
  float getZPositionAlongWire() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;



  /// Set the ID of the wire that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the wire that created this hit
  std::uint64_t& getCellID();

  /// Set the smeared distance of closest approach between the wire and the hit [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to smeared distance of closest approach between the wire and the hit [mm]
  float& getDistanceToWire();

  /// Set the smeared z position in the local wire coordinate system [mm]
  void setZPositionAlongWire(const float zPositionAlongWire);
  /// Get mutable reference to smeared z position in the local wire coordinate system [mm]
  float& getZPositionAlongWire();

  /// Set the time of the hit [ns].
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns].
  float& getTime();

  /// Set the energy deposited on the hit [GeV].
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited on the hit [GeV].
  float& getEDep();

  /// Set the error measured on eDep [GeV].
  void setEDepError(const float eDepError);
  /// Get mutable reference to error measured on eDep [GeV].
  float& getEDepError();

  /// Set the number of clusters associated to this hit
  void setClusterCount(const std::uint32_t clusterCount);
  /// Get mutable reference to number of clusters associated to this hit
  std::uint32_t& getClusterCount();






  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiLocalObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigiLocal& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigiLocal& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigiLocal& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigiLocal& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigiLocal& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigiLocal>;

  friend void swap(MutableDriftChamberDigiLocal& a, MutableDriftChamberDigiLocal& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiLocalObj
  explicit MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> m_obj{new DriftChamberDigiLocalObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiLocal& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigiLocal> {
  std::size_t operator()(const extension::MutableDriftChamberDigiLocal& obj) const {
    return std::hash<extension::DriftChamberDigiLocalObj*>{}(obj.m_obj.get());
  }
};


#endif
