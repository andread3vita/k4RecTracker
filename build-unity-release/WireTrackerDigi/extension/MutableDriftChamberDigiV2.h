// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigiV2_H
#define EXTENSION_MutableDriftChamberDigiV2_H

#include "extension/DriftChamberDigiV2Obj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigiV2.h"

#include "edm4hep/Vector3d.h"
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
class DriftChamberDigiV2Collection;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigiV2
 *  Digitized hit (before tracking) for Drift Chamber v2 (requires data extension).
 *  @author: A. Tolosa-Delgado, B. Francois, CERN
 */
class MutableDriftChamberDigiV2 {

  friend class DriftChamberDigiV2Collection;
  friend class DriftChamberDigiV2MutableCollectionIterator;
  friend class DriftChamberDigiV2;

public:
  using object_type = DriftChamberDigiV2;
  using collection_type = DriftChamberDigiV2Collection;

  /// default constructor
  MutableDriftChamberDigiV2() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigiV2(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const edm4hep::Vector3d& position, const edm4hep::Vector3d& directionSW, const float distanceToWire, const std::uint32_t nCluster);

  /// copy constructor
  MutableDriftChamberDigiV2(const MutableDriftChamberDigiV2& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigiV2(MutableDriftChamberDigiV2&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigiV2& operator=(MutableDriftChamberDigiV2 other) &; // Rebind this to other's internal object
  MutableDriftChamberDigiV2& operator=(MutableDriftChamberDigiV2 other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiV2 clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigiV2() = default;


public:

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited on the hit [GeV]
  float getEDep() const;

  /// Access the error measured on eDep [GeV]
  float getEDepError() const;

  /// Access the point on the sensitive wire (SW) which is closest to the simhit [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the direction of SW
  const edm4hep::Vector3d& getDirectionSW() const;

  /// Access the distance hit-wire [mm]
  float getDistanceToWire() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getNCluster() const;



  /// Set the ID of the sensor that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the sensor that created this hit
  std::uint64_t& getCellID();

  /// Set the type of the raw data hit
  void setType(const std::int32_t type);
  /// Get mutable reference to type of the raw data hit
  std::int32_t& getType();

  /// Set the quality bit flag of the hit
  void setQuality(const std::int32_t quality);
  /// Get mutable reference to quality bit flag of the hit
  std::int32_t& getQuality();

  /// Set the time of the hit [ns]
  void setTime(const float time);
  /// Get mutable reference to time of the hit [ns]
  float& getTime();

  /// Set the energy deposited on the hit [GeV]
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited on the hit [GeV]
  float& getEDep();

  /// Set the error measured on eDep [GeV]
  void setEDepError(const float eDepError);
  /// Get mutable reference to error measured on eDep [GeV]
  float& getEDepError();

  /// Set the point on the sensitive wire (SW) which is closest to the simhit [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to point on the sensitive wire (SW) which is closest to the simhit [mm]
  edm4hep::Vector3d& getPosition();

  /// Set the direction of SW
  void setDirectionSW(const edm4hep::Vector3d& directionSW);
  /// Get mutable reference to direction of SW
  edm4hep::Vector3d& getDirectionSW();

  /// Set the distance hit-wire [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to distance hit-wire [mm]
  float& getDistanceToWire();

  /// Set the number of clusters associated to this hit
  void setNCluster(const std::uint32_t nCluster);
  /// Get mutable reference to number of clusters associated to this hit
  std::uint32_t& getNCluster();



  void addToNElectrons(const std::uint16_t&);
  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;



  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiV2Obj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigiV2& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigiV2& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigiV2& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigiV2& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigiV2& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigiV2>;

  friend void swap(MutableDriftChamberDigiV2& a, MutableDriftChamberDigiV2& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiV2Obj
  explicit MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> m_obj{new DriftChamberDigiV2Obj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiV2& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigiV2> {
  std::size_t operator()(const extension::MutableDriftChamberDigiV2& obj) const {
    return std::hash<extension::DriftChamberDigiV2Obj*>{}(obj.m_obj.get());
  }
};


#endif
