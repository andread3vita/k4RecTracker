// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableSenseWireHit_H
#define EXTENSION_MutableSenseWireHit_H

#include "extension/SenseWireHitObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/SenseWireHit.h"

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
class SenseWireHitCollection;
class TrackerHit;
class TrackerHit_dev;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableSenseWireHit
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class MutableSenseWireHit {

  friend class SenseWireHitCollection;
  friend class SenseWireHitMutableCollectionIterator;
  friend class SenseWireHit;

public:
  using object_type = SenseWireHit;
  using collection_type = SenseWireHitCollection;

  /// default constructor
  MutableSenseWireHit() = default;

  /// Constructor initializing all members
  MutableSenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError);

  /// copy constructor
  MutableSenseWireHit(const MutableSenseWireHit& other) noexcept = default;

  /// move constructor
  MutableSenseWireHit(MutableSenseWireHit&& other) noexcept = default;

  /// copy-assignment operator
  MutableSenseWireHit& operator=(MutableSenseWireHit other) &; // Rebind this to other's internal object
  MutableSenseWireHit& operator=(MutableSenseWireHit other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHit clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableSenseWireHit() = default;


public:

  /// Access the ID of the sensor that created this hit
  std::uint64_t getCellID() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the time of the hit [ns]
  float getTime() const;

  /// Access the energy deposited by the hit [GeV]
  float getEDep() const;

  /// Access the error on eDep [GeV]
  float getEDepError() const;

  /// Access the angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float getWireStereoAngle() const;

  /// Access the azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  float getWireAzimuthalAngle() const;

  /// Access the point on the sense wire which is closest to the hit (center of the circle) [mm]
  const edm4hep::Vector3d& getPosition() const;

  /// Access the error on the hit position along the wire direction [mm]
  double getPositionAlongWireError() const;

  /// Access the distance between the hit and the wire (radius of the circle) [mm]
  float getDistanceToWire() const;

  /// Access the error on distanceToWire [mm]
  float getDistanceToWireError() const;



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

  /// Set the energy deposited by the hit [GeV]
  void setEDep(const float eDep);
  /// Get mutable reference to energy deposited by the hit [GeV]
  float& getEDep();

  /// Set the error on eDep [GeV]
  void setEDepError(const float eDepError);
  /// Get mutable reference to error on eDep [GeV]
  float& getEDepError();

  /// Set the angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  void setWireStereoAngle(const float wireStereoAngle);
  /// Get mutable reference to angle between the sense wire axis and the drift chamber axis (usually the z-axis) - use it together with wireAzimuthalAngle to get the wire direction
  float& getWireStereoAngle();

  /// Set the azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  void setWireAzimuthalAngle(const float wireAzimuthalAngle);
  /// Get mutable reference to azimuthal angle at the middle of the sense wire - use it together with wireStereoAngle to get the wire direction
  float& getWireAzimuthalAngle();

  /// Set the point on the sense wire which is closest to the hit (center of the circle) [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to point on the sense wire which is closest to the hit (center of the circle) [mm]
  edm4hep::Vector3d& getPosition();

  /// Set the error on the hit position along the wire direction [mm]
  void setPositionAlongWireError(const double positionAlongWireError);
  /// Get mutable reference to error on the hit position along the wire direction [mm]
  double& getPositionAlongWireError();

  /// Set the distance between the hit and the wire (radius of the circle) [mm]
  void setDistanceToWire(const float distanceToWire);
  /// Get mutable reference to distance between the hit and the wire (radius of the circle) [mm]
  float& getDistanceToWire();

  /// Set the error on distanceToWire [mm]
  void setDistanceToWireError(const float distanceToWireError);
  /// Get mutable reference to error on distanceToWire [mm]
  float& getDistanceToWireError();



  void addToNElectrons(const std::uint16_t&);
  std::size_t nElectrons_size() const;
  std::uint16_t getNElectrons(std::size_t) const;
  std::vector<std::uint16_t>::const_iterator nElectrons_begin() const;
  std::vector<std::uint16_t>::const_iterator nElectrons_end() const;
  podio::RelationRange<std::uint16_t> getNElectrons() const;

 /// Return the number of clusters associated to the hit
 auto getNClusters() const { return getNElectrons().size(); }
 


  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from SenseWireHitObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<SenseWireHitObj>{nullptr}; }

  constexpr bool operator==(const MutableSenseWireHit& other) const { return m_obj == other.m_obj; }
  bool operator==(const SenseWireHit& other) const;

  constexpr bool operator!=(const MutableSenseWireHit& other) const { return !(*this == other); }
  bool operator!=(const SenseWireHit& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableSenseWireHit& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableSenseWireHit>;

  friend void swap(MutableSenseWireHit& a, MutableSenseWireHit& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitObj
  explicit MutableSenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj);

  podio::utils::MaybeSharedPtr<SenseWireHitObj> m_obj{new SenseWireHitObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHit& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableSenseWireHit> {
  std::size_t operator()(const extension::MutableSenseWireHit& obj) const {
    return std::hash<extension::SenseWireHitObj*>{}(obj.m_obj.get());
  }
};


#endif
