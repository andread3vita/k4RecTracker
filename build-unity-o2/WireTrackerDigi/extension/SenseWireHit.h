// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_SenseWireHit_H
#define EXTENSION_SenseWireHit_H

#include "extension/SenseWireHitObj.h"

#include "edm4hep/Vector3d.h"
#include "podio/RelationRange.h"
#include <cstdint>
#include <vector>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
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


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::SenseWireHit& obj);
};

namespace extension {

class MutableSenseWireHit;
class SenseWireHitCollection;
class SenseWireHitCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class SenseWireHit
 *  Sense wire hit, before the hit position degeneracy is resolved. The circle representing possible positions is parametrized with its center, radius and normal vector (given by the wire direction).
 *  @author: EDM4hep authors
 */
class SenseWireHit {

  friend class MutableSenseWireHit;
  friend class SenseWireHitCollection;
  friend class extension::SenseWireHitCollectionData;
  friend class SenseWireHitCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const SenseWireHit & obj);

public:
  using mutable_type = MutableSenseWireHit;
  using collection_type = SenseWireHitCollection;

  /// default constructor
  SenseWireHit() = default;

  /// Constructor initializing all members
  SenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError);

  /// copy constructor
  SenseWireHit(const SenseWireHit& other) noexcept = default;

  /// move constructor
  SenseWireHit(SenseWireHit&& other) noexcept = default;

  /// copy-assignment operator
  SenseWireHit& operator=(SenseWireHit other) &; // Rebind this to other's internal object
  SenseWireHit& operator=(SenseWireHit other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableSenseWireHit clone(bool cloneRelations=true) const;

  /// destructor
  ~SenseWireHit() = default;

  /// converting constructor from mutable object
  SenseWireHit(const MutableSenseWireHit& other);

  static SenseWireHit makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::SenseWireHit";

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

  constexpr bool operator==(const SenseWireHit& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableSenseWireHit& other) const;

  constexpr bool operator!=(const SenseWireHit& other) const { return !(*this == other); }
  bool operator!=(const MutableSenseWireHit& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const SenseWireHit& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<SenseWireHit>;

  friend void swap(SenseWireHit& a, SenseWireHit& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing SenseWireHitObj
  explicit SenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj);
  SenseWireHit(SenseWireHitObj* obj);

  podio::utils::MaybeSharedPtr<SenseWireHitObj> m_obj{new SenseWireHitObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const SenseWireHit& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHit& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::SenseWireHit> {
  std::size_t operator()(const extension::SenseWireHit& obj) const {
    return std::hash<extension::SenseWireHitObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::SenseWireHit::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::SenseWireHit::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
