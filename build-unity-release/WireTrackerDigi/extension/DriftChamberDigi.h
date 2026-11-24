// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigi_H
#define EXTENSION_DriftChamberDigi_H

#include "extension/DriftChamberDigiObj.h"

#include "edm4hep/Vector3d.h"
#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"
#include "podio/detail/OrderKey.h"

#include <ostream>
#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiCollection;
class TrackerHit;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::DriftChamberDigi& obj);
};

namespace extension {

class MutableDriftChamberDigi;
class DriftChamberDigiCollection;
class DriftChamberDigiCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class DriftChamberDigi
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class DriftChamberDigi {

  friend class MutableDriftChamberDigi;
  friend class DriftChamberDigiCollection;
  friend class extension::DriftChamberDigiCollectionData;
  friend class DriftChamberDigiCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const DriftChamberDigi & obj);

public:
  using mutable_type = MutableDriftChamberDigi;
  using collection_type = DriftChamberDigiCollection;

  /// default constructor
  DriftChamberDigi() = default;

  /// Constructor initializing all members
  DriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position);

  /// copy constructor
  DriftChamberDigi(const DriftChamberDigi& other) noexcept = default;

  /// move constructor
  DriftChamberDigi(DriftChamberDigi&& other) noexcept = default;

  /// copy-assignment operator
  DriftChamberDigi& operator=(DriftChamberDigi other) &; // Rebind this to other's internal object
  DriftChamberDigi& operator=(DriftChamberDigi other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigi clone(bool cloneRelations=true) const;

  /// destructor
  ~DriftChamberDigi() = default;

  /// converting constructor from mutable object
  DriftChamberDigi(const MutableDriftChamberDigi& other);

  static DriftChamberDigi makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::DriftChamberDigi";

  /// Access the ID of the wire that created this hit
  std::uint64_t getCellID() const;

  /// Access the position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getLeftPosition() const;

  /// Access the position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  const edm4hep::Vector3d& getRightPosition() const;

  /// Access the time of the hit [ns].
  float getTime() const;

  /// Access the energy deposited on the hit [GeV].
  float getEDep() const;

  /// Access the error measured on eDep [GeV].
  float getEDepError() const;

  /// Access the number of clusters associated to this hit
  std::uint32_t getClusterCount() const;

  /// Access the type of the raw data hit
  std::int32_t getType() const;

  /// Access the quality bit flag of the hit
  std::int32_t getQuality() const;

  /// Access the hit position [mm]
  const edm4hep::Vector3d& getPosition() const;





  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}; }

  constexpr bool operator==(const DriftChamberDigi& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableDriftChamberDigi& other) const;

  constexpr bool operator!=(const DriftChamberDigi& other) const { return !(*this == other); }
  bool operator!=(const MutableDriftChamberDigi& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const DriftChamberDigi& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<DriftChamberDigi>;

  friend void swap(DriftChamberDigi& a, DriftChamberDigi& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiObj
  explicit DriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj);
  DriftChamberDigi(DriftChamberDigiObj* obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiObj> m_obj{new DriftChamberDigiObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigi& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigi& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::DriftChamberDigi> {
  std::size_t operator()(const extension::DriftChamberDigi& obj) const {
    return std::hash<extension::DriftChamberDigiObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::DriftChamberDigi::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigi::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
