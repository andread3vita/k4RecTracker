// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_MutableDriftChamberDigi_H
#define EXTENSION_MutableDriftChamberDigi_H

#include "extension/DriftChamberDigiObj.h"
// Make the immutable class available from its mutable version but not vice versa
#include "extension/DriftChamberDigi.h"

#include "edm4hep/Vector3d.h"
#include <cstdint>

#include "podio/utilities/MaybeSharedPtr.h"

#include <cstdint>

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json_fwd.hpp"
#endif

// forward declarations
namespace extension {
class DriftChamberDigiCollection;
class TrackerHit;
}


namespace extension {

/// \addtogroup extension_user_facing
/// @{

/** @class MutableDriftChamberDigi
 *  Drift chamber digitized hit (before tracking) in global coordinates. Assumes that the hits are radially in the middle of the cells
 *  @author: B. Francois, CERN
 */
class MutableDriftChamberDigi {

  friend class DriftChamberDigiCollection;
  friend class DriftChamberDigiMutableCollectionIterator;
  friend class DriftChamberDigi;

public:
  using object_type = DriftChamberDigi;
  using collection_type = DriftChamberDigiCollection;

  /// default constructor
  MutableDriftChamberDigi() = default;

  /// Constructor initializing all members
  MutableDriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position);

  /// copy constructor
  MutableDriftChamberDigi(const MutableDriftChamberDigi& other) noexcept = default;

  /// move constructor
  MutableDriftChamberDigi(MutableDriftChamberDigi&& other) noexcept = default;

  /// copy-assignment operator
  MutableDriftChamberDigi& operator=(MutableDriftChamberDigi other) &; // Rebind this to other's internal object
  MutableDriftChamberDigi& operator=(MutableDriftChamberDigi other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigi clone(bool cloneRelations=true) const;

  /// destructor
  ~MutableDriftChamberDigi() = default;


public:

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



  /// Set the ID of the wire that created this hit
  void setCellID(const std::uint64_t cellID);
  /// Get mutable reference to ID of the wire that created this hit
  std::uint64_t& getCellID();

  /// Set the position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  void setLeftPosition(const edm4hep::Vector3d& leftPosition);
  /// Get mutable reference to position of the hit assuming it was on the left side of the wire, radially in the middle of the cell [mm]
  edm4hep::Vector3d& getLeftPosition();

  /// Set the position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  void setRightPosition(const edm4hep::Vector3d& rightPosition);
  /// Get mutable reference to position of the hit assuming it was on the right side of the wire, radially in the middle of the cell [mm]
  edm4hep::Vector3d& getRightPosition();

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

  /// Set the type of the raw data hit
  void setType(const std::int32_t type);
  /// Get mutable reference to type of the raw data hit
  std::int32_t& getType();

  /// Set the quality bit flag of the hit
  void setQuality(const std::int32_t quality);
  /// Get mutable reference to quality bit flag of the hit
  std::int32_t& getQuality();

  /// Set the hit position [mm]
  void setPosition(const edm4hep::Vector3d& position);
  /// Get mutable reference to hit position [mm]
  edm4hep::Vector3d& getPosition();






  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiObj>{nullptr}; }

  constexpr bool operator==(const MutableDriftChamberDigi& other) const { return m_obj == other.m_obj; }
  bool operator==(const DriftChamberDigi& other) const;

  constexpr bool operator!=(const MutableDriftChamberDigi& other) const { return !(*this == other); }
  bool operator!=(const DriftChamberDigi& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const MutableDriftChamberDigi& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<MutableDriftChamberDigi>;

  friend void swap(MutableDriftChamberDigi& a, MutableDriftChamberDigi& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiObj
  explicit MutableDriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiObj> m_obj{new DriftChamberDigiObj{}, podio::utils::MarkOwned};
};

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigi& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::MutableDriftChamberDigi> {
  std::size_t operator()(const extension::MutableDriftChamberDigi& obj) const {
    return std::hash<extension::DriftChamberDigiObj*>{}(obj.m_obj.get());
  }
};


#endif
