// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#ifndef EXTENSION_DriftChamberDigiLocal_H
#define EXTENSION_DriftChamberDigiLocal_H

#include "extension/DriftChamberDigiLocalObj.h"

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
class DriftChamberDigiLocalCollection;
}


namespace podio::detail {
// Internal function used in less comparison operators of the datatypes and interface types
OrderKey getOrderKey(const extension::DriftChamberDigiLocal& obj);
};

namespace extension {

class MutableDriftChamberDigiLocal;
class DriftChamberDigiLocalCollection;
class DriftChamberDigiLocalCollectionData;

/// \addtogroup extension_user_facing
/// @{
/** @class DriftChamberDigiLocal
 *  Drift chamber digitized hit (before tracking) in local coordinates
 *  @author: B. Francois, CERN
 */
class DriftChamberDigiLocal {

  friend class MutableDriftChamberDigiLocal;
  friend class DriftChamberDigiLocalCollection;
  friend class extension::DriftChamberDigiLocalCollectionData;
  friend class DriftChamberDigiLocalCollectionIterator;
  friend podio::detail::OrderKey podio::detail::getOrderKey(const DriftChamberDigiLocal & obj);

public:
  using mutable_type = MutableDriftChamberDigiLocal;
  using collection_type = DriftChamberDigiLocalCollection;

  /// default constructor
  DriftChamberDigiLocal() = default;

  /// Constructor initializing all members
  DriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount);

  /// copy constructor
  DriftChamberDigiLocal(const DriftChamberDigiLocal& other) noexcept = default;

  /// move constructor
  DriftChamberDigiLocal(DriftChamberDigiLocal&& other) noexcept = default;

  /// copy-assignment operator
  DriftChamberDigiLocal& operator=(DriftChamberDigiLocal other) &; // Rebind this to other's internal object
  DriftChamberDigiLocal& operator=(DriftChamberDigiLocal other) && = delete; // Prevent rebinding temporary as the changes wouldn't persist

  /// create a mutable deep-copy of the object with identical relations
  /// if cloneRelations=false, the relations are not cloned and will be empty
  MutableDriftChamberDigiLocal clone(bool cloneRelations=true) const;

  /// destructor
  ~DriftChamberDigiLocal() = default;

  /// converting constructor from mutable object
  DriftChamberDigiLocal(const MutableDriftChamberDigiLocal& other);

  static DriftChamberDigiLocal makeEmpty();

public:

  static constexpr std::string_view typeName = "extension::DriftChamberDigiLocal";

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





  /// check whether the object is actually available
  bool isAvailable() const;
  /// disconnect from DriftChamberDigiLocalObj instance
  void unlink() { m_obj = podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>{nullptr}; }

  constexpr bool operator==(const DriftChamberDigiLocal& other) const { return m_obj == other.m_obj; }
  bool operator==(const MutableDriftChamberDigiLocal& other) const;

  constexpr bool operator!=(const DriftChamberDigiLocal& other) const { return !(*this == other); }
  bool operator!=(const MutableDriftChamberDigiLocal& other) const { return !(*this == other); }

  // less comparison operator, so that objects can be e.g. stored in sets.
  bool operator<(const DriftChamberDigiLocal& other) const { return podio::detail::getOrderKey(*this) < podio::detail::getOrderKey(other); }

  podio::ObjectID id() const { return getObjectID(); }

  const podio::ObjectID getObjectID() const;

  friend std::hash<DriftChamberDigiLocal>;

  friend void swap(DriftChamberDigiLocal& a, DriftChamberDigiLocal& b) {
    using std::swap;
    swap(a.m_obj, b.m_obj); // swap out the internal pointers
  }

private:
  /// constructor from existing DriftChamberDigiLocalObj
  explicit DriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj);
  DriftChamberDigiLocal(DriftChamberDigiLocalObj* obj);

  podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> m_obj{new DriftChamberDigiLocalObj{}, podio::utils::MarkOwned};
};

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiLocal& value);

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiLocal& value);
#endif


/// @}
} // namespace extension



template<>
struct std::hash<extension::DriftChamberDigiLocal> {
  std::size_t operator()(const extension::DriftChamberDigiLocal& obj) const {
    return std::hash<extension::DriftChamberDigiLocalObj*>{}(obj.m_obj.get());
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
constexpr std::string_view extension::DriftChamberDigiLocal::typeName;
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated"
constexpr std::string_view extension::DriftChamberDigiLocal::typeName;
  #pragma GCC diagnostic pop
#endif


#endif
