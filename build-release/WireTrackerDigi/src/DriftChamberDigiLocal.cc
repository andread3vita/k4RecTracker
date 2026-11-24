// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigiLocal.h"
#include "extension/MutableDriftChamberDigiLocal.h"
#include "extension/DriftChamberDigiLocalObj.h"
#include "extension/DriftChamberDigiLocalData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



DriftChamberDigiLocal::DriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount) {
  m_obj->data.cellID = cellID;
  m_obj->data.distanceToWire = distanceToWire;
  m_obj->data.zPositionAlongWire = zPositionAlongWire;
  m_obj->data.time = time;
  m_obj->data.eDep = eDep;
  m_obj->data.eDepError = eDepError;
  m_obj->data.clusterCount = clusterCount;
}

DriftChamberDigiLocal& DriftChamberDigiLocal::operator=(DriftChamberDigiLocal other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigiLocal DriftChamberDigiLocal::clone(bool cloneRelations) const {
  auto tmp = new DriftChamberDigiLocalObj(podio::ObjectID{}, m_obj->data);
  if (cloneRelations) {
  }
  else {
  }
  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

DriftChamberDigiLocal::DriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj) : m_obj(std::move(obj)) {}

DriftChamberDigiLocal::DriftChamberDigiLocal(const MutableDriftChamberDigiLocal& other): DriftChamberDigiLocal(other.m_obj) {}

DriftChamberDigiLocal::DriftChamberDigiLocal(DriftChamberDigiLocalObj* obj) : m_obj(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj>(obj)) {}

DriftChamberDigiLocal DriftChamberDigiLocal::makeEmpty() {
  return {nullptr};
}

std::uint64_t DriftChamberDigiLocal::getCellID() const { return m_obj->data.cellID; }
float DriftChamberDigiLocal::getDistanceToWire() const { return m_obj->data.distanceToWire; }
float DriftChamberDigiLocal::getZPositionAlongWire() const { return m_obj->data.zPositionAlongWire; }
float DriftChamberDigiLocal::getTime() const { return m_obj->data.time; }
float DriftChamberDigiLocal::getEDep() const { return m_obj->data.eDep; }
float DriftChamberDigiLocal::getEDepError() const { return m_obj->data.eDepError; }
std::uint32_t DriftChamberDigiLocal::getClusterCount() const { return m_obj->data.clusterCount; }






bool DriftChamberDigiLocal::isAvailable() const {
  return m_obj;
}

const podio::ObjectID DriftChamberDigiLocal::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool DriftChamberDigiLocal::operator==(const MutableDriftChamberDigiLocal& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiLocal& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " cellID : " << value.getCellID() << '\n';
  o << " distanceToWire : " << value.getDistanceToWire() << '\n';
  o << " zPositionAlongWire : " << value.getZPositionAlongWire() << '\n';
  o << " time : " << value.getTime() << '\n';
  o << " eDep : " << value.getEDep() << '\n';
  o << " eDepError : " << value.getEDepError() << '\n';
  o << " clusterCount : " << value.getClusterCount() << '\n';



  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiLocal& value) {
  j = nlohmann::json{
    {"cellID", value.getCellID()}
    ,{"distanceToWire", value.getDistanceToWire()}
    ,{"zPositionAlongWire", value.getZPositionAlongWire()}
    ,{"time", value.getTime()}
    ,{"eDep", value.getEDep()}
    ,{"eDepError", value.getEDepError()}
    ,{"clusterCount", value.getClusterCount()}
  };


}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::DriftChamberDigiLocal& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
