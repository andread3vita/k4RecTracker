// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigi.h"
#include "extension/MutableDriftChamberDigi.h"
#include "extension/DriftChamberDigiObj.h"
#include "extension/DriftChamberDigiData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



DriftChamberDigi::DriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position) {
  m_obj->data.cellID = cellID;
  m_obj->data.leftPosition = leftPosition;
  m_obj->data.rightPosition = rightPosition;
  m_obj->data.time = time;
  m_obj->data.eDep = eDep;
  m_obj->data.eDepError = eDepError;
  m_obj->data.clusterCount = clusterCount;
  m_obj->data.type = type;
  m_obj->data.quality = quality;
  m_obj->data.position = position;
}

DriftChamberDigi& DriftChamberDigi::operator=(DriftChamberDigi other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigi DriftChamberDigi::clone(bool cloneRelations) const {
  auto tmp = new DriftChamberDigiObj(podio::ObjectID{}, m_obj->data);
  if (cloneRelations) {
  }
  else {
  }
  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

DriftChamberDigi::DriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj) : m_obj(std::move(obj)) {}

DriftChamberDigi::DriftChamberDigi(const MutableDriftChamberDigi& other): DriftChamberDigi(other.m_obj) {}

DriftChamberDigi::DriftChamberDigi(DriftChamberDigiObj* obj) : m_obj(podio::utils::MaybeSharedPtr<DriftChamberDigiObj>(obj)) {}

DriftChamberDigi DriftChamberDigi::makeEmpty() {
  return {nullptr};
}

std::uint64_t DriftChamberDigi::getCellID() const { return m_obj->data.cellID; }
const edm4hep::Vector3d& DriftChamberDigi::getLeftPosition() const { return m_obj->data.leftPosition; }
const edm4hep::Vector3d& DriftChamberDigi::getRightPosition() const { return m_obj->data.rightPosition; }
float DriftChamberDigi::getTime() const { return m_obj->data.time; }
float DriftChamberDigi::getEDep() const { return m_obj->data.eDep; }
float DriftChamberDigi::getEDepError() const { return m_obj->data.eDepError; }
std::uint32_t DriftChamberDigi::getClusterCount() const { return m_obj->data.clusterCount; }
std::int32_t DriftChamberDigi::getType() const { return m_obj->data.type; }
std::int32_t DriftChamberDigi::getQuality() const { return m_obj->data.quality; }
const edm4hep::Vector3d& DriftChamberDigi::getPosition() const { return m_obj->data.position; }






bool DriftChamberDigi::isAvailable() const {
  return m_obj;
}

const podio::ObjectID DriftChamberDigi::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool DriftChamberDigi::operator==(const MutableDriftChamberDigi& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const DriftChamberDigi& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " cellID : " << value.getCellID() << '\n';
  o << " leftPosition : " << value.getLeftPosition() << '\n';
  o << " rightPosition : " << value.getRightPosition() << '\n';
  o << " time : " << value.getTime() << '\n';
  o << " eDep : " << value.getEDep() << '\n';
  o << " eDepError : " << value.getEDepError() << '\n';
  o << " clusterCount : " << value.getClusterCount() << '\n';
  o << " type : " << value.getType() << '\n';
  o << " quality : " << value.getQuality() << '\n';
  o << " position : " << value.getPosition() << '\n';



  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigi& value) {
  j = nlohmann::json{
    {"cellID", value.getCellID()}
    ,{"leftPosition", value.getLeftPosition()}
    ,{"rightPosition", value.getRightPosition()}
    ,{"time", value.getTime()}
    ,{"eDep", value.getEDep()}
    ,{"eDepError", value.getEDepError()}
    ,{"clusterCount", value.getClusterCount()}
    ,{"type", value.getType()}
    ,{"quality", value.getQuality()}
    ,{"position", value.getPosition()}
  };


}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::DriftChamberDigi& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
