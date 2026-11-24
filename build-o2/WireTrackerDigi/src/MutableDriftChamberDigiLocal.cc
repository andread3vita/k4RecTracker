// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigiLocal.h"
#include "extension/MutableDriftChamberDigiLocal.h"
#include "extension/DriftChamberDigiLocalObj.h"
#include "extension/DriftChamberDigiLocalData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

namespace extension {



MutableDriftChamberDigiLocal::MutableDriftChamberDigiLocal(const std::uint64_t cellID, const float distanceToWire, const float zPositionAlongWire, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount) {
  m_obj->data.cellID = cellID;
  m_obj->data.distanceToWire = distanceToWire;
  m_obj->data.zPositionAlongWire = zPositionAlongWire;
  m_obj->data.time = time;
  m_obj->data.eDep = eDep;
  m_obj->data.eDepError = eDepError;
  m_obj->data.clusterCount = clusterCount;
}

MutableDriftChamberDigiLocal& MutableDriftChamberDigiLocal::operator=(MutableDriftChamberDigiLocal other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigiLocal MutableDriftChamberDigiLocal::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new DriftChamberDigiLocalObj(podio::ObjectID{}, m_obj->data);
    return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr(new DriftChamberDigiLocalObj(*m_obj), podio::utils::MarkOwned));
}

MutableDriftChamberDigiLocal::MutableDriftChamberDigiLocal(podio::utils::MaybeSharedPtr<DriftChamberDigiLocalObj> obj) : m_obj(std::move(obj)) {}

std::uint64_t MutableDriftChamberDigiLocal::getCellID() const { return m_obj->data.cellID; }
float MutableDriftChamberDigiLocal::getDistanceToWire() const { return m_obj->data.distanceToWire; }
float MutableDriftChamberDigiLocal::getZPositionAlongWire() const { return m_obj->data.zPositionAlongWire; }
float MutableDriftChamberDigiLocal::getTime() const { return m_obj->data.time; }
float MutableDriftChamberDigiLocal::getEDep() const { return m_obj->data.eDep; }
float MutableDriftChamberDigiLocal::getEDepError() const { return m_obj->data.eDepError; }
std::uint32_t MutableDriftChamberDigiLocal::getClusterCount() const { return m_obj->data.clusterCount; }


void MutableDriftChamberDigiLocal::setCellID(const std::uint64_t cellID) { m_obj->data.cellID = cellID; }
std::uint64_t& MutableDriftChamberDigiLocal::getCellID() { return m_obj->data.cellID; }
void MutableDriftChamberDigiLocal::setDistanceToWire(const float distanceToWire) { m_obj->data.distanceToWire = distanceToWire; }
float& MutableDriftChamberDigiLocal::getDistanceToWire() { return m_obj->data.distanceToWire; }
void MutableDriftChamberDigiLocal::setZPositionAlongWire(const float zPositionAlongWire) { m_obj->data.zPositionAlongWire = zPositionAlongWire; }
float& MutableDriftChamberDigiLocal::getZPositionAlongWire() { return m_obj->data.zPositionAlongWire; }
void MutableDriftChamberDigiLocal::setTime(const float time) { m_obj->data.time = time; }
float& MutableDriftChamberDigiLocal::getTime() { return m_obj->data.time; }
void MutableDriftChamberDigiLocal::setEDep(const float eDep) { m_obj->data.eDep = eDep; }
float& MutableDriftChamberDigiLocal::getEDep() { return m_obj->data.eDep; }
void MutableDriftChamberDigiLocal::setEDepError(const float eDepError) { m_obj->data.eDepError = eDepError; }
float& MutableDriftChamberDigiLocal::getEDepError() { return m_obj->data.eDepError; }
void MutableDriftChamberDigiLocal::setClusterCount(const std::uint32_t clusterCount) { m_obj->data.clusterCount = clusterCount; }
std::uint32_t& MutableDriftChamberDigiLocal::getClusterCount() { return m_obj->data.clusterCount; }







bool MutableDriftChamberDigiLocal::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableDriftChamberDigiLocal::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableDriftChamberDigiLocal::operator==(const DriftChamberDigiLocal& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiLocal& value) {
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

