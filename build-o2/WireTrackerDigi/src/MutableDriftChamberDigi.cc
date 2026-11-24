// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigi.h"
#include "extension/MutableDriftChamberDigi.h"
#include "extension/DriftChamberDigiObj.h"
#include "extension/DriftChamberDigiData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

namespace extension {



MutableDriftChamberDigi::MutableDriftChamberDigi(const std::uint64_t cellID, const edm4hep::Vector3d& leftPosition, const edm4hep::Vector3d& rightPosition, const float time, const float eDep, const float eDepError, const std::uint32_t clusterCount, const std::int32_t type, const std::int32_t quality, const edm4hep::Vector3d& position) {
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

MutableDriftChamberDigi& MutableDriftChamberDigi::operator=(MutableDriftChamberDigi other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigi MutableDriftChamberDigi::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new DriftChamberDigiObj(podio::ObjectID{}, m_obj->data);
    return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableDriftChamberDigi(podio::utils::MaybeSharedPtr(new DriftChamberDigiObj(*m_obj), podio::utils::MarkOwned));
}

MutableDriftChamberDigi::MutableDriftChamberDigi(podio::utils::MaybeSharedPtr<DriftChamberDigiObj> obj) : m_obj(std::move(obj)) {}

std::uint64_t MutableDriftChamberDigi::getCellID() const { return m_obj->data.cellID; }
const edm4hep::Vector3d& MutableDriftChamberDigi::getLeftPosition() const { return m_obj->data.leftPosition; }
const edm4hep::Vector3d& MutableDriftChamberDigi::getRightPosition() const { return m_obj->data.rightPosition; }
float MutableDriftChamberDigi::getTime() const { return m_obj->data.time; }
float MutableDriftChamberDigi::getEDep() const { return m_obj->data.eDep; }
float MutableDriftChamberDigi::getEDepError() const { return m_obj->data.eDepError; }
std::uint32_t MutableDriftChamberDigi::getClusterCount() const { return m_obj->data.clusterCount; }
std::int32_t MutableDriftChamberDigi::getType() const { return m_obj->data.type; }
std::int32_t MutableDriftChamberDigi::getQuality() const { return m_obj->data.quality; }
const edm4hep::Vector3d& MutableDriftChamberDigi::getPosition() const { return m_obj->data.position; }


void MutableDriftChamberDigi::setCellID(const std::uint64_t cellID) { m_obj->data.cellID = cellID; }
std::uint64_t& MutableDriftChamberDigi::getCellID() { return m_obj->data.cellID; }
void MutableDriftChamberDigi::setLeftPosition(const edm4hep::Vector3d& leftPosition) { m_obj->data.leftPosition = leftPosition; }
edm4hep::Vector3d& MutableDriftChamberDigi::getLeftPosition() { return m_obj->data.leftPosition; }
void MutableDriftChamberDigi::setRightPosition(const edm4hep::Vector3d& rightPosition) { m_obj->data.rightPosition = rightPosition; }
edm4hep::Vector3d& MutableDriftChamberDigi::getRightPosition() { return m_obj->data.rightPosition; }
void MutableDriftChamberDigi::setTime(const float time) { m_obj->data.time = time; }
float& MutableDriftChamberDigi::getTime() { return m_obj->data.time; }
void MutableDriftChamberDigi::setEDep(const float eDep) { m_obj->data.eDep = eDep; }
float& MutableDriftChamberDigi::getEDep() { return m_obj->data.eDep; }
void MutableDriftChamberDigi::setEDepError(const float eDepError) { m_obj->data.eDepError = eDepError; }
float& MutableDriftChamberDigi::getEDepError() { return m_obj->data.eDepError; }
void MutableDriftChamberDigi::setClusterCount(const std::uint32_t clusterCount) { m_obj->data.clusterCount = clusterCount; }
std::uint32_t& MutableDriftChamberDigi::getClusterCount() { return m_obj->data.clusterCount; }
void MutableDriftChamberDigi::setType(const std::int32_t type) { m_obj->data.type = type; }
std::int32_t& MutableDriftChamberDigi::getType() { return m_obj->data.type; }
void MutableDriftChamberDigi::setQuality(const std::int32_t quality) { m_obj->data.quality = quality; }
std::int32_t& MutableDriftChamberDigi::getQuality() { return m_obj->data.quality; }
void MutableDriftChamberDigi::setPosition(const edm4hep::Vector3d& position) { m_obj->data.position = position; }
edm4hep::Vector3d& MutableDriftChamberDigi::getPosition() { return m_obj->data.position; }







bool MutableDriftChamberDigi::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableDriftChamberDigi::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableDriftChamberDigi::operator==(const DriftChamberDigi& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigi& value) {
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

