// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigiV2.h"
#include "extension/MutableDriftChamberDigiV2.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2Data.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

namespace extension {



MutableDriftChamberDigiV2::MutableDriftChamberDigiV2(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const edm4hep::Vector3d& position, const edm4hep::Vector3d& directionSW, const float distanceToWire, const std::uint32_t nCluster) {
  m_obj->data.cellID = cellID;
  m_obj->data.type = type;
  m_obj->data.quality = quality;
  m_obj->data.time = time;
  m_obj->data.eDep = eDep;
  m_obj->data.eDepError = eDepError;
  m_obj->data.position = position;
  m_obj->data.directionSW = directionSW;
  m_obj->data.distanceToWire = distanceToWire;
  m_obj->data.nCluster = nCluster;
}

MutableDriftChamberDigiV2& MutableDriftChamberDigiV2::operator=(MutableDriftChamberDigiV2 other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigiV2 MutableDriftChamberDigiV2::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new DriftChamberDigiV2Obj(podio::ObjectID{}, m_obj->data);
    tmp->m_nElectrons = new std::vector<std::uint16_t>();
    tmp->data.nElectrons_begin = 0;
    tmp->data.nElectrons_end = 0;
    return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(new DriftChamberDigiV2Obj(*m_obj), podio::utils::MarkOwned));
}

MutableDriftChamberDigiV2::MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> obj) : m_obj(std::move(obj)) {}

std::uint64_t MutableDriftChamberDigiV2::getCellID() const { return m_obj->data.cellID; }
std::int32_t MutableDriftChamberDigiV2::getType() const { return m_obj->data.type; }
std::int32_t MutableDriftChamberDigiV2::getQuality() const { return m_obj->data.quality; }
float MutableDriftChamberDigiV2::getTime() const { return m_obj->data.time; }
float MutableDriftChamberDigiV2::getEDep() const { return m_obj->data.eDep; }
float MutableDriftChamberDigiV2::getEDepError() const { return m_obj->data.eDepError; }
const edm4hep::Vector3d& MutableDriftChamberDigiV2::getPosition() const { return m_obj->data.position; }
const edm4hep::Vector3d& MutableDriftChamberDigiV2::getDirectionSW() const { return m_obj->data.directionSW; }
float MutableDriftChamberDigiV2::getDistanceToWire() const { return m_obj->data.distanceToWire; }
std::uint32_t MutableDriftChamberDigiV2::getNCluster() const { return m_obj->data.nCluster; }


void MutableDriftChamberDigiV2::setCellID(const std::uint64_t cellID) { m_obj->data.cellID = cellID; }
std::uint64_t& MutableDriftChamberDigiV2::getCellID() { return m_obj->data.cellID; }
void MutableDriftChamberDigiV2::setType(const std::int32_t type) { m_obj->data.type = type; }
std::int32_t& MutableDriftChamberDigiV2::getType() { return m_obj->data.type; }
void MutableDriftChamberDigiV2::setQuality(const std::int32_t quality) { m_obj->data.quality = quality; }
std::int32_t& MutableDriftChamberDigiV2::getQuality() { return m_obj->data.quality; }
void MutableDriftChamberDigiV2::setTime(const float time) { m_obj->data.time = time; }
float& MutableDriftChamberDigiV2::getTime() { return m_obj->data.time; }
void MutableDriftChamberDigiV2::setEDep(const float eDep) { m_obj->data.eDep = eDep; }
float& MutableDriftChamberDigiV2::getEDep() { return m_obj->data.eDep; }
void MutableDriftChamberDigiV2::setEDepError(const float eDepError) { m_obj->data.eDepError = eDepError; }
float& MutableDriftChamberDigiV2::getEDepError() { return m_obj->data.eDepError; }
void MutableDriftChamberDigiV2::setPosition(const edm4hep::Vector3d& position) { m_obj->data.position = position; }
edm4hep::Vector3d& MutableDriftChamberDigiV2::getPosition() { return m_obj->data.position; }
void MutableDriftChamberDigiV2::setDirectionSW(const edm4hep::Vector3d& directionSW) { m_obj->data.directionSW = directionSW; }
edm4hep::Vector3d& MutableDriftChamberDigiV2::getDirectionSW() { return m_obj->data.directionSW; }
void MutableDriftChamberDigiV2::setDistanceToWire(const float distanceToWire) { m_obj->data.distanceToWire = distanceToWire; }
float& MutableDriftChamberDigiV2::getDistanceToWire() { return m_obj->data.distanceToWire; }
void MutableDriftChamberDigiV2::setNCluster(const std::uint32_t nCluster) { m_obj->data.nCluster = nCluster; }
std::uint32_t& MutableDriftChamberDigiV2::getNCluster() { return m_obj->data.nCluster; }


void MutableDriftChamberDigiV2::addToNElectrons(const std::uint16_t& component) {
  m_obj->m_nElectrons->push_back(component);
  m_obj->data.nElectrons_end++;
}

std::vector<std::uint16_t>::const_iterator MutableDriftChamberDigiV2::nElectrons_begin() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_begin);
  return ret_value;
}

std::vector<std::uint16_t>::const_iterator MutableDriftChamberDigiV2::nElectrons_end() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_end);
  return ret_value;
}

std::size_t MutableDriftChamberDigiV2::nElectrons_size() const {
  return m_obj->data.nElectrons_end - m_obj->data.nElectrons_begin;
}

std::uint16_t MutableDriftChamberDigiV2::getNElectrons(std::size_t index) const {
  if (nElectrons_size() > index) {
    return m_obj->m_nElectrons->at(m_obj->data.nElectrons_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<std::uint16_t> MutableDriftChamberDigiV2::getNElectrons() const {
  auto begin = m_obj->m_nElectrons->begin();
  std::advance(begin, m_obj->data.nElectrons_begin);
  auto end = m_obj->m_nElectrons->begin();
  std::advance(end, m_obj->data.nElectrons_end);
  return {begin, end};
}






bool MutableDriftChamberDigiV2::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableDriftChamberDigiV2::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableDriftChamberDigiV2::operator==(const DriftChamberDigiV2& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableDriftChamberDigiV2& value) {
  j = nlohmann::json{
    {"cellID", value.getCellID()}
    ,{"type", value.getType()}
    ,{"quality", value.getQuality()}
    ,{"time", value.getTime()}
    ,{"eDep", value.getEDep()}
    ,{"eDepError", value.getEDepError()}
    ,{"position", value.getPosition()}
    ,{"directionSW", value.getDirectionSW()}
    ,{"distanceToWire", value.getDistanceToWire()}
    ,{"nCluster", value.getNCluster()}
    ,{"nElectrons", value.getNElectrons()}
  };


}
#endif


} // namespace extension

