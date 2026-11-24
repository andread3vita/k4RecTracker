// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/SenseWireHit.h"
#include "extension/MutableSenseWireHit.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHitData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

namespace extension {



MutableSenseWireHit::MutableSenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError) {
  m_obj->data.cellID = cellID;
  m_obj->data.type = type;
  m_obj->data.quality = quality;
  m_obj->data.time = time;
  m_obj->data.eDep = eDep;
  m_obj->data.eDepError = eDepError;
  m_obj->data.wireStereoAngle = wireStereoAngle;
  m_obj->data.wireAzimuthalAngle = wireAzimuthalAngle;
  m_obj->data.position = position;
  m_obj->data.positionAlongWireError = positionAlongWireError;
  m_obj->data.distanceToWire = distanceToWire;
  m_obj->data.distanceToWireError = distanceToWireError;
}

MutableSenseWireHit& MutableSenseWireHit::operator=(MutableSenseWireHit other) & {
  swap(*this, other);
  return *this;
}

MutableSenseWireHit MutableSenseWireHit::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new SenseWireHitObj(podio::ObjectID{}, m_obj->data);
    tmp->m_nElectrons = new std::vector<std::uint16_t>();
    tmp->data.nElectrons_begin = 0;
    tmp->data.nElectrons_end = 0;
    return MutableSenseWireHit(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(new SenseWireHitObj(*m_obj), podio::utils::MarkOwned));
}

MutableSenseWireHit::MutableSenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj) : m_obj(std::move(obj)) {}

std::uint64_t MutableSenseWireHit::getCellID() const { return m_obj->data.cellID; }
std::int32_t MutableSenseWireHit::getType() const { return m_obj->data.type; }
std::int32_t MutableSenseWireHit::getQuality() const { return m_obj->data.quality; }
float MutableSenseWireHit::getTime() const { return m_obj->data.time; }
float MutableSenseWireHit::getEDep() const { return m_obj->data.eDep; }
float MutableSenseWireHit::getEDepError() const { return m_obj->data.eDepError; }
float MutableSenseWireHit::getWireStereoAngle() const { return m_obj->data.wireStereoAngle; }
float MutableSenseWireHit::getWireAzimuthalAngle() const { return m_obj->data.wireAzimuthalAngle; }
const edm4hep::Vector3d& MutableSenseWireHit::getPosition() const { return m_obj->data.position; }
double MutableSenseWireHit::getPositionAlongWireError() const { return m_obj->data.positionAlongWireError; }
float MutableSenseWireHit::getDistanceToWire() const { return m_obj->data.distanceToWire; }
float MutableSenseWireHit::getDistanceToWireError() const { return m_obj->data.distanceToWireError; }


void MutableSenseWireHit::setCellID(const std::uint64_t cellID) { m_obj->data.cellID = cellID; }
std::uint64_t& MutableSenseWireHit::getCellID() { return m_obj->data.cellID; }
void MutableSenseWireHit::setType(const std::int32_t type) { m_obj->data.type = type; }
std::int32_t& MutableSenseWireHit::getType() { return m_obj->data.type; }
void MutableSenseWireHit::setQuality(const std::int32_t quality) { m_obj->data.quality = quality; }
std::int32_t& MutableSenseWireHit::getQuality() { return m_obj->data.quality; }
void MutableSenseWireHit::setTime(const float time) { m_obj->data.time = time; }
float& MutableSenseWireHit::getTime() { return m_obj->data.time; }
void MutableSenseWireHit::setEDep(const float eDep) { m_obj->data.eDep = eDep; }
float& MutableSenseWireHit::getEDep() { return m_obj->data.eDep; }
void MutableSenseWireHit::setEDepError(const float eDepError) { m_obj->data.eDepError = eDepError; }
float& MutableSenseWireHit::getEDepError() { return m_obj->data.eDepError; }
void MutableSenseWireHit::setWireStereoAngle(const float wireStereoAngle) { m_obj->data.wireStereoAngle = wireStereoAngle; }
float& MutableSenseWireHit::getWireStereoAngle() { return m_obj->data.wireStereoAngle; }
void MutableSenseWireHit::setWireAzimuthalAngle(const float wireAzimuthalAngle) { m_obj->data.wireAzimuthalAngle = wireAzimuthalAngle; }
float& MutableSenseWireHit::getWireAzimuthalAngle() { return m_obj->data.wireAzimuthalAngle; }
void MutableSenseWireHit::setPosition(const edm4hep::Vector3d& position) { m_obj->data.position = position; }
edm4hep::Vector3d& MutableSenseWireHit::getPosition() { return m_obj->data.position; }
void MutableSenseWireHit::setPositionAlongWireError(const double positionAlongWireError) { m_obj->data.positionAlongWireError = positionAlongWireError; }
double& MutableSenseWireHit::getPositionAlongWireError() { return m_obj->data.positionAlongWireError; }
void MutableSenseWireHit::setDistanceToWire(const float distanceToWire) { m_obj->data.distanceToWire = distanceToWire; }
float& MutableSenseWireHit::getDistanceToWire() { return m_obj->data.distanceToWire; }
void MutableSenseWireHit::setDistanceToWireError(const float distanceToWireError) { m_obj->data.distanceToWireError = distanceToWireError; }
float& MutableSenseWireHit::getDistanceToWireError() { return m_obj->data.distanceToWireError; }


void MutableSenseWireHit::addToNElectrons(const std::uint16_t& component) {
  m_obj->m_nElectrons->push_back(component);
  m_obj->data.nElectrons_end++;
}

std::vector<std::uint16_t>::const_iterator MutableSenseWireHit::nElectrons_begin() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_begin);
  return ret_value;
}

std::vector<std::uint16_t>::const_iterator MutableSenseWireHit::nElectrons_end() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_end);
  return ret_value;
}

std::size_t MutableSenseWireHit::nElectrons_size() const {
  return m_obj->data.nElectrons_end - m_obj->data.nElectrons_begin;
}

std::uint16_t MutableSenseWireHit::getNElectrons(std::size_t index) const {
  if (nElectrons_size() > index) {
    return m_obj->m_nElectrons->at(m_obj->data.nElectrons_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<std::uint16_t> MutableSenseWireHit::getNElectrons() const {
  auto begin = m_obj->m_nElectrons->begin();
  std::advance(begin, m_obj->data.nElectrons_begin);
  auto end = m_obj->m_nElectrons->begin();
  std::advance(end, m_obj->data.nElectrons_end);
  return {begin, end};
}






bool MutableSenseWireHit::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableSenseWireHit::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableSenseWireHit::operator==(const SenseWireHit& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHit& value) {
  j = nlohmann::json{
    {"cellID", value.getCellID()}
    ,{"type", value.getType()}
    ,{"quality", value.getQuality()}
    ,{"time", value.getTime()}
    ,{"eDep", value.getEDep()}
    ,{"eDepError", value.getEDepError()}
    ,{"wireStereoAngle", value.getWireStereoAngle()}
    ,{"wireAzimuthalAngle", value.getWireAzimuthalAngle()}
    ,{"position", value.getPosition()}
    ,{"positionAlongWireError", value.getPositionAlongWireError()}
    ,{"distanceToWire", value.getDistanceToWire()}
    ,{"distanceToWireError", value.getDistanceToWireError()}
    ,{"nElectrons", value.getNElectrons()}
  };


}
#endif


} // namespace extension

