// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/SenseWireHit.h"
#include "extension/MutableSenseWireHit.h"
#include "extension/SenseWireHitObj.h"
#include "extension/SenseWireHitData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



SenseWireHit::SenseWireHit(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const float wireStereoAngle, const float wireAzimuthalAngle, const edm4hep::Vector3d& position, const double positionAlongWireError, const float distanceToWire, const float distanceToWireError) {
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

SenseWireHit& SenseWireHit::operator=(SenseWireHit other) & {
  swap(*this, other);
  return *this;
}

MutableSenseWireHit SenseWireHit::clone(bool cloneRelations) const {
  auto tmp = new SenseWireHitObj(podio::ObjectID{}, m_obj->data);
  tmp->m_nElectrons = new std::vector<std::uint16_t>();
  if (cloneRelations) {
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_nElectrons->reserve(m_obj->m_nElectrons->size());
    for (size_t i = m_obj->data.nElectrons_begin; i < m_obj->data.nElectrons_end; i++) {
      tmp->m_nElectrons->emplace_back((*m_obj->m_nElectrons)[i]);
    }
    tmp->data.nElectrons_begin = 0;
    tmp->data.nElectrons_end = tmp->m_nElectrons->size();
  }
  else {
    tmp->data.nElectrons_begin = 0;
    tmp->data.nElectrons_end = 0;
  }
  return MutableSenseWireHit(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

SenseWireHit::SenseWireHit(podio::utils::MaybeSharedPtr<SenseWireHitObj> obj) : m_obj(std::move(obj)) {}

SenseWireHit::SenseWireHit(const MutableSenseWireHit& other): SenseWireHit(other.m_obj) {}

SenseWireHit::SenseWireHit(SenseWireHitObj* obj) : m_obj(podio::utils::MaybeSharedPtr<SenseWireHitObj>(obj)) {}

SenseWireHit SenseWireHit::makeEmpty() {
  return {nullptr};
}

std::uint64_t SenseWireHit::getCellID() const { return m_obj->data.cellID; }
std::int32_t SenseWireHit::getType() const { return m_obj->data.type; }
std::int32_t SenseWireHit::getQuality() const { return m_obj->data.quality; }
float SenseWireHit::getTime() const { return m_obj->data.time; }
float SenseWireHit::getEDep() const { return m_obj->data.eDep; }
float SenseWireHit::getEDepError() const { return m_obj->data.eDepError; }
float SenseWireHit::getWireStereoAngle() const { return m_obj->data.wireStereoAngle; }
float SenseWireHit::getWireAzimuthalAngle() const { return m_obj->data.wireAzimuthalAngle; }
const edm4hep::Vector3d& SenseWireHit::getPosition() const { return m_obj->data.position; }
double SenseWireHit::getPositionAlongWireError() const { return m_obj->data.positionAlongWireError; }
float SenseWireHit::getDistanceToWire() const { return m_obj->data.distanceToWire; }
float SenseWireHit::getDistanceToWireError() const { return m_obj->data.distanceToWireError; }



std::vector<std::uint16_t>::const_iterator SenseWireHit::nElectrons_begin() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_begin);
  return ret_value;
}

std::vector<std::uint16_t>::const_iterator SenseWireHit::nElectrons_end() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_end);
  return ret_value;
}

std::size_t SenseWireHit::nElectrons_size() const {
  return m_obj->data.nElectrons_end - m_obj->data.nElectrons_begin;
}

std::uint16_t SenseWireHit::getNElectrons(std::size_t index) const {
  if (nElectrons_size() > index) {
    return m_obj->m_nElectrons->at(m_obj->data.nElectrons_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<std::uint16_t> SenseWireHit::getNElectrons() const {
  auto begin = m_obj->m_nElectrons->begin();
  std::advance(begin, m_obj->data.nElectrons_begin);
  auto end = m_obj->m_nElectrons->begin();
  std::advance(end, m_obj->data.nElectrons_end);
  return {begin, end};
}





bool SenseWireHit::isAvailable() const {
  return m_obj;
}

const podio::ObjectID SenseWireHit::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool SenseWireHit::operator==(const MutableSenseWireHit& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const SenseWireHit& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " cellID : " << value.getCellID() << '\n';
  o << " type : " << value.getType() << '\n';
  o << " quality : " << value.getQuality() << '\n';
  o << " time : " << value.getTime() << '\n';
  o << " eDep : " << value.getEDep() << '\n';
  o << " eDepError : " << value.getEDepError() << '\n';
  o << " wireStereoAngle : " << value.getWireStereoAngle() << '\n';
  o << " wireAzimuthalAngle : " << value.getWireAzimuthalAngle() << '\n';
  o << " position : " << value.getPosition() << '\n';
  o << " positionAlongWireError : " << value.getPositionAlongWireError() << '\n';
  o << " distanceToWire : " << value.getDistanceToWire() << '\n';
  o << " distanceToWireError : " << value.getDistanceToWireError() << '\n';


  o << " nElectrons : ";
  for (unsigned i = 0; i < value.nElectrons_size(); ++i) {
    o << value.getNElectrons(i) << " ";
  }
  o << '\n';

  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHit& value) {
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


podio::detail::OrderKey podio::detail::getOrderKey(const extension::SenseWireHit& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
