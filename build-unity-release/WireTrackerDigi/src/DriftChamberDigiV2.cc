// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/DriftChamberDigiV2.h"
#include "extension/MutableDriftChamberDigiV2.h"
#include "extension/DriftChamberDigiV2Obj.h"
#include "extension/DriftChamberDigiV2Data.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



DriftChamberDigiV2::DriftChamberDigiV2(const std::uint64_t cellID, const std::int32_t type, const std::int32_t quality, const float time, const float eDep, const float eDepError, const edm4hep::Vector3d& position, const edm4hep::Vector3d& directionSW, const float distanceToWire, const std::uint32_t nCluster) {
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

DriftChamberDigiV2& DriftChamberDigiV2::operator=(DriftChamberDigiV2 other) & {
  swap(*this, other);
  return *this;
}

MutableDriftChamberDigiV2 DriftChamberDigiV2::clone(bool cloneRelations) const {
  auto tmp = new DriftChamberDigiV2Obj(podio::ObjectID{}, m_obj->data);
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
  return MutableDriftChamberDigiV2(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

DriftChamberDigiV2::DriftChamberDigiV2(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj> obj) : m_obj(std::move(obj)) {}

DriftChamberDigiV2::DriftChamberDigiV2(const MutableDriftChamberDigiV2& other): DriftChamberDigiV2(other.m_obj) {}

DriftChamberDigiV2::DriftChamberDigiV2(DriftChamberDigiV2Obj* obj) : m_obj(podio::utils::MaybeSharedPtr<DriftChamberDigiV2Obj>(obj)) {}

DriftChamberDigiV2 DriftChamberDigiV2::makeEmpty() {
  return {nullptr};
}

std::uint64_t DriftChamberDigiV2::getCellID() const { return m_obj->data.cellID; }
std::int32_t DriftChamberDigiV2::getType() const { return m_obj->data.type; }
std::int32_t DriftChamberDigiV2::getQuality() const { return m_obj->data.quality; }
float DriftChamberDigiV2::getTime() const { return m_obj->data.time; }
float DriftChamberDigiV2::getEDep() const { return m_obj->data.eDep; }
float DriftChamberDigiV2::getEDepError() const { return m_obj->data.eDepError; }
const edm4hep::Vector3d& DriftChamberDigiV2::getPosition() const { return m_obj->data.position; }
const edm4hep::Vector3d& DriftChamberDigiV2::getDirectionSW() const { return m_obj->data.directionSW; }
float DriftChamberDigiV2::getDistanceToWire() const { return m_obj->data.distanceToWire; }
std::uint32_t DriftChamberDigiV2::getNCluster() const { return m_obj->data.nCluster; }



std::vector<std::uint16_t>::const_iterator DriftChamberDigiV2::nElectrons_begin() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_begin);
  return ret_value;
}

std::vector<std::uint16_t>::const_iterator DriftChamberDigiV2::nElectrons_end() const {
  auto ret_value = m_obj->m_nElectrons->begin();
  std::advance(ret_value, m_obj->data.nElectrons_end);
  return ret_value;
}

std::size_t DriftChamberDigiV2::nElectrons_size() const {
  return m_obj->data.nElectrons_end - m_obj->data.nElectrons_begin;
}

std::uint16_t DriftChamberDigiV2::getNElectrons(std::size_t index) const {
  if (nElectrons_size() > index) {
    return m_obj->m_nElectrons->at(m_obj->data.nElectrons_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<std::uint16_t> DriftChamberDigiV2::getNElectrons() const {
  auto begin = m_obj->m_nElectrons->begin();
  std::advance(begin, m_obj->data.nElectrons_begin);
  auto end = m_obj->m_nElectrons->begin();
  std::advance(end, m_obj->data.nElectrons_end);
  return {begin, end};
}





bool DriftChamberDigiV2::isAvailable() const {
  return m_obj;
}

const podio::ObjectID DriftChamberDigiV2::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool DriftChamberDigiV2::operator==(const MutableDriftChamberDigiV2& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const DriftChamberDigiV2& value) {
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
  o << " position : " << value.getPosition() << '\n';
  o << " directionSW : " << value.getDirectionSW() << '\n';
  o << " distanceToWire : " << value.getDistanceToWire() << '\n';
  o << " nCluster : " << value.getNCluster() << '\n';


  o << " nElectrons : ";
  for (unsigned i = 0; i < value.nElectrons_size(); ++i) {
    o << value.getNElectrons(i) << " ";
  }
  o << '\n';

  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const DriftChamberDigiV2& value) {
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


podio::detail::OrderKey podio::detail::getOrderKey(const extension::DriftChamberDigiV2& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
