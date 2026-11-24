// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiAssociation.h"
#include "extension/MutableMCRecoDriftChamberDigiAssociation.h"
#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/MCRecoDriftChamberDigiAssociationData.h"

#include "extension/DriftChamberDigi.h"
#include "edm4hep/SimTrackerHit.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

namespace extension {



MutableMCRecoDriftChamberDigiAssociation::MutableMCRecoDriftChamberDigiAssociation(const float weight) {
  m_obj->data.weight = weight;
}

MutableMCRecoDriftChamberDigiAssociation& MutableMCRecoDriftChamberDigiAssociation::operator=(MutableMCRecoDriftChamberDigiAssociation other) & {
  swap(*this, other);
  return *this;
}

MutableMCRecoDriftChamberDigiAssociation MutableMCRecoDriftChamberDigiAssociation::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new MCRecoDriftChamberDigiAssociationObj(podio::ObjectID{}, m_obj->data);
    return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(new MCRecoDriftChamberDigiAssociationObj(*m_obj), podio::utils::MarkOwned));
}

MutableMCRecoDriftChamberDigiAssociation::MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj) : m_obj(std::move(obj)) {}

float MutableMCRecoDriftChamberDigiAssociation::getWeight() const { return m_obj->data.weight; }

const extension::DriftChamberDigi MutableMCRecoDriftChamberDigiAssociation::getDigi() const {
  if (!m_obj->m_digi) {
    return extension::DriftChamberDigi::makeEmpty();
  }
  return extension::DriftChamberDigi(*(m_obj->m_digi));
}

const edm4hep::SimTrackerHit MutableMCRecoDriftChamberDigiAssociation::getSim() const {
  if (!m_obj->m_sim) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_sim));
}


void MutableMCRecoDriftChamberDigiAssociation::setWeight(const float weight) { m_obj->data.weight = weight; }
float& MutableMCRecoDriftChamberDigiAssociation::getWeight() { return m_obj->data.weight; }

void MutableMCRecoDriftChamberDigiAssociation::setDigi(const extension::DriftChamberDigi& value) {
  m_obj->m_digi = std::make_unique<extension::DriftChamberDigi>(value);
}

void MutableMCRecoDriftChamberDigiAssociation::setSim(const edm4hep::SimTrackerHit& value) {
  m_obj->m_sim = std::make_unique<edm4hep::SimTrackerHit>(value);
}







bool MutableMCRecoDriftChamberDigiAssociation::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableMCRecoDriftChamberDigiAssociation::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableMCRecoDriftChamberDigiAssociation::operator==(const MCRecoDriftChamberDigiAssociation& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiAssociation& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["digi"] = nlohmann::json{value.getDigi().id()};

    j["sim"] = nlohmann::json{value.getSim().id()};


}
#endif


} // namespace extension

