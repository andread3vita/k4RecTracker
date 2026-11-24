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

#include <ostream>

namespace extension {



MCRecoDriftChamberDigiAssociation::MCRecoDriftChamberDigiAssociation(const float weight) {
  m_obj->data.weight = weight;
}

MCRecoDriftChamberDigiAssociation& MCRecoDriftChamberDigiAssociation::operator=(MCRecoDriftChamberDigiAssociation other) & {
  swap(*this, other);
  return *this;
}

MutableMCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociation::clone(bool cloneRelations) const {
  auto tmp = new MCRecoDriftChamberDigiAssociationObj(podio::ObjectID{}, m_obj->data);
  if (cloneRelations) {
  if (m_obj->m_digi) {
    tmp->m_digi = std::make_unique<extension::DriftChamberDigi>((*m_obj->m_digi));
  }
  if (m_obj->m_sim) {
    tmp->m_sim = std::make_unique<edm4hep::SimTrackerHit>((*m_obj->m_sim));
  }
  }
  else {
  }
  return MutableMCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

MCRecoDriftChamberDigiAssociation::MCRecoDriftChamberDigiAssociation(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj> obj) : m_obj(std::move(obj)) {}

MCRecoDriftChamberDigiAssociation::MCRecoDriftChamberDigiAssociation(const MutableMCRecoDriftChamberDigiAssociation& other): MCRecoDriftChamberDigiAssociation(other.m_obj) {}

MCRecoDriftChamberDigiAssociation::MCRecoDriftChamberDigiAssociation(MCRecoDriftChamberDigiAssociationObj* obj) : m_obj(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiAssociationObj>(obj)) {}

MCRecoDriftChamberDigiAssociation MCRecoDriftChamberDigiAssociation::makeEmpty() {
  return {nullptr};
}

float MCRecoDriftChamberDigiAssociation::getWeight() const { return m_obj->data.weight; }

const extension::DriftChamberDigi MCRecoDriftChamberDigiAssociation::getDigi() const {
  if (!m_obj->m_digi) {
    return extension::DriftChamberDigi::makeEmpty();
  }
  return extension::DriftChamberDigi(*(m_obj->m_digi));
}

const edm4hep::SimTrackerHit MCRecoDriftChamberDigiAssociation::getSim() const {
  if (!m_obj->m_sim) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_sim));
}






bool MCRecoDriftChamberDigiAssociation::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MCRecoDriftChamberDigiAssociation::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MCRecoDriftChamberDigiAssociation::operator==(const MutableMCRecoDriftChamberDigiAssociation& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiAssociation& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " weight : " << value.getWeight() << '\n';

  o << " digi : " << value.getDigi().id() << '\n';
  o << " sim : " << value.getSim().id() << '\n';


  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiAssociation& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["digi"] = nlohmann::json{value.getDigi().id()};

    j["sim"] = nlohmann::json{value.getSim().id()};


}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::MCRecoDriftChamberDigiAssociation& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
