// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/SenseWireHitSimTrackerHitLink.h"
#include "extension/MutableSenseWireHitSimTrackerHitLink.h"
#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHitSimTrackerHitLinkData.h"

#include "extension/SenseWireHit.h"
#include "edm4hep/SimTrackerHit.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



SenseWireHitSimTrackerHitLink::SenseWireHitSimTrackerHitLink(const float weight) {
  m_obj->data.weight = weight;
}

SenseWireHitSimTrackerHitLink& SenseWireHitSimTrackerHitLink::operator=(SenseWireHitSimTrackerHitLink other) & {
  swap(*this, other);
  return *this;
}

MutableSenseWireHitSimTrackerHitLink SenseWireHitSimTrackerHitLink::clone(bool cloneRelations) const {
  auto tmp = new SenseWireHitSimTrackerHitLinkObj(podio::ObjectID{}, m_obj->data);
  if (cloneRelations) {
  if (m_obj->m_from) {
    tmp->m_from = std::make_unique<extension::SenseWireHit>((*m_obj->m_from));
  }
  if (m_obj->m_to) {
    tmp->m_to = std::make_unique<edm4hep::SimTrackerHit>((*m_obj->m_to));
  }
  }
  else {
  }
  return MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

SenseWireHitSimTrackerHitLink::SenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj) : m_obj(std::move(obj)) {}

SenseWireHitSimTrackerHitLink::SenseWireHitSimTrackerHitLink(const MutableSenseWireHitSimTrackerHitLink& other): SenseWireHitSimTrackerHitLink(other.m_obj) {}

SenseWireHitSimTrackerHitLink::SenseWireHitSimTrackerHitLink(SenseWireHitSimTrackerHitLinkObj* obj) : m_obj(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj>(obj)) {}

SenseWireHitSimTrackerHitLink SenseWireHitSimTrackerHitLink::makeEmpty() {
  return {nullptr};
}

float SenseWireHitSimTrackerHitLink::getWeight() const { return m_obj->data.weight; }

const extension::SenseWireHit SenseWireHitSimTrackerHitLink::getFrom() const {
  if (!m_obj->m_from) {
    return extension::SenseWireHit::makeEmpty();
  }
  return extension::SenseWireHit(*(m_obj->m_from));
}

const edm4hep::SimTrackerHit SenseWireHitSimTrackerHitLink::getTo() const {
  if (!m_obj->m_to) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_to));
}






bool SenseWireHitSimTrackerHitLink::isAvailable() const {
  return m_obj;
}

const podio::ObjectID SenseWireHitSimTrackerHitLink::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool SenseWireHitSimTrackerHitLink::operator==(const MutableSenseWireHitSimTrackerHitLink& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const SenseWireHitSimTrackerHitLink& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " weight : " << value.getWeight() << '\n';

  o << " from : " << value.getFrom().id() << '\n';
  o << " to : " << value.getTo().id() << '\n';


  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const SenseWireHitSimTrackerHitLink& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["from"] = nlohmann::json{value.getFrom().id()};

    j["to"] = nlohmann::json{value.getTo().id()};


}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::SenseWireHitSimTrackerHitLink& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
