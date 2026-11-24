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

namespace extension {



MutableSenseWireHitSimTrackerHitLink::MutableSenseWireHitSimTrackerHitLink(const float weight) {
  m_obj->data.weight = weight;
}

MutableSenseWireHitSimTrackerHitLink& MutableSenseWireHitSimTrackerHitLink::operator=(MutableSenseWireHitSimTrackerHitLink other) & {
  swap(*this, other);
  return *this;
}

MutableSenseWireHitSimTrackerHitLink MutableSenseWireHitSimTrackerHitLink::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new SenseWireHitSimTrackerHitLinkObj(podio::ObjectID{}, m_obj->data);
    return MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr(new SenseWireHitSimTrackerHitLinkObj(*m_obj), podio::utils::MarkOwned));
}

MutableSenseWireHitSimTrackerHitLink::MutableSenseWireHitSimTrackerHitLink(podio::utils::MaybeSharedPtr<SenseWireHitSimTrackerHitLinkObj> obj) : m_obj(std::move(obj)) {}

float MutableSenseWireHitSimTrackerHitLink::getWeight() const { return m_obj->data.weight; }

const extension::SenseWireHit MutableSenseWireHitSimTrackerHitLink::getFrom() const {
  if (!m_obj->m_from) {
    return extension::SenseWireHit::makeEmpty();
  }
  return extension::SenseWireHit(*(m_obj->m_from));
}

const edm4hep::SimTrackerHit MutableSenseWireHitSimTrackerHitLink::getTo() const {
  if (!m_obj->m_to) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_to));
}


void MutableSenseWireHitSimTrackerHitLink::setWeight(const float weight) { m_obj->data.weight = weight; }
float& MutableSenseWireHitSimTrackerHitLink::getWeight() { return m_obj->data.weight; }

void MutableSenseWireHitSimTrackerHitLink::setFrom(const extension::SenseWireHit& value) {
  m_obj->m_from = std::make_unique<extension::SenseWireHit>(value);
}

void MutableSenseWireHitSimTrackerHitLink::setTo(const edm4hep::SimTrackerHit& value) {
  m_obj->m_to = std::make_unique<edm4hep::SimTrackerHit>(value);
}







bool MutableSenseWireHitSimTrackerHitLink::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableSenseWireHitSimTrackerHitLink::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableSenseWireHitSimTrackerHitLink::operator==(const SenseWireHitSimTrackerHitLink& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableSenseWireHitSimTrackerHitLink& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["from"] = nlohmann::json{value.getFrom().id()};

    j["to"] = nlohmann::json{value.getTo().id()};


}
#endif


} // namespace extension

