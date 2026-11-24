// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/MCRecoDriftChamberDigiV2Association.h"
#include "extension/MutableMCRecoDriftChamberDigiV2Association.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/MCRecoDriftChamberDigiV2AssociationData.h"

#include "extension/DriftChamberDigiV2.h"
#include "edm4hep/SimTrackerHit.h"

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



MCRecoDriftChamberDigiV2Association::MCRecoDriftChamberDigiV2Association(const float weight) {
  m_obj->data.weight = weight;
}

MCRecoDriftChamberDigiV2Association& MCRecoDriftChamberDigiV2Association::operator=(MCRecoDriftChamberDigiV2Association other) & {
  swap(*this, other);
  return *this;
}

MutableMCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2Association::clone(bool cloneRelations) const {
  auto tmp = new MCRecoDriftChamberDigiV2AssociationObj(podio::ObjectID{}, m_obj->data);
  if (cloneRelations) {
  if (m_obj->m_digi) {
    tmp->m_digi = std::make_unique<extension::DriftChamberDigiV2>((*m_obj->m_digi));
  }
  if (m_obj->m_sim) {
    tmp->m_sim = std::make_unique<edm4hep::SimTrackerHit>((*m_obj->m_sim));
  }
  }
  else {
  }
  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

MCRecoDriftChamberDigiV2Association::MCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj) : m_obj(std::move(obj)) {}

MCRecoDriftChamberDigiV2Association::MCRecoDriftChamberDigiV2Association(const MutableMCRecoDriftChamberDigiV2Association& other): MCRecoDriftChamberDigiV2Association(other.m_obj) {}

MCRecoDriftChamberDigiV2Association::MCRecoDriftChamberDigiV2Association(MCRecoDriftChamberDigiV2AssociationObj* obj) : m_obj(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj>(obj)) {}

MCRecoDriftChamberDigiV2Association MCRecoDriftChamberDigiV2Association::makeEmpty() {
  return {nullptr};
}

float MCRecoDriftChamberDigiV2Association::getWeight() const { return m_obj->data.weight; }

const extension::DriftChamberDigiV2 MCRecoDriftChamberDigiV2Association::getDigi() const {
  if (!m_obj->m_digi) {
    return extension::DriftChamberDigiV2::makeEmpty();
  }
  return extension::DriftChamberDigiV2(*(m_obj->m_digi));
}

const edm4hep::SimTrackerHit MCRecoDriftChamberDigiV2Association::getSim() const {
  if (!m_obj->m_sim) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_sim));
}






bool MCRecoDriftChamberDigiV2Association::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MCRecoDriftChamberDigiV2Association::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MCRecoDriftChamberDigiV2Association::operator==(const MutableMCRecoDriftChamberDigiV2Association& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const MCRecoDriftChamberDigiV2Association& value) {
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
void to_json(nlohmann::json& j, const MCRecoDriftChamberDigiV2Association& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["digi"] = nlohmann::json{value.getDigi().id()};

    j["sim"] = nlohmann::json{value.getSim().id()};


}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::MCRecoDriftChamberDigiV2Association& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
