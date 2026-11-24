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

namespace extension {



MutableMCRecoDriftChamberDigiV2Association::MutableMCRecoDriftChamberDigiV2Association(const float weight) {
  m_obj->data.weight = weight;
}

MutableMCRecoDriftChamberDigiV2Association& MutableMCRecoDriftChamberDigiV2Association::operator=(MutableMCRecoDriftChamberDigiV2Association other) & {
  swap(*this, other);
  return *this;
}

MutableMCRecoDriftChamberDigiV2Association MutableMCRecoDriftChamberDigiV2Association::clone(bool cloneRelations) const {
  if (!cloneRelations) {
    auto tmp = new MCRecoDriftChamberDigiV2AssociationObj(podio::ObjectID{}, m_obj->data);
    return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
  }
  return MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr(new MCRecoDriftChamberDigiV2AssociationObj(*m_obj), podio::utils::MarkOwned));
}

MutableMCRecoDriftChamberDigiV2Association::MutableMCRecoDriftChamberDigiV2Association(podio::utils::MaybeSharedPtr<MCRecoDriftChamberDigiV2AssociationObj> obj) : m_obj(std::move(obj)) {}

float MutableMCRecoDriftChamberDigiV2Association::getWeight() const { return m_obj->data.weight; }

const extension::DriftChamberDigiV2 MutableMCRecoDriftChamberDigiV2Association::getDigi() const {
  if (!m_obj->m_digi) {
    return extension::DriftChamberDigiV2::makeEmpty();
  }
  return extension::DriftChamberDigiV2(*(m_obj->m_digi));
}

const edm4hep::SimTrackerHit MutableMCRecoDriftChamberDigiV2Association::getSim() const {
  if (!m_obj->m_sim) {
    return edm4hep::SimTrackerHit::makeEmpty();
  }
  return edm4hep::SimTrackerHit(*(m_obj->m_sim));
}


void MutableMCRecoDriftChamberDigiV2Association::setWeight(const float weight) { m_obj->data.weight = weight; }
float& MutableMCRecoDriftChamberDigiV2Association::getWeight() { return m_obj->data.weight; }

void MutableMCRecoDriftChamberDigiV2Association::setDigi(const extension::DriftChamberDigiV2& value) {
  m_obj->m_digi = std::make_unique<extension::DriftChamberDigiV2>(value);
}

void MutableMCRecoDriftChamberDigiV2Association::setSim(const edm4hep::SimTrackerHit& value) {
  m_obj->m_sim = std::make_unique<edm4hep::SimTrackerHit>(value);
}







bool MutableMCRecoDriftChamberDigiV2Association::isAvailable() const {
  return m_obj;
}

const podio::ObjectID MutableMCRecoDriftChamberDigiV2Association::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool MutableMCRecoDriftChamberDigiV2Association::operator==(const MCRecoDriftChamberDigiV2Association& other) const {
  return m_obj == other.m_obj;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const MutableMCRecoDriftChamberDigiV2Association& value) {
  j = nlohmann::json{
    {"weight", value.getWeight()}
  };

    j["digi"] = nlohmann::json{value.getDigi().id()};

    j["sim"] = nlohmann::json{value.getSim().id()};


}
#endif


} // namespace extension

