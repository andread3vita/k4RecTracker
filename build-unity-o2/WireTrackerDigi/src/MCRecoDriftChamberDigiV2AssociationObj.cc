// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/MCRecoDriftChamberDigiV2AssociationObj.h"
#include "extension/DriftChamberDigiV2.h"
#include "edm4hep/SimTrackerHit.h"
#include "podio/ObjectID.h"

namespace extension {

/// \cond PODIO_INTERNAL

MCRecoDriftChamberDigiV2AssociationObj::MCRecoDriftChamberDigiV2AssociationObj() :
  data(),
  m_digi(nullptr),
  m_sim(nullptr)
{  }

MCRecoDriftChamberDigiV2AssociationObj::MCRecoDriftChamberDigiV2AssociationObj(const podio::ObjectID& id_, const MCRecoDriftChamberDigiV2AssociationData& data_) :
  id(id_), data(data_)
{  }

MCRecoDriftChamberDigiV2AssociationObj::MCRecoDriftChamberDigiV2AssociationObj(const MCRecoDriftChamberDigiV2AssociationObj& other) :
  data(other.data),
  m_digi(nullptr),
  m_sim(nullptr)
{
  if (other.m_digi) {
    m_digi = std::make_unique<extension::DriftChamberDigiV2>(*(other.m_digi));
  }
  if (other.m_sim) {
    m_sim = std::make_unique<edm4hep::SimTrackerHit>(*(other.m_sim));
  }
}

MCRecoDriftChamberDigiV2AssociationObj::~MCRecoDriftChamberDigiV2AssociationObj() = default;
/// \endcond
} // namespace extension

