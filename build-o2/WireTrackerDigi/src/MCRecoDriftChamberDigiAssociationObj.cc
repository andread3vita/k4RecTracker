// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/MCRecoDriftChamberDigiAssociationObj.h"
#include "extension/DriftChamberDigi.h"
#include "edm4hep/SimTrackerHit.h"
#include "podio/ObjectID.h"

namespace extension {

/// \cond PODIO_INTERNAL

MCRecoDriftChamberDigiAssociationObj::MCRecoDriftChamberDigiAssociationObj() :
  data(),
  m_digi(nullptr),
  m_sim(nullptr)
{  }

MCRecoDriftChamberDigiAssociationObj::MCRecoDriftChamberDigiAssociationObj(const podio::ObjectID& id_, const MCRecoDriftChamberDigiAssociationData& data_) :
  id(id_), data(data_)
{  }

MCRecoDriftChamberDigiAssociationObj::MCRecoDriftChamberDigiAssociationObj(const MCRecoDriftChamberDigiAssociationObj& other) :
  data(other.data),
  m_digi(nullptr),
  m_sim(nullptr)
{
  if (other.m_digi) {
    m_digi = std::make_unique<extension::DriftChamberDigi>(*(other.m_digi));
  }
  if (other.m_sim) {
    m_sim = std::make_unique<edm4hep::SimTrackerHit>(*(other.m_sim));
  }
}

MCRecoDriftChamberDigiAssociationObj::~MCRecoDriftChamberDigiAssociationObj() = default;
/// \endcond
} // namespace extension

