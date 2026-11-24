// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/DriftChamberDigiV2Obj.h"
#include "podio/ObjectID.h"

namespace extension {

/// \cond PODIO_INTERNAL

DriftChamberDigiV2Obj::DriftChamberDigiV2Obj() :
  data(),
  m_nElectrons(new std::vector<std::uint16_t>())
{  }

DriftChamberDigiV2Obj::DriftChamberDigiV2Obj(const podio::ObjectID& id_, const DriftChamberDigiV2Data& data_) :
  id(id_), data(data_)
{  }

DriftChamberDigiV2Obj::DriftChamberDigiV2Obj(const DriftChamberDigiV2Obj& other) :
  data(other.data),
  m_nElectrons(new std::vector<std::uint16_t>(*(other.m_nElectrons)))
{
}

DriftChamberDigiV2Obj::~DriftChamberDigiV2Obj() {
  if (id.index == podio::ObjectID::untracked) {
    delete m_nElectrons;
  }
}
/// \endcond
} // namespace extension

