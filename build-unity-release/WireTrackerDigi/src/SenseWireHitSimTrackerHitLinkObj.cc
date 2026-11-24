// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/SenseWireHitSimTrackerHitLinkObj.h"
#include "extension/SenseWireHit.h"
#include "edm4hep/SimTrackerHit.h"
#include "podio/ObjectID.h"

namespace extension {

/// \cond PODIO_INTERNAL

SenseWireHitSimTrackerHitLinkObj::SenseWireHitSimTrackerHitLinkObj() :
  data(),
  m_from(nullptr),
  m_to(nullptr)
{  }

SenseWireHitSimTrackerHitLinkObj::SenseWireHitSimTrackerHitLinkObj(const podio::ObjectID& id_, const SenseWireHitSimTrackerHitLinkData& data_) :
  id(id_), data(data_)
{  }

SenseWireHitSimTrackerHitLinkObj::SenseWireHitSimTrackerHitLinkObj(const SenseWireHitSimTrackerHitLinkObj& other) :
  data(other.data),
  m_from(nullptr),
  m_to(nullptr)
{
  if (other.m_from) {
    m_from = std::make_unique<extension::SenseWireHit>(*(other.m_from));
  }
  if (other.m_to) {
    m_to = std::make_unique<edm4hep::SimTrackerHit>(*(other.m_to));
  }
}

SenseWireHitSimTrackerHitLinkObj::~SenseWireHitSimTrackerHitLinkObj() = default;
/// \endcond
} // namespace extension

