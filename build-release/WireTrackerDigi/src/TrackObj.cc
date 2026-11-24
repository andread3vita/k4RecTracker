// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

#include "extension/TrackObj.h"
#include "extension/Track.h"
#include "podio/ObjectID.h"

namespace extension {

/// \cond PODIO_INTERNAL

TrackObj::TrackObj() :
  data(),
  m_trackerHits(new std::vector<extension::TrackerHit>()),
  m_tracks(new std::vector<extension::Track>()),
  m_subdetectorHitNumbers(new std::vector<std::int32_t>()),
  m_trackStates(new std::vector<edm4hep::TrackState>()),
  m_dxQuantities(new std::vector<edm4hep::Quantity>())
{  }

TrackObj::TrackObj(const podio::ObjectID& id_, const TrackData& data_) :
  id(id_), data(data_)
{  }

TrackObj::TrackObj(const TrackObj& other) :
  data(other.data),
  m_trackerHits(new std::vector<extension::TrackerHit>(*(other.m_trackerHits))),
  m_tracks(new std::vector<extension::Track>(*(other.m_tracks))),
  m_subdetectorHitNumbers(new std::vector<std::int32_t>(*(other.m_subdetectorHitNumbers))),
  m_trackStates(new std::vector<edm4hep::TrackState>(*(other.m_trackStates))),
  m_dxQuantities(new std::vector<edm4hep::Quantity>(*(other.m_dxQuantities)))
{
}

TrackObj::~TrackObj() {
  if (id.index == podio::ObjectID::untracked) {
    delete m_trackerHits;
    delete m_tracks;
    delete m_subdetectorHitNumbers;
    delete m_trackStates;
    delete m_dxQuantities;
  }
}
/// \endcond
} // namespace extension

