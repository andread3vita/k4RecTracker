// AUTOMATICALLY GENERATED FILE - DO NOT EDIT

// datamodel specific includes
#include "extension/Track.h"
#include "extension/MutableTrack.h"
#include "extension/TrackObj.h"
#include "extension/TrackData.h"


#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
#include "nlohmann/json.hpp"
#endif

#include <ostream>

namespace extension {



Track::Track(const std::int32_t type, const float chi2, const std::int32_t ndf, const float dEdx, const float dEdxError, const float radiusOfInnermostHit) {
  m_obj->data.type = type;
  m_obj->data.chi2 = chi2;
  m_obj->data.ndf = ndf;
  m_obj->data.dEdx = dEdx;
  m_obj->data.dEdxError = dEdxError;
  m_obj->data.radiusOfInnermostHit = radiusOfInnermostHit;
}

Track& Track::operator=(Track other) & {
  swap(*this, other);
  return *this;
}

MutableTrack Track::clone(bool cloneRelations) const {
  auto tmp = new TrackObj(podio::ObjectID{}, m_obj->data);
  tmp->m_trackerHits = new std::vector<extension::TrackerHit>();
  tmp->m_tracks = new std::vector<extension::Track>();
  tmp->m_subdetectorHitNumbers = new std::vector<std::int32_t>();
  tmp->m_trackStates = new std::vector<edm4hep::TrackState>();
  tmp->m_dxQuantities = new std::vector<edm4hep::Quantity>();
  if (cloneRelations) {
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_trackerHits->reserve(m_obj->m_trackerHits->size());
    for (size_t i = m_obj->data.trackerHits_begin; i < m_obj->data.trackerHits_end; i++) {
      tmp->m_trackerHits->emplace_back((*m_obj->m_trackerHits)[i]);
    }
    tmp->data.trackerHits_begin = 0;
    tmp->data.trackerHits_end = tmp->m_trackerHits->size();
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_tracks->reserve(m_obj->m_tracks->size());
    for (size_t i = m_obj->data.tracks_begin; i < m_obj->data.tracks_end; i++) {
      tmp->m_tracks->emplace_back((*m_obj->m_tracks)[i]);
    }
    tmp->data.tracks_begin = 0;
    tmp->data.tracks_end = tmp->m_tracks->size();
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_subdetectorHitNumbers->reserve(m_obj->m_subdetectorHitNumbers->size());
    for (size_t i = m_obj->data.subdetectorHitNumbers_begin; i < m_obj->data.subdetectorHitNumbers_end; i++) {
      tmp->m_subdetectorHitNumbers->emplace_back((*m_obj->m_subdetectorHitNumbers)[i]);
    }
    tmp->data.subdetectorHitNumbers_begin = 0;
    tmp->data.subdetectorHitNumbers_end = tmp->m_subdetectorHitNumbers->size();
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_trackStates->reserve(m_obj->m_trackStates->size());
    for (size_t i = m_obj->data.trackStates_begin; i < m_obj->data.trackStates_end; i++) {
      tmp->m_trackStates->emplace_back((*m_obj->m_trackStates)[i]);
    }
    tmp->data.trackStates_begin = 0;
    tmp->data.trackStates_end = tmp->m_trackStates->size();
    // If the current object has been read from a file, then the object may only have a slice of the relation vector
    // so this slice has to be copied in case we want to modify it
    tmp->m_dxQuantities->reserve(m_obj->m_dxQuantities->size());
    for (size_t i = m_obj->data.dxQuantities_begin; i < m_obj->data.dxQuantities_end; i++) {
      tmp->m_dxQuantities->emplace_back((*m_obj->m_dxQuantities)[i]);
    }
    tmp->data.dxQuantities_begin = 0;
    tmp->data.dxQuantities_end = tmp->m_dxQuantities->size();
  }
  else {
    tmp->data.trackerHits_begin = 0;
    tmp->data.trackerHits_end = 0;
    tmp->data.tracks_begin = 0;
    tmp->data.tracks_end = 0;
    tmp->data.subdetectorHitNumbers_begin = 0;
    tmp->data.subdetectorHitNumbers_end = 0;
    tmp->data.trackStates_begin = 0;
    tmp->data.trackStates_end = 0;
    tmp->data.dxQuantities_begin = 0;
    tmp->data.dxQuantities_end = 0;
  }
  return MutableTrack(podio::utils::MaybeSharedPtr(tmp, podio::utils::MarkOwned));
}

Track::Track(podio::utils::MaybeSharedPtr<TrackObj> obj) : m_obj(std::move(obj)) {}

Track::Track(const MutableTrack& other): Track(other.m_obj) {}

Track::Track(TrackObj* obj) : m_obj(podio::utils::MaybeSharedPtr<TrackObj>(obj)) {}

Track Track::makeEmpty() {
  return {nullptr};
}

std::int32_t Track::getType() const { return m_obj->data.type; }
float Track::getChi2() const { return m_obj->data.chi2; }
std::int32_t Track::getNdf() const { return m_obj->data.ndf; }
float Track::getDEdx() const { return m_obj->data.dEdx; }
float Track::getDEdxError() const { return m_obj->data.dEdxError; }
float Track::getRadiusOfInnermostHit() const { return m_obj->data.radiusOfInnermostHit; }



std::vector<extension::TrackerHit>::const_iterator Track::trackerHits_begin() const {
  auto ret_value = m_obj->m_trackerHits->begin();
  std::advance(ret_value, m_obj->data.trackerHits_begin);
  return ret_value;
}

std::vector<extension::TrackerHit>::const_iterator Track::trackerHits_end() const {
  auto ret_value = m_obj->m_trackerHits->begin();
  std::advance(ret_value, m_obj->data.trackerHits_end);
  return ret_value;
}

std::size_t Track::trackerHits_size() const {
  return m_obj->data.trackerHits_end - m_obj->data.trackerHits_begin;
}

extension::TrackerHit Track::getTrackerHits(std::size_t index) const {
  if (trackerHits_size() > index) {
    return m_obj->m_trackerHits->at(m_obj->data.trackerHits_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<extension::TrackerHit> Track::getTrackerHits() const {
  auto begin = m_obj->m_trackerHits->begin();
  std::advance(begin, m_obj->data.trackerHits_begin);
  auto end = m_obj->m_trackerHits->begin();
  std::advance(end, m_obj->data.trackerHits_end);
  return {begin, end};
}


std::vector<extension::Track>::const_iterator Track::tracks_begin() const {
  auto ret_value = m_obj->m_tracks->begin();
  std::advance(ret_value, m_obj->data.tracks_begin);
  return ret_value;
}

std::vector<extension::Track>::const_iterator Track::tracks_end() const {
  auto ret_value = m_obj->m_tracks->begin();
  std::advance(ret_value, m_obj->data.tracks_end);
  return ret_value;
}

std::size_t Track::tracks_size() const {
  return m_obj->data.tracks_end - m_obj->data.tracks_begin;
}

extension::Track Track::getTracks(std::size_t index) const {
  if (tracks_size() > index) {
    return m_obj->m_tracks->at(m_obj->data.tracks_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<extension::Track> Track::getTracks() const {
  auto begin = m_obj->m_tracks->begin();
  std::advance(begin, m_obj->data.tracks_begin);
  auto end = m_obj->m_tracks->begin();
  std::advance(end, m_obj->data.tracks_end);
  return {begin, end};
}


std::vector<std::int32_t>::const_iterator Track::subdetectorHitNumbers_begin() const {
  auto ret_value = m_obj->m_subdetectorHitNumbers->begin();
  std::advance(ret_value, m_obj->data.subdetectorHitNumbers_begin);
  return ret_value;
}

std::vector<std::int32_t>::const_iterator Track::subdetectorHitNumbers_end() const {
  auto ret_value = m_obj->m_subdetectorHitNumbers->begin();
  std::advance(ret_value, m_obj->data.subdetectorHitNumbers_end);
  return ret_value;
}

std::size_t Track::subdetectorHitNumbers_size() const {
  return m_obj->data.subdetectorHitNumbers_end - m_obj->data.subdetectorHitNumbers_begin;
}

std::int32_t Track::getSubdetectorHitNumbers(std::size_t index) const {
  if (subdetectorHitNumbers_size() > index) {
    return m_obj->m_subdetectorHitNumbers->at(m_obj->data.subdetectorHitNumbers_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<std::int32_t> Track::getSubdetectorHitNumbers() const {
  auto begin = m_obj->m_subdetectorHitNumbers->begin();
  std::advance(begin, m_obj->data.subdetectorHitNumbers_begin);
  auto end = m_obj->m_subdetectorHitNumbers->begin();
  std::advance(end, m_obj->data.subdetectorHitNumbers_end);
  return {begin, end};
}


std::vector<edm4hep::TrackState>::const_iterator Track::trackStates_begin() const {
  auto ret_value = m_obj->m_trackStates->begin();
  std::advance(ret_value, m_obj->data.trackStates_begin);
  return ret_value;
}

std::vector<edm4hep::TrackState>::const_iterator Track::trackStates_end() const {
  auto ret_value = m_obj->m_trackStates->begin();
  std::advance(ret_value, m_obj->data.trackStates_end);
  return ret_value;
}

std::size_t Track::trackStates_size() const {
  return m_obj->data.trackStates_end - m_obj->data.trackStates_begin;
}

edm4hep::TrackState Track::getTrackStates(std::size_t index) const {
  if (trackStates_size() > index) {
    return m_obj->m_trackStates->at(m_obj->data.trackStates_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<edm4hep::TrackState> Track::getTrackStates() const {
  auto begin = m_obj->m_trackStates->begin();
  std::advance(begin, m_obj->data.trackStates_begin);
  auto end = m_obj->m_trackStates->begin();
  std::advance(end, m_obj->data.trackStates_end);
  return {begin, end};
}


std::vector<edm4hep::Quantity>::const_iterator Track::dxQuantities_begin() const {
  auto ret_value = m_obj->m_dxQuantities->begin();
  std::advance(ret_value, m_obj->data.dxQuantities_begin);
  return ret_value;
}

std::vector<edm4hep::Quantity>::const_iterator Track::dxQuantities_end() const {
  auto ret_value = m_obj->m_dxQuantities->begin();
  std::advance(ret_value, m_obj->data.dxQuantities_end);
  return ret_value;
}

std::size_t Track::dxQuantities_size() const {
  return m_obj->data.dxQuantities_end - m_obj->data.dxQuantities_begin;
}

edm4hep::Quantity Track::getDxQuantities(std::size_t index) const {
  if (dxQuantities_size() > index) {
    return m_obj->m_dxQuantities->at(m_obj->data.dxQuantities_begin + index);
  }
  throw std::out_of_range("index out of bounds for existing references");
}

podio::RelationRange<edm4hep::Quantity> Track::getDxQuantities() const {
  auto begin = m_obj->m_dxQuantities->begin();
  std::advance(begin, m_obj->data.dxQuantities_begin);
  auto end = m_obj->m_dxQuantities->begin();
  std::advance(end, m_obj->data.dxQuantities_end);
  return {begin, end};
}





bool Track::isAvailable() const {
  return m_obj;
}

const podio::ObjectID Track::getObjectID() const {
  if (m_obj) {
    return m_obj->id;
  }
  return podio::ObjectID{};
}

bool Track::operator==(const MutableTrack& other) const {
  return m_obj == other.m_obj;
}

std::ostream& operator<<(std::ostream& o, const Track& value) {
  if (!value.isAvailable()) {
    return o << "[not available]";
  }
  o << " id: " << value.id() << '\n';
  o << " type : " << value.getType() << '\n';
  o << " chi2 : " << value.getChi2() << '\n';
  o << " ndf : " << value.getNdf() << '\n';
  o << " dEdx : " << value.getDEdx() << '\n';
  o << " dEdxError : " << value.getDEdxError() << '\n';
  o << " radiusOfInnermostHit : " << value.getRadiusOfInnermostHit() << '\n';


  o << " trackerHits : ";
  for (unsigned i = 0; i < value.trackerHits_size(); ++i) {
    o << value.getTrackerHits(i) << " ";
  }
  o << '\n';
  o << " tracks : ";
  for (unsigned i = 0; i < value.tracks_size(); ++i) {
    o << value.getTracks(i).id() << " ";
  }
  o << '\n';
  o << " subdetectorHitNumbers : ";
  for (unsigned i = 0; i < value.subdetectorHitNumbers_size(); ++i) {
    o << value.getSubdetectorHitNumbers(i) << " ";
  }
  o << '\n';
  o << " trackStates : ";
  for (unsigned i = 0; i < value.trackStates_size(); ++i) {
    o << value.getTrackStates(i) << " ";
  }
  o << '\n';
  o << " dxQuantities : ";
  for (unsigned i = 0; i < value.dxQuantities_size(); ++i) {
    o << value.getDxQuantities(i) << " ";
  }
  o << '\n';

  return o;
}

#if defined(PODIO_JSON_OUTPUT) && !defined(__CLING__)
void to_json(nlohmann::json& j, const Track& value) {
  j = nlohmann::json{
    {"type", value.getType()}
    ,{"chi2", value.getChi2()}
    ,{"ndf", value.getNdf()}
    ,{"dEdx", value.getDEdx()}
    ,{"dEdxError", value.getDEdxError()}
    ,{"radiusOfInnermostHit", value.getRadiusOfInnermostHit()}
    ,{"subdetectorHitNumbers", value.getSubdetectorHitNumbers()}
    ,{"trackStates", value.getTrackStates()}
    ,{"dxQuantities", value.getDxQuantities()}
  };


  j["trackerHits"] = nlohmann::json::array();
  for (const auto& v : value.getTrackerHits()) {
    j["trackerHits"].emplace_back(v.id());
  }

  j["tracks"] = nlohmann::json::array();
  for (const auto& v : value.getTracks()) {
    j["tracks"].emplace_back(v.id());
  }

}
#endif


} // namespace extension


podio::detail::OrderKey podio::detail::getOrderKey(const extension::Track& obj) {
  return podio::detail::OrderKey{obj.m_obj.get()};
}
