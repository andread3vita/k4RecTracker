#pragma once

// GAUDI
#include "Gaudi/Property.h"
// #include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Algorithm.h"
// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHit3D.h"
#include "extension/TrackerHit3D.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
}  // namespace edm4hep
#endif

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "extension/TrackCollection.h"
#include "podio/UserDataCollection.h"

// CPP
#include <algorithm>
#include <map>

/** @class Evaluation_efficiency
 *
 *  
 *  
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2024-06
 *
 */

class Evaluation_efficiency : public Gaudi::Algorithm {
public:
  explicit Evaluation_efficiency(const std::string&, ISvcLocator*);
  virtual ~Evaluation_efficiency();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:

  
  // Input tracker hit collection name
  mutable DataHandle<extension::TrackCollection> m_input_tracks{"inputTracks", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::MCParticleCollection> m_input_MCparticles{"inputMCparticles", Gaudi::DataHandle::Reader,this};

  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_CDC_sim{"inputHits_CDC_sim", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXIB_sim{"inputHits_VTXIB_sim", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXD_sim{"inputHits_VTXD_sim", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXOB_sim{"inputHits_VTXOB_sim", Gaudi::DataHandle::Reader, this};

  // Output track collection name
  // mutable DataHandle<podio::UserDataCollection<std::map<int, double>>> m_output_trackEff{"outputTracksEff", Gaudi::DataHandle::Writer, this};
  // mutable DataHandle<podio::UserDataCollection<std::map<int, double>>> m_output_MCpur{"outputMCPur", Gaudi::DataHandle::Writer, this};

  mutable DataHandle <podio::UserDataCollection <int>> m_output_tracking{"outputTrackingEff", Gaudi::DataHandle::Writer, this};
};
