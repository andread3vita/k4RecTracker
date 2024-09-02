// #pragma once

// // GAUDI
// #include "Gaudi/Property.h"
// // #include "GaudiAlg/GaudiAlgorithm.h"
// #include "Gaudi/Algorithm.h"
// // K4FWCORE
// #include "k4FWCore/DataHandle.h"

// // EDM4HEP
// #include "edm4hep/SimTrackerHitCollection.h"
// #if __has_include("edm4hep/TrackerHit3DCollection.h")
// #include "edm4hep/TrackerHit3DCollection.h"
// #else
// #include "edm4hep/TrackerHit3D.h"
// #include "extension/TrackerHit3D.h"
// namespace edm4hep {
//   using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
// }  // namespace edm4hep
// #endif

// #include "extension/TrackCollection.h"
// #include "extension/TrackCollection.h"
// #include "extension/TrackerHit3DCollection.h"
// // EDM4HEP extension
// #include "extension/DriftChamberDigiCollection.h"
// #include "extension/DriftChamberDigiLocalCollection.h"
// #include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
// #include "podio/UserDataCollection.h"

// #include "onnxruntime_cxx_api.h"
// #include <memory> 
// #include <vector> 
// #include <string>   

// // GENFIT
// //#include "WireMeasurement.h"

// /** @class GGTF_tracking_dbscan_alg
//  *
//  *  
//  *  
//  *  @author Maria Dolores Garcia, Brieuc Francois, Andrea De Vita
//  *  @date   2024-07
//  *
//  */

// class GGTF_tracking_dbscan_alg : public Gaudi::Algorithm {
// public:
//   explicit GGTF_tracking_dbscan_alg(const std::string&, ISvcLocator*);
//   virtual ~GGTF_tracking_dbscan_alg();
//   /**  Initialize.
//    *   @return status code
//    */
//   virtual StatusCode initialize() final;
//   /**  Execute.
//    *   @return status code
//    */
//   virtual StatusCode execute(const EventContext&) const final;
//   /**  Finalize.
//    *   @return status code
//    */
//   virtual StatusCode finalize() final;

// private:
//   /// Pointer to the ONNX enviroment
//   std::unique_ptr<Ort::Env> fEnv;
//   /// Pointer to the ONNX inference session
//   std::unique_ptr<Ort::Session> fSession;
//   /// ONNX settings
//   Ort::SessionOptions fSessionOptions;
//   /// ONNX memory info
//   const OrtMemoryInfo* fInfo;
//   struct MemoryInfo;
//   /// the input names represent the names given to the model
//   /// when defining  the model's architecture (if applicable)
//   /// they can also be retrieved from model.summary()
//   std::vector<const char*> fInames;
//   std::vector<const char*> fOnames;

//   std::string modelPath={};
//   // Input tracker hit collection name
//   mutable DataHandle<extension::DriftChamberDigiCollection> m_input_hits_CDC{"inputHits_CDC", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXIB{"inputHits_VTXIB", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXD{"inputHits_VTXD", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXOB{"inputHits_VTXOB", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<extension::MCRecoDriftChamberDigiAssociation> m_input_Association_CDC{"inputAssociation_CDC_sim", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_CDC_sim{"inputHits_CDC_sim", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXIB_sim{"inputHits_VTXIB_sim", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXD_sim{"inputHits_VTXD_sim", Gaudi::DataHandle::Reader, this};
//   mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_hits_VTXOB_sim{"inputHits_VTXOB_sim", Gaudi::DataHandle::Reader, this};

//   // Output track collection name
//   mutable DataHandle<extension::TrackCollection> m_output_tracks{"outputTracks", Gaudi::DataHandle::Writer, this};
//   mutable DataHandle<extension::TrackerHit3DCollection> m_output_hits{"outputHits", Gaudi::DataHandle::Writer, this};
//   mutable DataHandle<podio::UserDataCollection <double>> m_clustering_space{"clustering_space", Gaudi::DataHandle::Writer, this};
//   mutable DataHandle<podio::UserDataCollection <double>> m_tracks_space{"clustering_space_tracks", Gaudi::DataHandle::Writer, this};
// };
