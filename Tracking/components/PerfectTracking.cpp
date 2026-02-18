/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Standard Library
#include <algorithm>
#include <cmath>
#include <cstdlib>    // For getenv
#include <filesystem> // For std::filesystem::path
#include <fstream>    // For std::ifstream
#include <iostream>
#include <iterator> // For std::istreambuf_iterator
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// ONNX & Torch
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"
#include <ATen/ATen.h>
#include <torch/torch.h>

// ROOT
#include "TVector3.h"

// === Gaudi Framework ===
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

// === k4FWCore / k4Interface ===
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// === EDM4HEP & PODIO ===
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackCollection.h"

// === Project-specific ===
#include "utils.h"

/** @struct PerfectTracking
 *
 *
 */

struct PerfectTracking final : k4FWCore::MultiTransformer< std::tuple<edm4hep::TrackCollection> (
                                
                                const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
                                const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
                                const edm4hep::MCParticleCollection&)> {

  PerfectTracking(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {

                            KeyValues("InputPlanarHitCollections", {"InputPlanarHitCollections"}),
                            KeyValues("InputWireHitCollections", {"InputWireHitCollections"}),
                            KeyValues("InputMCParticles", {"InputMCParticles"})

                         },
                         {

                            KeyValues("OutputPerfectTracks", {"OutputPerfectTracks"})

                         }) {}

  StatusCode initialize() override {

    

    return StatusCode::SUCCESS;

  }

  std::tuple<edm4hep::TrackCollection> 
    operator()(const std::vector<   const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& planarHitLinks,
                                    const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& wireHitLinks,
                                    const edm4hep::MCParticleCollection& mcParticles) const override {
    
    //////////////////////////////////
    ////////// PERFECT TRACKING //////
    //////////////////////////////////

    // Create a new TrackCollection for storing the output tracks
    edm4hep::TrackCollection outputTracks;
    // Loop over MCParticles to create perfect tracks
    for (const auto& mcParticle : mcParticles) {

        // if ( mcParticle.getGeneratorStatus() != 1) {
        //     continue; // Skip non-final-state particles
        // }

        auto edm4hep_track = outputTracks.create();
        auto mcParticleObjectId = mcParticle.getObjectID();

        // Loop over planar hit link collections
        for (const auto& planarHitLinkCollection : planarHitLinks) {
            for (const auto& hitLink : *planarHitLinkCollection) {

                auto simHit = hitLink.getTo();
                auto digiHit = hitLink.getFrom();

                if (simHit.getParticle().getObjectID() == mcParticleObjectId) {
                    edm4hep_track.addToTrackerHits(digiHit);
                }
            }
        }

        // Loop over wire hit link collections
        for (const auto& wireHitLinkCollection : wireHitLinks) {
            for (const auto& hitLink : *wireHitLinkCollection) {    
                
                auto simHit = hitLink.getTo();
                auto digiHit = hitLink.getFrom();

                if (simHit.getParticle().getObjectID() == mcParticleObjectId) {
                    edm4hep_track.addToTrackerHits(digiHit);
                }
            }   
        }

        edm4hep_track.setType(1);

    }

    // Return the output collections as a tuple
    return std::make_tuple(std::move(outputTracks));
  }

  StatusCode finalize() override {

    return StatusCode::SUCCESS;
  }

private:

  

};

DECLARE_COMPONENT(PerfectTracking)