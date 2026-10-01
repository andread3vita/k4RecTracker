/*
 * Copyright (c) 2014-2026 Key4hep-Project.
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

#include "edm4hep/MCParticle.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "k4FWCore/Transformer.h"
#include "podio/ObjectID.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using DigiToSimHitLinkCollections = std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>;

struct ParticleHitCount {
  edm4hep::MCParticle particle;
  std::size_t hitCount{};
};

} // namespace

/** @struct TracksToMCParticlesLinker
 *
 * @brief Associate each reconstructed track with the MCParticle responsible for most of its hits.
 *
 * The TrackerHit-to-SimTrackerHit links are used to obtain the MCParticle of each hit on an input track. Each digitized
 * hit has one associated SimTrackerHit and therefore casts one vote for that SimTrackerHit's MCParticle. The MCParticle
 * with the largest number of votes is linked to the track. The link weight is the fraction of all hits on the track
 * that voted for the selected particle. Tracks without any truth-associated hits do not produce a link.
 *
 * In case of equal hit counts, the first particle encountered while visiting the track hits and input link collections
 * is selected, which makes the result stable for a fixed input ordering.
 *
 * @author Andrea De Vita
 * @date   2026-02
 */
struct TracksToMCParticlesLinker final : k4FWCore::Transformer<edm4hep::TrackMCParticleLinkCollection(
                                             const edm4hep::TrackCollection&, const DigiToSimHitLinkCollections&)> {
  TracksToMCParticlesLinker(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(
            name, svcLoc,
            {KeyValues("TrackCollection", {"InputTracks"}), KeyValues("DigiToSimHitsLinks", {"DigiToSimHitsLinks"})},
            {KeyValues("LinksTracksMCParticles", {"TracksMCParticlesLinks"})}) {}

  edm4hep::TrackMCParticleLinkCollection
  operator()(const edm4hep::TrackCollection& inputTracks,
             const DigiToSimHitLinkCollections& digiToSimHitLinkCollections) const override {
    edm4hep::TrackMCParticleLinkCollection outputLinks;

    // Index the truth particles contributing to every digitized hit once per event. There may be several link
    // collections because the hits on a track can originate from different tracker subdetectors.
    std::unordered_map<podio::ObjectID, edm4hep::MCParticle> particleByHit;
    for (const auto* linkCollection : digiToSimHitLinkCollections) {
      if (linkCollection == nullptr) {
        continue;
      }

      for (const auto& hitLink : *linkCollection) {
        const auto digiHit = hitLink.getFrom();
        const auto simHit = hitLink.getTo();
        if (!digiHit.isAvailable() || !simHit.isAvailable()) {
          continue;
        }

        const auto particle = simHit.getParticle();
        if (!particle.isAvailable()) {
          continue;
        }

        particleByHit.emplace(digiHit.getObjectID(), particle);
      }
    }

    for (const auto& track : inputTracks) {
      const auto trackHits = track.getTrackerHits();
      if (trackHits.empty()) {
        continue;
      }

      std::vector<ParticleHitCount> particleHitCounts;
      for (const auto& hit : trackHits) {
        const auto hitParticle = particleByHit.find(hit.getObjectID());
        if (hitParticle == particleByHit.end()) {
          continue;
        }

        const auto particleID = hitParticle->second.getObjectID();
        const auto count =
            std::find_if(particleHitCounts.begin(), particleHitCounts.end(),
                         [&particleID](const auto& entry) { return entry.particle.getObjectID() == particleID; });
        if (count == particleHitCounts.end()) {
          particleHitCounts.push_back({hitParticle->second, 1});
        } else {
          ++count->hitCount;
        }
      }

      if (particleHitCounts.empty()) {
        continue;
      }

      const auto bestParticle =
          std::max_element(particleHitCounts.begin(), particleHitCounts.end(),
                           [](const auto& lhs, const auto& rhs) { return lhs.hitCount < rhs.hitCount; });

      const auto weight = static_cast<float>(bestParticle->hitCount) / static_cast<float>(trackHits.size());
      auto outputLink = edm4hep::MutableTrackMCParticleLink();
      outputLink.setFrom(track);
      outputLink.setTo(bestParticle->particle);
      outputLink.setWeight(weight);
      outputLinks.push_back(outputLink);
    }

    return outputLinks;
  }
};

DECLARE_COMPONENT(TracksToMCParticlesLinker)
