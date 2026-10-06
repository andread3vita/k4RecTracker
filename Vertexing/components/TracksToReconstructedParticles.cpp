/*
 * Copyright (c) 2020-2026 Key4hep-Project.
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

#include <Gaudi/Property.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <k4FWCore/Transformer.h>

#include <cmath>
#include <limits>
#include <string>

struct TracksToReconstructedParticles final
    : k4FWCore::Transformer<edm4hep::ReconstructedParticleCollection(const edm4hep::TrackCollection&)> {
  TracksToReconstructedParticles(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc, {KeyValue("InputTracks", {"InputTracks"})},
                    {KeyValue("OutputReconstructedParticles", {"ReconstructedParticles"})}) {}

  edm4hep::ReconstructedParticleCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::ReconstructedParticleCollection particles;
    std::size_t tracksWithoutUsableIPState = 0;

    // For B in tesla and omega in mm^-1, this gives momentum in GeV.
    constexpr double curvatureToMomentum = 2.99792458e-4;

    for (const auto& track : tracks) {
      auto particle = particles.create();
      particle.addToTracks(track);

      edm4hep::TrackState stateAtIP{};
      bool hasUsableIPState = false;
      float charge = 1.0f;
      for (const auto& state : track.getTrackStates()) {
        if (state.location == edm4hep::TrackState::AtIP) {
          stateAtIP = state;
          hasUsableIPState = std::abs(state.omega) > std::numeric_limits<float>::epsilon();
          break;
        }
      }

      if (hasUsableIPState) {
        const double magneticField = m_magneticField.value();
        charge = stateAtIP.omega * magneticField < 0.0 ? -1.0f : 1.0f;

        const double pt =
            curvatureToMomentum * std::abs(magneticField) / std::abs(static_cast<double>(stateAtIP.omega));
        const double px = pt * std::cos(stateAtIP.phi);
        const double py = pt * std::sin(stateAtIP.phi);
        const double pz = pt * stateAtIP.tanLambda;
        const double momentumMagnitude = std::sqrt(px * px + py * py + pz * pz);

        particle.setMomentum(
            {static_cast<float>(px), static_cast<float>(py), static_cast<float>(pz)});

        // Tracking does not determine a particle species. Use a massless
        // four-vector and do not assign a PDG hypothesis. This supplies the
        // momentum needed by secondary vertexing without claiming PID.
        particle.setEnergy(static_cast<float>(momentumMagnitude));
      } else {
        ++tracksWithoutUsableIPState;
      }

      particle.setCharge(charge);
      particle.setMass(0.0f);
      particle.setPDG(0);
    }

    if (tracksWithoutUsableIPState != 0) {
      warning() << tracksWithoutUsableIPState
                << " track(s) have no usable AtIP curvature; their reconstructed-particle momentum remains zero."
                << endmsg;
    }

    debug() << "Created " << particles.size() << " reconstructed particles from " << tracks.size() << " tracks"
            << endmsg;
    return particles;
  }

private:
  Gaudi::Property<double> m_magneticField{this, "MagneticField", 2.0, "Magnetic field along z [T]"};
};

DECLARE_COMPONENT(TracksToReconstructedParticles)
