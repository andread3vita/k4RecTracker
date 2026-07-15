// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackCollection.h"
#include "extension/VertexCollection.h"

// Vertexing kernel
#include "LinearizedHelixVertexFitter.h"

// C++
#include <array>
#include <vector>

/** @class SimpleVertexFinder
 *
 *  Gaudi transformer that reconstructs vertex candidates from a set of
 *  reconstructed tracks.
 *
 *  @author Andrea De Vita
 */
struct SimpleVertexFinder final
    : k4FWCore::Transformer<extension::VertexCollection(const edm4hep::TrackCollection&)> {

  SimpleVertexFinder(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,

                    {KeyValues("InputFittedTracks", {"InputFittedTracks"})},
                    {KeyValues("OutputVerticesCandidates", {"OutputVerticesCandidates"})}) {}

  StatusCode initialize() override { return StatusCode::SUCCESS; }

  extension::VertexCollection operator()(const edm4hep::TrackCollection& fittedTracks) const override {

    extension::VertexCollection VerticesCandidates;

    // Create a new vertex candidate and add all the track to it
    auto vertexCandidate = VerticesCandidates.create();
    for (const auto& track : fittedTracks) {

      vertexCandidate.addToTracks(track);
    }

    return VerticesCandidates;
  }

private:
};

DECLARE_COMPONENT(SimpleVertexFinder)