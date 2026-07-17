// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4hep
#include "edm4hep/TrackCollection.h"
#include "extension/VertexCollection.h"

// C++
#include <array>
#include <vector>

#include <cxxabi.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>

namespace {

std::string demangleTypeName(const char* mangledName) {
  int status = 0;

  std::unique_ptr<char, void (*)(void*)> demangledName{abi::__cxa_demangle(mangledName, nullptr, nullptr, &status),
                                                       std::free};

  if (status == 0 && demangledName) {
    return demangledName.get();
  }

  return mangledName;
}

template <typename T>
std::string typeName() {
  return demangleTypeName(typeid(T).name());
}

template <typename T>
std::string objectTypeName(const T& object) {
  return demangleTypeName(typeid(object).name());
}

} // namespace

/** @class SimpleVertexFinder
 *
 *  Gaudi transformer that reconstructs vertex candidates from a set of
 *  reconstructed tracks.
 *
 *  @author Andrea De Vita
 */
struct SimpleVertexFinder final : k4FWCore::Transformer<extension::VertexCollection(const edm4hep::TrackCollection&)> {

  SimpleVertexFinder(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,

                    {KeyValues("InputFittedTracks", {"InputFittedTracks"})},
                    {KeyValues("OutputVerticesCandidates", {"OutputVerticesCandidates"})}) {}

  StatusCode initialize() override { return StatusCode::SUCCESS; }

  extension::VertexCollection operator()(const edm4hep::TrackCollection& fittedTracks) const override {

    std::cout << "\n========== SimpleVertexFinder debug ==========\n";

    std::cout << "fittedTracks declared type: " << typeName<decltype(fittedTracks)>() << '\n';

    std::cout << "fittedTracks underlying type: " << typeName<std::remove_cvref_t<decltype(fittedTracks)>>() << '\n';

    std::cout << "fittedTracks runtime type: " << objectTypeName(fittedTracks) << '\n';

    std::cout << "fittedTracks size: " << fittedTracks.size() << '\n';

    extension::VertexCollection verticesCandidates;

    std::cout << "verticesCandidates declared type: " << typeName<decltype(verticesCandidates)>() << '\n';

    std::cout << "verticesCandidates runtime type: " << objectTypeName(verticesCandidates) << '\n';

    auto vertexCandidate = verticesCandidates.create();

    std::cout << "vertexCandidate declared type: " << typeName<decltype(vertexCandidate)>() << '\n';

    std::cout << "vertexCandidate runtime type: " << objectTypeName(vertexCandidate) << '\n';

    std::size_t trackIndex = 0;

    for (const auto& track : fittedTracks) {
      std::cout << "\nTrack " << trackIndex << ":\n";

      std::cout << "  track declared type: " << typeName<decltype(track)>() << '\n';

      std::cout << "  track underlying type: " << typeName<std::remove_cvref_t<decltype(track)>>() << '\n';

      std::cout << "  track runtime type: " << objectTypeName(track) << '\n';

      vertexCandidate.addToTracks(track);

      std::cout << "  Added track to vertex candidate\n";

      ++trackIndex;
    }

    std::cout << "\nNumber of tracks added: " << trackIndex << '\n';

    std::cout << "Number of vertex candidates: " << verticesCandidates.size() << '\n';

    std::cout << "Return type: " << typeName<extension::VertexCollection>() << '\n';

    std::cout << "==============================================\n\n";

    return verticesCandidates;
  }

private:
};

DECLARE_COMPONENT(SimpleVertexFinder)