set(PackageName Tracking)

project(${PackageName})

find_package(GenFit REQUIRED)
find_package(MarlinUtil REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
find_package(PandoraSDK REQUIRED)
FOREACH( pkg  MarlinUtil PandoraSDK )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/components/*.cpp
    ${PROJECT_SOURCE_DIR}/src/genfit_interfaces/*.cpp
    ${PROJECT_SOURCE_DIR}/src/onnx_helpers/*.cpp
)


file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.hpp
  ${PROJECT_SOURCE_DIR}/include/genfit_interfaces/*.hpp
  ${PROJECT_SOURCE_DIR}/include/onnx_helpers/*.hpp
)

set(CMAKE_REQUIRED_LIBRARIES DD4hep::DDRec)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  k4FWCore::k4FWCore
  extensionDict
  OnnxRuntime
  torch
  GenFit::genfit2
  DD4hep::DDRec
  DD4hep::DDCore
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/genfit_interfaces>
  $<INSTALL_INTERFACE:include>
  ${MarlinUtil_INCLUDE_DIRS}
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)


SET(test_name "test_CreateInputForTrackTests")
ADD_TEST(NAME ${test_name} COMMAND bash -c "ddsim --enableGun --gun.distribution uniform --gun.energy '10*GeV' --gun.particle e- --numberOfEvents 100 --outputFile ddsim_output_edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml")
set_test_env(${test_name})

SET(test_name "test_TracksFromGenParticles")
ADD_TEST(NAME ${test_name} COMMAND k4run test/runTracksFromGenParticles.py)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES DEPENDS test_CreateInputForTrackTests)

#endif()
