set(PackageName Vertexing)

project(${PackageName})

set(VERTEX_EXTENSION_MODEL vertexExtension)

###############################################################################
# Vertex EDM extension
###############################################################################

PODIO_GENERATE_DATAMODEL(
  ${VERTEX_EXTENSION_MODEL}
  dataFormatExtension/vertex.yaml
  vertex_ext_headers
  vertex_ext_sources

  UPSTREAM_EDM
    edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml

  IO_BACKEND_HANDLERS
    ${PODIO_IO_HANDLERS}

  OUTPUT_FOLDER
    ${CMAKE_CURRENT_BINARY_DIR}
)

PODIO_ADD_DATAMODEL_CORE_LIB(
  ${VERTEX_EXTENSION_MODEL}
  "${vertex_ext_headers}"
  "${vertex_ext_sources}"

  OUTPUT_FOLDER
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(
  ${VERTEX_EXTENSION_MODEL}
  PUBLIC
    EDM4HEP::edm4hep
)

# Generated headers are placed in:
#
#   ${CMAKE_CURRENT_BINARY_DIR}/extension/
#
# Adding ${CMAKE_CURRENT_BINARY_DIR} as an include root allows:
#
#   #include "extension/VertexCollection.h"
#
target_include_directories(
  ${VERTEX_EXTENSION_MODEL}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(
  extension::vertex
  ALIAS
  ${VERTEX_EXTENSION_MODEL}
)

install(
  TARGETS
    ${VERTEX_EXTENSION_MODEL}

  EXPORT
    ${CMAKE_PROJECT_NAME}Targets

  RUNTIME DESTINATION
    "${CMAKE_INSTALL_BINDIR}"
  COMPONENT bin

  LIBRARY DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT shlib

  ARCHIVE DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT dev

  PUBLIC_HEADER DESTINATION
    "${CMAKE_INSTALL_INCLUDEDIR}/extension"
  COMPONENT dev
)

install(
  FILES
    dataFormatExtension/vertex.yaml

  DESTINATION
    "${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/extension"

  COMPONENT dev
)

###############################################################################
# Vertexing Gaudi components
###############################################################################

file(
  GLOB vertexing_sources
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/components/*.cpp"
)

file(
  GLOB vertexing_headers
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

gaudi_add_module(
  ${PackageName}

  SOURCES
    ${vertexing_sources}

  LINK
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    ${VERTEX_EXTENSION_MODEL}
)

target_include_directories(
  ${PackageName}

  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}>

  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# Ensure the PODIO-generated files exist before SimpleVertexFinder.cpp
# and the other Gaudi component sources are compiled.
add_dependencies(
  ${PackageName}
  ${VERTEX_EXTENSION_MODEL}
)

set_target_properties(
  ${PackageName}
  PROPERTIES
    PUBLIC_HEADER "${vertexing_headers}"
)

install(
  TARGETS
    ${PackageName}

  EXPORT
    ${CMAKE_PROJECT_NAME}Targets

  RUNTIME DESTINATION
    "${CMAKE_INSTALL_BINDIR}"
  COMPONENT bin

  LIBRARY DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT shlib

  ARCHIVE DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"
  COMPONENT dev

  PUBLIC_HEADER DESTINATION
    "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"
  COMPONENT dev
)

###############################################################################
# Vertexing tests
###############################################################################

file(
  GLOB vertexing_test_scripts
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/test/testSimpleVertexFinder/*.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/test/testSimpleVertexFinder/*.sh"
)

file(
  COPY
    ${vertexing_test_scripts}

  DESTINATION
    "${CMAKE_CURRENT_BINARY_DIR}/test/testSimpleVertexFinder"
)

install(
  FILES
    ${vertexing_test_scripts}

  DESTINATION
    "test/testSimpleVertexFinder"
)

set(test_name "test_vertexFinder")

add_test(
  NAME
    ${test_name}

  COMMAND
    sh test_vertexFinder.sh
)

set_test_env(${test_name})

set_tests_properties(
  ${test_name}
  PROPERTIES
    WORKING_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/test/testSimpleVertexFinder"
)