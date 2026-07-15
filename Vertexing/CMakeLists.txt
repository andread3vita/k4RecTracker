set(PackageName Vertexing)
project(${PackageName})

set(VERTEX_EXTENSION_MODEL vertexExtension)

# Build the vertex extension data model.
PODIO_GENERATE_DATAMODEL(
  ${VERTEX_EXTENSION_MODEL}
  dataFormatExtension/vertex.yaml
  vertex_ext_headers
  vertex_ext_sources
  UPSTREAM_EDM edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR}
)

PODIO_ADD_DATAMODEL_CORE_LIB(
  ${VERTEX_EXTENSION_MODEL}
  "${vertex_ext_headers}"
  "${vertex_ext_sources}"
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(
  ${VERTEX_EXTENSION_MODEL}
  PUBLIC
    EDM4HEP::edm4hep
)

# Provide an extension:: namespace for the CMake target.
add_library(
  extension::vertex
  ALIAS ${VERTEX_EXTENSION_MODEL}
)

list(
  APPEND VERTEX_EXTENSION_INSTALL_LIBS
  ${VERTEX_EXTENSION_MODEL}
)

install(
  TARGETS ${VERTEX_EXTENSION_INSTALL_LIBS}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
          COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT dev
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/extension"
                COMPONENT dev
)

install(
  FILES dataFormatExtension/vertex.yaml
  DESTINATION "${CMAKE_INSTALL_PREFIX}/extension"
  COMPONENT dev
)

# Build the Vertexing Gaudi components.
file(
  GLOB sources
  CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/components/*.cpp
)

file(
  GLOB headers
  CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/include/*.h
)

gaudi_add_module(
  ${PackageName}
  SOURCES
    ${sources}
  LINK
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    ${VERTEX_EXTENSION_MODEL}
)

# The generated headers are located under:
#
#   ${CMAKE_CURRENT_BINARY_DIR}/extension/
#
# This path must be PRIVATE because it is only valid in the build tree and
# must not appear in the exported target's INTERFACE_INCLUDE_DIRECTORIES.
target_include_directories(
  ${PackageName}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}>
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# Ensure the generated datamodel exists before compiling the Gaudi module.
add_dependencies(
  ${PackageName}
  ${VERTEX_EXTENSION_MODEL}
)

set_target_properties(
  ${PackageName}
  PROPERTIES
    PUBLIC_HEADER "${headers}"
)

install(
  TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
          COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT shlib
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT dev
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"
                COMPONENT dev
)

# Vertexing tests.
file(
  GLOB scripts
  CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/test/testSimpleVertexFinder/*.py
  ${PROJECT_SOURCE_DIR}/test/testSimpleVertexFinder/*.sh
)

file(
  COPY ${scripts}
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test/testSimpleVertexFinder
)

install(
  FILES ${scripts}
  DESTINATION test/testSimpleVertexFinder
)

set(test_name "test_vertexFinder")

add_test(
  NAME ${test_name}
  COMMAND sh test_vertexFinder.sh
)

set_test_env(${test_name})

set_tests_properties(
  ${test_name}
  PROPERTIES
    WORKING_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/test/testSimpleVertexFinder"
)