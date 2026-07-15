set(PackageName Vertexing)

project(${PackageName} LANGUAGES CXX)

set(VERTEX_EXTENSION_MODEL extensionVertexing)
set(VERTEX_EXTENSION_DICT VertexingExtensionDict)

PODIO_GENERATE_DATAMODEL(
  ${VERTEX_EXTENSION_MODEL}
  dataFormatExtension/vertex.yaml
  ext_headers
  ext_sources

  UPSTREAM_EDM
    edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml

  IO_BACKEND_HANDLERS
    ${PODIO_IO_HANDLERS}

  OUTPUT_FOLDER
    ${CMAKE_CURRENT_BINARY_DIR}
)

PODIO_ADD_DATAMODEL_CORE_LIB(
  ${VERTEX_EXTENSION_MODEL}
  "${ext_headers}"
  "${ext_sources}"

  OUTPUT_FOLDER
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(
  ${VERTEX_EXTENSION_MODEL}
  PUBLIC
    EDM4HEP::edm4hep
)

PODIO_ADD_ROOT_IO_DICT(
  ${VERTEX_EXTENSION_DICT}
  ${VERTEX_EXTENSION_MODEL}
  "${ext_headers}"
  src/selection.xml

  OUTPUT_FOLDER
    ${CMAKE_CURRENT_BINARY_DIR}
)

add_library(
  extensionVertexing::extension
  ALIAS
  ${VERTEX_EXTENSION_MODEL}
)

add_library(
  extensionVertexing::extensionDict
  ALIAS
  ${VERTEX_EXTENSION_DICT}
)

set(
  VERTEX_COMPATIBILITY_INCLUDE_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/extension"
)

file(
  MAKE_DIRECTORY
  "${VERTEX_COMPATIBILITY_INCLUDE_DIR}"
)

set(vertex_compatibility_headers)

foreach(generated_header IN LISTS ext_headers)
  get_filename_component(
    generated_header_name
    "${generated_header}"
    NAME
  )

  set(
    compatibility_header
    "${VERTEX_COMPATIBILITY_INCLUDE_DIR}/${generated_header_name}"
  )

  file(
    WRITE
    "${compatibility_header}"
    "#pragma once\n"
    "#include \"${VERTEX_EXTENSION_MODEL}/${generated_header_name}\"\n"
  )

  list(
    APPEND
    vertex_compatibility_headers
    "${compatibility_header}"
  )
endforeach()

file(
  GLOB vertexing_sources
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
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
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    ROOT::Core
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
    ${VERTEX_EXTENSION_MODEL}
    ${VERTEX_EXTENSION_DICT}
)

target_include_directories(
  ${PackageName}

  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>

  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

add_dependencies(
  ${PackageName}
  ${VERTEX_EXTENSION_MODEL}
  ${VERTEX_EXTENSION_DICT}
)

if(vertexing_headers)
  set_target_properties(
    ${PackageName}
    PROPERTIES
      PUBLIC_HEADER "${vertexing_headers}"
  )
endif()

set(
  EXTENSION_INSTALL_LIBS
  ${VERTEX_EXTENSION_MODEL}
  ${VERTEX_EXTENSION_DICT}
)

install(
  TARGETS
    ${EXTENSION_INSTALL_LIBS}

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
    "${CMAKE_INSTALL_INCLUDEDIR}/${VERTEX_EXTENSION_MODEL}"
  COMPONENT dev
)


install(
  FILES
    ${vertex_compatibility_headers}

  DESTINATION
    "${CMAKE_INSTALL_INCLUDEDIR}/extension"

  COMPONENT dev
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${VERTEX_EXTENSION_DICT}Dict.rootmap"

  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"

  COMPONENT dev
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/lib${VERTEX_EXTENSION_DICT}_rdict.pcm"

  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}"

  COMPONENT dev
)

install(
  FILES
    dataFormatExtension/vertex.yaml

  DESTINATION
    "${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/extensionVertexing"

  COMPONENT dev
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
    "${CMAKE_INSTALL_INCLUDEDIR}/${PackageName}"
  COMPONENT dev
)

file(
  GLOB vertexing_test_scripts
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/test/testSimpleVertexFinder/*.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/test/testSimpleVertexFinder/*.sh"
)

file(
  MAKE_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/test/testSimpleVertexFinder"
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

  COMPONENT dev
)

set(test_name test_vertexFinder)

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