set(PackageName Vertexing)
project(${PackageName})

set(VERTEX_EXTENSION_MODEL vertexExtension)

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
  PUBLIC EDM4HEP::edm4hep
)

PODIO_ADD_ROOT_IO_DICT(
  vertexExtensionDict
  ${VERTEX_EXTENSION_MODEL}
  "${vertex_ext_headers}"
  src/selection.xml
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR}
)

add_library(
  extension::vertexExtension
  ALIAS ${VERTEX_EXTENSION_MODEL}
)

add_library(
  extension::vertexExtensionDict
  ALIAS vertexExtensionDict
)

install(
  TARGETS
    ${VERTEX_EXTENSION_MODEL}
    vertexExtensionDict
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/extension"
)

install(
  FILES dataFormatExtension/vertex.yaml
  DESTINATION "${CMAKE_INSTALL_DATADIR}/k4RecTracker/extensions"
)