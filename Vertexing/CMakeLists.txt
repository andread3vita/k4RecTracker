set(PackageName Vertexing)

project(${PackageName})

# Build the extension data model and link it against the upstream model
PODIO_GENERATE_DATAMODEL(extension dataFormatExtension/Vertex.yaml ext_headers ext_sources
  UPSTREAM_EDM edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
PODIO_ADD_DATAMODEL_CORE_LIB(VertexingExtension "${ext_headers}" "${ext_sources}"
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(VertexingExtension PUBLIC EDM4HEP::edm4hep)
PODIO_ADD_ROOT_IO_DICT(VertexingExtensionDict VertexingExtension "${ext_headers}" src/selection.xml
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
#PODIO_ADD_SIO_IO_BLOCKS(VertexingExtension "${ext_headers}" "${ext_sources}")
add_library(extension::extensionDict ALIAS VertexingExtensionDict)
list(APPEND EXTENSION_INSTALL_LIBS VertexingExtension VertexingExtensionDict)

install(TARGETS ${EXTENSION_INSTALL_LIBS}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/extension"
  COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/VertexingExtensionDictDict.rootmap"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

install(FILES
  dataFormatExtension/Vertex.yaml
  DESTINATION "${CMAKE_INSTALL_PREFIX}/extension" COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/libVertexingExtensionDict_rdict.pcm"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

# Sources: the kernel (src/) and the Gaudi components (components/)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
  ${PROJECT_SOURCE_DIR}/components/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    k4FWCore::k4FWCore
    ROOT::Core
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
    extensionDict
    VertexingExtensionDict
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)

SET(test_name "test_vertexFinder")
ADD_TEST(NAME ${test_name} COMMAND sh +x test_vertexFinder.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Vertexing/test/testDeterministicAnnealingVertexFindereter")


SET(test_name "test_vertexFitter")
ADD_TEST(NAME ${test_name} COMMAND sh +x test_vertexFitter.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES DEPENDS test_vertexFinder)
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Vertexing/test/testLeastSquaresVertexFitter")
