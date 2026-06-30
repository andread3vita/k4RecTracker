set(PackageName Vertexing)

project(${PackageName})

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

SET(test_name "test_vertexFitter")
ADD_TEST(NAME ${test_name} COMMAND sh +x test_vertexFitter.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Vertexing/test/testLeastSquaresVertexFitter")
