set(PackageName Tracking)

project(${PackageName})

find_package(GenFit REQUIRED)

# Collect source files
file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/measurement/*.cpp
)

# Collect header files
file(GLOB headers
    ${PROJECT_SOURCE_DIR}/include/*.hpp
    ${PROJECT_SOURCE_DIR}/include/measurement/*.hpp
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

# Ensure proper linking
target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/measurement>
  $<INSTALL_INTERFACE:include>
  ${GENFIT_INCLUDE_DIR}
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

# Install scripts
file(GLOB scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PackageName}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)
