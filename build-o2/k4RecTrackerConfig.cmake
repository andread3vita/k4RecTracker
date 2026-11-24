set(k4RecTracker_VERSION )


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was Project_Config_Template.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# modify if other PATH_VARS are used
set_and_check(k4RecTracker_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")
set_and_check(k4RecTracker_LIB_DIR "${PACKAGE_PREFIX_DIR}/lib")

include(CMakeFindDependencyMacro)
# modify to reflect dependencies
find_dependency(ROOT COMPONENTS RIO Tree)
find_dependency(Gaudi)
find_dependency(EDM4HEP)
find_dependency(k4FWCore)

include("${CMAKE_CURRENT_LIST_DIR}/k4RecTrackerTargets.cmake")

# Adapt to existing components
# check_required_components(k4RecTrackerPlugins)
