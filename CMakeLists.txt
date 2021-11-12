#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(musrSim)

#----------------------------------------------------------------------------
# KAMIL CHANGES:
#include_directories ("/home/download/Root_5.24_source/root/include")
# Load some basic macros which are needed later on
include(FindROOT.cmake)

set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR})
include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES${ROOT_LIBRARY_DIR})
link_directories( ${LINK_DIRECTORIES})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(musrSim musrSim.cc ${sources} ${headers})
##FIND_PACKAGE( ROOT )
#SET(ROOT_CONFIG_EXECUTABLE /home/download/Root_5.24_source/root/bin/root-config)
## ask root-config for the library varaibles
#EXEC_PROGRAM( ${ROOT_CONFIG_EXECUTABLE}
# ARGS "--noldflags --noauxlibs --libs"
# OUTPUT_VARIABLE root_flags )
##SET(ROOT_LIBRARIES ${root_flags})
#EchoString("Kamil ${root_flags}  Kamil")
#target_link_libraries(musrSim ${Geant4_LIBRARIES} ${root_flags})
target_link_libraries(musrSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  run/101.mac
  run/102.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(musrSim DEPENDS musrSim)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS musrSim DESTINATION bin)


