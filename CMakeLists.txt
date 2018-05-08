#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(eicdirc)

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
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Generate dictionary for root
#
execute_process(COMMAND rootcint -f PrtDict.cxx -c ../src/PrtHit.h ../src/PrtEvent.h ../src/PrtLutNode.h ../src/PrtTrackInfo.h  ../src/PrtPhotonInfo.h ../src/PrtAmbiguityInfo.h ../src/PrtLinkDef.h)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/src
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/src/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#set (CMAKE_CXX_FLAGS "-std=gnu++11")
add_executable(eicdirc PrtDirc.cxx ${sources} PrtDict.cxx ${headers})
target_link_libraries(eicdirc -lGLU ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lSpectrum -lHistPainter)
add_library(prtdirclib SHARED ../src/PrtHit.cxx ../src/PrtEvent.cxx ../src/PrtLutNode.cxx ../src/PrtTrackInfo.cxx  ../src/PrtPhotonInfo.cxx ../src/PrtAmbiguityInfo.cxx PrtDict.cxx ../src/PrtHit.h ../src/PrtEvent.h ../src/PrtLutNode.h ../src/PrtTrackInfo.h  ../src/PrtPhotonInfo.h ../src/PrtAmbiguityInfo.h)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Prt. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Prt_SCRIPTS
    PrtDirc.out
    PrtDirc.mac
    gui.mac
    icons.mac
    run.png
    vis.mac
  )

foreach(_script ${Prt_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS eicdirc DESTINATION bin)

