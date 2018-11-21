#----------------------------------------------------------------------------
# Setup the project
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TOT_InjAnalysis)

#----------------------------------------------------------------------------
find_package( ROOT REQUIRED ) 
FOREACH( pkg ROOT)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

#----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${EXTERN_BOOST_INCLUDE_PREFIX})
link_directories( ${link_directories} ${EXTERN_BOOST_LIB_PREFIX}  )


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable 
#
add_executable(totCalibration totCalibration.cxx ${sources} ${headers})
add_executable(plotter plotter.cxx ${sources} ${headers})
add_executable(testFitter testFitter.cxx ${sources} ${headers})
add_executable(calibration calibration.cxx ${sources} ${headers})
#add_executable(pulseDelay pulseDelay.cxx ${sources} ${headers})
#add_executable(pulseShapeFit pulseShapeFit.cxx ${sources} ${headers})
#add_executable(pulseShapeFitJune2018 pulseShapeFitJune2018.cxx ${sources} ${headers})
add_executable(generateCalibrationFile generateCalibrationFile.cxx ${sources} ${headers})
target_link_libraries(totCalibration boost_program_options boost_iostreams)
target_link_libraries(plotter boost_program_options)
target_link_libraries(testFitter boost_program_options)
target_link_libraries(calibration boost_program_options)
#target_link_libraries(pulseDelay boost_program_options)
#target_link_libraries(pulseShapeFit boost_program_options)
#target_link_libraries(pulseShapeFitJune2018 boost_program_options)
target_link_libraries(generateCalibrationFile boost_program_options)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS totCalibration plotter testFitter calibration pulseDelay pulseShapeFit pulseShapeFitJune2018 generateCalibrationFile DESTINATION ${PROJECT_SOURCE_DIR}/bin)
