# version, package name and cmake version
SET(PACKAGE_VERSION "2.3")
SET(PACKAGE_NAME "HOP")
SET(gengetopt "gengetopt")

PROJECT(hopping C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

ADD_DEFINITIONS(
    -DPKG_VERSION=\"${PACKAGE_VERSION}\"
    -DPKG_NAME=\"${PACKAGE_NAME}\"
)

# determine, whether we want a static binary
SET(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")

# determine the build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "RELEASE")
ENDIF(NOT CMAKE_BUILD_TYPE)

# Some information output
MESSAGE(STATUS "Build and link statically: ${STATIC_LINKING}")
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# libraries #############################################

# do we want static libraries?
IF(STATIC_LINKING)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    FIND_LIBRARY(MATH_LIB NAMES libm.a)
ELSE(STATIC_LINKING)
    FIND_LIBRARY(MATH_LIB NAMES libm.so)
ENDIF(STATIC_LINKING)

# look for GSL
FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
SET(LIBS ${LIBS} ${GSL_LIBRARIES})
SET(LIBS ${LIBS} ${MATH_LIB})

FIND_PACKAGE( OpenMP REQUIRED)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# compilation #############################################

ADD_CUSTOM_COMMAND(
    COMMAND ${gengetopt} --conf-parser --show-required -l  
    --set-version=${PACKAGE_VERSION} 
    --set-package=${PACKAGE_NAME} 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/options.ggo
    OUTPUT cmdline.h cmdline.c
)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(INTERNAL_FILES 
    hop.c hop.h
    mc.c 
    mc_init.c 
    mc_hopping.c 
    mc_analyze.c 
    output.c 
    params.c
    be.c
    helper.c
    analytics.c
)

SET(EXTERNAL_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmdline.c 
    ${CMAKE_CURRENT_BINARY_DIR}/cmdline.h
    mgmres.c mgmres.h
)

# CFLAGS
SET_SOURCE_FILES_PROPERTIES(
    ${INTERNAL_FILES}
    PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -std=c99"
)
SET(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math")
SET(CMAKE_C_FLAGS_DEBUG "-g")
SET(CMAKE_C_FLAGS_PROFILE "-g -pg")

ADD_EXECUTABLE( hop 
    ${INTERNAL_FILES}
    ${EXTERNAL_FILES}
)

# linking ################################################

TARGET_LINK_LIBRARIES(hop ${LIBS})

SET_TARGET_PROPERTIES (hop PROPERTIES LINKER_FLAGS ${OpenMP_C_FLAGS})

IF(STATIC_LINKING)
    SET(CMAKE_EXE_LINKER_FLAGS "-static")
    SET_TARGET_PROPERTIES(hop PROPERTIES 
        LINK_SEARCH_END_STATIC 1)
ENDIF(STATIC_LINKING)


# installation #############################################

INSTALL(TARGETS hop DESTINATION bin)