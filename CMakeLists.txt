# version, package name and cmake version
set(PACKAGE_VERSION "2.4")
set(PACKAGE_NAME "HOP")
set(gengetopt "gengetopt")

project(hopping C)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

add_definitions(
    -DPKG_VERSION=\"${PACKAGE_VERSION}\"
    -DPKG_NAME=\"${PACKAGE_NAME}\"
)

# determine, whether we want a static binary
set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")

# determine the build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

# Some information output
message(STATUS "Build and link statically: ${STATIC_LINKING}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# libraries #############################################

# do we want static libraries?
if(STATIC_LINKING)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    find_library(MATH_LIB NAMES libm.a)
else(STATIC_LINKING)
    find_library(MATH_LIB NAMES libm.so)
endif(STATIC_LINKING)

# look for GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${GSL_LIBRARIES})
set(LIBS ${LIBS} ${MATH_LIB})

find_package(Lis)
if(LIS_FOUND)
    include_directories(${LIS_INCLUDE_DIR})
    set(LIBS ${LIBS} ${LIS_LIBRARIES})
    add_definitions(-DWITH_LIS)
endif(LIS_FOUND)

find_package( OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# compilation #############################################

add_custom_command(
    COMMAND ${gengetopt} --conf-parser --show-required -l  
    --set-version=${PACKAGE_VERSION} 
    --set-package=${PACKAGE_NAME} 
    -i ${CMAKE_CURRENT_SOURCE_DIR}/options.ggo
    OUTPUT cmdline.h cmdline.c
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(INTERNAL_FILES 
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

set(EXTERNAL_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmdline.c 
    ${CMAKE_CURRENT_BINARY_DIR}/cmdline.h
    mgmres.c mgmres.h
)

# CFLAGS
set_source_files_properties(
    ${INTERNAL_FILES}
    PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -std=c99"
)
set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_PROFILE "-g -pg")

add_executable( hop 
    ${INTERNAL_FILES}
    ${EXTERNAL_FILES}
)

# linking ################################################

target_link_libraries(hop ${LIBS})

set_target_properties(hop PROPERTIES LINKER_FLAGS ${OpenMP_C_FLAGS})

if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set_target_properties(hop PROPERTIES 
        LINK_SEARCH_END_STATIC 1)
endif(STATIC_LINKING)


# installation #############################################

install(TARGETS hop DESTINATION bin)