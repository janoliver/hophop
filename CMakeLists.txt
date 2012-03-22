# version
SET(PACKAGE_VERSION "1.2.2")
SET(PACKAGE_NAME "HOP")

# optional CFLAGS
SET(CMAKE_C_FLAGS "-O3 -ffast-math -march=native -fopenmp")

# gengetopt
SET(gengetopt "gengetopt")

# compilation -----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(hopping C)

ADD_CUSTOM_COMMAND(
        COMMAND ${gengetopt} --conf-parser --set-version=${PACKAGE_VERSION} --set-package=${PACKAGE_NAME} --show-required -l  -i options.ggo
        OUTPUT cmdline.h cmdline.c
)

CONFIGURE_FILE(
        "res/config.h.in"
        "config.h"
)

CONFIGURE_FILE(
        "res/options.ggo.in"
        "options.ggo"
)

SET_DIRECTORY_PROPERTIES(
        PROPERTIES 
        ADDITIONAL_MAKE_CLEAN_FILES "cmdline.h;cmdline.c;config.h;options.ggo"
)

ADD_EXECUTABLE(hop 
                   hop.c hop.h
                   mc.c 
                   mc_init.c 
                   mc_hopping.c 
                   mc_analyze.c 
                   output.c 
                   params.c
                   be.c
                   analytic.c
                   helper.c
                   cmdline.c cmdline.h
                   mgmres.c
)

TARGET_LINK_LIBRARIES(hop 
                          m 
                          gsl gslcblas 
                          gomp
)

INSTALL(TARGETS hop DESTINATION bin)
