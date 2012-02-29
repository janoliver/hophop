# version
SET(PACKAGE_VERSION 1.2)
SET(PACKAGE_NAME "hop")

# optional CFLAGS
SET(CMAKE_C_FLAGS "-Wall -O3 -ffast-math -march=native -fopenmp")

# gengetopt
SET(gengetopt "gengetopt")

# generate current time
execute_process(
        COMMAND date +%Y-%m-%d_%H:%M
        OUTPUT_VARIABLE BUILD_TIME
)

# compilation -----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(hopping C)

ADD_CUSTOM_COMMAND(
        COMMAND ${gengetopt} --conf-parser --set-version=${PACKAGE_VERSION} --set-package=${PACKAGE_NAME} --show-required -l  -i options.ggo
        OUTPUT cmdline.h cmdline.c
)

SET_DIRECTORY_PROPERTIES(
        PROPERTIES 
        ADDITIONAL_MAKE_CLEAN_FILES "cmdline.h;cmdline.c"
)

ADD_EXECUTABLE(hop 
                   hop.c hop.h
                   mc.c 
                   mc_init.c 
                   mc_hopping.c 
                   mc_analyze.c 
                   output.c 
                   cmdline.c cmdline.h
)

ADD_CUSTOM_TARGET(doc ALL
      COMMAND rst2html README.rst > doc/README.html)

TARGET_LINK_LIBRARIES(hop 
                          m 
                          gsl gslcblas 
                          gomp
)

#SET(CMAKE_CXX_FLAGS "${cflags} -fopenmp")

INSTALL(TARGETS hop DESTINATION bin)
