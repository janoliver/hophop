# version
SET(PACKAGE_VERSION 1.2)

# optional CFLAGS
SET(cflags "-Wall -O3 -ffast-math -march=native")

# gengetopt
#SET(gengetopt "/import/vtthomes/oelerich/software/bin/gengetopt")
SET(gengetopt "gengetopt")

# generate current time
execute_process(
        COMMAND date +%Y-%m-%d_%H:%M
        OUTPUT_VARIABLE BUILD_TIME
)

# compilation -----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(hopping C)

CONFIGURE_FILE(options.ggo.in options.ggo @ONLY) 

ADD_CUSTOM_COMMAND(
        COMMAND ${gengetopt} --conf-parser -i options.ggo
        OUTPUT cmdline.h cmdline.c
)

SET_DIRECTORY_PROPERTIES(
        PROPERTIES 
        ADDITIONAL_MAKE_CLEAN_FILES cmdline.h;cmdline.c;options.ggo
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

TARGET_LINK_LIBRARIES(hop 
                          m 
                          gsl gslcblas 
                          gomp
)

SET(CMAKE_CXX_FLAGS "${cflags} -fopenmp")

INSTALL(TARGETS hop DESTINATION bin)
