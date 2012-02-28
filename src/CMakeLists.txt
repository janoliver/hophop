# version
SET(PACKAGE_VERSION 1.2)

# optional CFLAGS
SET(cflags "-Wall -O3 -ffast-math -march=native")

# gengetopt
SET(gengetopt "gengetopt")
SET(rst2html "rst2html")
SET(rst2man "rst2man")
SET(gzip "gzip")

# generate current time
execute_process(
        COMMAND date +%Y-%m-%d_%H:%M
        OUTPUT_VARIABLE BUILD_TIME
)

# compilation -----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(hopping C)

ADD_CUSTOM_COMMAND(
        COMMAND ${gengetopt} --conf-parser -i ../res/options.ggo
        OUTPUT cmdline.h cmdline.c
)

ADD_CUSTOM_TARGET(doc ALL
        COMMAND ${rst2html} ../res/README.rst ../doc/README.html
        COMMAND ${rst2man} ../res/README.rst ../doc/hop.1
        COMMAND ${gzip} ../doc/hop.1
)

SET_DIRECTORY_PROPERTIES(
        PROPERTIES 
        ADDITIONAL_MAKE_CLEAN_FILES cmdline.h;cmdline.c;../doc/README.html
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
INSTALL(FILES ../doc/hop.1.gz DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man/man1/)
