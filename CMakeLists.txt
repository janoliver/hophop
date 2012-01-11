cmake_minimum_required (VERSION 2.6)
project (hopping C)

# version
set (PACKAGE_VERSION 1.1)

# optional CFLAGS
set(cflags "-Wall -O3 -ffast-math -march=native")

# gengetopt
set(gengetopt "/import/vtthomes/oelerich/software/bin/gengetopt")

CONFIGURE_FILE(options.ggo.in options.ggo @ONLY) 
add_custom_command(TARGET hop PRE_BUILD SOURCE options.ggo COMMAND ${gengetopt})
add_executable(hop hop.c mc.c mc_init.c mc_hopping.c mc_analyze.c output.c cmdline.c cmdline.h hop.h)
target_link_libraries(hop m gsl gslcblas gomp)
set(CMAKE_CXX_FLAGS "${cflags} -fopenmp")

install (TARGETS hop DESTINATION bin)
