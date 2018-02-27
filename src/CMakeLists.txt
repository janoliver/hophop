#
# hophop: Charge transport simulations in disordered systems
#
# Copyright (c) 2012-2018 Jan Oliver Oelerich <jan.oliver.oelerich@physik.uni-marburg.de>
# Copyright (c) 2012-2018 Disordered Many-Particle Physics Group, Philipps-Universität Marburg, Germany
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

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
        cli/cmdline.c
        cli/cmdline.h
        ${CMAKE_SOURCE_DIR}/3rdparty/mgmres/mgmres.c
        ${CMAKE_SOURCE_DIR}/3rdparty/mgmres/mgmres.h
        )

# CFLAGS
set_source_files_properties(
        ${INTERNAL_FILES}
        PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -std=c99"
)
set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_PROFILE "-g -pg")

add_executable(hophop
        ${INTERNAL_FILES}
        ${EXTERNAL_FILES}
        )

# linking ################################################

target_link_libraries(hophop ${LIBS})

set_target_properties(hophop PROPERTIES LINKER_FLAGS ${OpenMP_C_FLAGS})

install(TARGETS hophop DESTINATION bin)