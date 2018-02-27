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

find_program(GENGETOPT_EXECUTABLE gengetopt)

add_custom_target(gen-cli
        COMMAND ${GENGETOPT_EXECUTABLE} --conf-parser --show-required -l
        --set-version=${PACKAGE_VERSION}
        --set-package=${PACKAGE_NAME}
        -i ${CMAKE_SOURCE_DIR}/src/cli/options.ggo
        --output-dir=${CMAKE_SOURCE_DIR}/src/cli/
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating CLI header and c file"
        OUTPUT ${CMAKE_SOURCE_DIR}/src/cli/cmdline.h ${CMAKE_SOURCE_DIR}/src/cli/cmdline.c
        VERBATIM
        )
