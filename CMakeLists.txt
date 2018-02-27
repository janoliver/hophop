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

set(PACKAGE_VERSION "2.4")
set(PACKAGE_NAME "hophop")
set(PACKAGE_DESCRIPTION "A fast hopping charge transport simulation")
set(PACKAGE_AUTHOR "Jan Oliver Oelerich")
set(PACKAGE_AUTHOR_EMAIL "jan.oliver.oelerich@physik.uni-marburg.de")

project(hopping C)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

add_definitions(
    -DPKG_VERSION=\"${PACKAGE_VERSION}\"
    -DPKG_NAME=\"${PACKAGE_NAME}\"
)

# determine the build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

# Some information output
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# libraries #############################################

# math
set(LIBS ${LIBS} m)

# look for GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

find_package(Lis)
if(LIS_FOUND)
    include_directories(${LIS_INCLUDE_DIR})
    set(LIBS ${LIBS} ${LIS_LIBRARIES})
    add_definitions(-DWITH_LIS)
endif(LIS_FOUND)

find_package( OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_subdirectory(src)
add_subdirectory(src/cli)
add_subdirectory(docs)