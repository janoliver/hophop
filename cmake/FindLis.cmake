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

# - Try to find LIS
# Once done, this will define
#
#  LIS_FOUND
#  LIS_INCLUDE_DIRS
#  LIS_LIBRARIES
#
#  Environment variable LIS_ROOT_DIR can be set to give hints


IF (LIS_ROOT_DIR)

    FIND_LIBRARY(
            LIS_LIBRARY
            NAMES "lis"
            PATHS ${LIS_ROOT_DIR}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )

    FIND_PATH(
            LIS_INCLUDE_DIR
            NAMES "lis.h"
            PATHS ${LIS_ROOT_DIR}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )

ELSE (LIS_ROOT_DIR)

    FIND_LIBRARY(
            LIS_LIBRARY
            NAMES "lis"
    )

    FIND_PATH(
            LIS_INCLUDE_DIR
            NAMES "lis.h"
    )

ENDIF (LIS_ROOT_DIR)


set(LIS_LIBRARIES ${LIS_LIBRARY})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIS DEFAULT_MSG LIS_LIBRARIES LIS_INCLUDE_DIR)
