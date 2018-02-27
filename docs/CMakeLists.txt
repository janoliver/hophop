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

add_custom_target(docs-manual
        COMMAND sphinx-build -c ${CMAKE_CURRENT_SOURCE_DIR} -b html ${PROJECT_SOURCE_DIR}/docs manual
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Build html documentation"
        VERBATIM)

add_custom_target(docs-manual-tex
        COMMAND sphinx-build -c ${CMAKE_CURRENT_SOURCE_DIR} -b latex ${PROJECT_SOURCE_DIR}/docs manual-tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Build html documentation"
        VERBATIM)
