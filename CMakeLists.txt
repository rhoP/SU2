cmake_minimum_required(VERSION 3.14)

project(
  SU2
  VERSION 7.0
  LANGUAGES CXX C)
set(SU2_VERSION_NAME "Falcon" CACHE INTERNAL "")

set(CMAKE_CXX_STANDARD 14)


set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
include(cmake/ConfigureExternals.cmake)



set(CMAKE_MODULE_PATH "/home/rohit/libtorch/share/cmake/Torch;${CMAKE_MODULE_PATH}")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

su2_option(SU2_BUILD_CFD SETUP ON OFF BOOL "Build the SU2_CFD executable")


add_subdirectory(Common)


if(SU2_BUILD_CFD)
  add_subdirectory(SU2_CFD)
endif()

message(
  STATUS
    "
-------------------------------------------------------------------------
|    ___ _   _ ___                                                      |
|   / __| | | |_  )   Release ${SU2_VERSION} '${SU2_VERSION_NAME}'                            |
|   \\__ \\ |_| |/ /                                                      |
|   |___/\\___//___|   Suite                                             |
|                                                                       |
-------------------------------------------------------------------------
| The current SU2 release has been coordinated by the                   |
| SU2 International Developers Society <www.su2devsociety.org>          |
| with selected contributions from the open-source community.           |
-------------------------------------------------------------------------
| The main research teams contributing to the current release are:      |
| - Prof. Juan J. Alonso's group at Stanford University.                |
| - Prof. Piero Colonna's group at Delft University of Technology.      |
| - Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology. |
| - Prof. Alberto Guardone's group at Polytechnic University of Milan.  |
| - Prof. Rafael Palacios' group at Imperial College London.            |
| - Prof. Vincent Terrapon's group at the University of Liege.          |
| - Prof. Edwin van der Weide's group at the University of Twente.      |
| - Lab. of New Concepts in Aeronautics at Tech. Inst. of Aeronautics.  |
-------------------------------------------------------------------------
| Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,       |
|                      Tim Albring, and the SU2 contributors.           |
|                                                                       |
| SU2 is free software; you can redistribute it and/or                  |
| modify it under the terms of the GNU Lesser General Public            |
| License as published by the Free Software Foundation; either          |
| version 2.1 of the License, or (at your option) any later version.    |
|                                                                       |
| SU2 is distributed in the hope that it will be useful,                |
| but WITHOUT ANY WARRANTY; without even the implied warranty of        |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |
| Lesser General Public License for more details.                       |
|                                                                       |
| You should have received a copy of the GNU Lesser General Public      |
| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |
-------------------------------------------------------------------------

Build Configuration Summary:

    Source code location: ${CMAKE_SOURCE_DIR}
    Install location:     ${CMAKE_INSTALL_PREFIX}
    Version:              ${SU2_VERSION}
    C++ Compiler:         ${CMAKE_CXX_COMPILER}
    C Compiler:           ${CMAKE_C_COMPILER}
    MPI support:          ${SU2_ENABLE_MPI}
    Metis support:        ${SU2_ENABLE_METIS}
    Parmetis support:     ${SU2_ENABLE_PARMETIS}
    TecIO support:        ${SU2_ENABLE_TECIO}
    CGNS support:         ${SU2_ENABLE_CGNS}
    Mutation++ support:   ${SU2_ENABLE_MUTATIONPP}
    MKL support:          ${SU2_ENABLE_MKL}
    Datatype support:
        double            ${SU2_BUILD_NORMAL}
        codi_reverse      ${SU2_BUILD_REVERSE}
        codi_forward      ${SU2_BUILD_DIRECTDIFF}

    Build SU2_CFD:        ${SU2_BUILD_CFD}
    Build SU2_DOT:        ${SU2_BUILD_DOT}
    Build SU2_MSH:        ${SU2_BUILD_MSH}
    Build SU2_DEF:        ${SU2_BUILD_DEF}
    Build SU2_SOL:        ${SU2_BUILD_SOL}
    Build SU2_GEO:        ${SU2_BUILD_GEO}
    Build Py Wrapper:     ${SU2_BUILD_PY_WRAPPER}

Please be sure to add the \$SU2_HOME and \$SU2_RUN environment variables,
and update your \$PATH (and \$PYTHONPATH if applicable) with \$SU2_RUN.

Based on the input to this configuration, add these lines to your .bashrc file:

export SU2_RUN=\"${CMAKE_INSTALL_PREFIX}/bin\"
export SU2_HOME=\"${CMAKE_SOURCE_DIR}\"
export PATH=\$PATH:\$SU2_RUN
export PYTHONPATH=\$PYTHONPATH:\$SU2_RUN
")
