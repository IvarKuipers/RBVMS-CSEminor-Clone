#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.17)
project(RBVMS LANGUAGES CXX C CUDA)

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    coefficients.cpp
    evolution.cpp
    monitor.cpp
    precon.cpp
    weakform.cpp
    printInfo.cpp)

set(HEADERS
    coefficients.hpp
    evolution.hpp
    monitor.hpp
    precon.hpp
    weakform.hpp
    buildInfo.hpp)

#-------------------------------------------------------------------------------
# Create header
#-------------------------------------------------------------------------------
get_filename_component(tmp1 ${MFEM_DIR} REALPATH)
get_filename_component(tmp2 ${PROJECT_BINARY_DIR}/.. REALPATH)

if (${tmp1} STREQUAL ${tmp2})
  message(STATUS "Report RBVMS, MFEM and hypre in build info.")
  add_custom_command(OUTPUT buildInfo.hpp
                     PRE_BUILD
                     COMMAND /bin/bash
                             "${PROJECT_SOURCE_DIR}/../config/buildInfo.sh"
                             "${PROJECT_SOURCE_DIR}" "mfem" "hypre"
                     COMMENT "Getting git info")
else()
  message(STATUS "External MFEM and hypre libraries used.")
  message(STATUS "Only report RBVMS in build info.")
  add_custom_command(OUTPUT buildInfo.hpp
                     PRE_BUILD
                     COMMAND /bin/bash
                             "${PROJECT_SOURCE_DIR}/../config/buildInfo.sh"
                             "${PROJECT_SOURCE_DIR}"
                     COMMENT "Getting git info")
endif()

#-------------------------------------------------------------------------------
# Add packages
#-------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../config)

#-------------------------------------------------------------------------------
# Add MPI
#-------------------------------------------------------------------------------
find_package(MPI REQUIRED)

set(MPICXX mpicxx)
set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPICXX})

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}")

message(STATUS "Found mfem: in ${MFEM_DIR} (version ${MFEM_VERSION})")
message(STATUS "MFEM_CXX_FLAGS = ${MFEM_CXX_FLAGS}")
message(STATUS "MFEM_INCLUDE_DIRS = ${MFEM_INCLUDE_DIRS}")
message(STATUS "MFEM_LIBRARIES = ${MFEM_LIBRARIES}")

include_directories(${MFEM_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Add CUDA
#-------------------------------------------------------------------------------
find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbvms rbvms.cpp ${SOURCES} ${HEADERS})

target_link_libraries(rbvms PRIVATE ${MFEM_LIBRARIES} ${CUDA_LIBRARIES})
install(TARGETS rbvms DESTINATION bin)
