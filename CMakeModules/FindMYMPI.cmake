# taken from https://github.com/libigl/eigen/blob/master/cmake/FindMPI.cmake
# The part for finding suffixes is removed since it found the wrong sufixes
# - Find the MPI library
#
# Usage:
#   find_package(MPI [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   MPI_FOUND               ... true if fftw is found on the system
#   MPI_LIB                 ... full path to fftw library
#   MPI_INCLUDE             ... fftw include directory
#
# The following variables will be checked by the function
#   MPI_USE_STATIC_LIBS    ... if true, only static libraries are found
#   MPI_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#


#If environment variable MPI_DIR is specified, it has same effect as MPI_ROOT

#set( MPI_ROOT "/opt/packages/openmpi/openmpi-2.0.1/" )
set( MPI_ROOT $ENV{MPI_ROOT} )

if(NOT MPI_ROOT)
set(MPI_ROOT $ENV{MPIHOME})
endif()

if(NOT MPI_ROOT)
set(MPI_ROOT $ENV{MPI_DIR})
endif()
if(NOT MPI_ROOT)
set(MPI_ROOT $ENV{MPI_HOME})
endif()
if(NOT MPI_ROOT)
# stampede
set(MPI_ROOT $ENV{I_MPI_ROOT})
endif()
if(NOT MPI_ROOT)
# Titan
set(MPI_ROOT $ENV{MPICH_DIR})
endif()

message("MPI : " ${MPI_ROOT})

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT MPI_ROOT )
  pkg_check_modules( PKG_MPI QUIET "mpi" )
endif()

#Search for dynamic libs (.so)
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )


find_library(MPI_LIB
    NAMES "mpich" "mpi" 
    PATHS ${MPI_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
)
# if dynamic lib does not exist look for static
if(NOT MPI_LIB)

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )

find_library(MPI_LIB
    NAMES "mpich" "mpi" 
    PATHS ${MPI_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
)

endif()

#message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})

#find includes
  find_path(
    MPI_INCLUDE
    NAMES "mpi.h"
    PATHS ${MPI_ROOT}
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )


#message(" inside cmake ..... " $MPI_INCLUDE} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI DEFAULT_MSG
                                  MPI_INCLUDE MPI_LIB)

mark_as_advanced(MPI_INCLUDE MPI_LIB )

