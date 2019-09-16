#
# Find Zoltan include directories and libraries
#
# ZOLTAN_INCLUDE             - list of include paths to find netcdf.h
# ZOLTAN_LIB                 - list of libraries to link against when using Zoltan
# ZOLTAN_FOUND               - Do not attempt to use Zoltan if "no", "0", or undefined.


#set(ZOLTAN_PREFIX "" CACHE STRING "zoltan install directory")

#message("zoltan prefix: " ${ZOLTAN_PREFIX})

#find_path(ZOLTAN_INCLUDE_DIR zoltan.h PATHS "${ZOLTAN_PREFIX}/include")

#find_library(ZOLTAN_LIBRARY zoltan PATHS "${ZOLTAN_PREFIX}/lib")

#message("zoltan include: " ${ZOLTAN_INCLUDE_DIR})


#set (ZOLTAN_ROOT "" CACHE PATH "Path to search for Zoltan header and library files" )
# above lined modified by Dr_J
#hard coded directory, Zoltan is a library, so no executables are installed 
# PATH variable (it does not have a /bin directory (executable)) so I had to hardcode this 
#set(ZOLTAN_ROOT  "/cm/shared/apps/zoltan/cuda75/openmpi-2.0.1/gcc-4.8.1/3.8.3/" )

#set(ZOLTAN_ROOT  "/opt/packages/zoltan/3.83/" )

set(ZOLTAN_ROOT  $ENV{ZOLTANROOT} )

#for Bridges PSC
if(NOT ZOLTAN_ROOT)
set(ZOLTAN_ROOT $ENV{ZOLTAN_DIR})
endif()

#for Stampede TACC
if(NOT ZOLTAN_ROOT)
set(ZOLTAN_ROOT $ENV{TACC_ZOLTAN_DIR})
endif()

#for Titan ORNL
if(NOT ZOLTAN_ROOT)
set(ZOLTAN_ROOT $ENV{CRAY_TRILINOS_PREFIX_DIR})
endif()

#for COMET SDSC
if(NOT ZOLTAN_ROOT)
set(ZOLTAN_ROOT $ENV{TRILINOSHOME})
endif()


message("zoltan_dir : " ${ZOLTAN_ROOT})

set (ZOLTAN_FOUND NO CACHE INTERNAL "Found Zoltan components successfully." )

# look for the shared (dynamic) library first
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )

find_library( 
  ZOLTAN_LIB   
  NAMES "zoltan"  
  PATHS ${ZOLTAN_ROOT}
  PATH_SUFFIXES "lib" "lib64"
  NO_DEFAULT_PATH
)

# if dynamic lib (.so) is not available look for the static lib (.a)

if(NOT ZOLTAN_LIB)

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )

find_library( 
  ZOLTAN_LIB   
  NAMES "zoltan"  
  PATHS ${ZOLTAN_ROOT}
  PATH_SUFFIXES "lib" "lib64"
  NO_DEFAULT_PATH
)

endif()

# add include path 
find_path( 
  ZOLTAN_INCLUDE 
  NAMES "zoltan.h"
  PATHS ${ZOLTAN_ROOT} 
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH
)


message("zoltan_lib : " ${ZOLTAN_LIB})
message("zoltan_include : " ${ZOLTAN_INCLUDE})


include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (ZOLTAN_ROOT ZOLTAN_INCLUDE ZOLTAN_LIB)
#find_package_handle_standard_args (ZOLTAN "Zoltan not found, check environment variables ZOLTAN_ROOT" ZOLTAN_INCLUDE ZOLTAN_LIB)
find_package_handle_standard_args (ZOLTAN DEFAULT_MSG ZOLTAN_INCLUDE ZOLTAN_LIB)

mark_as_advanced(
  ZOLTAN_ROOT
  ZOLTAN_INCLUDE
  ZOLTAN_LIB
)


