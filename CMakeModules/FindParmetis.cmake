# - Try to find parmetis
# Once done this will define
#  PARMETIS_FOUND - System has PARMETIS
#  PARMETIS_INCLUDE - The PARMETIS include directories
#  PARMETIS_LIB - The libraries needed to use PARMETIS
#  PARMETIS_DEFINITIONS - Compiler switches required for using PARMETIS

#set(PARMETIS_ROOT "/opt/packages/parmetis/gnu_openmpi/4.0.3/")
set(PARMETIS_ROOT $ENV{PARMETISROOT})


if(NOT PARMETIS_ROOT)
set(PARMETIS_ROOT $ENV{PARMETISHOME})
endif()

#bridges
if(NOT PARMETIS_ROOT)
set(PARMETIS_ROOT $ENV{PARMETIS_DIR})
endif()

#stampede
if(NOT PARMETIS_ROOT)
set(PARMETIS_ROOT $ENV{TACC_PARMETIS_DIR})
endif()

#titan
if(NOT PARMETIS_ROOT)
set(PARMETIS_ROOT $ENV{CRAY_TPSL_PREFIX_DIR})
endif()

message("parmetis prefix:" ${PARMETIS_ROOT})

message("parmetis root:" ${PARMETIS_ROOT})

# search for dynamic lib first
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )

find_library(
  PARMETIS_LIB 
  NAMES "parmetis"  
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES "lib" "lib64"
  NO_DEFAULT_PATH
)

# if dynamic lib does not exist look for static
if(NOT PARMETIS_LIB)

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )

find_library(
  PARMETIS_LIB 
  NAMES "parmetis" 
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES "lib" "lib64"
  NO_DEFAULT_PATH
)
endif()

find_path(
  PARMETIS_INCLUDE 
  NAMES "parmetis.h" 
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH
)

if(NOT EXISTS "${PARMETIS_INCLUDE}")
  message(FATAL_ERROR "parmetis include dir not found")
endif()


if(NOT EXISTS "${PARMETIS_LIB}")
  message(FATAL_ERROR "parmetis library not found")
endif()


find_library(
  METIS_LIB 
  NAMES metis 
  PATHS ${PARMETIS_ROOT}
  PATH_SUFFIXES "lib" "lib64"
  NO_DEFAULT_PATH
)
if(NOT EXISTS "${METIS_LIB}")
  message(FATAL_ERROR "metis library not found")
endif()

set(PARMETIS_LIBS ${PARMETIS_LIB} ${METIS_LIB})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PARMETIS  DEFAULT_MSG
    PARMETIS_LIBS METIS_LIB PARMETIS_INCLUDE)

    message (STATUS "---     ParMETIS Configuration ::")
    message (STATUS "        INCLUDES  : ${PARMETIS_INCLUDE}")
    message (STATUS "        LIBRARIES : ${PARMETIS_LIB}")


mark_as_advanced(PARMETIS_INCLUDE PARMETIS_LIBS )
