

##set( HDF5_ROOT "/pylon5/ac561fp/hasbesta/packages/hdf5-1.8.17/" )

#set(HDF5 "OFF")

set(HDF5_ROOT $ENV{HDF5_ROOT})

if(NOT HDF5_ROOT)
set(HDF5_ROOT $ENV{HDF5HOME})
endif()

if(NOT HDF5_ROOT)
set(HDF5_ROOT $ENV{HDF5_DIR})
endif()

if(NOT HDF5_ROOT)
set(HDF5_ROOT $ENV{TACC_HDF5_DIR})
endif()


message("HDF5_ROOT  :" ${HDF5_ROOT})
# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT HDF5_ROOT )
  pkg_check_modules( PKG_FFTW QUIET "hdf5" )
endif()

#Search for dynamic libs (.so)
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )

message("Suffixes"  ${CMAKE_FIND_LIBRARY_SUFFIXES})


if( HDF5_ROOT )

 find_library(
    HDF5_LIB
    NAMES "hdf5" "hdf5-parallel" "hdf5-shared"
    PATHS ${HDF5_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
)

# if dynamic lib (.so) is not available look for the static lib (.a)
if(NOT HDF5_LIB)

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )

 find_library(
    HDF5_LIB
    NAMES "hdf5" "hdf5-parallel" "hdf5-shared"
    PATHS ${HDF5_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
)


endif()




#find includes
  find_path(
    HDF5_INCLUDE
    NAMES "hdf5.h"
    PATHS ${HDF5_ROOT} 
    PATH_SUFFIXES "include"
    NO_DEFAULT_PATH
  )

endif( HDF5_ROOT )


message("HDF5 INCLUDE  :" ${HDF5_INCLUDE})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_INCLUDE HDF5_LIB)

mark_as_advanced(HDF5_C_INCLUDE_DIR HDF5_hdf5_LIBRARY_REALEASE)



