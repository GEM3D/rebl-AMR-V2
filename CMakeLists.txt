cmake_minimum_required(VERSION 2.8.11)

project(amrGem)

#set(dir "/home/jaber/GEM_AMR/obj")

#disable in source build
#set(CMAKE_DISABLE_SOURCE_CHANGES ON)
#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)


message(STATUS "cmake is building on a ${CMAKE_SYSTEM_NAME} system")
set(dir ${PROJECT_SOURCE_DIR})
message("project source dir" ${dir})

set(CMAKE_FILES_DIRECTORY ${dir}/build)
message("Filed Directory" ${CMAKE_FILES_DIRECTORY})

#define some colors (Unix-based systems)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()


# get the number of physical procesors 

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${dir}/bin/)
#set(CMAKE_BUILD_DIRECTORY ${dir}/build)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/build)

message("project binary dir" ${EXECUTABLE_OUTPUT_PATH})
message("CMAKE SOURCE DIR:" ${CMAKE_SOURCE_DIR})
message("CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})


# Out-of-Source build is required by amrGem, this if statement prevents
# in-source-builds 

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()


message("cmake Home Dir:" ${CMAKE_HOME_DIRECTORY})
message("cmake Build Dir:" ${CMAKE_BUILD_DIRECTORY})
message("cmake Binary Dir:" ${CMAKE_BINARY_DIR})
message("cmake source dir:" ${CMAKE_SOURCE_DIR})

############################################################################
#
#                    SECTION: FINDING-PACKAGE
#
#        1) MPI and HDF5 have native FindXXX.cmake's provided by cmake
#        2) For Zoltan and ParMETIS one needs to specify FindXXX.cmake 
#        3) These non-native to cmake modules are included in CMakeModules
#
############################################################################

#set the minimum required versions for external packages

set(MPI_MIN_VERSION "2.0.1")
set(HDF5_MIN_VERSION "1.8.16")
set(ZOLTAN_MIN_VERSION "3.8.3")
set(PARMETIS_MIN_VERSION "4.0.3")

# tell the user the recommended package versions

message("")
message("${Green}**************************************************************************${ColourReset}")
message("")
message("${White}            suggested version for MPI:        ${ColourReset} " ${Red}${MPI_MIN_VERSION}${ColourReset})
message("${White}            suggested version for HDF5:       ${ColourReset} " ${Red}${HDF5_MIN_VERSION}${ColourReset})
message("${White}            suggested version for ZOLTAN:     ${ColourReset} " ${Red}${ZOLTAN_MIN_VERSION}${ColourReset})
message("${White}            suggested version for PARMETIS:   ${ColourReset} " ${Red}${PARMETIS_MIN_VERSION}${ColourReset})
message("")
message("${Green}**************************************************************************${ColourReset}")

# specify the directory including non-native FindXXX.cmake 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
#message( " cmake moduel path:" ${CMAKE_MODULE_PATH})
#message("GEM_AMR Source" ${GEM_AMR_SOURCE_DIR})

# find MPI
find_package(MYMPI ${MPI_MIN_VERSION} EXACT REQUIRED)
if(MPI_FOUND)
include_directories(${MPI_INCLUDE})
message( "${Red} MPI Include  : ${ColourReset}" ${MPI_INCLUDE})
message( "${Red} MPI Library  : ${ColourReset}" ${MPI_LIB} )
#message("${Green} 1. MPI Library : Found ${ColourReset}")
endif()

#find HDF5
find_package(MYHDF5 ${HDF5_MIN_VERSION} REQUIRED)
if(HDF5_FOUND)
include_directories(${HDF5_INCLUDE})
message( "${Red} HDF5 Library Found: ${ColourReset}" ${HDF5_INCLUDE})
#message("${Green} 2. HDF5 Library : Found ${ColourReset}")
endif()

#include(ProcessorCount)
#ProcessorCount(N)
#message("number of processors: "  ${N})
#if(NOT N EQUAL 0)
#  set(CTEST_BUILD_FLAGS -j${N})
#  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
#  set(CMAKE_MAKE_PROGRAM "${CMAKE_MAKE_PROGRAM} -j ${N}")	
#endif()
#message("cmake make program" ${CMAKE_MAKE_PROGRAM})

#add_library(zoltan /usr/local/zoltan/lib)
#link_directories(/usr/local/zoltan/lib)

find_package(Parmetis ${PARMETIS_MIN_VERSION} EXACT REQUIRED)
if(PARMETIS_FOUND)
message("${Green} 3. ParMETIS Library : Found: ${PARMETIS_LIBS} ${ColourReset}")
#message("version" ${PACKAGE_FIND_VERSION})
endif()

find_package(Zoltan ${ZOLTAN_MIN_VERSION} REQUIRED)
if(ZOLTAN_FOUND)
include_directories(${ZOLTAN_INCLUDE})
message("${Red} ZOLTAN Library Found: ${ColourReset}" ${ZOLTAN_INCLUDE})
message("${Green} 4. ZOLTAN Library : Found${ColourReset}")
endif()
#message( "ZOLTAN lib:  "  ${ZOLTAN_LIBRARIES})
#message( "ZOLTAN Include:  "  ${ZOLTAN_INCLUDE_DIR})
#message( "ZOLTAN Include:  "  ${ZOLTAN_INCLUDE_PATH})
#include_directories(/usr/local/zoltan/include)
#include_directories(/usr/local/zoltan/lib)

#include_directories(/usr/local/hdf5/1.18.15/include)
#message( "MPI" ${HDF5_C_LIBRARIES})

find_program(VALGRIND "valgrind")
if(VALGRIND)
message("${Red} Valgrind Found ${ColourReset}")
endif()

############################################################################
#
#                    SECTION: BUILD AND INSTALL
#
#        1) Build all the *.cpp and *.c files in "/src" directory
#        2) Link with external package libraries 
#        3) Install the final executable in "/bin"
#
############################################################################

include_directories(${dir}/src/include)
# dont forget to modify the address of the source files
file(GLOB OBJ "${dir}/src/*.cpp" "${dir}/src/*.c")
#set(CMAKE_BUILD_FILES_DIRECTORY ${dir})
#set(CMAKE_BUILD_DIRECTORY ${dir}/bin)      

set(EXECUTABLE_OUTPUT_PATH ${dir}/bin)

#add_definitions(-Wall -Wno-unused-variable  -std=c++0x)

# Define all header files for linking, this is required such that one doesnt
# have to hardcode the directories for I/O operations and hence ensure
# portability

set(HEADER_FILES "${dir}/src/include/typedefs.h"
"${dir}/src/include/parse_stl.h"
"${dir}/src/include/tree.h" "${dir}/src/include/forest.h"
"${dir}/src/include/communicate.h" "${dir}/src/include/definitions.h"
"${dir}/src/include/datatype.h" "${dir}/src/include/phdf5.h"
"${dir}/src/include/templatePhdf5.h"  
"${dir}/src/include/templateForest.h"
"${dir}/src/include/params.h"
) 

#set (CMAKE_CXX_STANDARD 11)
#set_property(TARGET tgt PROPERTY CXX_STANDARD 11)

if(UNIX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++1y -g  ")
#    SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -Wall -std=gnu++0x" )
SET(CMAKE_EXEC_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS}   -std=c++1y " )
endif()


#message( "obj" ${OBJ})
message( "*.h files" ${HEADER_FILES})

# Generate Object Files
add_executable(amrGem ${OBJ} ${HEADER_FILES}) 


#link_directories(/usr/local/zoltan/lib /usr/local/zoltan/include)
#set_target_properties(objects PROPERTIES LINKER_LANGUAGE Cxx)

target_link_libraries(amrGem  ${MPI_LIB} ${HDF5_LIB} ${ZOLTAN_LIB} ${PARMETIS_LIBS}  )

#install the final executable in bin directory

#install(TARGETS amrGem DESTINATION ${dir/}bin)

############################################################################
#
#                    SECTION: ARCHIEVING
#
#	 1) it gets the date and appends that to the project name
#	 2) compresses the "/src" in a TGZ format
# 	 3) Puts the resulting file in the "/archieve" directory
#
############################################################################

MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
         string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\2.\\1.\\3"
${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

TODAY(DATE)
message(STATUS "Compilation date = ${DATE}")


# the following two lines specify the file type for binary and source
# generation separately

SET(CPACK_GENERATOR "TGZ")
SET(CPACK_SOURCE_GENERATOR "TGZ")

#INCLUDE(CPack)
# Do not rely on CPACK_xxx variables, they are not set till INCLUDE(CPack) is
# invoked, so there is no way to set them and as cpack to name the file as you
# like, I hard coded the project name and the date as a file_name to be used
# by cpack, this works fine 

#message("package name: " ${CPACK_PACKAGE_NAME})

set(CPACK_PACKAGE_FILE_NAME "amrGem-${DATE}")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "amrGem-${DATE}")

# set the directory you would like to be compressed, For under development
# projects it is normally the src directory that includes *.cpp and *.h files

SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "${dir}/src;/")

SET(CPACK_PACKAGE_VENDOR "BSU")

# If you like to add any other files to archieve, need to use install file
# command, now in addition to source I am compresing the CMakeList.txt along
# with itimply installing it to the source directory. since src is the primary
# directory to archieve, make install command will copy this to src directory
# though, use the  following line if you would like to include
# CMakeList.txt in your archieve otherwise comment it out, one way to 
# get rid of this is to use an alias to remove that file from src after
# packing

install(FILES ${dir}/CMakeLists.txt DESTINATION ${dir}/src)

message("file name:" ${CPACK_PACKAGE_FILE_NAME})

#install(TARGETS  amrGem ARCHIVE  DESTINATION archieve)

# This variable (CPACK_OUTPUT_FILE_PREFIX)  will put the result of amrGem.date.tar.gz into the folder
# named archieve, this is only available in  in CMake versions >= 2.8.3 

set(CPACK_OUTPUT_FILE_PREFIX "${dir}/archieve")

#message("PACKAGE Install Dir:" ${CPACK_PACKAGE_INSTALL_DIRECTORY})

INCLUDE(CPack)

############################################################################
#
#                          SECTION: CTEST
#
############################################################################


#enable_testing()


#message("mpiexec" ${MPIEXEC})

#set(MPIRUN  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}  4) 


#add_test(NAME bunny COMMAND ${MPIRUN} "${dir}/bin/amrGem" "${dir}/input/bunny.stl")


# to invoke the ctest type; note that run.txt is GEM_AMR direcory
# ctest  <../run.txt -V

#set(valgrind ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}  4  ${VALGRIND} --tool=memcheck
#--leak-check=yes --show-reachable=yes --gen-suppressions=yes
#--collect-jumps=no  --track-fds=yes "${dir}/bin/{amrGem}")
#message("valgrind command" ${valgrind})


#add_test(NAME memcheck COMMAND valgrind "${dir}/bin/amrGem" "${dir}/input/porsche.stl")


#set(CTEST_MEMORYCHECK_COMMAND "/cm/shared/apps/valgrind/gcc/64/3.11.0/")
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--log-file=vlogs --leak-check=full --xml=yes")

# to run the test "ctest <../run.txt"; run.txt is a directory to specify the
# inputs for it  

#
#  Doxygen Integrated to Generate Documentation for the Code
#  Specify the INPUT parameter in Doxyfile as the src and include files
#  You need to tell doxygen where to look for files by giving the directory
#

#option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)

if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)

  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it
correctly")
  endif()

  add_custom_target (doc ALL
                     ${DOXYGEN_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
                     )

endif()


