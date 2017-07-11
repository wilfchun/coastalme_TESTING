#########################################################################################
#
# CMakeLists.txt for CoastalME (Coastal Modelling Environment)
#
#########################################################################################
# Set CMake required version
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)


#########################################################################################
# Detect the build environment
if (UNIX)
   if (UNIX AND NOT APPLE AND NOT CYGWIN)
      set(LINUX TRUE)
      message("Creating Linux-Unix Makefile for CoastalME")
   endif (UNIX AND NOT APPLE AND NOT CYGWIN)

   if (CYGWIN)
      message("Creating Cygwin Makefile for CoastalME")
   endif (CYGWIN)

   if (APPLE)
      message(FATAL_ERROR "CoastalME cannot be built on OSX... yet")
   endif (APPLE)
endif (UNIX)

if (WIN32)
   message(FATAL_ERROR "CoastalME cannot be built on Windows, using the Visual Studio compiler... yet")
endif (WIN32)

if (MINGW)
   message(FATAL_ERROR "CoastyalME cannot be built on Windows, using the MinGW compiler... yet")
endif (MINGW)

message("")

#########################################################################################
# Set up project-specific defaults
project(CME)

set(CMAKE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")      # Assume that we are in the src folder
set(CMAKE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..")
#set(CMAKE_INSTALL_DIR "$ENV{HOME}/coast/CoastalME")      # May be changed by the user

file(GLOB CME_SOURCE_FILES *.cpp)
set(CME_EXECUTABLE cme)

# Set the path for our own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#########################################################################################
# Detect the compiler and set compiler-specific flags etc.
message("Checking compiler")

# We are using c++
enable_language(CXX)

# We need c++11, first set the c++11 flag (note that later on we check the compiler for c++11 capability)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Set the CoastalME CPU flag TODO Need to change this for different CPU architectures
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dx86")

# Use 64-bit libraries
set(FIND_LIBRARY_USE_LIB64_PATHS true)

if (CMAKE_COMPILER_IS_GNUCC)
   message("C++ compiler is ${CMAKE_CXX_COMPILER_ID} g++ ${CMAKE_CXX_COMPILER_VERSION}")

   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -W -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wredundant-decls")

   # For Debug builds under gcc
   if (CMAKE_BUILD_TYPE MATCHES Debug)
      if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.8")
         # Undefined behaviour sanitizer introduced with gcc 4.9, see http://developerblog.redhat.com/2014/10/16/gcc-undefined-behavior-sanitizer-ubsan/
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fsanitize=shift -fsanitize=integer-divide-by-zero -fsanitize=unreachable -fsanitize=null -fsanitize=return -fsanitize=address")
      endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.8")

      # TODO I think these were introduced in a later gcc version, need to check
      if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
         # Undefined behaviour sanitizer introduced with gcc 4.9, see http://developerblog.redhat.com/2014/10/16/gcc-undefined-behavior-sanitizer-ubsan/
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=vla-bound -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fsanitize=bounds -fsanitize=alignment -fsanitize=object-size -fsanitize=vptr")
         # TODO not yet working
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmpx -fcheck-pointer-bounds -fsanitize-blacklist=blacklist.txt")
      endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")

#      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
   endif (CMAKE_BUILD_TYPE MATCHES Debug)

   # Finally set some gcc-specific linker flags
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed -Wl,--no-undefined")
endif(CMAKE_COMPILER_IS_GNUCC)

if (MSVC)
   message("Compiler is Visual Studio")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif (MSVC)

message("")

#########################################################################################
# If not specified, assume a Debug build
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug)
endif (NOT CMAKE_BUILD_TYPE)

# Sort out build types
if (CMAKE_BUILD_TYPE MATCHES Debug)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Debug CACHE PATH "Directory for Debug builds" FORCE)
   # For Debug builds on any platform, add the model's "internal" debug flag (used for extra logging) etc.
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
endif (CMAKE_BUILD_TYPE MATCHES Debug)

if (CMAKE_BUILD_TYPE MATCHES Release)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Release CACHE PATH "Directory for Release builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES Release)

if (CMAKE_BUILD_TYPE MATCHES gcov)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 --coverage")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/gcov CACHE PATH "Directory for gcov builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES gcov)

if (CMAKE_BUILD_TYPE MATCHES Callgrind)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Callgrind CACHE PATH "Directory for Callgrind builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES Callgrind)


#########################################################################################
# Find libraries
message("Finding libraries")

# The model cannot be built without GDAL
find_package(GDAL REQUIRED)
#message(STATUS "GDAL_LIBRARIES=${GDAL_LIBRARIES}")
#message(STATUS "GDAL_INCLUDE_DIRS=${GDAL_INCLUDE_DIRS}")
set(LIBS ${LIBS} ${GDAL_LIBRARIES})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GDAL_INCLUDE_DIRS})
#message(STATUS "LIBS=${LIBS}")
#message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")

# Added by DFM, stolen from https://github.com/qgis/QGIS/blob/master/cmake/FindGDAL.cmake
set(GDAL_CONFIG_PREFER_PATH "$ENV{GDAL_HOME}/bin" CACHE STRING "preferred path to GDAL (gdal_config)")
set(GDAL_CONFIG_PREFER_FWTOOLS_PATH "$ENV{FWTOOLS_HOME}/bin_safe" CACHE STRING "preferred path to GDAL (gdal_config) from FWTools")
find_program(GDAL_CONFIG gdal-config
   ${GDAL_CONFIG_PREFER_PATH}
   ${GDAL_CONFIG_PREFER_FWTOOLS_PATH}
   /usr/local/bin/
   /usr/bin/
)
#message("GDAL_CONFIG=${GDAL_CONFIG}")

if (GDAL_CONFIG)
   # Extract the GDAL version
   exec_program(${GDAL_CONFIG}
      ARGS --version
      OUTPUT_VARIABLE GDAL_VERSION
   )
   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" GDAL_VERSION_MAJOR "${GDAL_VERSION}")
   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" GDAL_VERSION_MINOR "${GDAL_VERSION}")

   message (STATUS "GDAL is version ${GDAL_VERSION}")
#message("GDAL_VERSION_MAJOR=${GDAL_VERSION_MAJOR}")
#message("GDAL_VERSION_MINOR=${GDAL_VERSION_MINOR}")

   # Now check the GDAL version: CoastalME can only be built with GDAL 2.1 or above
   if (GDAL_VERSION_MAJOR LESS 2)
      message(FATAL_ERROR "CoastalME requires GDAL 2.1 or later. You have version ${GDAL_VERSION}. Please upgrade")
   endif (GDAL_VERSION_MAJOR LESS 2)

   if (GDAL_VERSION_MAJOR GREATER 2)
      if (GDAL_VERSION_MINOR LESS 1)
         message(FATAL_ERROR "CoastalME requires GDAL 2.1 or later. You have version ${GDAL_VERSION}. Please upgrade")
      endif (GDAL_VERSION_MINOR LESS 1)
   endif (GDAL_VERSION_MAJOR GREATER 2)
else (GDAL_CONFIG)
   message(FATAL_ERROR "${GDAL_CONFIG} not found, cannot determine GDAL version. Check your GDAL installation")
endif (GDAL_CONFIG)


#########################################################################################
# The important bits
include_directories(SYSTEM ${CMAKE_INCLUDE_PATH})
include_directories(${CME_SOURCE_DIR})
add_executable(${CME_EXECUTABLE} ${CME_SOURCE_FILES})

install(TARGETS ${CME_EXECUTABLE} RUNTIME DESTINATION ${CMAKE_INSTALL_DIR})

# Check to see if the compiler supports c++11
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD_REQUIRED ON)


#########################################################################################
# Sort out CShore

# The cshore library needs gfortran on Linux TODO What about Windows?
if (UNIX)
   find_package(GFortranLibs REQUIRED)
   include_directories(${GFORTRAN_INCLUDE_DIR})
   target_link_libraries(${CME_EXECUTABLE} ${LIBS} ${LIBGFORTRAN_LIBRARIES})
   target_link_libraries(${CME_EXECUTABLE} ${LIBS} ${LIBQUADMATH_LIBRARIES})

   if (UNIX AND NOT APPLE AND NOT CYGWIN)
      if (CMAKE_BUILD_TYPE MATCHES Debug)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
      else (CMAKE_BUILD_TYPE MATCHES Debug)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
      endif (CMAKE_BUILD_TYPE MATCHES Debug)
   endif (UNIX AND NOT APPLE AND NOT CYGWIN)

   if (CYGWIN)
      if (CMAKE_BUILD_TYPE MATCHES Debug)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
      else (CMAKE_BUILD_TYPE MATCHES Debug)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
      endif (CMAKE_BUILD_TYPE MATCHES Debug)
   endif (CYGWIN)

   if (APPLE)
      # TODO
   endif (APPLE)

   # Link with the cshore library
   target_link_libraries(${CME_EXECUTABLE} ${LIBS} ${CME_SOURCE_DIR}/lib/libcshore.a)
endif (UNIX)

if (WIN32)
   # TODO
endif (WIN32)


#########################################################################################
# Tell the user what we have found
message("")
message("Writing Makefile")
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_DIR = ${CMAKE_INSTALL_DIR}")
# message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
# message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
# message(STATUS "LIBS=${LIBS}")
# message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
# message(STATUS "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")


#########################################################################################
# Some extra messages
if (CMAKE_BUILD_TYPE MATCHES Callgrind)
   message("")
   message("When the build has finished, use valgrind/callgrid as follows:")
   message("")
   message("valgrind --tool=callgrind ./cme")
   message("")
   message("Then run:")
   message("callgrind_annotate --auto=yes callgrind.out.XXXX > ./profile/callgrind/callgrind.txt")
   message("where XXXX is the number of the callgrind.out.XXXX that was produced by valgrind")
   message("")
   message("Then look at callgrind.txt")
   message("")

endif (CMAKE_BUILD_TYPE MATCHES Callgrind)

if (CMAKE_BUILD_TYPE MATCHES gcov)
   message("")
   message("When the build has finished, use gcov/lcov as follows:")
   message("")
   message("./cme")
   message("lcov --capture --directory ./src/CMakeFiles/cme.dir/ --output-file ./profile/lcov_output/coverage.info")
   message("cd ./profile/lcov_output")
   message("genhtml coverage.info")
   message("")
   message("Then look at index.html in your browser")
   message("")
endif (CMAKE_BUILD_TYPE MATCHES gcov)

#########################################################################################
