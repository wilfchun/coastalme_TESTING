#########################################################################################
#
# CMakeLists.txt for CoastalME (Coastal Modelling Environment)
#
#########################################################################################
# Set CMake defaults
cmake_minimum_required(VERSION 2.8.12)


#########################################################################################
# Detect the environment, and set up compiler-independent stuff
if (UNIX)
   message("Creating Linux-Unix Makefile for CoastalME")
endif (UNIX)
if (WIN32)
   message("Creating Windows Makefile for CoastalME")
endif(WIN32)
if (CYGWIN)
   message("Creating Cygwin Makefile for CoastalME")
   message("NOTE: CShore is not yet set up for Cygwin")
endif (CYGWIN)
if (MINGW)
   message("Creating MinGW Makefile for CoastalME")
   message("NOTE: CShore is not yet set up for MinGW")
endif (MINGW)

# Set up some c++ and architecture defaults
enable_language(CXX)
set(FIND_LIBRARY_USE_LIB64_PATHS true)

# Set up project-specific defaults
project(CME)
file(GLOB CME_SOURCE_FILES *.cpp)
set(CME_EXECUTABLE cme)
#set(CME_INSTALL_DIR "$ENV{HOME}/coast/CoastalME")          # May be changed by the user
set(CME_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..")               # May be changed by the user

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug)                              # If not specified, assume a Debug build
endif (NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE MATCHES Debug)
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Debug CACHE PATH "Directory for Debug builds" FORCE)
   # For Debug builds on any platform, add the model's "internal" debug flag (used for extra logging) etc.
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
endif (CMAKE_BUILD_TYPE MATCHES Debug)

if (CMAKE_BUILD_TYPE MATCHES Release)
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Release CACHE PATH "Directory for Release builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES Release)

if (CMAKE_BUILD_TYPE MATCHES gcov)
#   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 --coverage")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/gcov CACHE PATH "Directory for gcov builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES gcov)

if (CMAKE_BUILD_TYPE MATCHES Callgrind)
#   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Callgrind CACHE PATH "Directory for Callgrind builds" FORCE)
endif (CMAKE_BUILD_TYPE MATCHES Callgrind)

# Set up compiler-independent flags and build directories for the various build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dx86")          # TODO Need to change this for different CPU architectures


#########################################################################################
# Detect the compiler and set compiler-specific flags etc. TODO This isn't finished
message("Checking compiler")

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
set(GDAL_CONFIG_PREFER_PATH "$ENV{GDAL_HOME}/bin" CACHE STRING "preferred path to GDAL (gdal-config)")
set(GDAL_CONFIG_PREFER_FWTOOLS_PATH "$ENV{FWTOOLS_HOME}/bin_safe" CACHE STRING "preferred path to GDAL (gdal-config) from FWTools")
find_program(GDAL_CONFIG gdal-config
   ${GDAL_CONFIG_PREFER_PATH}
   ${GDAL_CONFIG_PREFER_FWTOOLS_PATH}
   /usr/local/bin/
   /usr/bin/
)
#message("DBG GDAL_CONFIG ${GDAL_CONFIG}")

if (GDAL_CONFIG)
   # Extract the GDAL version
   exec_program(${GDAL_CONFIG}
      ARGS --version
      OUTPUT_VARIABLE GDAL_VERSION
   )
   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" GDAL_VERSION_MAJOR "${GDAL_VERSION}")
   string(REGEX REPLACE "([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" GDAL_VERSION_MINOR "${GDAL_VERSION}")

   message ("GDAL is version ${GDAL_VERSION}")
#message("DBG GDAL_VERSION_MAJOR ${GDAL_VERSION_MAJOR}")
#message("DBG GDAL_VERSION_MINOR ${GDAL_VERSION_MINOR}")

   # Now check the GDAL version, and copy the appropriate files into the source folder
   if (GDAL_VERSION_MAJOR GREATER 1)
      if (UNIX OR CYGWIN)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_vector_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_vector.cpp")
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_utils_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_utils.cpp")
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_raster_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_raster.cpp")
      else (UNIX OR CYGWIN)
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_vector_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_vector.cpp")
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_utils_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_utils.cpp")
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_raster_gdal_version_2.cpp" "${CMAKE_SOURCE_DIR}/gis_raster.cpp")
      endif (UNIX OR CYGWIN)
   else (GDAL_VERSION_MAJOR GREATER 1)
      if (UNIX OR CYGWIN)
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_vector_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_vector.cpp")
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_utils_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_utils.cpp")
         execute_process(COMMAND cp -f "${CMAKE_SOURCE_DIR}/gdal_versions/gis_raster_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_raster.cpp")
      else (UNIX OR CYGWIN)
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_vector_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_vector.cpp")
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_utils_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_utils.cpp")
         execute_process(COMMAND copy "${CMAKE_SOURCE_DIR}/gdal_versions/gis_raster_gdal_version_1.cpp" "${CMAKE_SOURCE_DIR}/gis_raster.cpp")
      endif (UNIX OR CYGWIN)
   endif (GDAL_VERSION_MAJOR GREATER 1)

endif (GDAL_CONFIG)

#########################################################################################
# Sort out CShore for Linux and Windows
set(CSHORE_DIR "${CME_INSTALL_DIR}/cshore")
#message(STATUS "CSHORE_DIR=${CSHORE_DIR}")
if (WIN32)
   execute_process(COMMAND copy "${CSHORE_DIR}/cshore_usace_win.out" "${CSHORE_DIR}/cshore.out")
endif (WIN32)
if (UNIX)
   execute_process(COMMAND cp -f "${CSHORE_DIR}/cshore_usace_linux.out" "${CSHORE_DIR}/cshore.out")
endif (UNIX)


#########################################################################################
# Finally, the important bits
message("Writing Makefile")
include_directories(SYSTEM ${CMAKE_INCLUDE_PATH})
include_directories(${CME_SOURCE_DIR})
add_executable(${CME_EXECUTABLE} ${CME_SOURCE_FILES})
target_link_libraries(${CME_EXECUTABLE} ${LIBS})
install(TARGETS ${CME_EXECUTABLE} RUNTIME DESTINATION ${CME_INSTALL_DIR})

# Tell the user what has happened
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
#message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
#message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
#message(STATUS "LIBS=${LIBS}")
message(STATUS "CME_INSTALL_DIR=${CME_INSTALL_DIR}")
#message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")

# Some extra messages
if (CMAKE_BUILD_TYPE MATCHES Callgrind)
   message("")
   message("When the build has finished, use valgrind/callgrid as follows:")
   message("")
   message("valgrind --tool=callgrind ./cme")
   message("")
   message("Then run:")
   message("callgrind_annotate --auto=yes callgrind.out.XXXX > ./profile/callgrind/callgrind.txt")
   message("where XXXX is the number of the callgrind.out.XXXX that was produced by vallgrind")
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
