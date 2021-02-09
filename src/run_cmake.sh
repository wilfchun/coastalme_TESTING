#!/bin/sh

# Change this to select build type

buildtype=DEBUG
#buildtype=RELEASE
#buildtype=RELWITHDEBINFO     # Not yet implemented in CMakeLists.txt
#buildtype=MINSIZEREL         # Not yet implemented in CMakeLists.txt
#buildtype=GCOV
#buildtype=VALGRIND

# Change this to select the CShore library type
#cshorelibrary=STATIC
cshorelibrary=SHARED

# Change this to select CShore input/output method
#cshoreinout=FILE
cshoreinout=ARG
#cshoreinout=BOTH

echo "CoastalME: starting CMake for Linux (using gcc, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout .

echo
echo "Finished CMake for Linux (using gcc, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
