#!/bin/sh

# Change this to select build type
buildtype=Debug
#buildtype=Release
#buildtype=RelWithDebInfo
#buildtype=MinSizeRel
#buildtype=gcov
#buildtype=Callgrind

# Change this to select CShore input/output method
#cshoreinout=FILE
cshoreinout=ARG
#cshoreinout=BOTH

echo "CoastalME: starting CMake for Linux (using gcc, $buildtype build, CShore input/output method=$cshoreinout)"

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_INOUT=$cshoreinout .

echo
echo "Finished CMake for Linux (using gcc, $buildtype build, CShore input/output method=$cshoreinout)"
