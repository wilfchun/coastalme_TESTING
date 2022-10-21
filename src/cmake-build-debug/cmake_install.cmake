# Install script for directory: /home/wilfc/github/coastalme_TESTING/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/sbin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme"
         RPATH "$ORIGIN/src/lib:$$ORIGIN/src/lib:/home/wilfc/github/coastalme_TESTING/src/lib:/usr/local/lib")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/wilfc/github/coastalme_TESTING/src/../cme")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/wilfc/github/coastalme_TESTING/src/.." TYPE EXECUTABLE FILES "/home/wilfc/github/coastalme_TESTING/src/cmake-build-debug/Debug/cme")
  if(EXISTS "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme"
         OLD_RPATH "/home/wilfc/github/coastalme_TESTING/src/lib::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "$ORIGIN/src/lib:$$ORIGIN/src/lib:/home/wilfc/github/coastalme_TESTING/src/lib:/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/sbin/strip" "$ENV{DESTDIR}/home/wilfc/github/coastalme_TESTING/src/../cme")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/wilfc/github/coastalme_TESTING/src/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
