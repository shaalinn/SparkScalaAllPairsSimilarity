# Install script for directory: /home/shalin/Downloads/l2ap/build/GKlib

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/shalin/Downloads/l2ap/build/GKlib/libGKlib.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkblas.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkutils.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_defs.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_arch.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/GKlib.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkmemory.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/ms_stdint.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_macros.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkrandom.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gkregex.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/ms_inttypes.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_types.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkpqueue.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_proto.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mksort.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_mkpqueue2.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_externs.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_struct.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/gk_getopt.h"
    "/home/shalin/Downloads/l2ap/build/GKlib/ms_stat.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/shalin/Downloads/l2ap/build/GKlib/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/shalin/Downloads/l2ap/build/GKlib/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
