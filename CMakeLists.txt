cmake_minimum_required(VERSION 3.10)
project(seqGraph VERSION 1.0.1 DESCRIPTION "seq graph")

set(CMAKE_CXX_STANDARD 14)
set(DEFAULT_BUILD_TYPE DEBUG)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
include(GNUInstallDirs)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(SOURCE_FILES src/Edge.cpp
        src/EndPoint.cpp
        src/Exception.cpp
        src/Graph.cpp
        src/Junction.cpp
        src/Vertex.cpp
        src/Weight.cpp
        src/matching.cpp include/cxxopts.hpp)

add_executable(matching main.cpp ${SOURCE_FILES} include/util.h)
