cmake_minimum_required(VERSION 3.19)
project(seqMap)

set(CMAKE_CXX_STANDARD 14)

#add_library(seqMap SHARED library.cpp library.h)

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

aux_source_directory(./src  SEQMAP_SRC)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(test test.cpp ${SEQMAP_SRC} ./include)
#target_link_libraries(countSupport -lCbcSolver -lCbc -lOsi -lOsiCbc -lCgl -lOsiClp -lClp -lCoinUtils -llapack -lbz2 -lz -lhts)
