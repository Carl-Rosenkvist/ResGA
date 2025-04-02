cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SMASHER)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


# Tell cmake where to find modules (also append ../../cmake for testing)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# We can assume that SMASH_DIR is set here as it is required by FindSMASH.cmake
add_definitions("-DSMASH_TOP_LEVEL_DIR=\"$ENV{SMASH}\"")

set(SOURCES
src/main.cc
src/Decays.cc
src/MonteCarloEmulator.cc
)
set(HEADERS
    include/Decays.h
    include/MonteCarloEmulator.h
)



add_executable(run ${SOURCES} ${HEADERS})

target_include_directories(run PRIVATE include)
# Set the relevant generic compiler flags (optimisation + warnings) Try to be close to what SMASH
# uses
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -pthread"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -std=c++17")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
find_package(SMASH)
if(SMASH_FOUND)
    message(STATUS "SMASH libraries = ${SMASH_LIBRARIES}")
    include_directories(${SMASH_INCLUDE_DIR})
    target_link_libraries(run ${SMASH_LIBRARIES})
else()
    message(FATAL_ERROR "SMASH libraries not found!")
endif()
