cmake_minimum_required(VERSION 2.6)
project(xshooter C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_GSL ON)
set(NO_LAPACK ON)

find_package(phypp REQUIRED)

# Include dependencies
message(STATUS ${PHYPP_INCLUDE_DIRS})
include_directories(${PHYPP_INCLUDE_DIRS})

# Build tools
add_executable(rename rename.cpp)
target_link_libraries(rename ${PHYPP_LIBRARIES})
install(TARGETS rename DESTINATION .)

add_executable(correct_dark correct_dark.cpp)
target_link_libraries(correct_dark ${PHYPP_LIBRARIES})
install(TARGETS correct_dark DESTINATION .)

add_executable(get_telluric get_telluric.cpp)
target_link_libraries(get_telluric ${PHYPP_LIBRARIES})
install(TARGETS get_telluric DESTINATION .)
