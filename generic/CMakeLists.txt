cmake_minimum_required(VERSION 2.6)
project(generic C CXX)

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
add_executable(extract2d extract2d.cpp)
target_link_libraries(extract2d ${PHYPP_LIBRARIES})
install(TARGETS extract2d DESTINATION .)

add_executable(get_fluxes get_fluxes.cpp)
target_link_libraries(get_fluxes ${PHYPP_LIBRARIES})
install(TARGETS get_fluxes DESTINATION .)

add_executable(flux_rescale flux_rescale.cpp)
target_link_libraries(flux_rescale ${PHYPP_LIBRARIES})
install(TARGETS flux_rescale DESTINATION .)

add_executable(stack1d stack1d.cpp)
target_link_libraries(stack1d ${PHYPP_LIBRARIES})
install(TARGETS stack1d DESTINATION .)

add_executable(stack2d stack2d.cpp)
target_link_libraries(stack2d ${PHYPP_LIBRARIES})
install(TARGETS stack2d DESTINATION .)

add_executable(fit_star_spectrum fit_star_spectrum.cpp)
target_link_libraries(fit_star_spectrum ${PHYPP_LIBRARIES})
install(TARGETS fit_star_spectrum DESTINATION .)
