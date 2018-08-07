cmake_minimum_required(VERSION 2.6)
project(kmos C CXX)

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

add_executable(reduce reduce.cpp)
target_link_libraries(reduce ${PHYPP_LIBRARIES})
install(TARGETS reduce DESTINATION .)

add_executable(fill_nan fill_nan.cpp)
target_link_libraries(fill_nan ${PHYPP_LIBRARIES})
install(TARGETS fill_nan DESTINATION .)

add_executable(extract_ifu extract_ifu.cpp)
target_link_libraries(extract_ifu ${PHYPP_LIBRARIES})
install(TARGETS extract_ifu DESTINATION .)

add_executable(make_profile make_profile.cpp)
target_link_libraries(make_profile ${PHYPP_LIBRARIES})
install(TARGETS make_profile DESTINATION .)

add_executable(make_shifts make_shifts.cpp)
target_link_libraries(make_shifts ${PHYPP_LIBRARIES})
install(TARGETS make_shifts DESTINATION .)

add_executable(median_sub median_sub.cpp)
target_link_libraries(median_sub ${PHYPP_LIBRARIES})
install(TARGETS median_sub DESTINATION .)

add_executable(merge_cs_nn merge_cs_nn.cpp)
target_link_libraries(merge_cs_nn ${PHYPP_LIBRARIES})
install(TARGETS merge_cs_nn DESTINATION .)

add_executable(stack_cubes stack_cubes.cpp)
target_link_libraries(stack_cubes ${PHYPP_LIBRARIES})
install(TARGETS stack_cubes DESTINATION .)

add_executable(cube2slit cube2slit.cpp)
target_link_libraries(cube2slit ${PHYPP_LIBRARIES})
install(TARGETS cube2slit DESTINATION .)

add_executable(contsub contsub.cpp)
target_link_libraries(contsub ${PHYPP_LIBRARIES})
install(TARGETS contsub DESTINATION .)

add_executable(cdetect cdetect.cpp)
target_link_libraries(cdetect ${PHYPP_LIBRARIES})
install(TARGETS cdetect DESTINATION .)

add_executable(clinefit clinefit.cpp)
target_link_libraries(clinefit ${PHYPP_LIBRARIES})
install(TARGETS clinefit DESTINATION .)

add_executable(multispecfit multispecfit.cpp)
target_link_libraries(multispecfit ${PHYPP_LIBRARIES})
install(TARGETS multispecfit DESTINATION .)

add_executable(multispecfit_sof multispecfit_sof.cpp)
target_link_libraries(multispecfit_sof ${PHYPP_LIBRARIES})
install(TARGETS multispecfit_sof DESTINATION .)

