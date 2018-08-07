cmake_minimum_required(VERSION 2.6)
project(speclib-superbuild C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected (-DCMAKE_BUILD_TYPE=...)")
    message(STATUS "Default to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(ExternalProject)

# Set default install dir to "here"
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "No install directory selected (-DCMAKE_INSTALL_PREFIX=...)")
    message(STATUS "Default to ${PROJECT_SOURCE_DIR}/bin")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "default install path" FORCE)
endif()

# Configure and build cfitsio
ExternalProject_Add(cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    GIT_REPOSITORY https://github.com/HEASARC/cfitsio.git
    GIT_TAG ff6bdb0c2f891ae6709f6bb2378e24e005b7c92b
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DUSE_PTHREADS=ON
        -DBUILD_TOOLS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build WCSlib
ExternalProject_Add(wcslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/wcslib
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/wcslib
    URL ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/external/wcslib/configure
        --prefix=${PROJECT_SOURCE_DIR}/external
        --without-cfitsio --without-pgplot --disable-fortran --disable-utils
    BUILD_COMMAND make install && rm ${PROJECT_SOURCE_DIR}/external/lib/libwcs.so
    BUILD_IN_SOURCE 1)

# Configure and build phy++
ExternalProject_Add(phy++
    DEPENDS cfitsio wcslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/phypp
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/phypp
    GIT_REPOSITORY https://github.com/cschreib/phypp.git
    GIT_TAG 45f9a633792aa4ef6782e974180b8945a3fa94a7
    CMAKE_ARGS
        -DNO_REFLECTION=ON
        -DNO_FFTW=ON
        -DNO_LAPACK=ON
        -DNO_GSL=ON
        -DPHYPP_INPLACE_BUILD=ON
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build sub-projects
ExternalProject_Add(generic
    DEPENDS phy++
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/generic
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DPHYPP_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/generic
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/phypp/cmake
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

ExternalProject_Add(xshooter
    DEPENDS phy++
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/xshooter
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DPHYPP_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/xshooter
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/phypp/cmake
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

ExternalProject_Add(mosfire
    DEPENDS phy++
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/mosfire
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DPHYPP_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/mosfire
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/phypp/cmake
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

ExternalProject_Add(kmos
    DEPENDS phy++
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/kmos
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DPHYPP_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/kmos
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/phypp/cmake
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
