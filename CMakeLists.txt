# Copyright © 2023 Zen Shawn. All rights reserved.
# Created at 12:49:24, July 10, 2023

cmake_minimum_required(VERSION 3.14)

option(BUILD_SHARED_LIBS "build shared libs" OFF)
option(BUILD_EXAMPLES "build examples" ON)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "build type, Release/Debug/MinSizeRel/RelWithDebInfo")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(iso-contour VERSION 0.1.0 DESCRIPTION "Introduction about iso-contour")

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra)
endif()

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/out" CACHE PATH "install prefix")
include(GNUInstallDirs)

find_package(spdlog CONFIG REQUIRED)
find_package(CGAL CONFIG REQUIRED)

add_subdirectory(src)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

include(WriteBasicConfigVersionFile)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
)
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
)