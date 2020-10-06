cmake_minimum_required(VERSION 3.12)


# MSVC runtime library flags are selected by an abstraction.
cmake_policy(SET CMP0091 NEW)

project(DiscregridAll)

option(DISCREGRID_USE_MSVC_STATIC_RUNTIME "Use MSVC static runtime library" OFF)

include(GNUInstallDirs)

if(WIN32)
    if(DISCREGRID_USE_MSVC_STATIC_RUNTIME)
        message(STATUS "Using MSVC static runtime")
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

        set(CMAKE_INSTALL_UCRT_LIBRARIES FALSE)
        set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
    else()
        set(CMAKE_INSTALL_UCRT_LIBRARIES TRUE)
        set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
    endif()

    set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP FALSE)
    set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
    include(InstallRequiredSystemLibraries)
endif()


list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Require C++11 compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if((CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT))
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Prefix prepended to install directories" FORCE)
    message(STATUS "Setting default CMAKE_INSTALL_PREFIX to: ${CMAKE_INSTALL_PREFIX}")
endif()

option(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" OFF)
if (BUILD_AS_SHARED_LIBS)
	set(BUILD_AS_SHARED_LIBS ON)
	if(WIN32)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON CACHE BOOL "Export all symbols")
	endif(WIN32)
endif (BUILD_AS_SHARED_LIBS)

set(CMAKE_DEBUG_POSTFIX "_d")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(discregrid)

option(BUILD_CMD_EXECUTABLE "Build command line executable" ON)
if(BUILD_CMD_EXECUTABLE)
	add_subdirectory(cmd)
endif(BUILD_CMD_EXECUTABLE)

