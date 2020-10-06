cmake_minimum_required(VERSION 3.12)

project(Discregrid)

set(HEADERS
	include/Discregrid/discrete_grid.hpp
	include/Discregrid/cubic_lagrange_discrete_grid.hpp
)

set(HEADERS_ACCELERATION

	include/Discregrid/acceleration/bounding_sphere_hierarchy.hpp
	include/Discregrid/acceleration/bounding_sphere.hpp
	include/Discregrid/acceleration/kd_tree.hpp
	include/Discregrid/acceleration/kd_tree.inl
)

set(HEADERS_MESH
	include/Discregrid/mesh/triangle_mesh.hpp
	include/Discregrid/mesh/entity_containers.hpp
	include/Discregrid/mesh/entity_iterators.hpp
	include/Discregrid/mesh/halfedge.hpp
)

set(HEADERS_GEOMETRY
	include/Discregrid/geometry/mesh_distance.hpp

)

set(HEADERS_GEOMETRY_PRIVATE
	src/geometry/point_triangle_distance.hpp
)

set(HEADERS_UTILITY
	include/Discregrid/utility/serialize.hpp
	include/Discregrid/utility/lru_cache.hpp

)

set(HEADERS_UTILITY_PRIVATE
	src/utility/timing.hpp
	src/utility/spinlock.hpp
)

set(SOURCES
	src/discrete_grid.cpp
	src/cubic_lagrange_discrete_grid.cpp
)

set(SOURCES_DATA
)

set(SOURCES_ACCELERATION
	src/acceleration/bounding_sphere_hierarchy.cpp
)

set(SOURCES_MESH
	src/mesh/entity_containers.cpp
	src/mesh/entity_iterators.cpp
	src/mesh/triangle_mesh.cpp
)

set(SOURCES_GEOMETRY
	src/geometry/mesh_distance.cpp
	src/geometry/point_triangle_distance.cpp
)

set(SOURCES_UTILITY
	src/utility/timing.cpp
)

macro(SOURCEGROUP name)
	string(TOLOWER ${name} name_lower)
	string(SUBSTRING ${name_lower} 0 1 FIRST_LETTER)
	string(TOUPPER ${FIRST_LETTER} FIRST_LETTER)
	string(REGEX REPLACE "^.(.*)" "${FIRST_LETTER}\\1" dname "${name_lower}")
	source_group("Source Files ${dname}" FILES ${SOURCES_${name}})
	source_group("Header Files ${dname}" FILES ${HEADERS_${name}})
endmacro()
SOURCEGROUP(MESH)
SOURCEGROUP(ACCELERATION)
SOURCEGROUP(DATA)
SOURCEGROUP(GEOMETRY)
SOURCEGROUP(GEOMETRY_PRIVATE)
SOURCEGROUP(UTILITY)
SOURCEGROUP(UTILITY_PRIVATE)

set(DISCREGRID_SOURCE_FILES
	${HEADERS}
	${SOURCES}
	${HEADERS_ACCELERATION}
	${SOURCES_ACCELERATION}
	${SOURCES_DATA}
	${HEADERS_MESH}
	${SOURCES_MESH}
	${HEADERS_GEOMETRY}
    ${HEADERS_GEOMETRY_PRIVATE}
	${SOURCES_GEOMETRY}
	${HEADERS_UTILITY}
    ${HEADERS_UTILITY_PRIVATE}
	${SOURCES_UTILITY}
)

set(DISCREGRID_PUBLIC_HEADERS
    ${HEADERS}
    ${HEADERS_ACCELERATION}
    ${HEADERS_MESH}
    ${HEADERS_GEOMETRY}
    ${HEADERS_UTILITY}
)

# Eigen library.
find_package(Eigen3 REQUIRED)

# OpenMP support.
find_package(OpenMP REQUIRED)

if(BUILD_AS_SHARED_LIBS)
	add_library(Discregrid SHARED ${DISCREGRID_SOURCE_FILES})
else()
	add_library(Discregrid ${DISCREGRID_SOURCE_FILES})
endif()

target_link_libraries(Discregrid PUBLIC
    Eigen3::Eigen
    OpenMP::OpenMP_CXX)

target_include_directories(Discregrid
    PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Discregrid>
)

target_include_directories(Discregrid
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

include(GNUInstallDirs)

install(TARGETS Discregrid
        EXPORT DiscregridTargets
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        RUNTIME DESTINATION  "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION  "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION  "${CMAKE_INSTALL_LIBDIR}"
)
foreach ( file ${DISCREGRID_PUBLIC_HEADERS} )
    get_filename_component( dir ${file} DIRECTORY )
    # use CMAKE_INSTALL_PREFIX instead of CMAKE_INSTALL_INCLUDEDIR since
    # each file has already a parent directory named "include"
    install( FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX}/${dir} )
endforeach()


include(CMakePackageConfigHelpers)

configure_package_config_file(DiscregridConfig.cmake.in
    ${CMAKE_BINARY_DIR}/Discregrid/DiscregridConfig.cmake
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Discregrid"
)

install(FILES ${CMAKE_BINARY_DIR}/Discregrid/DiscregridConfig.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Discregrid")

install(EXPORT DiscregridTargets
        NAMESPACE "Discregrid::"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Discregrid"
)

# Add an alias library for the projects which consumes Discregrid library directly
# ie without resorting to the DiscregridTargets file
add_library(Discregrid::Discregrid ALIAS Discregrid)

