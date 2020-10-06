cmake_minimum_required(VERSION 3.12)

find_package( Discregrid REQUIRED )

add_executable(GenerateDensityMap
	main.cpp
	gauss_quadrature.hpp
	gauss_quadrature.cpp
	sph_kernel.hpp
)

target_link_libraries(GenerateDensityMap
    PUBLIC Discregrid::Discregrid
)

target_include_directories(GenerateDensityMap
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../extern>)

target_compile_definitions(GenerateDensityMap PRIVATE -D_USE_MATH_DEFINES)

set_target_properties(GenerateDensityMap PROPERTIES FOLDER Cmd)

install(TARGETS GenerateDensityMap
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

