cmake_minimum_required(VERSION 3.10)

# Declare configuration file to embed supplied mesh files.
set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/resources/")
configure_file("../resource_path.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/resource_path.hpp")

find_package(Discregrid REQUIRED)

add_executable(GenerateSDF
	main.cpp
	resource_path.hpp
)

target_link_libraries(GenerateSDF
    PUBLIC Discregrid::Discregrid
)
target_include_directories(GenerateSDF
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../extern>)

set_target_properties(GenerateSDF PROPERTIES FOLDER Cmd)

install(TARGETS GenerateSDF
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

