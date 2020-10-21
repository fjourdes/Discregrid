cmake_minimum_required(VERSION 3.10)

find_package( Discregrid REQUIRED )

add_executable(DiscreteFieldToBitmap
	main.cpp
	bmp_file.hpp
	bmp_file.cpp
)

target_link_libraries(DiscreteFieldToBitmap PUBLIC Discregrid::Discregrid )

target_include_directories(DiscreteFieldToBitmap
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../extern>)

target_compile_options(DiscreteFieldToBitmap
    PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:-Wno-multichar>
    $<$<CXX_COMPILER_ID:MSVC>:/wd4996>
)

target_compile_definitions(DiscreteFieldToBitmap
    PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:-DNOMINMAX>
)

set_target_properties(DiscreteFieldToBitmap PROPERTIES FOLDER Cmd)

install(TARGETS DiscreteFieldToBitmap
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

