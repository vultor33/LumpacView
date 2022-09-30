cmake_minimum_required ( VERSION 3.3 )
project(stereoisomer)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#headers
include_directories(src)

#cpps
file(GLOB CPP_FILES src/*.cpp src/*/*.cpp)
add_executable(stereoisomer ${CPP_FILES})

target_link_libraries(stereoisomer ${BUILD_LIBS})

