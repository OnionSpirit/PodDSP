cmake_minimum_required(VERSION 3.19)
project(poddsp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

include(FetchContent)

file(GLOB SRC source/*.cpp)
file(GLOB HDR include/*.h)

add_library(${PROJECT_NAME} STATIC ${SRC} ${HDR})

FetchContent_Declare(liquid
        GIT_REPOSITORY https://github.com/jgaeddert/liquid-dsp.git
        GIT_TAG master
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/liquid-dsp)
FetchContent_MakeAvailable(liquid)


target_link_libraries(${PROJECT_NAME} PRIVATE liquid)
