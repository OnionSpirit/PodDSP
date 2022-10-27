cmake_minimum_required(VERSION 3.19)
project(PodDSP_test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(GTEST_CATCH_EXCEPTIONS 0)
set(GTEST_BREAK_ON_FAILURE 1)

include(FetchContent)

FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG main
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/google/googletest)
FetchContent_MakeAvailable(googletest)

add_executable(${PROJECT_NAME}
        main.cpp
        All_Tests.cpp
        )

#
add_subdirectory(vssdsp)
target_include_directories(${PROJECT_NAME} PRIVATE vssdsp/include)
target_link_libraries(${PROJECT_NAME} PRIVATE vssdsp)
#

target_link_libraries(${PROJECT_NAME} PRIVATE gtest)
