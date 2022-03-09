# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.19.8)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(liquid-populate NONE)

include(ExternalProject)
ExternalProject_Add(liquid-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "git://github.com/jgaeddert/liquid-dsp.git" "GIT_TAG" "master"
                    SOURCE_DIR          "/home/moskalev-i/CLionProjects/TestF/PodDSP/liquid-dsp"
                    BINARY_DIR          "/home/moskalev-i/CLionProjects/TestF/PodDSP/cmake-build-debug/_deps/liquid-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


