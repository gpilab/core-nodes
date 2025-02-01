# Instructions to build the project
# =================================
# Prerequisites:
# 1. CMake version 3.2 to 3.18
# 2. Python with development headers
# 3. pybind11 library

# Steps to build:
# 1. Navigate to the project directory:
#    cd path_to_core_directory/core/spiral
#
# 2. Create a build directory:
#    mkdir build && cd build
#
# 3. Run CMake to configure the project:
#    cmake ..
#
# 4. Build the project:
#    cmake --build .
#
# 5. The built Python module will be copied to one folder up from the build directory.
# Author: Guru Krishnamoorthy
# Date: 2024-06-01

cmake_minimum_required(VERSION 3.2...3.18)

project(WPGEN VERSION 1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 11)
if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64e")
        set(CMAKE_OSX_ARCHITECTURES arm64 arm64e)
    else()
        set(CMAKE_OSX_ARCHITECTURES x86_64)
    endif()
endif()

# Set visibility settings
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# Find pybind11 package
# =============================
# Locate Dependencies
# =============================

# Locate Python
find_package(Python REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Found Python: ${Python_EXECUTABLE} (${Python_VERSION})")

# Locate Pybind11
find_package(pybind11 REQUIRED HINTS 
    ${Python_SITEARCH} # This variable points to the Python site-packages for the current Python version
)

# Create a library target for wpgen
add_library(wpgen STATIC cpp/WPGen.cpp cpp/WPCompose.cpp cpp/SpiralOrder.cpp)
set_target_properties(wpgen PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set the source files for the Python module
set(SOURCES cpp/binder.cpp)

# Create the Python module target
pybind11_add_module(whirledpeas MODULE ${SOURCES})
target_link_libraries(whirledpeas PRIVATE wpgen)

# Copy the built library to one folder up
add_custom_command(TARGET whirledpeas POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:whirledpeas> ${CMAKE_CURRENT_BINARY_DIR}/../)