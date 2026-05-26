# ============================
# BASIC DEFINITIONS
# ============================

cmake_minimum_required(VERSION 3.10)
project(QuarkoniumTMD CXX)

set(QuarkoniumTMD_VERSION 1.0.0)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ============================
# DEPENDENCIES
# ============================

find_package(apfelxx REQUIRED)

# ============================
# MAIN LIBRARY
# ============================

file(GLOB_RECURSE source_files src/*.cc)

if(source_files)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
    add_library(QuarkoniumTMD SHARED ${source_files})

    target_include_directories(QuarkoniumTMD
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
            $<INSTALL_INTERFACE:include>
    )

    target_compile_options(QuarkoniumTMD PRIVATE -Wall -Wextra -O3)

    target_link_libraries(QuarkoniumTMD PUBLIC apfelxx::apfelxx)
endif()

# ============================
# PROGRAMS
# ============================

add_subdirectory(run)
