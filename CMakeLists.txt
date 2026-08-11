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

# apfelxx (no CMake config file shipped, so find headers/lib manually)
find_path(APFELXX_INCLUDE_DIR apfel/apfelxx.h REQUIRED)
find_library(APFELXX_LIBRARY apfelxx REQUIRED)
add_library(apfelxx::apfelxx INTERFACE IMPORTED GLOBAL)
set_target_properties(apfelxx::apfelxx PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${APFELXX_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${APFELXX_LIBRARY}"
)

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
   execute_process(COMMAND ${LHAPDF_CONFIG} --cflags OUTPUT_VARIABLE LHAPDF_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND ${LHAPDF_CONFIG} --libs OUTPUT_VARIABLE LHAPDF_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
   set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
else(LHAPDF_CONFIG)
  message(FATAL_ERROR "LHAPDF not found!")
endif(LHAPDF_CONFIG)

# Wrap the lhapdf-config output in a proper INTERFACE target too,
# so it can be linked/propagated the same way as apfelxx::apfelxx
add_library(LHAPDF::LHAPDF INTERFACE IMPORTED GLOBAL)
separate_arguments(LHAPDF_CXX_FLAGS_LIST UNIX_COMMAND "${LHAPDF_CXX_FLAGS}")
separate_arguments(LHAPDF_LIBRARIES_LIST UNIX_COMMAND "${LHAPDF_LIBRARIES}")
set_target_properties(LHAPDF::LHAPDF PROPERTIES
    INTERFACE_COMPILE_OPTIONS "${LHAPDF_CXX_FLAGS_LIST}"
    INTERFACE_LINK_LIBRARIES "${LHAPDF_LIBRARIES_LIST}"
)

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

    target_link_libraries(QuarkoniumTMD PUBLIC apfelxx::apfelxx LHAPDF::LHAPDF)
endif()

# ============================
# PROGRAMS
# ============================

add_subdirectory(run)
