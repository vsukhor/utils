cmake_minimum_required(VERSION 3.15)

project(utils
    VERSION 1.3.0
    DESCRIPTION "A cooection of misc code for common use in other repos."
    LANGUAGES CXX)

add_compile_options(-O3 -msse4 -fPIC)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

add_library(utils STATIC
            utils/biochemical/pdb.cpp)

target_include_directories(utils PUBLIC utils/biochemical
                                        utils/arrays
                                        utils/graph
                                        utils/common
                                        utils/config
                                        utils/random
                                        utils/volumetric_data)

target_compile_features(utils PRIVATE cxx_std_20)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#    option(EXPORT_COMPILE_COMMANDS "Compile commands are visible" ON)
    set(CMAKE_VERBOSE_MAKEFILE OFF)

    # Docs with Doxygen
    option(BUILD_DOC "Build documentation" ON)
    if(BUILD_DOC)
        add_subdirectory(docs)
    endif()

    option(BUILD_TESTS "Build test files" ON)
    if(BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()

endif()
