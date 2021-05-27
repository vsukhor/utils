cmake_minimum_required(VERSION 3.15)

project(utils
    VERSION 1.1.0
    DESCRIPTION "Utility library."
    LANGUAGES CXX)

    add_compile_options(-O3 -msse4 -fPIC)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

    add_library(utils STATIC
                utils/biochemical/pdb.cpp
                utils/common/threads.cpp
                utils/common/misc.cpp
                utils/common/msgr.cpp
                utils/common/exceptions.cpp)

    target_compile_features(utils PRIVATE cxx_std_17)

    target_include_directories(utils PUBLIC utils/biochemical
                                            utils/arrays
                                            utils/graph
                                            utils/common
                                            utils/config
                                            utils/random
                                            utils/volumetric_data)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#    option(EXPORT_COMPILE_COMMANDS "Compile commands are visible" ON)
    set(CMAKE_VERBOSE_MAKEFILE ON)

    # Docs with Doxygen
    option(BUILD_DOC "Build documentation" ON)
    if(BUILD_DOC)
        find_package(Doxygen)
        if(Doxygen_FOUND)
            # set input and output files
            set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/doxyfile.in)
            set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

            # request to configure the file
            configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
            message("Doxygen build started: ${CMAKE_CURRENT_SOURCE_DIR}")

            # note the option ALL which allows to build the docs together with
            # the application
            add_custom_target(doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM)
        else()
            message(STATUS "Doxygen not found, not building docs")
        endif()
    endif()
endif()

# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR UTILS_BUILD_TESTING) AND BUILD_TESTING)
    add_subdirectory(tests)
endif()


