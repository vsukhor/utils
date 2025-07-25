include(FetchContent)
set(FETCHCONTENT_QUIET off)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

include(GoogleTest)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

add_definitions(-DUSE_UTILS_ASSERT)

add_executable(unittests
	test_matrix33.cpp
	test_matrix44.cpp
	test_misc.cpp
	test_msgr.cpp
	test_parameter.cpp
	test_reader.cpp
)

include_directories(unittests PUBLIC
	../utils/biochemical
	../utils/arrays
	../utils/graph
	../utils/common
	../utils/config
	../utils/random
	../utils/volumetric_data
)

target_link_libraries(unittests PRIVATE $<TARGET_FILE:utils>)
target_link_libraries(unittests PRIVATE gtest gmock)
target_link_libraries(unittests PRIVATE gtest_main)

gtest_discover_tests(unittests)

add_custom_command(
    TARGET unittests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:unittests> ${CMAKE_CURRENT_SOURCE_DIR}
)
