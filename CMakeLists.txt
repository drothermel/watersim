cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(final)

###############################################
# GTest - http://code.google.com/p/googletest/
###############################################
add_subdirectory(gtest170)
enable_testing(true)

###############################################
# Water
###############################################
include_directories(src gtest170/include gtest170)
add_executable(unit-tests unit-tests.cpp
							Water.cpp)

target_link_libraries(unit-tests gtest gtest_main)
add_test(
	NAME unit-tests 
	COMMAND unit-tests)

set(GCC_COMPILE_FLAGS "-lm")
add_definitions($(GCC_COMPILE_FLAGS))