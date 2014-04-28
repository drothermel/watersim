cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(final)

add_executable(water_main Main.cpp
				Water.cpp)

set(GCC_COMPILE_FLAGS "-lm")
add_definitions($(GCC_COMPILE_FLAGS))
