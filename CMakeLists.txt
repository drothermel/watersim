#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(final)

#This tells CMake to fib.c and name it fibonacci
add_executable(water_main Main.cpp 
					Water.cpp)

set(GCC_COMPILE_FLAGS "-lm")
add_definitions($(GCC_COMPILE_FLAGS))