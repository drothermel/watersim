######### Unit Testing ##########

Currently I am building gtest inside of my project which seems incredibly
unnecessary, but after spending an hour trying to figure out how to just
link the libraries in my make file, I'm finally giving up and just using
the way that works.

(For the future when I try to fix this again:
	I have placed the library files [libgtest.a and libgtest_main.a]
	in the /usr/lib/ directory and tried both have CMake find the
	package and to directly include /usr/lib/ as a link_directory and
	then to link_libraries(unit-tests gtest gtest_main) which is 
	supposedly looking in link_directories for libgtest.a and
	libgtest_main.a.  The reason I'm thinking this isn't happening
	is that the gtest/gtest.h include statement in the unittest.cpp
	file will not compile because it cannot be found.  No idea.
)
