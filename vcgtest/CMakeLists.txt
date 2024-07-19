cmake_minimum_required(VERSION 3.13)
project ( vcgtest )

include_directories(${CMAKE_SOURCE_DIR}/vcglib)

add_executable(vcgtest main.cpp)

target_link_libraries(
	vcgtest
	PUBLIC
		vcglib
	)
