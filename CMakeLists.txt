cmake_minimum_required(VERSION 3.6)
project(PCM_R_tree)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(pcm-r-tree ${SOURCE_FILES})