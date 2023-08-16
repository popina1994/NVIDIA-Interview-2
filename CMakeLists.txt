cmake_minimum_required(VERSION 3.25)
project(NVIDIA-Interview-2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_STANDARD 20)

add_executable(untitled main.cpp)

#add_subdirectory(Google-Test)
