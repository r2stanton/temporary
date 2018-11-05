cmake_minimum_required(VERSION 3.9)
project(temporary)

set(CMAKE_CXX_STANDARD 11)

add_executable(temporary main.cpp Ray.cpp Ray.h Coords.cpp Coords.h derivs.cpp)