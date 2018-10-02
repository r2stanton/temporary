cmake_minimum_required(VERSION 3.9)
project(lensing)

set(CMAKE_CXX_STANDARD 11)

add_executable(lensing main.cpp Ray.cpp Ray.h Coords.cpp Coords.h)