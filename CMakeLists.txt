cmake_minimum_required(VERSION 3.14)
project(try)

set(CMAKE_CXX_STANDARD 14)

add_executable(try main.cpp functions.cpp constants.h functions.h)