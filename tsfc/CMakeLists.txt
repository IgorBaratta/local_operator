cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

add_executable(${PROJECT_NAME} main.cpp)

find_package(MPI 3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)

message(${CMAKE_CXX_FLAGS})