cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

add_executable(${PROJECT_NAME} matrix.cpp ${CMAKE_CURRENT_SOURCE_DIR}/problem.c)
add_executable(${PROJECT_NAME}-mf mf.cpp ${CMAKE_CURRENT_SOURCE_DIR}/problem.c)

set(UFC_INCLUDE_DIRS $ENV{UFC_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${UFC_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME}-mf SYSTEM PUBLIC ${UFC_INCLUDE_DIRS})

message(${CMAKE_CXX_FLAGS})
message(${CMAKE_C_FLAGS})