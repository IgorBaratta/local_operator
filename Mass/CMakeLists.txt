cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

add_executable(${PROJECT_NAME} main.cpp)
add_compile_definitions(PRECISION=${PRECISION})
add_compile_definitions(BATCH_SIZE=${BATCH_SIZE})

if(${DEGREE} GREATER_EQUAL 1 AND ${DEGREE} LESS_EQUAL 8)
    add_compile_definitions(DEGREE=${DEGREE})
else()
    message(FATAL_ERROR "message to display")
endif()


if(${BATCH_SIZE}>1)
    add_compile_definitions(USE_VECTOR_EXTENSIONS=True)
endif()

find_package(MPI 3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)

message(${CMAKE_CXX_FLAGS})