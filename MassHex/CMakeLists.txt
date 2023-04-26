cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors -g -Wall")

add_executable(${PROJECT_NAME} main.cpp)

if (DEFINED OPTIMIZE_SUM_FACTORIZATION)
    message( STATUS "OPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE_SUM_FACTORIZATION}")
    add_compile_definitions(OPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE_SUM_FACTORIZATION})
else()
    message( FATAL_ERROR "OPTIMIZE_SUM_FACTORIZATION should be defined at compitle time.")
endif()

if (DEFINED PRECISION)
    message( STATUS "PRECISION=${PRECISION}")
    add_compile_definitions(PRECISION=${PRECISION})
else()
    message( FATAL_ERROR "PRECISION should be defined at compitle time.")
endif()

if (DEFINED BATCH_SIZE)
    message( STATUS "BATCH_SIZE=${BATCH_SIZE}")
    add_compile_definitions(BATCH_SIZE=${BATCH_SIZE})
else()
    message( FATAL_ERROR "BATCH_SIZE should be defined at compitle time.")
endif()

if(${DEGREE} GREATER_EQUAL 1 AND ${DEGREE} LESS_EQUAL 15)
    message( STATUS "DEGREE=${DEGREE}")
    add_compile_definitions(DEGREE=${DEGREE})
else()
    message(FATAL_ERROR "Degree >= 1 and Degree <= 15")
endif()

find_package(MPI 3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)

message(${CMAKE_CXX_FLAGS})