cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -DNDEBUG")

add_executable(${PROJECT_NAME} main.cpp)

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

if(${DEGREE} GREATER_EQUAL 1 AND ${DEGREE} LESS_EQUAL 8)
    add_compile_definitions(DEGREE=${DEGREE})
else()
    message(FATAL_ERROR "message to display")
endif()


find_package(MPI 3 REQUIRED)
find_package(likwid REQUIRED)

set(BLA_VENDOR FLAME)
find_package(BLAS REQUIRED)


target_include_directories(${PROJECT_NAME} PUBLIC $ENV{BLIS_DIR}/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${LIKWID_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX likwid::likwid BLAS::BLAS)

message(${CMAKE_CXX_FLAGS})