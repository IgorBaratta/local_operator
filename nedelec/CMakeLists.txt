cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME benchmark)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-DNDEBUG ${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_C_FLAGS "-DNDEBUG ${CMAKE_C_FLAGS} -Wall")

include(FindPackageHandleStandardArgs)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

find_package(Basix REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)

add_executable(${PROJECT_NAME} main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/nedelec.c)

set(UFC_INCLUDE_DIRS $ENV{UFC_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${UFC_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${xtensor_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${xtensor_blas_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${xtensor_blas_LIBRARY})
target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})

message(${CMAKE_CXX_FLAGS})
message(${CMAKE_C_FLAGS})