cmake_minimum_required(VERSION 3.13)
project(SSE_PTHREADS_MPI C)

set(CMAKE_C_STANDARD 99)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
SET(CMAKE_C_COMPILER mpicc)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(SSE_PTHREADS_MPI main.c)