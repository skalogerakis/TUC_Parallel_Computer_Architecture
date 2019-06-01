cmake_minimum_required(VERSION 3.13)
project(SSE_PTHREADS C)

set(CMAKE_C_STANDARD 99)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(SSE_PTHREADS main.c)