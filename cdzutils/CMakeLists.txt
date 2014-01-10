
cmake_minimum_required (VERSION 2.6)
set( CMAKE_VERBOSE_MAKEFILE off )

project (cdzutils)

add_definitions(-std=c99)

include_directories(include)
include_directories(src)

add_library(cdzutils 
    SHARED
    src/cdzutils.c 
    include/cdzutils.h
    )

add_executable(cdzutils-test 
    src/test.c 
    include/cdzutils.h
    )

target_link_libraries(cdzutils-test cdzutils)

