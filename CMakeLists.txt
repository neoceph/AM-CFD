cmake_minimum_required(VERSION 3.10)

set(Project AM-CFD)
set(ProjectLib AM-CFD_lib)

project(${Project})

# set C++ compiler flags 
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

# check if openMP compilers are available and link them
find_package(OpenMP 4.0 REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

enable_testing()

include_directories("${CMAKE_SOURCE_DIR}/src/headers")

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(libs/googletest)