set(BINARY ${CMAKE_PROJECT_NAME})

# set C++ compiler flags 
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

# check if openMP compilers are available and link them
find_package(OpenMP 4.0 REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)

file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/sources/*.cpp")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/src/headers/*.h")

add_executable(${BINARY} main.cpp)

add_library(${BINARY}_lib ${SOURCES} ${HEADERS})
message(STATUS "SOURCES = ${HEADERS}")
target_link_libraries(${BINARY} ${BINARY}_lib)