set(SubProject ${CMAKE_PROJECT_NAME}_lib) # setting the library name

# specifying the project name
project(${SubProject})

# specifying the directory where header files are located
# include_directories("/headers")

# selecting all the files with .cpp extension and assigning them to a variable 'SOURCES'
file(GLOB SOURCES "sources/*.*")
file(GLOB HEADERS "headers/*.*")

# adding to the library all the files that are selected and assigned to the variable name 'SOURCES'
add_library(
    ${SubProject} # Name of the project that is extracted from the root directory CMakeLists.txt
    ${SOURCES} # 
    ${HEADERS}
    )

install(TARGETS ${SubProject} DESTINATION lib)