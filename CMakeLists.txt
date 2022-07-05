cmake_minimum_required(VERSION 3.21)

# Read the JSON file. JSON file contains all information regarding project titles and named variables including developers.
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/CMakeSetup.json CMakeSetup_STRING)

# extracting project name and devleoper names from the json file.
string(JSON Project_NAME GET ${CMakeSetup_STRING} 0 ProjectName)
string(JSON Project_DEVELOPER GET ${CMakeSetup_STRING} 0 Developer)

# CPack and git versioning stuffs

# Appends the cmake/modules path to MAKE_MODULE_PATH variable.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Building a static library would require uncommenting the following two variable setting.
# set(BUILD_SHARED_LIBS OFF)
# set(CMAKE_EXE_LINKER_FLAGS "-static")


# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# Make a version file containing the current version from git.
if(EXISTS "${PROJECT_SOURCE_DIR}/versionForSourcePackage.cmake")
  include("${PROJECT_SOURCE_DIR}/versionForSourcePackage.cmake")
else()
  include(GetGitRevisionDescription)
  git_describe(VERSION --tags --dirty=-dirty)
endif()

#parse the version information into pieces.
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# setting project and project library variable names.
set(Project ${Project_NAME})
set(ProjectLib ${Project}_lib)

set(CMAKE_INSTALL_PREFIX "/${Project_NAME}_${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# settting the project name.
project(${Project} C CXX)

# specifying the C++ standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# check if openMP compilers are available and link them
find_package(OpenMP 4.0 REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# intializing and enabling test.
include(CTest)
enable_testing()

# specifying the header directory
include_directories(
    ${ProjectLib} "${CMAKE_SOURCE_DIR}/src/headers"
    )

# setting the directory of the source for doxygen documentation
set(PROJECT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src")

# specifying the directory where all the sources are located.
add_subdirectory(
    src
    )

add_subdirectory(
    test
    )

# adding git submodules chekout from the source

# Download all the submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# update submodules as needed
    option(GIT_SBUMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STSTUS "Submodule update")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT
            )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# check all the submodules
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extLibs/googletest/CMakeLists.txt")
    message(FATAL_ERROR "The googletest submodules was not downloaded! GIT_SUBMODULE is either turned off or failed. Please update submodules")
endif()

# setting up sphinx-build with readthedocs theme
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/usermanual)    
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CMakeSetup.json ${CMAKE_CURRENT_BINARY_DIR}/CMakeSetup.json
    COPYONLY
    )

message("Copying Sphinx setup files")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc/readthedocs/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/usermanual/)

message("Sphinx Document generation started")
add_custom_command(
    OUTPUT ../html/index.html
    COMMAND make html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/usermanual
    VERBATIM
    )

message("Sphinx PDF documentation generation started")
add_custom_target(
    HTML_MANUAL ALL
    COMMAND make latexpdf
    DEPENDS ../html/index.html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/usermanual
    )

# setting up doxygen and code documenation html and pdf
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # creating a /build/doc directory for generating doxyfile outputs in the separate build directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/api/latex)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)

    message(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    message(${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT})
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
    add_custom_target(
        docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
        )

    message("Doxygen Document generation started")
    add_custom_command(
        OUTPUT refman.tex
        COMMAND make docs
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        VERBATIM
        )

    message("Doxygen PDF generation started")
    add_custom_target(
        API_PDF ALL
        COMMAND make pdf
        DEPENDS refman.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/api/latex
        )

else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation!")
endif (DOXYGEN_FOUND)


# adding the googletest library separately as including the subdirectory with the previous command did not work!
add_subdirectory(
    extLibs/googletest EXCLUDE_FROM_ALL
    )

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")

# The following will transfer the version from git to the source package.
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${VERSION_SHORT}-Source")

# set(CPACK_INSTALL_COMMANDS "${CMAKE_COMMAND} -E make_directory \
#     ${PROJECT_BINARY_DIR}/_CPack_Packages/Linux-Source/${CPACK_SOURCE_GENERATOR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}/"
#     "${CMAKE_COMMAND} -E copy \
#     ${PROJECT_BINARY_DIR}/versionForSourcePackage.cmake \
#     ${PROJECT_BINARY_DIR}/_CPack_Packages/Linux-Source/${CPACK_SOURCE_GENERATOR}/${CPACK_SOURCE_PACKAGE_FILE_NAME}/")

# Exclude the build and .git directory from the source package.
set(CPACK_SOURCE_IGNORE_FILES "${PROJECT_SOURCE_DIR}/.git/;${PROJECT_BINARY_DIR}/;${CPACK_SOURCE_IGNORE_FILES}")


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)

configure_file(
   ${PROJECT_SOURCE_DIR}/versionForSourcePackage.cmake.in
   ${PROJECT_BINARY_DIR}/versionForSourcePackage.cmake @ONLY)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

set(source_files "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

#Add the version_file to the project build.
# generating executable program by specifying the location of the file
add_executable(${Project} src/main.cpp ${version_file})
set_target_properties(${Project} PROPERTIES LINK_FLAGS "-Wl,-rpath,./")

target_link_libraries(${Project} ${ProjectLib} -static-libgcc -static-libstdc++)


# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# specifying install rules for dependency libraries
# target_link_libraries(${Project} -static-libgcc -static-libstdc++)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/dependencies/lib DESTINATION .)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)


# specifying install rules for documentation 
string(TOLOWER ${Project_NAME} User_Guide)
include(InstallRequiredSystemLibraries)
install(
    TARGETS ${Project}
    )
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/api/html DESTINATION help/api)
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/usermanual/documentation/html DESTINATION help/manual)
install(FILES ${PROJECT_BINARY_DIR}/doc/api/latex/refman.pdf DESTINATION help/api RENAME "Api Documentation.pdf")
install(FILES ${PROJECT_BINARY_DIR}/doc/usermanual/documentation/latex/${User_Guide}.pdf DESTINATION help/manual RENAME "User Guide.pdf")

# CPACK ruls
string(TIMESTAMP YEAR_MONTH "%Y_%m")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
# set(CPACK_GENERATOR "RPM;DEB;TGZ;TXZ;TZ;STGZ;NSIS;NSIS64") # NSIS 3.0 and above is not available for Ubuntu 18.04
set(CPACK_GENERATOR "STGZ")
set(CPACK_PACKAGE_NAME ${Project})
set(CPACK_PACKAGE_RELEASE ${YEAR_MONTH})
set(CPACK_PACKAGE_CONTACT ${Project_DEVELOPER})
set(CPACK_PACKAGE_VENDOR "My Company")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION_SHORT}-${CPACK_PACKAGE_RELEASE}-${CMAKE_SYSTEM_PROCESSOR}")
include(CPack)